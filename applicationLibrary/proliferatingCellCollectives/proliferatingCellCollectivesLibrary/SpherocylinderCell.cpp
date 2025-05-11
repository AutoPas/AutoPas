/**
 * @file SpherocylinderCell.cpp
 * @author Manuel Lerchner
 * @date 11/05/2025
 */

#include "SpherocylinderCell.h"

#include <algorithm>
#include <random>
#include <sstream>

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Quaternion.h"

using autopas::utils::ArrayMath::cross;
using autopas::utils::ArrayMath::dot;
using autopas::utils::ArrayMath::L2Norm;
using autopas::utils::ArrayMath::operator-;
using autopas::utils::ArrayMath::operator+;
using autopas::utils::ArrayMath::operator*;

namespace pccLib {

SpherocylinderCell::SpherocylinderCell(const std::array<double, 3> &position,
                                       const std::array<double, 3> &linear_velocity,
                                       const std::array<double, 3> &angular_velocity,
                                       const std::array<double, 4> &quaternion, double length0,
                                       unsigned long moleculeId)
    : mdLib::MultisiteMoleculeLJ(position, linear_velocity, quaternion, angular_velocity, moleculeId, 0),
      _length0(length0),
      _length(length0),
      _diameter(0.5 * length0),
      _maxLength(2.0 * length0) {}

void SpherocylinderCell::grow(double dt, double tao, double lamb) {
  double growth = (_length / tao) * std::exp(-lamb * _stress);
  _length += growth * dt;
  if (_length > _maxLength) _length = _maxLength;
}

std::optional<SpherocylinderCell> SpherocylinderCell::divide() {
  if (_length < _maxLength) {
    return std::nullopt;
  }

  std::array<double, 3> dir = getDirectionVector();
  std::array<double, 3> newCenterLeft = _r - dir * (0.25 * _length);
  std::array<double, 3> newCenterRight = _r + dir * (0.25 * _length);

  static thread_local std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<double> angle_dist(-M_PI / 64.0, M_PI / 64.0);
  double angle = angle_dist(gen);

  std::array<double, 4> dq = {std::cos(angle), 0.0, 0.0, std::sin(angle)};
  std::array<double, 4> newOrientationLeft = autopas::utils::quaternion::qMul(getQuaternion(), dq);
  std::array<double, 4> newOrientationRight = autopas::utils::quaternion::qMul(getQuaternion(), dq);

  SpherocylinderCell new_cell(newCenterRight, _v, getAngularVel(), getQuaternion(), _length0, _id + 1);

  _r = newCenterLeft;
  _v = _v;
  _angularVel = _angularVel;
  setQuaternion(newOrientationLeft);
  _length = _length0;
  _id = _id + 1;

  return new_cell;
}

std::array<double, 3> SpherocylinderCell::getDirectionVector() const {
  static const std::vector<std::array<double, 3>> defaultDirection = {{1.0, 0.0, 0.0}};
  return autopas::utils::quaternion::rotateVectorOfPositions(getQuaternion(), defaultDirection)[0];
}

std::pair<std::array<double, 3>, std::array<double, 3>> SpherocylinderCell::getEndpoints() const {
  double halfLength = _length / 2.0;
  std::array<double, 3> dir = getDirectionVector();
  std::array<double, 3> end1 = _r + dir * (halfLength - _diameter / 2.0);
  std::array<double, 3> end2 = _r - dir * (halfLength - _diameter / 2.0);
  return {end1, end2};
}

std::optional<std::tuple<double, std::array<double, 3>, std::array<double, 3>>> SpherocylinderCell::getCollisionInfo(
    const SpherocylinderCell &other) const {
  auto [p1, p2] = getEndpoints();
  auto [q1, q2] = other.getEndpoints();

  auto [closestP1, closestP2, minDist] =
      minimumDistanceBetweenLineSegments(p1, p2, q1, q2, true, true, true, true, true);

  auto overlap = (_diameter + other._diameter) / 2.0 - minDist;
  if (overlap <= 0) {
    return std::nullopt;
  }

  auto normal = closestP1 - closestP2;
  if (L2Norm(normal) > 0) {
    normal = normal * (1.0 / L2Norm(normal));
  } else {
    normal = {1.0, 0.0, 0.0};
  }

  auto contactPoint = (closestP1 + closestP2) * 0.5;
  return std::make_tuple(overlap, normal, contactPoint);
}

std::string SpherocylinderCell::toString() const {
  std::ostringstream oss;
  oss << "SpherocylinderCell{";
  oss << "id=" << _id;
  oss << ", r=[" << _r[0] << ", " << _r[1] << ", " << _r[2] << "]";
  oss << ", v=[" << _v[0] << ", " << _v[1] << ", " << _v[2] << "]";
  oss << ", f=[" << _f[0] << ", " << _f[1] << ", " << _f[2] << "]";
  oss << ", torque=[" << _torque[0] << ", " << _torque[1] << ", " << _torque[2] << "]";
  oss << ", length=" << _length;
  oss << ", diameter=" << _diameter;
  oss << ", stress=" << _stress;
  oss << "}";
  return oss.str();
}

std::tuple<std::array<double, 3>, std::array<double, 3>, double> minimumDistanceBetweenLineSegments(
    const std::array<double, 3> &a0, const std::array<double, 3> &a1, const std::array<double, 3> &b0,
    const std::array<double, 3> &b1, bool clampAll, bool clampA0, bool clampA1, bool clampB0, bool clampB1) {
  if (clampAll) {
    clampA0 = clampA1 = clampB0 = clampB1 = true;
  }

  std::array<double, 3> A = a1 - a0;
  std::array<double, 3> B = b1 - b0;
  double magA = L2Norm(A);
  double magB = L2Norm(B);

  std::array<double, 3> _A = (magA > 0) ? (A * (1.0 / magA)) : A;
  std::array<double, 3> _B = (magB > 0) ? (B * (1.0 / magB)) : B;

  std::array<double, 3> cross = autopas::utils::ArrayMath::cross(_A, _B);
  double denom = std::pow(L2Norm(cross), 2);

  // If lines are nearly parallel (cross product magnitude squared < 1e-10),
  // handle as a special case to avoid numerical instability
  if (denom < 1e-10) {
    double d0 = dot(_A, b0 - a0);

    if (clampA0 || clampA1 || clampB0 || clampB1) {
      double d1 = dot(_A, b1 - a0);

      if (d0 <= 0 && d1 <= 0) {
        if (clampA0 && clampB1) {
          if (std::abs(d0) < std::abs(d1)) return {a0, b0, L2Norm(a0 - b0)};
          return {a0, b1, L2Norm(a0 - b1)};
        }
      } else if (d0 >= magA && d1 >= magA) {
        if (clampA1 && clampB0) {
          if (std::abs(d0) < std::abs(d1)) return {a1, b0, L2Norm(a1 - b0)};
          return {a1, b1, L2Norm(a1 - b1)};
        }
      }
    }
    std::array<double, 3> diff = (_A * d0) + a0 - b0;
    return {a0, b0, L2Norm(diff)};
  }

  std::array<double, 3> t = b0 - a0;
  double c1 = dot(_A, _A);
  double c2 = dot(_A, _B);
  double c3 = dot(_B, _B);
  double c4 = dot(_A, t);
  double c5 = dot(_B, t);

  double den = c1 * c3 - c2 * c2;
  double t0, t1;
  if (std::abs(den) < 1e-10) {
    t0 = 0;
    t1 = (c3 > 1e-10) ? dot(_B, a0 - b0) / c3 : 0;
  } else {
    t0 = (c2 * c5 - c3 * c4) / den;
    t1 = (c1 * c5 - c2 * c4) / den;
  }

  std::array<double, 3> pA = a0 + (_A * t0);
  std::array<double, 3> pB = b0 + (_B * t1);

  // Skip clamping if no constraints are specified
  if (!(clampA0 || clampA1 || clampB0 || clampB1)) {
    return {pA, pB, L2Norm(pA - pB)};
  }

  // First pass: Clamp points directly to endpoints
  if (clampA0 && t0 < 0) pA = a0;
  if (clampA1 && t0 > magA) pA = a1;
  if (clampB0 && t1 < 0) pB = b0;
  if (clampB1 && t1 > magB) pB = b1;

  // Second pass: Project point B onto segment B after clamping A
  if ((clampA0 && t0 < 0) || (clampA1 && t0 > magA)) {
    double dotB = dot(_B, pA - b0);
    dotB = std::clamp(dotB, clampB0 ? 0.0 : -INFINITY, clampB1 ? magB : INFINITY);
    pB = b0 + (_B * dotB);
  }

  // Third pass: Project point A onto segment A after clamping B
  if ((clampB0 && t1 < 0) || (clampB1 && t1 > magB)) {
    double dotA = dot(_A, pB - a0);
    dotA = std::clamp(dotA, clampA0 ? 0.0 : -INFINITY, clampA1 ? magA : INFINITY);
    pA = a0 + (_A * dotA);
  }

  return {pA, pB, L2Norm(pA - pB)};
}

}  // namespace pccLib
