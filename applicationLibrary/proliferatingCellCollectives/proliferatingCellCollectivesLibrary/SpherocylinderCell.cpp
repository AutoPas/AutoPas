/**
 * @file SpherocylinderCell.cpp
 * @author Manuel Lerchner
 * @date 11/05/2025
 */

#include "SpherocylinderCell.h"

#include <algorithm>
#include <random>
#include <sstream>

#include "DCPQuery.h"
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
                                       const std::array<double, 4> &quaternion, double length0, double stress,
                                       unsigned long moleculeId)
    : mdLib::MultisiteMoleculeLJ(position, linear_velocity, quaternion, angular_velocity, moleculeId, 0),
      _length0(length0),
      _length(length0),
      _stress(stress),
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
  double angle = angle_dist(gen) / (1 + _stress);

  std::array<double, 4> dqLeft = {std::cos(angle), 0.0, 0.0, std::sin(angle)};
  std::array<double, 4> dqRight = {std::cos(-angle), 0.0, 0.0, std::sin(-angle)};
  std::array<double, 4> newOrientationLeft = autopas::utils::quaternion::qMul(getQuaternion(), dqLeft);
  std::array<double, 4> newOrientationRight = autopas::utils::quaternion::qMul(getQuaternion(), dqRight);

  SpherocylinderCell right_cell(newCenterRight, _v, _angularVel, newOrientationRight, _length0, _stress, 2 * _id + 2);
  SpherocylinderCell left_cell(newCenterLeft, _v, _angularVel, newOrientationLeft, _length0, _stress, 2 * _id + 3);

  *this = std::move(left_cell);

  return right_cell;
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

std::optional<std::tuple<double, std::array<double, 3>, std::array<double, 3>, std::array<double, 3>>>
SpherocylinderCell::getCollisionInfo(const SpherocylinderCell &other) const {
  // First rough check to see if the cells are too far apart to collide
  double diffsq = dot(_r - other._r, _r - other._r);
  if (diffsq > std::pow(_length + other._length, 2)) {
    return std::nullopt;
  }

  // Otherwise, perform the spherocylinder-spherocylinder collision check
  auto [p1, p2] = getEndpoints();
  auto [q1, q2] = other.getEndpoints();

  static DCPQuery DistSegSeg3;

  double s, t = 0;
  std::array<double, 3> closestP1, closestP2;
  auto minDist = DistSegSeg3(p1, p2, q1, q2, closestP1, closestP2, s, t);

  auto overlap = (_diameter + other._diameter) / 2.0 - minDist;

  if (overlap <= 0) {
    return std::nullopt;
  }

  auto normal = closestP1 - closestP2;
  double norm = L2Norm(normal);
  if (norm > 0) {
    normal = normal * (1.0 / norm);
  } else {
    normal = {1.0, 0.0, 0.0};
  }
  return std::make_tuple(overlap, normal, closestP1, closestP2);
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

}  // namespace pccLib
