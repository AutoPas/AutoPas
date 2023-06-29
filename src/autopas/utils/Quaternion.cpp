/**
 * @file Quaternion.cpp
 * @author S. Newcome
 * @date 17/08/2022
 */

#include "Quaternion.h"

#include "ExceptionHandler.h"

namespace autopas::utils::quaternion {

std::vector<std::array<double, 3>> rotateVectorOfPositions(const std::array<double, 4> &q,
                                                           const std::vector<std::array<double, 3>> &positionVector) {
  const auto ww = q[0] * q[0];
  const auto wx = q[0] * q[1];
  const auto wy = q[0] * q[2];
  const auto wz = q[0] * q[3];
  const auto xx = q[1] * q[1];
  const auto xy = q[1] * q[2];
  const auto xz = q[1] * q[3];
  const auto yy = q[2] * q[2];
  const auto yz = q[2] * q[3];
  const auto zz = q[3] * q[3];

  const auto r00 = ww + xx - yy - zz;
  const auto r01 = 2. * (xy - wz);
  const auto r02 = 2. * (xz + wy);
  const auto r10 = 2. * (xy + wz);
  const auto r11 = ww - xx + yy - zz;
  const auto r12 = 2. * (yz - wx);
  const auto r20 = 2. * (xz - wy);
  const auto r21 = 2. * (yz + wx);
  const auto r22 = ww - xx - yy + zz;

  std::vector<std::array<double, 3>> rotatedPositions(positionVector.size());

  for (size_t i = 0; i < positionVector.size(); ++i) {
    rotatedPositions[i][0] = r00 * positionVector[i][0] + r01 * positionVector[i][1] + r02 * positionVector[i][2];
    rotatedPositions[i][1] = r10 * positionVector[i][0] + r11 * positionVector[i][1] + r12 * positionVector[i][2];
    rotatedPositions[i][2] = r20 * positionVector[i][0] + r21 * positionVector[i][1] + r22 * positionVector[i][2];
  }

  return rotatedPositions;
}

std::array<double, 3> rotatePosition(const std::array<double, 4> &q, const std::array<double, 3> &pos) {
  // Alan Watt and Mark Watt (1992) Advanced Animation and Rendering Techniques: Theory and Practice
  const auto ww = q[0] * q[0];
  const auto wx = q[0] * q[1];
  const auto wy = q[0] * q[2];
  const auto wz = q[0] * q[3];
  const auto xx = q[1] * q[1];
  const auto xy = q[1] * q[2];
  const auto xz = q[1] * q[3];
  const auto yy = q[2] * q[2];
  const auto yz = q[2] * q[3];
  const auto zz = q[3] * q[3];

  const auto r00 = ww + xx - yy - zz;
  const auto r01 = 2. * (xy - wz);
  const auto r02 = 2. * (xz + wy);
  const auto r10 = 2. * (xy + wz);
  const auto r11 = ww - xx + yy - zz;
  const auto r12 = 2. * (yz - wx);
  const auto r20 = 2. * (xz - wy);
  const auto r21 = 2. * (yz + wx);
  const auto r22 = ww - xx - yy + zz;

  // rotated position
  return {
      r00 * pos[0] + r01 * pos[1] + r02 * pos[2],
      r10 * pos[0] + r11 * pos[1] + r12 * pos[2],
      r20 * pos[0] + r21 * pos[1] + r22 * pos[2],
  };
}

std::array<double, 3> rotatePositionBackwards(const std::array<double, 4> &q, const std::array<double, 3> &pos) {
  return rotatePosition([q]() -> std::array<double, 4> { return {q[0], -q[1], -q[2], -q[3]}; }(), pos);
}

std::array<double, 4> qMul(const std::array<double, 4> &q1, const std::array<double, 4> &q2) {
  return {q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3],
          q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2],
          q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1],
          q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0]};
}

std::array<double, 4> qMul(const std::array<double, 4> &q, const std::array<double, 3> &v) {
  return qMul(q, [v]() -> std::array<double, 4> { return {0, v[0], v[1], v[2]}; }());
}

std::array<double, 4> qMul(const std::array<double, 3> &v, const std::array<double, 4> &q) {
  return qMul([v]() -> std::array<double, 4> { return {0, v[0], v[1], v[2]}; }(), q);
}

std::array<double, 4> qConjugate(const std::array<double, 4> &q) { return {q[0], -q[1], -q[2], -q[3]}; }

std::array<double, 3> convertQuaternionTo3DVec(const std::array<double, 4> &q) {
  if (q[0] > 1e-13 or q[0] < -1e-13) {
    autopas::utils::ExceptionHandler::exception(
        "Calling convertQuaternionTo3DVec on a quaternion with non-zero scalar part!");
  }
  return {q[1], q[2], q[3]};
}

std::array<double, 4> qMirror(const std::array<double, 4> &q, const int &dimensionNormalToMirror) {
  if (dimensionNormalToMirror == 0) {
    return {q[0], q[1], -q[2], -q[3]};
  } else if (dimensionNormalToMirror == 1) {
    return {q[0], -q[1], q[2], -q[3]};
  } else if (dimensionNormalToMirror == 2) {
    return {q[0], -q[1], -q[2], q[3]};
  } else {
    autopas::utils::ExceptionHandler::exception("Calling qMirror with dimensionNormalToMirror not 0, 1, or 2!");
    return q;
  }
}

}  // namespace autopas::utils::quaternion
