/**
* @file Quaternion.cpp
* @author S. Newcome
* @date 17/08/2022
*/

#include "Quaternion.h"

namespace autopas::utils::quaternion {

/**
 * Compute rotational matrix from quaternion
 * @param q Quaternion
 * @return rotational matrix
 */
std::array<double,9> calculateRotationalMatrix(const std::array<double,4> q) {
  const auto ww = q[0]*q[0];
  const auto wx = q[0]*q[1];
  const auto wy = q[0]*q[2];
  const auto wz = q[0]*q[3];
  const auto xx = q[1]*q[1];
  const auto xy = q[1]*q[2];
  const auto xz = q[1]*q[3];
  const auto yy = q[2]*q[2];
  const auto yz = q[2]*q[3];
  const auto zz = q[3]*q[3];

  std::array<double,9> rotMatrix{};
  rotMatrix[0] = ww+xx-yy-zz;
  rotMatrix[1] = 2.*(xy-wz);
  rotMatrix[2] = 2.*(xy-wz);
  rotMatrix[3] = 2.*(xy+wz);
  rotMatrix[4] = ww-xx+yy-zz;
  rotMatrix[5] = 2.*(yz-wx);
  rotMatrix[6] = 2.*(xz-wy);
  rotMatrix[7] = 2.*(yz+wx);
  rotMatrix[8] = ww-xx-yy+zz;

  return rotMatrix;
}

/**
 * Rotates a std::vector of 3D positions
 * @param q Quaternion defining rotation
 * @param positionVector std::vector of arrays of 3 doubles, defining positions
 * @return rotated positions
 * @note we could instead pre-compute the rotation matrix -
 */
std::vector<std::array<double,3>> rotateVectorOfPositions(const std::array<double,4> q, const std::vector<std::array<double,3>> positionVector) {
  const auto ww = q[0]*q[0];
  const auto wx = q[0]*q[1];
  const auto wy = q[0]*q[2];
  const auto wz = q[0]*q[3];
  const auto xx = q[1]*q[1];
  const auto xy = q[1]*q[2];
  const auto xz = q[1]*q[3];
  const auto yy = q[2]*q[2];
  const auto yz = q[2]*q[3];
  const auto zz = q[3]*q[3];

  const auto r00 = ww+xx-yy-zz;
  const auto r01 = 2.*(xy-wz);
  const auto r02 = 2.*(xz+wy);
  const auto r10 = 2.*(xy+wz);
  const auto r11 = ww-xx+yy-zz;
  const auto r12 = 2.*(yz-wx);
  const auto r20 = 2.*(xz-wy);
  const auto r21 = 2.*(yz+wx);
  const auto r22 = ww-xx-yy+zz;

  std::vector<std::array<double,3>> rotatedPositions(positionVector.size());

  for (size_t i=0; i < positionVector.size(); ++i) {
    rotatedPositions[i][0] = r00 * positionVector[i][0] + r01 * positionVector[i][1] + r02 * positionVector[i][2];
    rotatedPositions[i][1] = r10 * positionVector[i][0] + r11 * positionVector[i][1] + r12 * positionVector[i][2];
    rotatedPositions[i][2] = r20 * positionVector[i][0] + r21 * positionVector[i][1] + r22 * positionVector[i][2];
  }

  return rotatedPositions;

}

/**
 * Rotates a single 3D position
 * @param q Quaternion defining rotation
 * @param position array of 3 doubles, defining position
 * @return rotated position
 */
std::array<double,3> rotatePosition(const std::array<double,4> q, const std::array<double,3> pos) {
  // todo investigate the more efficient version discussed in wikipedia page https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
  // Alan Watt and Mark Watt (1992) Advanced Animation and Rendering Techniques: Theory and Practice
  const auto ww = q[0]*q[0];
  const auto wx = q[0]*q[1];
  const auto wy = q[0]*q[2];
  const auto wz = q[0]*q[3];
  const auto xx = q[1]*q[1];
  const auto xy = q[1]*q[2];
  const auto xz = q[1]*q[3];
  const auto yy = q[2]*q[2];
  const auto yz = q[2]*q[3];
  const auto zz = q[3]*q[3];

  const auto r00 = ww+xx-yy-zz;
  const auto r01 = 2.*(xy-wz);
  const auto r02 = 2.*(xz+wy);
  const auto r10 = 2.*(xy+wz);
  const auto r11 = ww-xx+yy-zz;
  const auto r12 = 2.*(yz-wx);
  const auto r20 = 2.*(xz-wy);
  const auto r21 = 2.*(yz+wx);
  const auto r22 = ww-xx-yy+zz;

  std::array<double,3> rotatedPosition;

  rotatedPosition[0] = r00 * pos[0] + r01 * pos[1] + r02 * pos[2];
  rotatedPosition[1] = r10 * pos[0] + r11 * pos[1] + r12 * pos[2];
  rotatedPosition[2] = r20 * pos[0] + r21 * pos[1] + r22 * pos[2];

  return rotatedPosition;
}

/**
 * Rotates a single 3D position backwards
 * @param q Quaternion defining rotation
 * @param position array of 3 doubles, defining position
 * @return rotated position
 */
std::array<double,3> rotatePositionBackwards(const std::array<double,4> q, const std::array<double,3> pos) {
  const auto ww =  q[0]*q[0];
  const auto wx = -q[0]*q[1];
  const auto wy = -q[0]*q[2];
  const auto wz = -q[0]*q[3];
  const auto xx =  q[1]*q[1];
  const auto xy =  q[1]*q[2];
  const auto xz =  q[1]*q[3];
  const auto yy =  q[2]*q[2];
  const auto yz =  q[2]*q[3];
  const auto zz =  q[3]*q[3];

  const auto r00 = ww+xx-yy-zz;
  const auto r01 = 2.*(xy-wz);
  const auto r02 = 2.*(xy-wz);
  const auto r10 = 2.*(xy+wz);
  const auto r11 = ww-xx+yy-zz;
  const auto r12 = 2.*(yz-wx);
  const auto r20 = 2.*(xz-wy);
  const auto r21 = 2.*(yz+wx);
  const auto r22 = ww-xx-yy+zz;

  std::array<double,3> rotatedPosition;

  rotatedPosition[0] = r00 * pos[0] + r01 * pos[1] + r02 * pos[2];
  rotatedPosition[1] = r10 * pos[0] + r11 * pos[1] + r12 * pos[2];
  rotatedPosition[2] = r20 * pos[0] + r21 * pos[1] + r22 * pos[2];

  return rotatedPosition;
}

/**
 * Quaternion multiplication
 * @param q1 quaternion 1
 * @param q2 quaternion 2
 * @return q1 times q2
 */
std::array<double,4> qMul(const std::array<double,4> q1, const std::array<double,4> q2) {
  return { q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3],
          q1[0]*q2[1] + q2[0]*q1[1] + q1[2]*q2[3] - q1[3]*q2[2],
          q1[0]*q2[2] + q2[0]*q1[2] + q1[3]*q2[1] - q1[1]*q2[3],
          q1[0]*q2[3] + q2[0]*q1[3] + q1[1]*q2[2] - q1[2]*q2[1] };
}

/**
 * Quaternion multiplication (converting v into a quaternion (0,v))
 * @param q quaternion
 * @param v 3D-vector
 * @return q times (0,v)
 */
std::array<double,4> qMul(const std::array<double,4> q, const std::array<double,3> v) {
  return { - q[1]*v[0] - q[2]*v[1] - q[3]*v[2],
          q[0]*v[0] + q[2]*v[2] - q[3]*v[1],
          q[0]*v[1] + q[3]*v[0] - q[1]*v[2],
          q[0]*v[2] + q[1]*v[1] - q[2]*v[0] };
}

/**
 * Quaternion multiplication (converting v into a quaternion (0,v))
 * @param v 3D-vector
 * @param q quaternion
 * @return (0,v) times q
 */
std::array<double,4> qMul(const std::array<double,3> v, const std::array<double,4> q) {
  return {- v[0]*q[1] - v[1]*q[2] - v[2]*q[3],
          q[0]*v[0] + v[1]*q[3] - v[2]*q[2],
          q[0]*v[1] + v[2]*q[1] - v[0]*q[3],
          q[0]*v[2] + v[0]*q[2] - v[1]*q[1] };
}

/**
 * Quaternion conjugation.
 * @param q quaternion
 * @return conjugated quaternion
 */
std::array<double,4> qConjugate(const std::array<double,4> q) {
  return {q[0],-q[1],-q[2],-q[3]};
}

}
