/**
 * @file Quaternion.h
 * @author S. Newcome
 * @date 17/02/2022
*/

#include <array>
#include <vector>

/**
 * Array utils specifically for handling quaternions
 */

namespace autopas::utils::quaternion {

/**
 * Compute rotational matrix from quaternion
 * @param q Quaternion
 * @return rotational matrix
 */
std::array<std::array<double,3>,3> calculateRotationalMatrix(const std::array<double,4> q) {
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

  std::array<std::array<double,3>,3> rotMatrix{};
  rotMatrix[0][0] = ww+xx-yy-zz;
  rotMatrix[0][1] = 2.*(xy-wz);
  rotMatrix[0][2] = 2.*(xy-wz);
  rotMatrix[1][0] = 2.*(xy+wz);
  rotMatrix[1][1] = ww-xx+yy-zz;
  rotMatrix[1][2] = 2.*(yz-wx);
  rotMatrix[2][0] = 2.*(xz-wy);
  rotMatrix[2][1] = 2.*(yz+wx);
  rotMatrix[2][2] = ww-xx-yy+zz;

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
  const auto r02 = 2.*(xy-wz);
  const auto r10 = 2.*(xy+wz);
  const auto r11 = ww-xx+yy-zz;
  const auto r12 = 2.*(yz-wx);
  const auto r20 = 2.*(xz-wy);
  const auto r21 = 2.*(yz+wx);
  const auto r22 = ww-xx-yy+zz;

  std::vector<std::array<double,3>> rotatedPositions(positionVector.size());

  for (int i=0; i<positionVector.size(); ++i) {
    rotatedPositions[i][0] = r00 * positionVector[i][0] + r01 * positionVector[i][1] + r02 * positionVector[i][2];
    rotatedPositions[i][1] = r10 * positionVector[i][0] + r11 * positionVector[i][1] + r12 * positionVector[i][2];
    rotatedPositions[i][2] = r20 * positionVector[i][0] + r21 * positionVector[i][1] + r22 * positionVector[i][2];
  }

  return rotatedPositions;

}

}