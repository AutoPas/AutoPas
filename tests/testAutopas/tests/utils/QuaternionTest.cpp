/**
 * @file QuaternionTest.cpp
 * @author S. Newcome
 * @date 16/08/2022
 */

#include <gtest/gtest.h>

#include "QuaternionTest.h"

using namespace autopas;

#define PI 3.14159265359

const std::array<char,3> axes{'x','y','z'};

/**
 * Returns normalized quaternion from direction r and angle theta.
 * @param r
 * @param theta
 * @return
 */
std::array<double, 4> returnNormalizedQuaternion(std::array<double, 3> r, double theta) {
  const auto cosThetaDiv2 = std::cos(theta * 0.5);
  const auto sinThetaDiv2 = std::sin(theta * 0.5);

  const std::array<double, 4> quaternion = {cosThetaDiv2, r[0] * sinThetaDiv2, r[1] * sinThetaDiv2, r[2] * sinThetaDiv2};

  const auto norm = std::sqrt(utils::ArrayMath::dot(quaternion, quaternion));

  return utils::ArrayMath::mulScalar(quaternion, 1/norm);
}

std::array<double, 3> returnRotationInAxes(std::array<double,3> pos, int unrotatedAxis, double theta) {
  int axisA = 0;
  int axisB = 0;
  if (unrotatedAxis == 0) {
    axisA = 1;
    axisB = 2;
  } else if (unrotatedAxis == 1) {
    axisA = 2;
    axisB = 0;
  } else {
    axisA = 0;
    axisB = 1;
  }
  std::array<double, 3> newPos;
  newPos[axisA] =  pos[axisA] * std::cos(theta) - pos[axisB] * std::sin(theta);
  newPos[axisB] =  pos[axisA] * std::sin(theta) + pos[axisB] * std::cos(theta);
  newPos[unrotatedAxis] = pos[unrotatedAxis];

  return newPos;
}

/**
 * Tests rotation by defining directionvector to be rotated about and angle of rotation.
 * From this a quaternion is definined and used to rotate vector.
 *
 * For a given vector of positions, the rotation is tested in all axes plus a direction with compenents in mujltiple axes
 * In each case, no, quarter, half, three quarter rotations are tested
 */
TEST(QuaternionTest, testRotatePosition) {
  const std::vector<std::array<double, 3>> dirVec = {{1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}};
  const std::vector<double> thetaVec = {0., PI / 2, PI, 3 * PI/2};

  const std::vector<std::array<double, 3>> posVec = {{0.,0.,0.},{1.,0.,0.},{0.5,1.,1.}, {-4,-1,0}};

  for (int axis = 0; axis < 3; ++axis) {
    for (auto theta : thetaVec) {
      for (auto pos : posVec) {
        const auto expectedRotatedPos = returnRotationInAxes(pos, axis, theta);
        const auto q = returnNormalizedQuaternion(dirVec[axis], theta);
        const auto rotatedPos = utils::quaternion::rotatePosition(q, pos);

        for (int i = 0; i < 3; ++i) {
          ASSERT_NEAR(expectedRotatedPos[i], rotatedPos[i], 1e-13)
              << "pos = {" << pos[0] << ", " << pos[1] << ", " << pos[2] << "}: Error in " << axes[i]
              << "-axis for rotation with axis " << axes[axis] << " fixed and theta = " << theta;
        }
      }
    }
  }
}