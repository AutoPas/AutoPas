/**
 * @file Quaternion.h
 * @author S. Newcome
 * @date 17/02/2022
*/

#pragma once

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
std::array<double,9> calculateRotationalMatrix(const std::array<double,4> q);

/**
 * Rotates a std::vector of 3D positions
 * @param q Quaternion defining rotation
 * @param positionVector std::vector of arrays of 3 doubles, defining positions
 * @return rotated positions
 * @note we could instead pre-compute the rotation matrix -
 */
std::vector<std::array<double,3>> rotateVectorOfPositions(const std::array<double,4> q, const std::vector<std::array<double,3>> positionVector);

/**
 * Rotates a single 3D position
 * @param q Quaternion defining rotation
 * @param position array of 3 doubles, defining position
 * @return rotated position
 */
std::array<double,3> rotatePosition(const std::array<double,4> q, const std::array<double,3> pos);

/**
 * Rotates a single 3D position backwards
 * @param q Quaternion defining rotation
 * @param position array of 3 doubles, defining position
 * @return rotated position
 */
std::array<double,3> rotatePositionBackwards(const std::array<double,4> q, const std::array<double,3> pos);

/**
 * Quaternion multiplication
 * @param q1 quaternion 1
 * @param q2 quaternion 2
 * @return q1 times q2
 */
std::array<double,4> qMul(const std::array<double,4> q1, const std::array<double,4> q2);

/**
 * Quaternion multiplication (converting v into a quaternion (0,v))
 * @param q quaternion
 * @param v 3D-vector
 * @return q times (0,v)
 */
std::array<double,4> qMul(const std::array<double,4> q, const std::array<double,3> v);

/**
 * Quaternion multiplication (converting v into a quaternion (0,v))
 * @param v 3D-vector
 * @param q quaternion
 * @return (0,v) times q
 */
std::array<double,4> qMul(const std::array<double,3> v, const std::array<double,4> q);

/**
 * Quaternion conjugation.
 * @param q quaternion
 * @return conjugated quaternion
 */
std::array<double,4> qConjugate(const std::array<double,4> q);

/**
 * Convert quaternion to 3d-vec
 * @param q quaternion
 * @return Quaternion without scalar part (i.e. 0th element)
 */
std::array<double, 3> convertQuaternionTo3DVec(const std::array<double, 4> q);
}