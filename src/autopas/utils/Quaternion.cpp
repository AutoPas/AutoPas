/**
 * @file Quaternion.cpp
 * @author S. Newcome
 * @date 17/08/2022
 */

#include "Quaternion.h"
#include "../../src/autopas/utils/ArrayMath.h"

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
  return rotatePosition({q[0], -q[1], -q[2], -q[3]}, pos);
}

std::array<double, 4> qMul(const std::array<double, 4> &q1, const std::array<double, 4> &q2) {
  return {q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3],
          q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2],
          q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1],
          q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0]};
}

std::array<double, 4> qMul(const std::array<double, 4> &q, const std::array<double, 3> &v) {
  return qMul(q, {0, v[0], v[1], v[2]});
}

std::array<double, 4> qMul(const std::array<double, 3> &v, const std::array<double, 4> &q) {
  return qMul({0, v[0], v[1], v[2]}, q);
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

std::array<double, 4> normalize(const std::array<double, 4>& q){
  double norm = std::sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  return {q[0]/norm, q[1]/norm, q[2]/norm, q[3]/norm};
}

/**
 * Returns quaternion describing rotation between the vectors u and v
 * @param u Vector1
 * @param v Vector2
 * @return Quaternion describing the rotation to get from u to v
 */
std::array<double, 4> getRotationBetweenVectors(std::array<double, 3> u, std::array<double,3> v){
  using autopas::utils::ArrayMath::normalize;
  using autopas::utils::quaternion::normalize;
  using autopas::utils::ArrayMath::dot;
  using autopas::utils::ArrayMath::L2Norm;
  using autopas::utils::ArrayMath::cross;

  //https://stackoverflow.com/questions/1171849/finding-quaternion-representing-the-rotation-from-one-vector-to-another
  //@todo (johnny) discuss with am how to get rid of this function or put it into util
  if(u[0]-v[0] < 0.001 && u[1]-v[1] < 0.001 && u[2]-v[2] < 0.001){  //if u == v   //this case will also be chosen if the site position is the same as the center of mass
    return {1, 0, 0, 0};  //no rotation
  }else if(u[0]+v[0] < 0.001 && u[1]+v[1] < 0.001 && u[2]+v[2] < 0.001) {    //if u == -v
    std::array<double, 3> orthogonal =  normalize(getAnyOrthogonalVector(u));
    return {0, orthogonal[0], orthogonal[1], orthogonal[2]};
  }

  float k_cos_theta = dot(u, v);
  float k = sqrt(L2Norm(u) * L2Norm(v));

  std::array<double, 3> orthogonal = cross(u,v);
  return normalize({k_cos_theta + k, orthogonal[0], orthogonal[1], orthogonal[2]});
}

std::array<double, 3> getAnyOrthogonalVector(std::array<double, 3> v)
{
  using autopas::utils::ArrayMath::cross;
  using autopas::utils::ArrayMath::L2Norm;

  const std::array<double, 3> hopefully_not_parallel_to_v{v[0] + 1.5, 0, 0};
  const std::array<double, 3> hopefully_orthogonal = cross(v, hopefully_not_parallel_to_v);
  if(L2Norm(hopefully_orthogonal) < 0.01){  //if hopefully_not_parallel was parallel (l2 norm might be a little overkill for that)
    return cross(v, {v[0], v[1]+1.5, v[2]});  //you can't be parallel to v + (1.5, 0,0) and v+(0,1.5,0) at the same time
  }
  return hopefully_orthogonal;  //at this point the vector is truly orthogonal
}

}  // namespace autopas::utils::quaternion
