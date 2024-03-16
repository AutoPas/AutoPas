//
// Created by jan on 1/26/24.
//

#ifndef AUTOPAS_ATLOOKUPTABLE_H
#define AUTOPAS_ATLOOKUPTABLE_H

#include <cmath>
#include <vector>

#include "LookUpTableTypes.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/logging/Logger.h"

namespace ForceLookUpTable {

template<PositionType positionType, IntervalType intervalType, InterpolationType interpolationType, typename floatType = double, typename intType = unsigned long>
class ATLookUpTable {};

template <IntervalType intervalType, InterpolationType interpolationType, typename floatType, typename intType>
class ATLookUpTable<relative, intervalType, interpolationType, floatType, intType> {

  using Entry = std::pair<std::array<std::array<floatType, 3>, 3>, floatType>;

 public:
  ATLookUpTable() {
    //AutoPasLog(DEBUG, "Default constructor called.");
  }

  // list: cutoffSquared, nu, ... (numberOfPoints)
  // Extremely unreadable and user-error-prone

  ATLookUpTable(std::initializer_list<floatType> args) {
//    std::cout << "LUT created.\n";
    if (args.size() < 3) {  // Fail gracefully
      throw autopas::utils::ExceptionHandler::AutoPasException("Not enough arguments for ATLookUpTable creation");
    }

    cutoffSquared = args.begin()[0];
    nu = args.begin()[1];
    if constexpr (intervalType == evenSpacing) {
      if (args.size() != 3) {  // Fail
        return;
      }
      numberOfPoints = static_cast<intType>(args.begin()[2]);
      if (numberOfPoints == 0)
        throw autopas::utils::ExceptionHandler::AutoPasException("At least one point needed for LUT.");
      pointDistance = cutoffSquared / numberOfPoints;
      fillTableEvenSpacing();
    }
  }

  Entry retrieveValue(const std::array<double, 3>& positionI, const std::array<double, 3>& positionJ, const std::array<double, 3>& positionK, floatType distSquaredIJ, floatType distSquaredJK, floatType distSquaredKI) {
    AutoPasLog(DEBUG, "Retrieved value from AT-LUT");
    if constexpr (interpolationType == nextNeighbor) {
      return getNextNeighbor(positionI[0], positionI[1], positionI[2], positionJ[0], positionJ[1], positionJ[2], positionK[0], positionK[1], positionK[2], distSquaredIJ, distSquaredJK, distSquaredKI);
    }
  }

 private:
  std::vector<Entry> lut; // Pair of potential energy and Triplets of forces for newton3, each of which consist of 3 values fpr the 3 dimensions
  intType numberOfPoints; // For even spacing
  floatType pointDistance; // For even spacing
  floatType cutoffSquared;

  // Temporary until we know what we are doing
  floatType nu;

  // Vector Index

  size_t getIndexNoP(size_t a, size_t b, size_t c) {
    AutoPasLog(DEBUG, "For {} {} {} return index {}", a, b, c, (a * (a+1) * (a+2)) / 6 + (b * (b + 1)) / 2 + c);
    return (a * (a+1) * (a+2)) / 6 + (b * (b + 1)) / 2 + c;
  }

  // Fill functions

  void fillTableEvenSpacing () {
    // std::cout << "Building Table\n";
    // std::cout << "Number of points: " << numberOfPoints << " Point distance: " << pointDistance << "\n";
    floatType j1, k1, k2;

    for (floatType distA = pointDistance / 2; distA < cutoffSquared; distA += pointDistance) {
      for (floatType distB = pointDistance / 2; distB <= distA; distB += pointDistance) {
        for (floatType distC = pointDistance / 2; distC <= distB; distC += pointDistance) {
          floatType cX, cY;
          // Use roots because distABC are technically distanceSquared
          floatType rootA, rootB, rootC;
          rootA = std::sqrt(distA);
          rootB = std::sqrt(distB);
          rootC = std::sqrt(distC);
          // This code is now slower for no reason, but more readable, so it will remain until everything has been validated to work
          cX = (rootB * rootB - rootC * rootC + rootA * rootA) / (2 * rootA);
          cY = std::sqrt((rootB * rootB) - (cX * cX));
          Entry val = ATFunctor(0, 0, 0, rootA, 0, 0, cX, cY, 0);
          lut.push_back(val);
        }
      }
    }

//    auto first = 0;
//    auto total = 0;
//    auto nan = 0;
//    auto last = 0;
//    for (auto i=0; i<lut.size(); i++) {
//      if (!std::isnan(lut.at(i).second)) {
//        if (first == 0)
//          first = i;
//        last = i;
//        total++;
//      }
//      else {
//        nan++;
//      }
//    }
//    std::cout << "First: " << first << " Last: " << last << " Total: " << total << " Nan: " << nan << " Size: "<< lut.size() << "\n";
  };

  // Interpolation functions

  Entry getNextNeighbor(floatType i1, floatType i2, floatType i3, floatType j1, floatType j2, floatType j3, floatType k1, floatType k2, floatType k3, floatType distSquaredIJ, floatType distSquaredJK, floatType distSquaredKI) {
    AutoPasLog(DEBUG, "Input was {} {} {} | {} {} {} | {} {} {} | dist: IJ {} JK {} KI {}", i1, i2, i3, j1, j2, j3, k1, k2, k3, distSquaredIJ, distSquaredJK, distSquaredKI);

    /*
     * 1. Round distances
     * 2. Get forces from lut
     * 3. Sort points from longest to shortest
     * 4. Translate triangle so the longest point is at origin and get entry from LUT
     * 5. Find rotation of triangle
     * 6. Rotate force vectors
     * 7. Return
     */

    if constexpr (intervalType == evenSpacing) {
      floatType IJround = std::floor(distSquaredIJ / pointDistance);
      floatType JKround = std::floor(distSquaredJK / pointDistance);
      floatType KIround = std::floor(distSquaredKI / pointDistance);

      if (IJround == numberOfPoints)
        IJround--;
      if (JKround == numberOfPoints)
        JKround--;
      if (KIround == numberOfPoints)
        KIround--;

      AutoPasLog(DEBUG, "IJround: {}    JKround: {}    KIround: {}", IJround, JKround, KIround);

      // Sort points from longest to shortest

      auto I = 0;
      auto J = 0;
      auto K = 0;

      IJround >= JKround ? K++ : I++;
      JKround > KIround ? I++ : J++;
      KIround > IJround ? J++ : K++;

      AutoPasLog(DEBUG, "Points: I: {}    J: {}   K: {}", I, J, K);

      Entry forces;
      Entry ret;
      Entry final;

      if (I == 0) {
        if (J == 1) {
          ret = rotate(i1, i2, i3, j1, j2, j3, k1, k2, k3, getIndexNoP(IJround, KIround, JKround));
          final = ret;
        }
        else {
          ret = rotate(i1, i2, i3, k1, k2, k3, j1, j2, j3, getIndexNoP(KIround, IJround, JKround));
          final = {{ret.first[0], ret.first[2], ret.first[1]}, ret.second};
        }
      }
      else if (J == 0) {
        if (I == 1) {
          ret = rotate(j1, j2, j3, i1, i2, i3, k1, k2, k3, getIndexNoP(IJround, JKround, KIround));
          final = {{ret.first[1], ret.first[0], ret.first[2]}, ret.second};
        }
        else {
          ret = rotate(j1, j2, j3, k1, k2, k3, i1, i2, i3, getIndexNoP(JKround, IJround, KIround));
          final = {{ret.first[2], ret.first[0], ret.first[1]}, ret.second};
        }
      }
      else {
        if (I == 1) {
          ret = rotate(k1, k2, k3, i1, i2, i3, j1, j2, j3, getIndexNoP(KIround, JKround, IJround));
          final = {{ret.first[1], ret.first[2], ret.first[0]}, ret.second};
        }
        else {
          ret = rotate(k1, k2, k3, j1, j2, j3, i1, i2, i3, getIndexNoP(JKround, KIround, IJround));
          final = {{ret.first[2], ret.first[1], ret.first[0]}, ret.second};
        }
      }
      AutoPasLog(DEBUG, [i1, i2, i3, j1, j2, j3, k1, k2, k3, this]() -> std::string {
        Entry compareEntry = ATFunctor(i1, i2, i3, j1, j2, j3, k1, k2, k3);
        std::array<std::array<floatType, 3>, 3> compareNormalized = {norm3(compareEntry.first[0][0], compareEntry.first[0][1], compareEntry.first[0][2]), norm3(compareEntry.first[1][0], compareEntry.first[1][1], compareEntry.first[1][2]), norm3(compareEntry.first[2][0], compareEntry.first[2][1], compareEntry.first[2][2])};
        return "Perfect normalized: " + std::to_string(compareNormalized[0][0]) + " " + std::to_string(compareNormalized[0][1]) + " " + std::to_string(compareNormalized[0][2]) + " | " + std::to_string(compareNormalized[1][0]) + " " + std::to_string(compareNormalized[1][1]) + " " + std::to_string(compareNormalized[1][2]) + " | " + std::to_string(compareNormalized[2][0]) + " " + std::to_string(compareNormalized[2][1]) + " " + std::to_string(compareNormalized[2][2]);
      }());
      AutoPasLog(DEBUG, [&final, this]() -> std::string {
        std::array<std::array<floatType, 3>, 3> finalNormalized = {norm3(final.first[0][0], final.first[0][1], final.first[0][2]), norm3(final.first[1][0], final.first[1][1], final.first[1][2]), norm3(final.first[2][0], final.first[2][1], final.first[2][2])};
        return "Return  normalized: " + std::to_string(finalNormalized[0][0]) + " " + std::to_string(finalNormalized[0][1]) + " " + std::to_string(finalNormalized[0][2]) + " | " + std::to_string(finalNormalized[1][0]) + " " + std::to_string(finalNormalized[1][1]) + " " + std::to_string(finalNormalized[1][2]) + " | " + std::to_string(finalNormalized[2][0]) + " " + std::to_string(finalNormalized[2][1]) + " " + std::to_string(finalNormalized[2][2]);
      }());
      AutoPasLog(DEBUG, "Perfect value: {}", ATFunctor(i1, i2, i3, j1, j2, j3, k1, k2, k3));
      AutoPasLog(DEBUG, "Return  value: {}", final);
      AutoPasLog(DEBUG, "Relative Error: {}", relError(final, ATFunctor(i1, i2, i3, j1, j2, j3, k1, k2, k3)));

      return final;
      }
  }

  // Functor stub

  Entry ATFunctor(floatType i1, floatType i2, floatType i3, floatType j1, floatType j2, floatType j3, floatType k1, floatType k2, floatType k3) {
    using namespace autopas::utils::ArrayMath::literals;

    const std::array<floatType, 3> displacementIJ{j1 - i1, j2 - i2, j3 - i3};
    const std::array<floatType, 3> displacementJK{k1 - j1, k2 - j2, k3 - j3};
    const std::array<floatType, 3> displacementKI{i1 - k1, i2 - k2, i3 - k3};

    const double distSquaredIJ = autopas::utils::ArrayMath::dot(displacementIJ, displacementIJ);
    const double distSquaredJK = autopas::utils::ArrayMath::dot(displacementJK, displacementJK);
    const double distSquaredKI = autopas::utils::ArrayMath::dot(displacementKI, displacementKI);

    const double IJDotKI = autopas::utils::ArrayMath::dot(displacementIJ, displacementKI);
    const double IJDotJK = autopas::utils::ArrayMath::dot(displacementIJ, displacementJK);
    const double JKDotKI = autopas::utils::ArrayMath::dot(displacementJK, displacementKI);
    const double allDotProducts = IJDotKI * IJDotJK * JKDotKI;

    const double allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
    const double allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
    const double factor = 3.0 * nu / allDistsTo5;

    const auto forceIDirectionJK = displacementJK * IJDotKI * (IJDotJK - JKDotKI);  // return from LUT (* factor)
    const auto forceIDirectionIJ =
        displacementIJ * (IJDotJK * JKDotKI - distSquaredJK * distSquaredKI +
                          5.0 * allDotProducts / distSquaredIJ);  // return from LUT (* factor)
    const auto forceIDirectionKI =
        displacementKI * (-IJDotJK * JKDotKI + distSquaredIJ * distSquaredJK -
                          5.0 * allDotProducts / distSquaredKI);  // return from LUT (* factor)

    auto forceI = (forceIDirectionJK + forceIDirectionIJ + forceIDirectionKI) * factor;

    auto forceJ = forceI;
    auto forceK = forceI;

    // newton3
    const auto forceJDirectionKI = displacementKI * IJDotJK * (JKDotKI - IJDotKI);
    const auto forceJDirectionIJ = displacementIJ * (-IJDotKI * JKDotKI + distSquaredJK * distSquaredKI - 5.0 * allDotProducts / distSquaredIJ);
    const auto forceJDirectionJK = displacementJK * (IJDotKI * JKDotKI - distSquaredIJ * distSquaredKI + 5.0 * allDotProducts / distSquaredJK);
    forceJ = (forceJDirectionKI + forceJDirectionIJ + forceJDirectionJK) * factor;
    forceK = (forceI + forceJ) * (-1.0);
    const auto potentialEnergy = factor * (allDistsSquared - 3.0 * allDotProducts) / 9.0;

    return std::make_pair(std::array{forceI, forceJ, forceK}, potentialEnergy);
  }

  // Error functions for debugging

  Entry absError(Entry e, Entry acc) {
    using namespace autopas::utils::ArrayMath::literals;
    return std::make_pair(std::array{e.first[0] - acc.first[0], e.first[1] - acc.first[1], e.first[2] - acc.first[2]} ,std::abs(e.second - acc.second));
  }

  Entry relError(Entry e, Entry acc) {
    using namespace autopas::utils::ArrayMath::literals;
    return std::make_pair(std::array{(e.first[0] - acc.first[0]) / acc.first[0], (e.first[1] - acc.first[1]) / acc.first[1], (e.first[2] - acc.first[2]) / acc.first[2]},std::abs(e.second - acc.second));
  }

  std::array<floatType, 3> norm3(floatType x1, floatType x2, floatType x3) {
    floatType div = std::sqrt(x1 * x1 + x2 * x2 + x3 * x3);
    return {x1 / div, x2 / div, x3 / div};
  }

  std::array<floatType, 4> quaternionMultiply(std::array<floatType, 4> v1, std::array<floatType, 4> v2) {
    std::array<floatType, 4> ret;
    ret[0] = v1[0]*v2[0] - v1[1]*v2[1] - v1[2]*v2[2] - v1[3]*v2[3];
    ret[1] = v1[0]*v2[1] + v1[1]*v2[0] - v1[2]*v2[3] + v1[3]*v2[2];
    ret[2] = v1[0]*v2[2] + v1[1]*v2[3] + v1[2]*v2[0] - v1[3]*v2[1];
    ret[3] = v1[0]*v2[3] - v1[1]*v2[2] + v1[2]*v2[1] + v1[3]*v2[0];

    return ret;
  }

  Entry rotate(floatType a1, floatType a2, floatType a3, floatType b1, floatType b2, floatType b3, floatType c1, floatType c2, floatType c3, size_t index) {

    using namespace autopas::utils::ArrayMath;

    std::array<floatType, 4> rot1Quaternion;
    std::array<floatType, 4> rot1InverseQuaternion;
    std::array<floatType, 4>rot2Quaternion;
    std::array<floatType, 4> rot2InverseQuaternion;

    b1 -= a1;
    b2 -= a2;
    b3 -= a3;
    c1 -= a1;
    c2 -= a2;
    c2 -= a3;

    AutoPasLog(DEBUG, "B: {} {} {}    C: {} {} {}", b1, b2, b3, c1, c2, c3);

    auto targetB = norm3(b1, b2, b3);
    auto targetC = norm3(c1, c2, c3);
    std::array<floatType, 3> sourceB = {1, 0, 0};
    AutoPasLog(DEBUG, "targetB normalized: {}", targetB);
    AutoPasLog(DEBUG, "targetC normalized: {}", targetC);

    // Use B and C being (1,0,0) and (0,1,0) because that way they are already unit vectors but the angle calculations stay the same

    // Find quaternion that rotates target B to (1,0,0)
    // auto norm = std::sqrt(((1 + targetB[0]) * (1 + targetB[0])) + (targetB[2] * targetB[2]) + (targetB[1] * targetB[1]));
    //AutoPasLog(DEBUG, "Quaternion norm is {}", norm);
    // Second entry is always zero because it only depends on elements from the cross product that contain B2 and B3, which are 0
//    rot1Quaternion = {(1 + targetB[0]) / norm, 0, targetB[2] / norm, -targetB[1] / norm};
//    rot1InverseQuaternion = {rot1Quaternion[0], 0, -rot1Quaternion[2], -rot1Quaternion[3]};
    auto cross = autopas::utils::ArrayMath::cross(targetB, sourceB);
    std::array<floatType, 4> targetBQuat = {0, targetB[0], targetB[1], targetB[2]};
    std::array<floatType, 4> sourceBQuat = {0, sourceB[0], sourceB[1], sourceB[2]};
    rot1Quaternion = {(1 + dot(targetBQuat, sourceBQuat)), cross[0], cross[1], cross[2]};
    rot1Quaternion = normalize(rot1Quaternion);
    rot1InverseQuaternion = {rot1Quaternion[0], -rot1Quaternion[1], -rot1Quaternion[2], -rot1Quaternion[3]};

    // Rotate targetB for debugging purposes
    // TODO: Rotations will be optimized in the way that the quaternion representing a vector has a real part of 0, so any calculation involving quatTargetB[0] is 0 and can be removed
    std::array<floatType, 4> tempQuat = {};
    AutoPasLog(DEBUG, "{}", [targetBQuat, rot1InverseQuaternion, rot1Quaternion, this]() -> std::string {
      std::array<floatType, 4> res = quaternionMultiply(quaternionMultiply(rot1InverseQuaternion, targetBQuat), rot1Quaternion);
      return "TargetB after rotation should be x 0 0 is " + std::to_string(res[1]) + " " + std::to_string(res[2]) + " " + std::to_string(res[3]);
    }());

    // Rotate targetC
    // invQuat * quatTargetC
    // TODO: Abuse that rot1InverseQuaternion[1] ist also 0
//    tempQuat = {
//        -rot1InverseQuaternion[1]*targetC[0] - rot1InverseQuaternion[2]*targetC[1] - rot1InverseQuaternion[3]*targetC[2],
//        rot1InverseQuaternion[0]*targetC[0] - rot1InverseQuaternion[2]*targetC[2] + rot1InverseQuaternion[3]*targetC[1],
//        rot1InverseQuaternion[0]*targetC[1] + rot1InverseQuaternion[1]*targetC[2] - rot1InverseQuaternion[3]*targetC[0],
//        rot1InverseQuaternion[0]*targetC[2] - rot1InverseQuaternion[1]*targetC[1] + rot1InverseQuaternion[2]*targetC[0]
//    };
    std::array<floatType, 4> targetCQuat = {0, targetC[0], targetC[1], targetC[2]};
    tempQuat = quaternionMultiply(rot1InverseQuaternion, targetCQuat);
    // tempQuat * quat
    tempQuat = quaternionMultiply(tempQuat, rot1Quaternion);
    targetC = {tempQuat[1], tempQuat[2], tempQuat[3]};
    AutoPasLog(DEBUG, "TargetC after first rotation is {}", targetC);

    // Find 2-D transformation that rotates C onto targetC
    // Use C = (tagetC[0], sqrt(target[1]^2 + target[2]^2), 0) to ensure rotation around the xy-axis at the cost of some operations.
    // TODO: Try optimizing the 2D rotation

    auto C1 = std::sqrt(targetC[1] * targetC[1] + targetC[2] * targetC[2]);
    std::array<floatType, 3> sourceC = {targetC[0], C1, 0};

    // norm = std::sqrt(((1 + targetC[0] * targetC[0] + targetC[1] * C1) * (1 + targetC[0] * targetC[0] + targetC[1] * C1)) + ((C1*targetC[2]) * (C1*targetC[2])) + ((targetC[0] * targetC[2]) * (targetC[0] * targetC[2])) + ((targetC[0] * targetC[1] - C1 * targetC[0]) * (targetC[0] * targetC[1] - C1 * targetC[0])));
    // rot2Quaternion = {(1 + targetC[0] * targetC[0] + targetC[1] * C1) / norm, (C1*targetC[2]) / norm, -(targetC[0] * targetC[2]) / norm, (targetC[0] * targetC[1] - C1 * targetC[0]) / norm};
    cross = autopas::utils::ArrayMath::cross(sourceC, targetC);
    targetCQuat = {0, targetC[0], targetC[1], targetC[2]};
    std::array<floatType, 4> sourceCQuat = {0, sourceC[0], sourceC[1], sourceC[2]};
    rot2Quaternion = {(1 + dot(targetBQuat, sourceBQuat)), cross[0], cross[1], cross[2]};
    rot2Quaternion = normalize(rot2Quaternion);
    rot2InverseQuaternion = {rot2Quaternion[0], -rot2Quaternion[1], -rot2Quaternion[2], -rot2Quaternion[3]};

    // Rotate C for debugging purposes
    AutoPasLog(DEBUG, "quatC is {}", quaternionMultiply(quaternionMultiply(rot2InverseQuaternion, sourceCQuat), rot2Quaternion));

    // Initialize forceQuaternion and rotate one force after the other
    Entry forces = lut[index];
    Entry ret{};

    // invQuat2 * tempQuat
    // TODO: Godbolt
    for (auto i=0; i<3; i++) {
      tempQuat = {0.0, forces.first[i][0], forces.first[i][1], forces.first[i][2]};
      AutoPasLog(DEBUG, "Force quat {} before rotation", tempQuat);
//      tempQuat = {-rot2InverseQuaternion[1] * tempQuat[1] - rot2InverseQuaternion[2] * tempQuat[2] - rot2InverseQuaternion[3] * tempQuat[3],
//                  rot2InverseQuaternion[0] * tempQuat[1] - rot2InverseQuaternion[2] * tempQuat[3] + rot2InverseQuaternion[3] * tempQuat[2],
//                  rot2InverseQuaternion[0] * tempQuat[2] + rot2InverseQuaternion[1] * tempQuat[3] - rot2InverseQuaternion[3] * tempQuat[1],
//                  rot2InverseQuaternion[0] * tempQuat[3] - rot2InverseQuaternion[1] * tempQuat[2] + rot2InverseQuaternion[2] * tempQuat[1]
//      };
//      AutoPasLog(DEBUG, "After first half rotation {}: {}", i, tempQuat);
      tempQuat = quaternionMultiply(rot1InverseQuaternion, tempQuat);
      // tempQuat * Quat2
      tempQuat = quaternionMultiply(tempQuat, rot2Quaternion);

      // tempQuat now
      AutoPasLog(DEBUG, "After first rotation {}: {}", i, tempQuat);
      tempQuat[0] = 0;

      // Rotate force onto original position
      // quat1 * force

//      tempQuat = {-rot1Quaternion[1] * tempQuat[1] - rot1Quaternion[2] * tempQuat[2] - rot1Quaternion[3] * tempQuat[3],
//                  rot1Quaternion[0] * tempQuat[1] - rot1Quaternion[2] * tempQuat[3] + rot1Quaternion[3] * tempQuat[2],
//                  rot1Quaternion[0] * tempQuat[2] + rot1Quaternion[1] * tempQuat[3] - rot1Quaternion[3] * tempQuat[1],
//                  rot1Quaternion[0] * tempQuat[3] - rot1Quaternion[1] * tempQuat[2] + rot1Quaternion[2] * tempQuat[1]
//      };
      tempQuat = quaternionMultiply(rot1Quaternion, tempQuat);
      AutoPasLog(DEBUG, "After second half rotation {}: {}", i, tempQuat);
      // force * invQuat1
      tempQuat = quaternionMultiply(tempQuat, rot1InverseQuaternion);

      AutoPasLog(DEBUG, "After second rotation {}: {}", i, tempQuat);

      // Forces should be rotated now
      ret.first[i][0] = tempQuat[1];
      ret.first[i][1] = tempQuat[2];
      ret.first[i][2] = tempQuat[3];
    }
    ret.second = forces.second;
    return ret;
  }
};

}


#endif  // AUTOPAS_ATLOOKUPTABLE_H
