//
// Created by jan on 1/26/24.
//

#ifndef AUTOPAS_ATLOOKUPTABLE_H
#define AUTOPAS_ATLOOKUPTABLE_H

#include <cmath>
#include <vector>

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/logging/Logger.h"
#include "LookUpTableTypes.h"

namespace ForceLookUpTable {

template <IntervalType intervalType, InterpolationType interpolationType, typename floatType = double, typename intType = unsigned long>
class ATLookUpTable {

  using Entry = std::pair<std::array<std::array<floatType, 3>, 3>, floatType>;

 public:
  ATLookUpTable() {
    //AutoPasLog(DEBUG, "Default constructor called.");
  }

  // list: cutoffSquared, nu, ... (numberOfPoints)
  // Extremely unreadable and user-error-prone

  ATLookUpTable(std::initializer_list<floatType> args) {
    std::cout << "LUT created.\n";
    if (args.size() < 3) {  // Fail gracefully
      //AutoPasLog(CRITICAL, "Args only has {} elements, but needs at least 3.", args.size());
      return;
    }
    std::cout << "Size of entry: " << sizeof(Entry) << "\n";
//    auto f = ATFunctor(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0);
//    auto t = ATFunctor(0., 0., 0., 1.0, 1.0, 1.0, 2.0, 2.0, 2.0);
//    std::cout << f.second << " " << t.second << " : " << f.first[0][0] << " " << t.first[0][0] << " | " << f.first[0][1] << " " << t.first[0][1] << " | " << f.first[0][2] << " " << t.first[0][2] << " | " << f.first[1][0] << " " << t.first[1][0] << " | " << f.first[1][1] << " " << t.first[1][1] << " | " << f.first[1][2] << " " << t.first[1][2] << " | " << f.first[2][0] << " " << t.first[2][0] << " | " << f.first[2][1] << " " << t.first[2][1] << " | " << f.first[2][2] << " " << t.first[2][2] << "\n";
//    std::cout << (f.second == t.second && f.first[0][0] == t.first[0][0] && f.first[0][1] == t.first[0][1] && f.first[0][2] == t.first[0][2] && f.first[1][0] == t.first[1][0] && f.first[1][1] == t.first[1][1] && f.first[1][2] == t.first[1][2] && f.first[2][0] == t.first[2][0] && f.first[2][1] == t.first[2][1] && f.first[2][2] == t.first[2][2]) << "\n";
//    f = ATFunctor(5.0, 2.0, -4.0, 4.0, 2.0, 0.0, 3.0, 1.0, 6.0);
//    t = ATFunctor(0.0, 0.0, 0.0, -1.0, 0.0, 4.0, -2.0, -1.0, 10.0);
//    std::cout << f.second << " " << t.second << " : " << f.first[0][0] << " " << t.first[0][0] << " | " << f.first[0][1] << " " << t.first[0][1] << " | " << f.first[0][2] << " " << t.first[0][2] << " | " << f.first[1][0] << " " << t.first[1][0] << " | " << f.first[1][1] << " " << t.first[1][1] << " | " << f.first[1][2] << " " << t.first[1][2] << " | " << f.first[2][0] << " " << t.first[2][0] << " | " << f.first[2][1] << " " << t.first[2][1] << " | " << f.first[2][2] << " " << t.first[2][2] << "\n";
//    std::cout << (f.second == t.second && f.first[0][0] == t.first[0][0] && f.first[0][1] == t.first[0][1] && f.first[0][2] == t.first[0][2] && f.first[1][0] == t.first[1][0] && f.first[1][1] == t.first[1][1] && f.first[1][2] == t.first[1][2] && f.first[2][0] == t.first[2][0] && f.first[2][1] == t.first[2][1] && f.first[2][2] == t.first[2][2]) << "\n";
//    f = ATFunctor(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0);
//    t = ATFunctor(1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0);
//    std::cout << f.second << " " << t.second << " : " << f.first[0][0] << " " << t.first[0][0] << " | " << f.first[0][1] << " " << t.first[0][1] << " | " << f.first[0][2] << " " << t.first[0][2] << " | " << f.first[1][0] << " " << t.first[1][0] << " | " << f.first[1][1] << " " << t.first[1][1] << " | " << f.first[1][2] << " " << t.first[1][2] << " | " << f.first[2][0] << " " << t.first[2][0] << " | " << f.first[2][1] << " " << t.first[2][1] << " | " << f.first[2][2] << " " << t.first[2][2] << "\n";
//    std::cout << (f.second == t.second && f.first[0][0] == t.first[0][0] && f.first[0][1] == t.first[0][1] && f.first[0][2] == t.first[0][2] && f.first[1][0] == t.first[1][0] && f.first[1][1] == t.first[1][1] && f.first[1][2] == t.first[1][2] && f.first[2][0] == t.first[2][0] && f.first[2][1] == t.first[2][1] && f.first[2][2] == t.first[2][2]) << "\n";
//    f =  ATFunctor(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0);
//    t = ATFunctor(0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0);
//    std::cout << f.second << " " << t.second << " : " << f.first[0][0] << " " << t.first[0][0] << " | " << f.first[0][1] << " " << t.first[0][1] << " | " << f.first[0][2] << " " << t.first[0][2] << " | " << f.first[1][0] << " " << t.first[1][0] << " | " << f.first[1][1] << " " << t.first[1][1] << " | " << f.first[1][2] << " " << t.first[1][2] << " | " << f.first[2][0] << " " << t.first[2][0] << " | " << f.first[2][1] << " " << t.first[2][1] << " | " << f.first[2][2] << " " << t.first[2][2] << "\n";
//    std::cout << (f.second == t.second && f.first[0][0] == t.first[0][0] && f.first[0][1] == t.first[0][1] && f.first[0][2] == t.first[0][2] && f.first[1][0] == t.first[1][0] && f.first[1][1] == t.first[1][1] && f.first[1][2] == t.first[1][2] && f.first[2][0] == t.first[2][0] && f.first[2][1] == t.first[2][1] && f.first[2][2] == t.first[2][2]) << "\n";
//    f =  ATFunctor(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0);
//    t = ATFunctor(1.0, 0.0, 1.0, 2.0, 1.0, 2.0, 3.0, 2.0, 3.0);
//    std::cout << f.second << " " << t.second << " : " << f.first[0][0] << " " << t.first[0][0] << " | " << f.first[0][1] << " " << t.first[0][1] << " | " << f.first[0][2] << " " << t.first[0][2] << " | " << f.first[1][0] << " " << t.first[1][0] << " | " << f.first[1][1] << " " << t.first[1][1] << " | " << f.first[1][2] << " " << t.first[1][2] << " | " << f.first[2][0] << " " << t.first[2][0] << " | " << f.first[2][1] << " " << t.first[2][1] << " | " << f.first[2][2] << " " << t.first[2][2] << "\n";
//    std::cout << (f.second == t.second && f.first[0][0] == t.first[0][0] && f.first[0][1] == t.first[0][1] && f.first[0][2] == t.first[0][2] && f.first[1][0] == t.first[1][0] && f.first[1][1] == t.first[1][1] && f.first[1][2] == t.first[1][2] && f.first[2][0] == t.first[2][0] && f.first[2][1] == t.first[2][1] && f.first[2][2] == t.first[2][2]) << "\n";
//    f =  ATFunctor(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0);
//    t = ATFunctor(1.0, 1.0, 0.0, 2.0, 2.0, 1.0, 3.0, 3.0, 2.0);
//    std::cout << f.second << " " << t.second << " : " << f.first[0][0] << " " << t.first[0][0] << " | " << f.first[0][1] << " " << t.first[0][1] << " | " << f.first[0][2] << " " << t.first[0][2] << " | " << f.first[1][0] << " " << t.first[1][0] << " | " << f.first[1][1] << " " << t.first[1][1] << " | " << f.first[1][2] << " " << t.first[1][2] << " | " << f.first[2][0] << " " << t.first[2][0] << " | " << f.first[2][1] << " " << t.first[2][1] << " | " << f.first[2][2] << " " << t.first[2][2] << "\n";
//    std::cout << (f.second == t.second && f.first[0][0] == t.first[0][0] && f.first[0][1] == t.first[0][1] && f.first[0][2] == t.first[0][2] && f.first[1][0] == t.first[1][0] && f.first[1][1] == t.first[1][1] && f.first[1][2] == t.first[1][2] && f.first[2][0] == t.first[2][0] && f.first[2][1] == t.first[2][1] && f.first[2][2] == t.first[2][2]) << "\n";


    cutoffSquared = args.begin()[0];
    nu = args.begin()[1];
    if constexpr (intervalType == evenSpacing) {
      if (args.size() != 3) {  // Fail
        //AutoPasLog(CRITICAL, "Args has {} elements, but needs 3 for even spacing.", args.size());
        return;
      }
      numberOfPoints = static_cast<intType>(args.begin()[2]);
      if (numberOfPoints == 0)
        throw autopas::utils::ExceptionHandler::AutoPasException("At least one point needed for LUT.");
      pointDistance = cutoffSquared / numberOfPoints;
      fillTableEvenSpacing();
    }
  }

  Entry retrieveValue(const std::array<double, 3>& displacementIJ, const std::array<double, 3>& displacementJK, const std::array<double, 3>& displacementKI, floatType distSquaredIJ, floatType distSquaredJK, floatType distSquaredKI) {
    AutoPasLog(DEBUG, "Retrieved value from AT-LUT");
    if constexpr (interpolationType == nextNeighbor) {
      return getNextNeighbor(displacementIJ[0], displacementIJ[1], displacementIJ[2], displacementJK[0], displacementJK[1], displacementJK[2], displacementKI[0], displacementKI[1], displacementKI[2], distSquaredIJ, distSquaredJK, distSquaredKI);
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

  int getIndexNoP(size_t a, size_t b, size_t c) {
    // Is there something better?!
    size_t index;
    // optimize this
    for (auto ia = 1; ia <= a; ia++) {
      index += (ia * (ia + 1)) / 2;
    }
    index += (b * (b + 1)) / 2;
    index += c;
    AutoPasLog(DEBUG, "For {} {} {} return index {}", a, b, c, index);
    return lut[index];
  }



  // Fill functions

  void fillTableEvenSpacing () {
    std::cout << "Building Table\n";
    std::cout << "Number of points: " << numberOfPoints << " Point distance: " << pointDistance << "\n";
    floatType j1, k1, k2;

    for (floatType distA = pointDistance / 2; distA < cutoffSquared; distA += pointDistance) {
      for (floatType distB = pointDistance / 2; distB <= distA; distB += pointDistance) {
        for (floatType distC = pointDistance / 2; distC <= distB; distC += pointDistance) {
          floatType cX, cY;
          cX = (distB * distB - distC * distC + distA * distA) / (2 * distA);
          cY = std::sqrt(distB * distB - ((distB * distB - distC * distC + distA * distA) * (distB * distB - distC * distC + distA * distA)) / (4 * distA * distA));
          Entry val = ATFunctor(0, 0, 0, distA, 0, 0, cX, cY, 0);
          lut.push_back(val);
        }
      }
    }
    auto first = 0;
    auto total = 0;
    auto nan = 0;
    auto last = 0;
    for (auto i=0; i<lut.size(); i++) {
      if (!std::isnan(lut.at(i).second)) {
        if (first == 0)
          first = i;
        last = i;
        total++;
      }
      else {
        nan++;
      }
    }
    std::cout << "First: " << first << " Last: " << last << " Total: " << total << " Nan: " << nan << " Size: "<< lut.size() << "\n";
  };

  // Interpolation functions

  Entry getNextNeighbor(floatType i1, floatType i2, floatType i3, floatType j1, floatType j2, floatType j3, floatType k1, floatType k2, floatType k3, floatType distSquaredIJ, floatType distSquaredJK, floatType distSquaredKI) {
    //using namespace autopas::utils::ArrayMath::literals;
    //auto accurate = ATFunctor(i1, i2, i3, j1, j2, j3, k1, k2, k3);
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


    //static Entry totalRelError;
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

      // Sort points from longest to shortest

      auto I = 0;
      auto J = 0;
      auto K = 0;

      IJround >= JKround ? K++ : I++;
      JKround > KIround ? I++ : J++;
      KIround > IJround ? J++ : K++;

      Entry forces;

      // 4. Translate triangle so the longest point is at origin and get entry from LUT

      std::tuple<floatType, floatType, floatType, floatType> rot1Quaternion;
      std::tuple<floatType, floatType, floatType, floatType> rot1InverseQuaternion;
      std::tuple<floatType, floatType, floatType, floatType> rot2Quaternion;
      std::tuple<floatType, floatType, floatType, floatType> rot2InverseQuaternion;

      if (I == 0) {
        j1 -= i1;
        k1 -= i1;
        j2 -= i2;
        k2 -= i2;
        j3 -= i3;
        k3 -= i3;
        i1 = 0;
        i2 = 0;
        i3 = 0;

        // Use B and C being (1,0,0) and (0,1,0) because that way they are already unit vectors but the angle calculations stay the same

        if (J == 1) {
          forces = lut[getIndexNoP(IJround, KIround, JKround)];


          // Check if J = B or J = -B
          if (j2 == 0 && j3 == 0) {
            if (j1 > 0) {
              // Same
            }
            else {
              // Opposite -> Flip
            }
          }
          else {
            auto targetB = norm(j1, j2, j3);
            auto targetC = norm(k1, k2, k3);

            // Find quaternion that rotates target B to (1,0,0)
            auto norm = std::sqrt((1 + targetB[0] * 1 + targetB[0]) + (targetB[2] * targetB[2]) + (targetB[1] * targetB[1]));
            rot1Quaternion = {1 + targetB[0] / norm, 0, targetB[2] / norm, -targetB[1] / norm};
            rot1InverseQuaternion = {rot1Quaternion[0], 0, -rot1Quaternion[2], -rot1Quaternion[3]};

            // Rotate targetB for debugging purposes
            // Rotations are optimized in the way that the quaternion representing a vector has a real part of 0, so any calculation involving targetB[0] is 0 and can be removed
            // invQuat * targetB
            std::tuple<floatType, floatType, floatType, floatType> tempQuat = {
                -rot1InverseQuaternion[1]*targetB[1] - rot1InverseQuaternion[2]*targetB[2] - rot1InverseQuaternion[3]*targetB[3],
                rot1InverseQuaternion[0]*targetB[1] - rot1InverseQuaternion[2]*targetB[3] + rot1InverseQuaternion[3]*targetB[2],
                rot1InverseQuaternion[0]*targetB[2] + rot1InverseQuaternion[1]*targetB[3] - rot1InverseQuaternion[3]*targetB[1],
                rot1InverseQuaternion[0]*targetB[3] - rot1InverseQuaternion[1]*targetB[2] + rot1InverseQuaternion[2]*targetB[1]
            };
            // tempQuat * quat
            tempQuat = {
                tempQuat[0]*rot1Quaternion[0] - tempQuat[1]*rot1Quaternion[1] - tempQuat[2]*rot1Quaternion[2] - tempQuat[3]*rot1Quaternion[3],
                tempQuat[0]*rot1Quaternion[1] + tempQuat[1]*rot1Quaternion[0] - tempQuat[2]*rot1Quaternion[3] + tempQuat[3]*rot1Quaternion[2],
                tempQuat[0]*rot1Quaternion[2] + tempQuat[1]*rot1Quaternion[3] + tempQuat[2]*rot1Quaternion[0] - tempQuat[3]*rot1Quaternion[1],
                tempQuat[0]*rot1Quaternion[3] - tempQuat[1]*rot1Quaternion[2] + tempQuat[2]*rot1Quaternion[1] + tempQuat[3]*rot1Quaternion[0]
            };
            targetB = {tempQuat[1], tempQuat[2], tempQuat[3]};
            AutoPasLog(DEBUG, "TargetB is {}", targetB);

            // Rotate targetC
            // invQuat * targetC
            tempQuat = {
                -rot1InverseQuaternion[1]*targetC[1] - rot1InverseQuaternion[2]*targetC[2] - rot1InverseQuaternion[3]*targetC[3],
                rot1InverseQuaternion[0]*targetC[1] - rot1InverseQuaternion[2]*targetC[3] + rot1InverseQuaternion[3]*targetC[2],
                rot1InverseQuaternion[0]*targetC[2] + rot1InverseQuaternion[1]*targetC[3] - rot1InverseQuaternion[3]*targetC[1],
                rot1InverseQuaternion[0]*targetC[3] - rot1InverseQuaternion[1]*targetC[2] + rot1InverseQuaternion[2]*targetC[1]
            };
            // tempQuat * quat
            tempQuat = {
                tempQuat[0]*rot1Quaternion[0] - tempQuat[1]*rot1Quaternion[1] - tempQuat[2]*rot1Quaternion[2] - tempQuat[3]*rot1Quaternion[3],
                tempQuat[0]*rot1Quaternion[1] + tempQuat[1]*rot1Quaternion[0] - tempQuat[2]*rot1Quaternion[3] + tempQuat[3]*rot1Quaternion[2],
                tempQuat[0]*rot1Quaternion[2] + tempQuat[1]*rot1Quaternion[3] + tempQuat[2]*rot1Quaternion[0] - tempQuat[3]*rot1Quaternion[1],
                tempQuat[0]*rot1Quaternion[3] - tempQuat[1]*rot1Quaternion[2] + tempQuat[2]*rot1Quaternion[1] + tempQuat[3]*rot1Quaternion[0]
            };
            targetC = {tempQuat[1], tempQuat[2], tempQuat[3]};
            AutoPasLog(DEBUG, "TargetC is {}", targetC);

            // Find quaternion that rotates C onto targetC

            // Rotate forces by the inverse of the B rotation






          }

        }
        else {
          forces = lut[getIndexNoP(KIround, IJround, JKround)];
        }
      }
      else if (J == 0) {
        i1 -= j1;
        k1 -= j1;
        i2 -= j2;
        k2 -= j2;
        i3 -= j3;
        k3 -= j3;
        j1 = 0;
        j2 = 0;
        j3 = 0;

        if (I == 1)
          forces = lut[getIndexNoP(IJround, JKround, KIround)];
        else
          forces = lut[getIndexNoP(JKround, IJround, KIround)];
      }
      else {
        i1 -= k1;
        j1 -= k1;
        i2 -= k2;
        j2 -= k2;
        i3 -= k3;
        j3 -= k3;
        k1 = 0;
        k2 = 0;
        k3 = 0;

        if (I == 1)
          forces = lut[getIndexNoP(KIround, JKround, IJround)];
        else
          forces = lut[getIndexNoP(JKround, KIround, IJround)];
      }


      // TODO: fin


      // How slow is std::floor?
//      auto relErr = relError(ret, accurate);
//      totalRelError.first[0] += relErr.first[0];
//      totalRelError.first[1] += relErr.first[1];
//      totalRelError.first[2] += relErr.first[2];
//      totalRelError.second += relErr.second;
//      AutoPasLog(DEBUG, "Total rel Error: {}", totalRelError);
      //AutoPasLog(DEBUG, "Return {} instead of {}\nAbs: {} Rel: {}", ret, accurate, absError(ret, accurate), relError(ret, accurate));
      //AutoPasLog(DEBUG, "Used {} {} {} | {} {} {}", (pointDistance / 2) - 2*cutoffSquared + ij1 * pointDistance, (pointDistance / 2) - 2*cutoffSquared + ij2 * pointDistance, (pointDistance / 2) - 2*cutoffSquared + ij3 * pointDistance, (pointDistance / 2) - 2*cutoffSquared + ik1 * pointDistance, (pointDistance / 2) - 2*cutoffSquared + ik2 * pointDistance, (pointDistance / 2) - 2*cutoffSquared + ik3 * pointDistance);
      return forces;
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

    // factor also needs to be included in LUT
    auto forceI = (forceIDirectionJK + forceIDirectionIJ + forceIDirectionKI) * factor;

    // Easy part, what else do we need?
    // -> forceI, J, K
    // Directions?

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

  inline std::tuple<floatType, floatType, floatType> norm(floatType x1, floatType x2, floatType x3) {
    floatType div = std::sqrt(x1 * x1 + x2 * x2 + x3 * x3);
    return {x1 / div, x2 / div, x3 / div};
  }

  inline std::tuple<floatType, floatType, floatType, floatType> quaternionMultiply(std::tuple<floatType, floatType, floatType, floatType> v1, std::tuple<floatType, floatType, floatType, floatType> v2) {
    std::tuple<floatType, floatType, floatType, floatType> ret;
    ret[0] = v1[0]*v2[0] - v1[1]*v2[1] - v1[2]*v2[2] - v1[3]*v2[3];
    ret[1] = v1[0]*v2[1] + v1[1]*v2[0] - v1[2]*v2[3] + v1[3]*v2[2];
    ret[2] = v1[0]*v2[2] + v1[1]*v2[3] + v1[2]*v2[0] - v1[3]*v2[1];
    ret[3] = v1[0]*v2[3] - v1[1]*v2[2] + v1[2]*v2[1] + v1[3]*v2[0];

    return ret;
  }


};

}


#endif  // AUTOPAS_ATLOOKUPTABLE_H
