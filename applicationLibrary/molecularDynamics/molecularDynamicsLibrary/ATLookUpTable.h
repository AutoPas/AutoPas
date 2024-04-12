/**
 * @file ATLookUpTable.h
 * @author J. Hampe
 * @date 26.1.2024
 */

#ifndef AUTOPAS_ATLOOKUPTABLE_H
#define AUTOPAS_ATLOOKUPTABLE_H

#include <cmath>
#include <vector>

#include "LookUpTableTypes.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/logging/Logger.h"
#include "autopas/utils/Timer.h"

namespace ForceLookUpTable {

template <PositionType positionType, IntervalType intervalType, InterpolationType interpolationType,
          typename floatType = double, typename intType = unsigned long>
class ATLookUpTable {};

/**
 * The Look-Up Table that only stores the possible distances between the particles and rotates the direction of the
 * forces on every request
 * @tparam intervalType How the stored support points are laid out, currently only an even spacing is supported
 * @tparam interpolationType How the values between the support points are computed, currently only a jump to the next
 * lowest neighbor is supported
 * @tparam floatType
 * @tparam intType
 */
template <IntervalType intervalType, InterpolationType interpolationType, typename floatType, typename intType>
class ATLookUpTable<relative, intervalType, interpolationType, floatType, intType> {
  using Entry = std::pair<std::array<std::array<floatType, 3>, 3>, floatType>;

 public:
  ATLookUpTable() {
    // AutoPasLog(DEBUG, "Default constructor called.");
  }

//  ATLookUpTable<relative, intervalType, interpolationType, floatType, intType>& operator=(ATLookUpTable<relative, intervalType, interpolationType, floatType, intType>&& other) noexcept {
//    if (this == &other)
//      return *this;
//
//    LUTtimer = std::move(other.LUTtimer);
//    lut = other.lut;
//    numberOfPoints = other.numberOfPoints;
//    pointDistance = other.pointDistance;
//    cutoffSquared = other.cutoffSquared;
//    nu = other.nu;
//    return *this;
//  }

  std::vector<autopas::utils::Timer>* LUTtimers;

  // list: cutoffSquared, nu, ... (numberOfPoints)
  // Extremely unreadable and user-error-prone

  /**
   * Constructor of the look-up table
   * Takes an initializer list of floats because it has to be agnostic to the different parameters needed for each
   * combination of intervalType and interpolationType
   * @param args Initializer list that (for evenSpacing and nextNeighbor) takes the form of {cutoffSquared, nu,
   * numberOfPoints}
   */
  ATLookUpTable(std::initializer_list<floatType> args, std::vector<autopas::utils::Timer>* timer) {
    LUTtimers = timer;
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
  /**
   * Retrieves a value from the lut
   * @param positionI Position of particle I
   * @param positionJ Position of particle J
   * @param positionK Position of particle K
   * @param distSquaredIJ Squared distance between I and J
   * @param distSquaredJK Squared distance between J and K
   * @param distSquaredKI Squared distance between K and I
   * @return An Entry containing the force vector for every particle, in the correct orientation
   */
  Entry retrieveValue(const std::array<double, 3> &positionI, const std::array<double, 3> &positionJ,
                      const std::array<double, 3> &positionK, floatType distSquaredIJ, floatType distSquaredJK,
                      floatType distSquaredKI) {
    AutoPasLog(DEBUG, "Retrieved value from AT-LUT");
    if constexpr (interpolationType == nextNeighbor) {
      return getNextNeighbor(positionI[0], positionI[1], positionI[2], positionJ[0], positionJ[1], positionJ[2],
                             positionK[0], positionK[1], positionK[2], distSquaredIJ, distSquaredJK, distSquaredKI);
    }
  }

 private:
  std::vector<Entry> lut;   // Pair of potential energy and Triplets of forces for newton3, each of which consist of 3
                            // values for the 3 dimensions
  intType numberOfPoints;   // For even spacing
  floatType pointDistance;  // For even spacing
  floatType cutoffSquared;

  // Temporary until we know what we are doing
  floatType nu;

  // Vector Index

  /**
   * From three 1D indices get the actual 1D index in the lookup table.
   * @param a
   * @param b
   * @param c
   * @return
   */
  size_t getLUTIndex(size_t a, size_t b, size_t c) {
    AutoPasLog(DEBUG, "For {} {} {} return index {}", a, b, c, (a * (a + 1) * (a + 2)) / 6 + (b * (b + 1)) / 2 + c);
    //TODO think about eliminating divisions
    return (a * (a + 1) * (a + 2)) / 6 + (b * (b + 1)) / 2 + c;
  }

  // Fill functions

  void fillTableEvenSpacing() {
    // std::cout << "Building Table\n";
    // std::cout << "Number of points: " << numberOfPoints << " Point distance: " << pointDistance << "\n";
    floatType j1, k1, k2;

    // TODO: explain why pointDistance / 2
    for (floatType distA = pointDistance / 2; distA < cutoffSquared; distA += pointDistance) {
      for (floatType distB = pointDistance / 2; distB <= distA; distB += pointDistance) {
        for (floatType distC = pointDistance / 2; distC <= distB; distC += pointDistance) {
          floatType cX, cY;
          // Use roots because distABC are technically distanceSquared
          const floatType rootA = std::sqrt(distA);
          // Circle equation
          cX = (distB - distC + distA) / (2 * rootA);
          cY = std::sqrt((distB) - (cX * cX));
          lut.push_back(ATFunctor(0, 0, 0, rootA, 0, 0, cX, cY, 0));
        }
      }
    }

    auto first = 0;
    auto total = 0;
    auto nan = 0;
    auto zero = 0;
    auto inf = 0;
    auto last = 0;
    for (auto i=0; i<lut.size(); i++) {
      if (std::isnan(lut.at(i).second)) {
        nan++;
        continue;
      }
      if (std::isinf(lut.at(i).second)) {
        inf++;
        continue;
      }
      if (lut.at(i).second == 0.0) {
        zero++;
        continue;
      }
      total++;
    }
    std::cout << "Size: "<< lut.size() << " Total: " << total << " Nan: " << nan << " Inf: " << inf << " Zero: " << zero << "\n";
  };

  // Interpolation functions

  /**
   * TODO: take std::arrays instead of individual values
   */
  Entry getNextNeighbor(floatType i1, floatType i2, floatType i3, floatType j1, floatType j2, floatType j3,
                        floatType k1, floatType k2, floatType k3, floatType distSquaredIJ, floatType distSquaredJK,
                        floatType distSquaredKI) {
    //(*LUTtimers)[0].start(); // Timer 2; Timer
    AutoPasLog(CRITICAL, "Input was {} {} {} | {} {} {} | {} {} {} | dist: IJ {} JK {} KI {}", i1, i2, i3, j1, j2, j3, k1,
               k2, k3, distSquaredIJ, distSquaredJK, distSquaredKI);

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
      // calculate 1D LUT index per dimension
      const auto calc1DLUTIndex = [&] (const auto distSquared) {
        // Edge case: any distance is exactly the longest distance in the table (=cutoff)
        if (distSquared >= cutoffSquared) {
          return static_cast<floatType>(numberOfPoints - 1);
        } else {
          return std::floor(distSquared / pointDistance);
        }
      };
      const auto IJ1DIndex = calc1DLUTIndex(distSquaredIJ);
      const auto JK1DIndex = calc1DLUTIndex(distSquaredJK);
      const auto KI1DIndex = calc1DLUTIndex(distSquaredKI);

      AutoPasLog(DEBUG, "IJround: {}    JKround: {}    KIround: {}", IJ1DIndex, JK1DIndex, KI1DIndex);

      // Rank the points by distance.
      // 0 touches the two longest sides
      // 1 touches the longest and shortest side
      // 2 touches the medium and shortest side
      uint8_t distRankI = 0;
      uint8_t distRankJ = 0;
      uint8_t distRankK = 0;

      IJ1DIndex >= JK1DIndex ? distRankK++ : distRankI++;
      JK1DIndex > KI1DIndex ? distRankI++ : distRankJ++;
      KI1DIndex > IJ1DIndex ? distRankJ++ : distRankK++;

      AutoPasLog(DEBUG, "Points: distRankI: {}    distRankJ: {}    distRankK: {}", distRankI, distRankJ, distRankK);

      Entry forces;
      Entry ret;
      Entry final;

      if (distRankI == 0) {
        if (distRankJ == 1) {
          // I < J < K
          ret = rotate(i1, i2, i3, j1, j2, j3, k1, k2, k3, getLUTIndex(IJ1DIndex, KI1DIndex, JK1DIndex));
          final = ret;
        } else {
          // I < K < J
          ret = rotate(i1, i2, i3, k1, k2, k3, j1, j2, j3, getLUTIndex(KI1DIndex, IJ1DIndex, JK1DIndex));
          final = {{ret.first[0], ret.first[2], ret.first[1]}, ret.second};
        }
      } else if (distRankJ == 0) {
        if (distRankI == 1) {
          // J < I < K
          ret = rotate(j1, j2, j3, i1, i2, i3, k1, k2, k3, getLUTIndex(IJ1DIndex, JK1DIndex, KI1DIndex));
          final = {{ret.first[1], ret.first[0], ret.first[2]}, ret.second};
        } else {
          // J < K < I
          ret = rotate(j1, j2, j3, k1, k2, k3, i1, i2, i3, getLUTIndex(JK1DIndex, IJ1DIndex, KI1DIndex));
          final = {{ret.first[2], ret.first[0], ret.first[1]}, ret.second};
        }
      } else {
        if (distRankI == 1) {
          // K < I < J
          ret = rotate(k1, k2, k3, i1, i2, i3, j1, j2, j3, getLUTIndex(KI1DIndex, JK1DIndex, IJ1DIndex));
          final = {{ret.first[1], ret.first[2], ret.first[0]}, ret.second};
        } else {
          // K < J < I
          ret = rotate(k1, k2, k3, j1, j2, j3, i1, i2, i3, getLUTIndex(JK1DIndex, KI1DIndex, IJ1DIndex));
          final = {{ret.first[2], ret.first[1], ret.first[0]}, ret.second};
        }
      }
      AutoPasLog(CRITICAL, [&]() -> std::string {
        Entry compareEntry = ATFunctor(i1, i2, i3, j1, j2, j3, k1, k2, k3);
        std::array<std::array<floatType, 3>, 3> compareNormalized = {
            norm3(compareEntry.first[0][0], compareEntry.first[0][1], compareEntry.first[0][2]),
            norm3(compareEntry.first[1][0], compareEntry.first[1][1], compareEntry.first[1][2]),
            norm3(compareEntry.first[2][0], compareEntry.first[2][1], compareEntry.first[2][2])};
        std::array<std::array<floatType, 3>, 3> finalNormalized = {
            norm3(final.first[0][0], final.first[0][1], final.first[0][2]),
            norm3(final.first[1][0], final.first[1][1], final.first[1][2]),
            norm3(final.first[2][0], final.first[2][1], final.first[2][2])};
        return "Diffs:\nIJ: " + std::to_string(i1-j1) + " " + std::to_string(i2-j2) + " " + std::to_string(i3-j3) + " | \n"
               + "IK: " + std::to_string(i1-k1) + " " + std::to_string(i2-k2) + " " + std::to_string(i3-k3) + " | \n"
               + "JK: " + std::to_string(j1-k1) + " " + std::to_string(j2-k2) + " " + std::to_string(j3-k3) + "\n"
               + "Result\n"
               + "    " + std::to_string(finalNormalized[0][0] - compareNormalized[0][0])
               + " " + std::to_string(finalNormalized[0][1] - compareNormalized[0][1])
               + " " + std::to_string(finalNormalized[0][2] - compareNormalized[0][2]) + "\n"
               + "    " + std::to_string(finalNormalized[1][0] - compareNormalized[1][0])
               + " " + std::to_string(finalNormalized[1][1] - compareNormalized[1][1])
               + " " + std::to_string(finalNormalized[1][2] - compareNormalized[1][2]) + "\n"
               + "    " + std::to_string(finalNormalized[2][0] - compareNormalized[2][0])
               + " " + std::to_string(finalNormalized[2][1] - compareNormalized[2][1])
               + " " + std::to_string(finalNormalized[2][2] - compareNormalized[2][2]);
      }());
      AutoPasLog(CRITICAL, [&]() -> std::string {
        Entry compareEntry = ATFunctor(i1, i2, i3, j1, j2, j3, k1, k2, k3);
        std::array<std::array<floatType, 3>, 3> compareNormalized = {
            norm3(compareEntry.first[0][0], compareEntry.first[0][1], compareEntry.first[0][2]),
            norm3(compareEntry.first[1][0], compareEntry.first[1][1], compareEntry.first[1][2]),
            norm3(compareEntry.first[2][0], compareEntry.first[2][1], compareEntry.first[2][2])};
        return "Perfect normalized: " + std::to_string(compareNormalized[0][0]) + " " +
               std::to_string(compareNormalized[0][1]) + " " + std::to_string(compareNormalized[0][2]) + " | " +
               std::to_string(compareNormalized[1][0]) + " " + std::to_string(compareNormalized[1][1]) + " " +
               std::to_string(compareNormalized[1][2]) + " | " + std::to_string(compareNormalized[2][0]) + " " +
               std::to_string(compareNormalized[2][1]) + " " + std::to_string(compareNormalized[2][2]);
      }());
      AutoPasLog(CRITICAL, [&]() -> std::string {
        std::array<std::array<floatType, 3>, 3> finalNormalized = {
            norm3(final.first[0][0], final.first[0][1], final.first[0][2]),
            norm3(final.first[1][0], final.first[1][1], final.first[1][2]),
            norm3(final.first[2][0], final.first[2][1], final.first[2][2])};
        return "Return normalized: " + std::to_string(finalNormalized[0][0]) + " " +
               std::to_string(finalNormalized[0][1]) + " " + std::to_string(finalNormalized[0][2]) + " | " +
               std::to_string(finalNormalized[1][0]) + " " + std::to_string(finalNormalized[1][1]) + " " +
               std::to_string(finalNormalized[1][2]) + " | " + std::to_string(finalNormalized[2][0]) + " " +
               std::to_string(finalNormalized[2][1]) + " " + std::to_string(finalNormalized[2][2]);
      }());
      // TODO: move calculations into logger macros
      using namespace autopas::utils::ArrayMath;
      const auto expectedEntry = ATFunctor(i1, i2, i3, j1, j2, j3, k1, k2, k3);
      std::array<double ,3> expectedVecLengths = {
          L2Norm(expectedEntry.first[0]),
          L2Norm(expectedEntry.first[1]),
          L2Norm(expectedEntry.first[2]),
      };
      std::array<double ,3> retVecLengths = {
          L2Norm(final.first[0]),
          L2Norm(final.first[1]),
          L2Norm(final.first[2]),
      };
      AutoPasLog(DEBUG, "\n"
                 "Perfect value: {}\n"
                 "Return  value: {}\n"
                 "Absolute Error: {}\n"
                 "Relative Error: {}\n"
                 "Perfect value vector lengths: {}\n"
                 "Return  value vector lengths: {}\n",
                 expectedEntry,
                 final,
                 absError(final, expectedEntry),
                 relError(final, expectedEntry),
                 expectedVecLengths,
                 retVecLengths
      );


      //(*LUTtimers)[0].stop(); // Timer 2 stop
      return final;
    }
  }

  // Functor stub

  Entry ATFunctor(floatType i1, floatType i2, floatType i3, floatType j1, floatType j2, floatType j3, floatType k1,
                  floatType k2, floatType k3) {
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
    const auto forceJDirectionIJ =
        displacementIJ * (-IJDotKI * JKDotKI + distSquaredJK * distSquaredKI - 5.0 * allDotProducts / distSquaredIJ);
    const auto forceJDirectionJK =
        displacementJK * (IJDotKI * JKDotKI - distSquaredIJ * distSquaredKI + 5.0 * allDotProducts / distSquaredJK);
    forceJ = (forceJDirectionKI + forceJDirectionIJ + forceJDirectionJK) * factor;
    forceK = (forceI + forceJ) * (-1.0);
    const auto potentialEnergy = factor * (allDistsSquared - 3.0 * allDotProducts) / 9.0;

    return std::make_pair(std::array{forceI, forceJ, forceK}, potentialEnergy);
  }

  // Error functions for debugging

  Entry absError(Entry e, Entry acc) {
    using namespace autopas::utils::ArrayMath::literals;
    return std::make_pair(std::array{e.first[0] - acc.first[0], e.first[1] - acc.first[1], e.first[2] - acc.first[2]},
                          std::abs(e.second - acc.second));
  }

  Entry relError(Entry e, Entry acc) {
    using namespace autopas::utils::ArrayMath::literals;
    return std::make_pair(
        std::array{(e.first[0] - acc.first[0]) / acc.first[0], (e.first[1] - acc.first[1]) / acc.first[1],
                   (e.first[2] - acc.first[2]) / acc.first[2]},
        std::abs(e.second - acc.second));
  }

  // TODO replace this with autopas::utils::ArrayMath::normalize
  std::array<floatType, 3> norm3(floatType x1, floatType x2, floatType x3) {
    floatType div =
        std::sqrt(x1 * x1 + x2 * x2 +
                  x3 * x3);  // sqrt(x1*x1 + x2*x2 + x3*x3 + 3*eps*x1*x1 + 3*eps*x2*x2 + 2*eps*x3*x3) * (1+eps)
    return {x1 / div, x2 / div,
            x3 / div};  // { (x1 / sqrt(x1*x1 + x2*x2 + x3*x3 + 3*eps*x1*x1 + 3*eps*x2*x2 + 2*eps*x3*x3) * (1+eps)) *
                        // (1+eps) , (x2 / sqrt(x1*x1 + x2*x2 + x3*x3 + 3*eps*x1*x1 + 3*eps*x2*x2 + 2*eps*x3*x3) *
                        // (1+eps)) * (1+eps), (x3 / sqrt(x1*x1 + x2*x2 + x3*x3 + 3*eps*x1*x1 + 3*eps*x2*x2 +
                        // 2*eps*x3*x3) * (1+eps)) * (1+eps) }
  }

  std::array<floatType, 4> quaternionMultiply(std::array<floatType, 4> v1, std::array<floatType, 4> v2) {
    //(*LUTtimers)[0].start(); // Timer 9 start
    std::array<floatType, 4> ret;
    ret[0] = v1[0] * v2[0] - v1[1] * v2[1] - v1[2] * v2[2] - v1[3] * v2[3];
    ret[1] = v1[0] * v2[1] + v1[1] * v2[0] - v1[2] * v2[3] + v1[3] * v2[2];
    ret[2] = v1[0] * v2[2] + v1[1] * v2[3] + v1[2] * v2[0] - v1[3] * v2[1];
    ret[3] = v1[0] * v2[3] - v1[1] * v2[2] + v1[2] * v2[1] + v1[3] * v2[0];
    //(*LUTtimers)[0].stop(); // Timer 9 stop
    return ret;
  }

  /**
   * Assumption: A touches the longest sides, B the longest and the shortest, C the medium and shortest.
   * TODO: take std::arrays instead of individual values
   * @param a1
   * @param a2
   * @param a3
   * @param b1
   * @param b2
   * @param b3
   * @param c1
   * @param c2
   * @param c3
   * @param index
   * @return
   */
  Entry rotate(floatType a1, floatType a2, floatType a3, floatType b1, floatType b2, floatType b3, floatType c1,
               floatType c2, floatType c3, size_t index) {
    using namespace autopas::utils::ArrayMath;

    // Move all coordinates so that A is at 0,0,0
    b1 -= a1;  // (b1-a1)*(1 + eps)
    b2 -= a2;  // (b2-a2)*(1 + eps)
    b3 -= a3;  // (b3-a3)*(1 + eps)
    c1 -= a1;  // (c1-a1)*(1 + eps)
    c2 -= a2;  // (c2-a2)*(1 + eps)
    c2 -= a3;  // (c3-a3)*(1 + eps)

    AutoPasLog(DEBUG, "B: {} {} {}    C: {} {} {}", b1, b2, b3, c1, c2, c3);

    (*LUTtimers)[1].start(); // Timer norm3 start
    const auto bNorm = norm3(b1, b2, b3);  // ((eps + 1)^3 (b1 - a1))/sqrt((eps + 1)^2 (3 eps (a1 - b1)^2 + (a1 - b1)^2 + 3 eps (a2 - b2)^2 + (a2 - b2)^2 + 2 eps (a3 - b3)^2 + (a3 - b3)^2))
    // Can maybe remove
    const auto cNorm = norm3(c1, c2, c3);  // ((eps + 1)^3 (c1 - a1))/sqrt((eps + 1)^2 (3 eps (a1 - c1)^2 + (a1 - c1)^2 + 3 eps (a2 - c2)^2 + (a2 - c2)^2 + 2 eps (a3 - c3)^2 + (a3 - c3)^2))
    (*LUTtimers)[1].stop(); // Timer norm3 stop
    const std::array<floatType, 3> xAxis = {1., 0., 0.};
    AutoPasLog(DEBUG, "targetB normalized: {}", bNorm);
    AutoPasLog(DEBUG, "targetC normalized: {}", cNorm);

    // Use B and C being (1,0,0) and (0,1,0) because that way they are already unit vectors but the angle calculations
    // stay the same

    // Find quaternion that rotates B to (1,0,0)
    (*LUTtimers)[2].start(); // Timer rot1Quat start
    const std::array<floatType, 3> crossB = cross(bNorm, xAxis);
    const std::array<floatType, 4> rot1Quaternion = {1. + dot(bNorm, xAxis), crossB[0], crossB[1], crossB[2]};
    const std::array<floatType, 4> rot1QuaternionNormalized = normalize(rot1Quaternion);
    const std::array<floatType, 4> rot1InverseQuaternionNormalized = {rot1QuaternionNormalized[0], -rot1QuaternionNormalized[1], -rot1QuaternionNormalized[2], -rot1QuaternionNormalized[3]};
    (*LUTtimers)[2].stop(); // Timer rot1Quat stop

    //    Sanity check: if we actually rotate B it should be at 1,0,0 because it is normalized
    AutoPasLog(DEBUG, "{}", [&]() -> std::string {
      const std::array<floatType, 4> targetBQuat = {0, bNorm[0], bNorm[1], bNorm[2]};
      const std::array<floatType, 4> sourceBQuat = {0, xAxis[0], xAxis[1], xAxis[2]};
      std::array<floatType, 4> res =
          quaternionMultiply(quaternionMultiply(rot1InverseQuaternionNormalized, targetBQuat), rot1QuaternionNormalized);
      return "TargetB after rotation should be x 0 0 is " + std::to_string(res[1]) + " " + std::to_string(res[2]) +
             " " + std::to_string(res[3]);
    }());

    // Rotate C with previously calculated rotation
    const std::array<floatType, 4> targetCQuat = {0., cNorm[0], cNorm[1], cNorm[2]};
    (*LUTtimers)[0].start(); // Timer 4 start
    const std::array<floatType, 4> targetCQuatRotated = quaternionMultiply(quaternionMultiply(rot1InverseQuaternionNormalized, targetCQuat), rot1QuaternionNormalized);
    (*LUTtimers)[0].stop(); // Timer 4 stop
    AutoPasLog(DEBUG, "TargetC after first rotation is {}", targetCQuatRotated);

    // Find 2-D transformation that rotates C onto targetC
    // TODO: Try optimizing the 2D rotation
    (*LUTtimers)[3].start(); // Timer rot2Quat start
    // Project C onto yz plane because we want to exactly rotate around the x-axis
    const std::array<floatType, 3> targetC2DNormed = norm3(0., targetCQuatRotated[2], targetCQuatRotated[3]);
    const std::array<floatType, 3> yAxis = {0., 1., 0.};

    const std::array<floatType, 3> crossC = cross(yAxis, targetC2DNormed);
    const std::array<floatType, 4> rot2Quaternion = {(1. + dot(yAxis, targetC2DNormed)), crossC[0], crossC[1], crossC[2]};
    const std::array<floatType, 4> rot2QuaternionNormalized = normalize(rot2Quaternion);
    const std::array<floatType, 4> rot2InverseQuaternionNormalized = {rot2QuaternionNormalized[0], -rot2QuaternionNormalized[1], -rot2QuaternionNormalized[2], -rot2QuaternionNormalized[3]};
    (*LUTtimers)[3].stop(); // Timer rot2Quat stop

    AutoPasLog(DEBUG, "SourceC: {} | TargetC: {}", yAxis, targetCQuatRotated);
    AutoPasLog(DEBUG, "C cross: {}", crossC);
    // Rotate C for debugging purposes
    AutoPasLog(DEBUG, "sourceCQuat rotated is {}",
               quaternionMultiply(quaternionMultiply(rot2InverseQuaternionNormalized, {0, yAxis[0], yAxis[1], yAxis[2]}), rot2QuaternionNormalized));
    AutoPasLog(DEBUG, "targetC rotated is {}", quaternionMultiply(quaternionMultiply(rot2QuaternionNormalized, {0., 0., targetC2DNormed[1], targetC2DNormed[2]}), rot2InverseQuaternionNormalized));

    // fetch values from table
    std::array<std::array<floatType, 3>, 3> forces{};
    floatType potentialEnergy{};
    std::tie(forces, potentialEnergy) = lut[index];

    AutoPasLog(DEBUG, "rot1Quat: {}", rot1QuaternionNormalized);
    AutoPasLog(DEBUG, "rot1InvQuat: {}", rot1InverseQuaternionNormalized);
    AutoPasLog(DEBUG, "rot2Quat: {}", rot2QuaternionNormalized);
    AutoPasLog(DEBUG, "rot2InvQuat: {}", rot2InverseQuaternionNormalized);

    // invQuat2 * tempQuat
    // TODO: Godbolt
    // Initialize forceQuaternion and rotate one force after the other
    Entry ret{};
    // For all force vectors from the table, undo both previous rotations
    for (size_t i = 0; i < 3; i++) {
      const std::array<floatType, 4> forceQuat = {0.0, forces[i][0], forces[i][1], forces[i][2]};
      AutoPasLog(DEBUG, "Force quat {} before rotation", forceQuat);
      (*LUTtimers)[0].start(); // Timer 4 start
      const std::array<floatType, 4> forceQuatFirstRot = quaternionMultiply(quaternionMultiply(rot2InverseQuaternionNormalized, forceQuat), rot2QuaternionNormalized);
      (*LUTtimers)[0].stop(); // Timer 4 stop

      // tempQuat now
      AutoPasLog(DEBUG, "After first rotation {}: {}", i, forceQuat);
      AutoPasLog(DEBUG, "Norm after first rotation {}: {}", i, normalize(forceQuat));
      //Could set quat[0] to 0 because it should be

      (*LUTtimers)[0].start(); // Timer 4 start
      // force * invQuat1
      const std::array<floatType, 4> forceQuatSecondRot = quaternionMultiply(quaternionMultiply(rot1QuaternionNormalized, forceQuatFirstRot), rot1InverseQuaternionNormalized);
      (*LUTtimers)[0].stop(); // Timer 4 stop

      AutoPasLog(DEBUG, "After second rotation {}: {}", i, forceQuatSecondRot);

      // Forces should be rotated now
      ret.first[i][0] = forceQuatSecondRot[1];
      ret.first[i][1] = forceQuatSecondRot[2];
      ret.first[i][2] = forceQuatSecondRot[3];
    }
    ret.second = potentialEnergy;
    //(*LUTtimers)[0].stop(); // Timer 3 stop; Timer 1 stop
    AutoPasLog(DEBUG, "After first rotation normalized:\n{}", [&]() -> std::string {
      std::array<floatType, 4> forces1Quat = {0.0, forces[0][0], forces[0][1], forces[0][2]};
      std::array<floatType, 4> forces2Quat = {0.0, forces[1][0], forces[1][1], forces[1][2]};
      std::array<floatType, 4> forces3Quat = {0.0, forces[2][0], forces[2][1], forces[2][2]};
      Entry perfect = ATFunctor(2.703957905911746, -0.0026572954522118756, 6.610864923489291, 3.3706888660320056, 0.3937881695307143, 7.713014487295848, 2.031792390399456, -1.60843618566604, 7.721119056770151);
      std::array<floatType, 4> perfect1Quat = {0.0, perfect.first[0][0], perfect.first[0][1], perfect.first[0][2]};
      std::array<floatType, 4> perfect2Quat = {0.0, perfect.first[1][0], perfect.first[1][1], perfect.first[1][2]};
      std::array<floatType, 4> perfect3Quat = {0.0, perfect.first[2][0], perfect.first[2][1], perfect.first[2][2]};
      forces1Quat = quaternionMultiply(quaternionMultiply(rot2InverseQuaternionNormalized, forces1Quat), rot2QuaternionNormalized);
      forces2Quat = quaternionMultiply(quaternionMultiply(rot2InverseQuaternionNormalized, forces2Quat), rot2QuaternionNormalized);
      forces3Quat = quaternionMultiply(quaternionMultiply(rot2InverseQuaternionNormalized, forces3Quat), rot2QuaternionNormalized);

      perfect1Quat = quaternionMultiply(quaternionMultiply(rot1InverseQuaternionNormalized, perfect1Quat), rot1QuaternionNormalized);
      perfect2Quat = quaternionMultiply(quaternionMultiply(rot1InverseQuaternionNormalized, perfect2Quat), rot1QuaternionNormalized);
      perfect3Quat = quaternionMultiply(quaternionMultiply(rot1InverseQuaternionNormalized, perfect3Quat), rot1QuaternionNormalized);
      return "Perfect is " + std::to_string(perfect.first[0][0]) + " " + std::to_string(perfect.first[0][1]) + " " + std::to_string(perfect.first[0][2]) + " | "
             + std::to_string(perfect.first[1][0]) + " " + std::to_string(perfect.first[1][1]) + " " + std::to_string(perfect.first[1][2]) + " | "
             + std::to_string(perfect.first[2][0]) + " " + std::to_string(perfect.first[2][1]) + " " + std::to_string(perfect.first[2][2]) + "\n"
             + "Forces: " + std::to_string(forces1Quat[0]) + " " + std::to_string(forces1Quat[1]) + " " + std::to_string(forces1Quat[2]) + " " + std::to_string(forces1Quat[3]) + " | "
                        + std::to_string(forces2Quat[0]) + " " + std::to_string(forces2Quat[1]) + " " + std::to_string(forces2Quat[2]) + " " + std::to_string(forces2Quat[3]) + " | "
                        + std::to_string(forces3Quat[0]) + " " + std::to_string(forces3Quat[1]) + " " + std::to_string(forces3Quat[2]) + " " + std::to_string(forces3Quat[3]) + "\nPerfect: "
                        + std::to_string(perfect1Quat[0]) + " " + std::to_string(perfect1Quat[1]) + " " + std::to_string(perfect1Quat[2]) + " " + std::to_string(perfect1Quat[3]) + " | "
                        + std::to_string(perfect2Quat[0]) + " " + std::to_string(perfect2Quat[1]) + " " + std::to_string(perfect2Quat[2]) + " " + std::to_string(perfect2Quat[3]) + " | "
                        + std::to_string(perfect3Quat[0]) + " " + std::to_string(perfect3Quat[1]) + " " + std::to_string(perfect3Quat[2]) + " " + std::to_string(perfect3Quat[3]);
    }());
    return ret;
  }
};

/**
 * The Look-Up table that stores the force values based on a number of possible arrangements of the particles to each
 * other in 3D space
 * @tparam intervalType How the stored support points are laid out, currently only an even spacing is supported
 * @tparam interpolationType How the values between the support points are computed, currently only a jump to the next
 * lowest neighbor is supported
 * @tparam floatType
 * @tparam intType
 */
template <IntervalType intervalType, InterpolationType interpolationType, typename floatType, typename intType>
class ATLookUpTable<absolute, intervalType, interpolationType, floatType, intType> {

 using Entry = std::pair<std::array<std::array<floatType, 3>, 3>, floatType>;

 public:
  ATLookUpTable() {
    // AutoPasLog(DEBUG, "Default constructor called.");
  }

  // list: cutoffSquared, nu, ... (numberOfPoints)
  // Extremely unreadable and user-error-prone

  /**
   * Constructor of the look-up table
   * Takes an initializer list of floats because it has to be agnostic to the different parameters needed for each
   * combination of intervalType and interpolationType
   * @param args Initializer list that (for evenSpacing and nextNeighbor) takes the form of {cutoffSquared, nu,
   * numberOfPoints}
   */
  ATLookUpTable(std::initializer_list<floatType> args) {
    std::cout << "LUT created.\n";
    if (args.size() < 3) {  // Fail gracefully
      // AutoPasLog(CRITICAL, "Args only has {} elements, but needs at least 3.", args.size());
      return;
    }
    // std::cout << "Size of entry: " << sizeof(Entry) << "\n";

    cutoffSquared = args.begin()[0];
    nu = args.begin()[1];
    if constexpr (intervalType == evenSpacing) {
      if (args.size() != 3) {  // Fail
        // AutoPasLog(CRITICAL, "Args has {} elements, but needs 3 for even spacing.", args.size());
        return;
      }
      numberOfPoints = static_cast<intType>(args.begin()[2]);
      if (numberOfPoints == 0)
        throw autopas::utils::ExceptionHandler::AutoPasException("At least one point needed for LUT.");
      pointDistance = (cutoffSquared * 4) / numberOfPoints;
      invPointDistance = 1 / pointDistance;
      two = numberOfPoints;
      three = two * numberOfPoints;
      four = three * numberOfPoints;
      five = four * numberOfPoints;
      six = five * numberOfPoints;
      fillTableEvenSpacing();
    }
  }

  /**
   * Retrieves a value from the lut
   * @param displacementIJ
   * @param displacementJK
   * @param displacementKI
   * @return An Entry containing the retrieved corresponding forces, as well as the global potential component
   */
  Entry retrieveValue(const std::array<double, 3> &displacementIJ, const std::array<double, 3> &displacementJK,
                      const std::array<double, 3> &displacementKI) {
    AutoPasLog(DEBUG, "Retrieved value from AT-LUT");
    if constexpr (interpolationType == nextNeighbor) {
      return getNextNeighbor(displacementIJ[0], displacementIJ[1], displacementIJ[2], displacementJK[0],
                             displacementJK[1], displacementJK[2], displacementKI[0], displacementKI[1],
                             displacementKI[2]);
    }
  }

 private:
  std::vector<Entry> lut;   // Pair of potential energy and Triplets of forces for newton3, each of which consist of 3
                            // values fpr the 3 dimensions
  intType numberOfPoints;   // For even spacing
  floatType pointDistance;  // For even spacing
  floatType invPointDistance;
  floatType cutoffSquared;
  size_t two;
  size_t three;
  size_t four;
  size_t five;
  size_t six;

  // Temporary until we know what we are doing
  floatType nu;

  // Vector Index

  int getIndexNoP(size_t j1, size_t j2, size_t j3, size_t k1, size_t k2, size_t k3) {
    return j1 + j2 * two + j3 * three + k1 * four + k2 * five + k3 * six;
  }

  // Fill functions

  void fillTableEvenSpacing() {
    std::cout << "Building Table\n";
    std::cout << "Number of points: " << numberOfPoints << " Point distance: " << pointDistance << " Inverse pointDistance: " << invPointDistance << "\n";
    floatType i1, i2, i3, j1, j2, j3, k1, k2, k3;
    j1 = j2 = j3 = k1 = k2 = k3 = (pointDistance / 2) - 2 * cutoffSquared;  // Why 2 * cutoffSquared?
    uint64_t dupes = 0;
    for (auto kc3 = 0; kc3 < numberOfPoints; kc3++) {
      for (auto kc2 = 0; kc2 < numberOfPoints; kc2++) {
        for (auto kc1 = 0; kc1 < numberOfPoints; kc1++) {
          for (auto jc3 = 0; jc3 < numberOfPoints; jc3++) {
            for (auto jc2 = 0; jc2 < numberOfPoints; jc2++) {
              for (auto jc1 = 0; jc1 < numberOfPoints; jc1++) {
                auto n = ATFunctor(0, 0, 0, j1, j2, j3, k1, k2, k3);
                lut.push_back(n);
                j1 += pointDistance;
              }
              j2 += pointDistance;
            }
            j3 += pointDistance;
          }
          k1 += pointDistance;
        }
        k2 += pointDistance;
      }
      std::cout << "Done with kc3=" << kc3 << " k3=" << k3 << "\n";
      k3 += pointDistance;
    }
    auto first = 0;
    auto total = 0;
    auto nan = 0;
    auto last = 0;
    for (auto i = 0; i < lut.size(); i++) {
      if (!std::isnan(lut.at(i).second)) {
        if (first == 0) first = i;
        last = i;
        total++;
      } else {
        nan++;
      }
    }
    std::cout << "First: " << first << " Last: " << last << " Total: " << total << " Nan: " << nan
              << " Duplicates: " << dupes << " Size: " << lut.size() << "\n";
  };

  // Interpolation functions

  Entry getNextNeighbor(floatType i1, floatType i2, floatType i3, floatType j1, floatType j2, floatType j3,
                        floatType k1, floatType k2, floatType k3) {
    // using namespace autopas::utils::ArrayMath::literals;
    // auto accurate = ATFunctor(i1, i2, i3, j1, j2, j3, k1, k2, k3);
    AutoPasLog(DEBUG, "Input was {} {} {} | {} {} {} | {} {} {}", i1, i2, i3, j1, j2, j3, k1, k2, k3);
    // static Entry totalRelError;
    if constexpr (intervalType == evenSpacing) {
      j1 -= i1;
      j2 -= i2;
      j3 -= i3;
      k1 -= i1;
      k2 -= i2;
      k3 -= i3;
      // AutoPasLog(DEBUG, "Subtracted: {} {} {} | {} {} {}", j1, j2, j3, k1, k2, k3);
      //  2 * cutoffSquared als member
      j1 += 2 * cutoffSquared;
      j2 += 2 * cutoffSquared;
      j3 += 2 * cutoffSquared;
      k1 += 2 * cutoffSquared;
      k2 += 2 * cutoffSquared;
      k3 += 2 * cutoffSquared;
      // / pointDistance als member
      size_t ij1 = std::floor(j1 * invPointDistance);
      size_t ij2 = std::floor(j2 * invPointDistance);
      size_t ij3 = std::floor(j3 * invPointDistance);
      size_t ik1 = std::floor(k1 * invPointDistance);
      size_t ik2 = std::floor(k2 * invPointDistance);
      size_t ik3 = std::floor(k3 * invPointDistance);
      if (ij1 == numberOfPoints) ij1--;
      if (ij2 == numberOfPoints) ij2--;
      if (ij3 == numberOfPoints) ij3--;
      if (ik1 == numberOfPoints) ik1--;
      if (ik2 == numberOfPoints) ik2--;
      if (ik3 == numberOfPoints) ik3--;

//      auto ret = lut[(
//          getIndexNoP(ij1, ij2, ij3, ik1, ik2, ik3))];  // How slow is std::floor?
//                                                        //      auto relErr = relError(ret, accurate);
//                                                              totalRelError.first[0] += relErr.first[0];
//                                                              totalRelError.first[1] += relErr.first[1];
//                                                              totalRelError.first[2] += relErr.first[2];
//                                                              totalRelError.second += relErr.second;
//                                                              AutoPasLog(DEBUG, "Total rel Error: {}", totalRelError);
//       AutoPasLog(DEBUG, "Return {} instead of {}\nAbs: {} Rel: {}", ret, accurate, absError(ret, accurate),
//       relError(ret, accurate)); AutoPasLog(DEBUG, "Used {} {} {} | {} {} {}", (pointDistance / 2) - 2*cutoffSquared +
//       ij1 * pointDistance, (pointDistance / 2) - 2*cutoffSquared + ij2 * pointDistance, (pointDistance / 2) -
//       2*cutoffSquared + ij3 * pointDistance, (pointDistance / 2) - 2*cutoffSquared + ik1 * pointDistance,
//       (pointDistance / 2) - 2*cutoffSquared + ik2 * pointDistance, (pointDistance / 2) - 2*cutoffSquared + ik3 *
//       pointDistance);
      return lut[getIndexNoP(ij1, ij2, ij3, ik1, ik2, ik3)];
    }
  }

  // Functor stub

  Entry ATFunctor(floatType i1, floatType i2, floatType i3, floatType j1, floatType j2, floatType j3, floatType k1,
                  floatType k2, floatType k3) {
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
    const auto forceJDirectionIJ =
        displacementIJ * (-IJDotKI * JKDotKI + distSquaredJK * distSquaredKI - 5.0 * allDotProducts / distSquaredIJ);
    const auto forceJDirectionJK =
        displacementJK * (IJDotKI * JKDotKI - distSquaredIJ * distSquaredKI + 5.0 * allDotProducts / distSquaredJK);
    forceJ = (forceJDirectionKI + forceJDirectionIJ + forceJDirectionJK) * factor;
    forceK = (forceI + forceJ) * (-1.0);
    const auto potentialEnergy = factor * (allDistsSquared - 3.0 * allDotProducts) / 9.0;

    return std::make_pair(std::array{forceI, forceJ, forceK}, potentialEnergy);
  }

  // Error functions for debugging

  Entry absError(Entry e, Entry acc) {
    using namespace autopas::utils::ArrayMath::literals;
    return std::make_pair(std::array{e.first[0] - acc.first[0], e.first[1] - acc.first[1], e.first[2] - acc.first[2]},
                          std::abs(e.second - acc.second));
  }

  Entry relError(Entry e, Entry acc) {
    using namespace autopas::utils::ArrayMath::literals;
    return std::make_pair(
        std::array{(e.first[0] - acc.first[0]) / acc.first[0], (e.first[1] - acc.first[1]) / acc.first[1],
                   (e.first[2] - acc.first[2]) / acc.first[2]},
        std::abs(e.second - acc.second));
  }
};

}  // namespace ForceLookUpTable

#endif  // AUTOPAS_ATLOOKUPTABLE_H
