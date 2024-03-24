//
// Created by jan on 1/1/24.
//

#ifndef AUTOPAS_LJLOOKUPTABLE_H
#define AUTOPAS_LJLOOKUPTABLE_H

#include <cmath>
#include <vector>

#include "LookUpTableTypes.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/logging/Logger.h"

namespace ForceLookUpTable {

/**
 * Look-Up Table for the Lennard-Jones Functor
 * @tparam intervalType How the stored support points are laid out, currently only an even spacing is supported
 * @tparam interpolationType How the values between the support points are computed, currently a jump to the next lowest
 * neighbor and linear interpolation is supported
 * @tparam floatType
 * @tparam intType
 */
template <IntervalType intervalType, InterpolationType interpolationType, typename floatType = double,
          typename intType = unsigned long>
class LJLookUpTable {
 public:
  LJLookUpTable() {
    // AutoPasLog(DEBUG, "Default constructor called.");
  }

  // list: cutoffSquared, sigmaSquared, epsilon24, ... (numberOfPoints)
  // Extremely unreadable and user-error-prone

  /**
   * Constructor of the look-up table
   * Takes an initializer list of floats because it has to be agnostic to the different parameters needed for each
   * combination of intervalType and interpolationType
   * @param args Initializer list that (for evenSpacing and nextNeighbor or linearInterpolation) takes the form of
   * {cutoffSquared, sigmaSquared, epsilon24, numberOfPoints}
   */
  LJLookUpTable(std::initializer_list<floatType> args) {
    // AutoPasLog(DEBUG, "LUT created.");
    if (args.size() < 3) {  // Fail gracefully
      // AutoPasLog(CRITICAL, "Args only has {} elements, but needs at least 3.", args.size());
      return;
    }
    cutoffSquared = args.begin()[0];
    sigmaSquared = args.begin()[1];
    epsilon24 = args.begin()[2];
    if constexpr (intervalType == evenSpacing) {
      if (args.size() != 4) {  // Fail
        // AutoPasLog(CRITICAL, "Args has {} elements, but needs 4 for even spacing.", args.size());
        return;
      }
      numberOfPoints = static_cast<intType>(args.begin()[3]);
      if (numberOfPoints == 0)
        throw autopas::utils::ExceptionHandler::AutoPasException("At least one point needed for LUT.");
      pointDistance = cutoffSquared / numberOfPoints;
      fillTableEvenSpacing();
    }
  }

  /**
   * Retrieves a value from the lut
   * @param distanceSquared
   * @return The stored/interpolated force
   */
  floatType retrieveValue(floatType distanceSquared) {
    AutoPasLog(DEBUG, "Retrieved Value.");
    if constexpr (interpolationType == nextNeighbor) {
      return getNextNeighbor(distanceSquared);
    }
    if constexpr (interpolationType == linear) {
      return getLinear(distanceSquared);
    }
  }

 private:
  std::vector<floatType> lut;  // If we need to save the points, just use every other spot for point and then data
  intType numberOfPoints;      // For even spacing
  floatType pointDistance;     // For even spacing
  floatType cutoffSquared;

  // Temporary until we know what we are doing
  floatType sigmaSquared;
  floatType epsilon24;

  // Fill functions

  void fillTableEvenSpacing() {
    for (auto i = 0; i < numberOfPoints; i++) {
      lut.push_back(LJFunctor((pointDistance / 2) + (i * pointDistance)));
    }
    // Using std::cout because Logger isn't initiated yet, when constructors are called.
    //      std::cout <<  "Table filled evenly spaced with distance " << pointDistance << "\n Content: ";
    //      for (auto i=0; i<numberOfPoints; i++) {
    //        std::cout << i << " : " << (pointDistance/2) + (i * pointDistance) << " : " << lut.at(i) << " | ";
    //      }
    //      std::cout << "\n";
  };

  // Interpolation functions

  floatType getNextNeighbor(floatType dr2) {
    if constexpr (intervalType == evenSpacing) {
      if (dr2 == cutoffSquared) {
        auto ret = lut.at(numberOfPoints - 1);
        AutoPasLog(DEBUG, "dr2 was cutoffSquared. Returning {}", ret);
        return ret;
      }
      // Compute perfect value for lowest dr2 where error is greatest
      if (dr2 < pointDistance / 2) {
        AutoPasLog(DEBUG, "dr2 is less than half pointDistance. Computing perfect value {}", LJFunctor(dr2));
        return LJFunctor(dr2);
      }
      auto ret = lut.at(std::floor(dr2 / pointDistance));  // How slow is std::floor?
      AutoPasLog(DEBUG, "Return {} instead of {}", ret, LJFunctor(dr2));
      static float totalError;
      AutoPasLog(DEBUG, "Error: {} | totalError {}", LJFunctor(dr2) - ret, [dr2, ret, this]() {
        totalError += abs(LJFunctor(dr2) - ret);
        return totalError;
      }());
      return ret;
    }
  }

  floatType getLinear(floatType dr2) {
    AutoPasLog(DEBUG, "Linear! Point distance: {} | Number of points {} | dr2 {}", pointDistance, numberOfPoints, dr2);
    if constexpr (intervalType == evenSpacing) {
      if (dr2 >= cutoffSquared - pointDistance / 2) {
        AutoPasLog(DEBUG, "dr2 {} is over last point {}. Returning {} Correct {}", dr2,
                   cutoffSquared - pointDistance / 2, lut.at(numberOfPoints - 1), LJFunctor(dr2));
        AutoPasLog(DEBUG, "Error {}", LJFunctor(dr2) - lut.at(numberOfPoints - 1));
        return lut.at(numberOfPoints - 1);
      }
      if (dr2 <= (pointDistance / 2)) {
        AutoPasLog(DEBUG, "dr2 {} was less or equal than half pointDistance {}. Returning perfect value {}", dr2,
                   pointDistance / 2, LJFunctor(dr2));
        return LJFunctor(dr2);
      }
      if (std::fmod(dr2, pointDistance) == 0) {
        AutoPasLog(DEBUG, "dr2 is multiple of pointDistance! Remainder is {}", std::fmod(dr2, pointDistance));
        return lut.at(dr2 / pointDistance);
      }
      floatType lowerX = std::floor((dr2 - pointDistance / 2) / pointDistance);
      floatType upperX = std::ceil((dr2 - pointDistance / 2) / pointDistance);
      floatType lowerY = lut.at(lowerX);
      floatType upperY = lut.at(upperX);
      lowerX = lowerX * pointDistance + pointDistance / 2;
      upperX = upperX * pointDistance + pointDistance / 2;
      auto ret =
          lowerY +
          (dr2 - lowerX) * (upperY - lowerY) /
              (upperX -
               lowerX);  // Content: 0 : 0.5 : 5760 | 1 : 1.5 : -1.93141 | 2 : 2.5 : -0.535757 | 3 : 3.5 : -0.152473
      static float totalError;
      AutoPasLog(DEBUG,
                 "lowerX: {} | upperX: {} | lowerY {} | upperY {}\n                            Input: {} | Return: {} "
                 "| Correct: {} | Error {} | totalError {}",
                 lowerX, upperX, lowerY, upperY, dr2, ret, LJFunctor(dr2), LJFunctor(dr2) - ret, [dr2, ret, this]() {
                   totalError += abs(LJFunctor(dr2) - ret);
                   return totalError;
                 }());
      return ret;
    }
  }

  // Functor stub

  floatType LJFunctor(floatType dr2) {
    floatType invdr2 = 1. / dr2;
    floatType lj6 = sigmaSquared * invdr2;
    lj6 = lj6 * lj6 * lj6;
    floatType lj12 = lj6 * lj6;
    floatType lj12m6 = lj12 - lj6;
    return epsilon24 * (lj12 + lj12m6) * invdr2;
  }
};

}  // namespace ForceLookUpTable

#endif  // AUTOPAS_LJLOOKUPTABLE_H
