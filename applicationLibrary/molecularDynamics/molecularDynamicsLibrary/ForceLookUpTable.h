//
// Created by jan on 1/1/24.
//

#ifndef AUTOPAS_FORCELOOKUPTABLE_H
#define AUTOPAS_FORCELOOKUPTABLE_H

#include <cmath>
#include <vector>

#include "autopas/utils/logging/Logger.h"

namespace ForceLookUpTable {
enum FunctorType { LJAoS };
enum IntervalType { evenSpacing };
enum InterpolationType { nextNeighbor };


template <IntervalType intervalType, InterpolationType interpolationType, FunctorType functorType, typename floatType = double, typename intType = unsigned long>
class ForceLookUpTable {
 public:

  ForceLookUpTable() {
      AutoPasLog(DEBUG, "Default constructor called.");
  }

  // list: cutoffSquared, sigmaSquared, epsilon24, ... (numberOfPoints)
  // Extremely unreadable and user-error-prone
  ForceLookUpTable(std::initializer_list<floatType> args) {
    AutoPasLog(DEBUG, "LUT created.");
    if (args.size() < 3) {  // Fail gracefully
      AutoPasLog(CRITICAL, "Args only has {} elements, but needs at least 3.", args.size());
      return;
    }
    if constexpr (intervalType == evenSpacing) {
      if (args.size() != 4) {  // Fail
        AutoPasLog(CRITICAL, "Args has {} elements, but needs 4 for even spacing.", args.size());
        return;
      }
      evenSpacingInitializer(static_cast<intType>(args.begin()[3]));
    }
  }


  floatType retrieveSingleValue(floatType distanceSquared) {
    if constexpr (interpolationType == nextNeighbor) {
      return nextNeighborSingle(distanceSquared);
    }
  }



 private:
  std::vector<floatType> lut; // If we need to save the points, just use every other spot for point and then data
  intType numberOfPoints; // For even spacing
  floatType pointDistance; // For even spacing
  floatType cutoffSquared;

  // Temporary until we know what we are doing
  floatType sigmaSquared;
  floatType epsilon24;

  // Initializers

  void evenSpacingInitializer(intType nOP) {
    if (nOP == 0) { // Fail gracefully
      AutoPasLog(CRITICAL, "Don't set number of points to 0 when using evenSpacing.");
      return;
    }
    numberOfPoints = nOP;
    pointDistance = sqrt(cutoffSquared) / numberOfPoints; // May break when floatType is not double. cutoffSquared obviously can't be negative.
    if constexpr (functorType == LJAoS) {
      fillTableLJEvenSpacing();
    }
  }

  // Fill functions

  void fillTableLJEvenSpacing () {
      for (auto i = 0; i<numberOfPoints; i++) {
        lut.push_back(LJFunctor((i+1) * pointDistance));
      }
  };

  // Interpolation functions

  floatType nextNeighborSingle(floatType dr2) {
    if constexpr (intervalType == evenSpacing) {
      if (dr2 == cutoffSquared)
        return lut.at(numberOfPoints-1);
      return lut.at(std::floor(dr2 / pointDistance)); // How slow is std::floor?
    }
  }

  // Functor stubs

  floatType LJFunctor(floatType dr2) {
    floatType invdr2 = 1. / dr2;
    floatType lj6 = sigmaSquared * invdr2;
    lj6 = lj6 * lj6 * lj6;
    floatType lj12 = lj6 * lj6;
    floatType lj12m6 = lj12 - lj6;
    return epsilon24 * (lj12 + lj12m6) * invdr2;
  }

};

}

#endif  // AUTOPAS_FORCELOOKUPTABLE_H
