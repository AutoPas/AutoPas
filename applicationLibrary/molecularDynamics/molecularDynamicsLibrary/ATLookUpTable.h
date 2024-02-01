//
// Created by jan on 1/26/24.
//

#ifndef AUTOPAS_ATLOOKUPTABLE_H
#define AUTOPAS_ATLOOKUPTABLE_H

#include <cmath>
#include <vector>

#include "autopas/utils/logging/Logger.h"
#include "LookUpTableTypes.h"

namespace ForceLookUpTable {

template <IntervalType intervalType, InterpolationType interpolationType, typename floatType = double, typename intType = unsigned long>
class ATLookUpTable {
 public:
  ATLookUpTable() {
    //AutoPasLog(DEBUG, "Default constructor called.");
  }

  // list: cutoffSquared, sigmaSquared, epsilon24, ... (numberOfPoints)
  // Extremely unreadable and user-error-prone

  ATLookUpTable(std::initializer_list<floatType> args) {
    //AutoPasLog(DEBUG, "LUT created.");
    if (args.size() < 3) {  // Fail gracefully
      //AutoPasLog(CRITICAL, "Args only has {} elements, but needs at least 3.", args.size());
      return;
    }
    cutoffSquared = args.begin()[0];
    sigmaSquared = args.begin()[1];
    epsilon24 = args.begin()[2];
    if constexpr (intervalType == evenSpacing) {
      if (args.size() != 4) {  // Fail
        //AutoPasLog(CRITICAL, "Args has {} elements, but needs 4 for even spacing.", args.size());
        return;
      }
      numberOfPoints = static_cast<intType>(args.begin()[3]);
      if (numberOfPoints == 0)
        throw autopas::utils::ExceptionHandler::AutoPasException("At least one point needed for LUT.");
      pointDistance = cutoffSquared / numberOfPoints;
      fillTableEvenSpacing();
    }
  }

  floatType retrieveValue(floatType distanceSquared) {
    if constexpr (interpolationType == nextNeighbor) {
      return getNextNeighbor(distanceSquared);
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

  // Fill functions

  void fillTableEvenSpacing () {

  };

  // Interpolation functions

  floatType getNextNeighbor(floatType dr2) {
    if constexpr (intervalType == evenSpacing) {
      if (dr2 == cutoffSquared)
        return lut.at(numberOfPoints-1);
      auto ret = lut.at(std::floor(dr2 / pointDistance)); // How slow is std::floor?
      auto accurate = LJFunctor(dr2);
      AutoPasLog(DEBUG, "Return {} instead of {}", ret, accurate);
      return ret;
    }
  }

  // Functor stub



};

}


#endif  // AUTOPAS_ATLOOKUPTABLE_H
