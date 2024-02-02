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

  // list: cutoffSquared, nu, ... (numberOfPoints)
  // Extremely unreadable and user-error-prone

  ATLookUpTable(std::initializer_list<floatType> args) {
    //AutoPasLog(DEBUG, "LUT created.");
    if (args.size() < 3) {  // Fail gracefully
      //AutoPasLog(CRITICAL, "Args only has {} elements, but needs at least 3.", args.size());
      return;
    }
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
      pointDistance = (cutoffSquared * 2) / numberOfPoints;
      fillTableEvenSpacing();
    }
  }

  std::pair<std::array<std::array<floatType, 3>, 3>, floatType> retrieveValue(std::array<double, 3>& displacementIJ, std::array<double, 3>& displacementJK, std::array<double, 3>& displacementKI) {
    if constexpr (interpolationType == nextNeighbor) {
      return getNextNeighbor(displacementIJ, displacementJK, displacementKI);
    }
  }

 private:
  std::vector<std::pair<std::array<std::array<floatType, 3>, 3>, floatType>> lut; // Pair of potential energy and Triplets of forces for newton3, each of which consist of 3 values fpr the 3 dimensions
  intType numberOfPoints; // For even spacing
  floatType pointDistance; // For even spacing
  floatType cutoffSquared;

  // Temporary until we know what we are doing
  floatType nu;

  // Vector Index

  int getIndexNoP(size_t i1, size_t i2, size_t i3, size_t j1, size_t j2, size_t j3, size_t k1, size_t k2, size_t k3) {
    static const size_t two = numberOfPoints;
    static const size_t three = two * numberOfPoints;
    static const size_t four = three * numberOfPoints;
    static const size_t five = four * numberOfPoints;
    static const size_t six = five * numberOfPoints;
    static const size_t seven = six * numberOfPoints;
    static const size_t eight = seven * numberOfPoints;
    static const size_t nine = eight * numberOfPoints;

    return i1 + i2 * two + i3 * three + j1 * four + j2 * five + j3 * six + k1 * seven + k2 * eight + k3 * nine;
  }



  // Fill functions

  void fillTableEvenSpacing () {
    floatType i1, i2, i3, j1, j2, j3, k1, k2, k3 = (pointDistance / 2) - cutoffSquared;
    for (auto ic1 = 0; ic1 < numberOfPoints; ic1++) {
      for (auto ic2 = 0; ic2 < numberOfPoints; ic2++) {
        for (auto ic3 = 0; ic3 < numberOfPoints; ic3++) {
          for (auto jc1 = 0; jc1 < numberOfPoints; jc1++) {
            for (auto jc2 = 0; jc2 < numberOfPoints; jc2++) {
              for (auto jc3 = 0; jc3 < numberOfPoints; jc3++) {
                for (auto kc1 = 0; kc1 < numberOfPoints; kc1++) {
                  for (auto kc2 = 0; kc2 < numberOfPoints; kc2++) {
                    for (auto kc3 = 0; kc3 < numberOfPoints; kc3++) {
                      lut.at(getIndexNoP(ic1, ic2, ic3, jc1, jc2, jc3, kc1, kc2, kc3)) = ATFunctor(i1, i2, i3, j1, j2, j3, k1, k2, k3);
                      kc3 += pointDistance;
                    }
                    kc2 += pointDistance;
                  }
                  kc1 += pointDistance;
                }
                jc3 += pointDistance;
              }
              jc2 += pointDistance;
            }
            jc1 += pointDistance;
          }
          ic3 += pointDistance;
        }
        ic2 += pointDistance;
      }
      ic1 += pointDistance;
    }
  };

  // Interpolation functions

  std::pair<std::array<std::array<floatType, 3>, 3>, floatType> getNextNeighbor(floatType i1, floatType i2, floatType i3, floatType j1, floatType j2, floatType j3, floatType k1, floatType k2, floatType k3) {
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

  std::pair<std::array<std::array<floatType, 3>, 3>, floatType> ATFunctor(floatType i1, floatType i2, floatType i3, floatType j1, floatType j2, floatType j3, floatType k1, floatType k2, floatType k3) {

    const auto displacementIJ = {j1 - i1, j2 - i2, j3 - i3};
    const auto displacementJK = {k1 - j1, k2 - j2, k3 - j3};
    const auto displacementKI = {i1 - k1, i2 - k2, i3 - k3};

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

    return std::make_pair({forceI, forceJ, forceK}, potentialEnergy);
  }


};

}


#endif  // AUTOPAS_ATLOOKUPTABLE_H
