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
//    auto f =  ATFunctor(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0);
//    auto t = ATFunctor(2.0, 1.0, 1.0, 3.0, 2.0, 2.0, 4.0, 3.0, 3.0);
//    std::cout << f.second << " " << t.second << " : " << f.first[0][0] << " " << t.first[0][0] << " | " << f.first[0][1] << " " << t.first[0][1] << " | " << f.first[0][2] << " " << t.first[0][2] << " | " << f.first[1][0] << " " << t.first[1][0] << " | " << f.first[1][1] << " " << t.first[1][1] << " | " << f.first[1][2] << " " << t.first[1][2] << " | " << f.first[2][0] << " " << t.first[2][0] << " | " << f.first[2][1] << " " << t.first[2][1] << " | " << f.first[2][2] << " " << t.first[2][2] << "\n";
//    std::cout << (f.second == t.second && f.first[0][0] == t.first[0][0] && f.first[0][1] == t.first[0][1] && f.first[0][2] == t.first[0][2] && f.first[1][0] == t.first[1][0] && f.first[1][1] == t.first[1][1] && f.first[1][2] == t.first[1][2] && f.first[2][0] == t.first[2][0] && f.first[2][1] == t.first[2][1] && f.first[2][2] == t.first[2][2]) << "\n";
//    f =  ATFunctor(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0);
//    t = ATFunctor(1.0, 2.0, 1.0, 2.0, 3.0, 2.0, 3.0, 4.0, 3.0);
//    std::cout << f.second << " " << t.second << " : " << f.first[0][0] << " " << t.first[0][0] << " | " << f.first[0][1] << " " << t.first[0][1] << " | " << f.first[0][2] << " " << t.first[0][2] << " | " << f.first[1][0] << " " << t.first[1][0] << " | " << f.first[1][1] << " " << t.first[1][1] << " | " << f.first[1][2] << " " << t.first[1][2] << " | " << f.first[2][0] << " " << t.first[2][0] << " | " << f.first[2][1] << " " << t.first[2][1] << " | " << f.first[2][2] << " " << t.first[2][2] << "\n";
//    std::cout << (f.second == t.second && f.first[0][0] == t.first[0][0] && f.first[0][1] == t.first[0][1] && f.first[0][2] == t.first[0][2] && f.first[1][0] == t.first[1][0] && f.first[1][1] == t.first[1][1] && f.first[1][2] == t.first[1][2] && f.first[2][0] == t.first[2][0] && f.first[2][1] == t.first[2][1] && f.first[2][2] == t.first[2][2]) << "\n";
//    f =  ATFunctor(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0);
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
      pointDistance = (cutoffSquared * 2) / numberOfPoints;
      fillTableEvenSpacing();
    }
  }

  std::pair<std::array<std::array<floatType, 3>, 3>, floatType> retrieveValue(const std::array<double, 3>& displacementIJ, const std::array<double, 3>& displacementJK, const std::array<double, 3>& displacementKI) {
    AutoPasLog(DEBUG, "Retrieved value from AT-LUT");
    if constexpr (interpolationType == nextNeighbor) {
      return getNextNeighbor(displacementIJ[0], displacementIJ[1], displacementIJ[2], displacementJK[0], displacementJK[1], displacementJK[2], displacementKI[0], displacementKI[1], displacementKI[2]);
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
    std::cout << "Building Table\n";
    std::cout << "Number of points: " << numberOfPoints << " Point distance: " << pointDistance << "\n";
    floatType i1, i2, i3, j1, j2, j3, k1, k2, k3;
    i1 = i2 = i3 = j1 = j2 = j3 = k1 = k2 = k3 = (pointDistance / 2) - cutoffSquared;
    uint64_t dupes = 0;
    for (auto kc3 = 0; kc3 < numberOfPoints; kc3++) {
      for (auto kc2 = 0; kc2 < numberOfPoints; kc2++) {
        for (auto kc1 = 0; kc1 < numberOfPoints; kc1++) {
          for (auto jc3 = 0; jc3 < numberOfPoints; jc3++) {
            for (auto jc2 = 0; jc2 < numberOfPoints; jc2++) {
              for (auto jc1 = 0; jc1 < numberOfPoints; jc1++) {
                for (auto ic3 = 0; ic3 < numberOfPoints; ic3++) {
                  for (auto ic2 = 0; ic2 < numberOfPoints; ic2++) {
                    for (auto ic1 = 0; ic1 < numberOfPoints; ic1++) {
                      auto n = ATFunctor(i1, i2, i3, j1, j2, j3, k1, k2, k3);
                      for (auto i = 0; i<lut.size(); i++) {
                        auto val = lut[i];
                        if (n.second == val.second && n.first[0][0] == val.first[0][0] && n.first[0][1] == val.first[0][1] && n.first[0][2] == val.first[0][2] && n.first[1][0] == val.first[1][0] && n.first[1][1] == val.first[1][1] && n.first[1][2] == val.first[1][2] && n.first[2][0] == val.first[2][0] && n.first[2][1] == val.first[2][1] && n.first[2][2] == val.first[2][2]) {
                          std::cout << "Duplicate: " << i1 << " " << i2 << " " << i3 << " | " << j1 << " " << j2 << " "
                                    << j3 << " | " << k1 << " " << k2 << " " << k3 << "\n";
                          dupes++;
                        }
                      }
                      lut.push_back(n);
                      i1 += pointDistance;
                    }
                    i2 += pointDistance;
                  }
                  i3 += pointDistance;
                }
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
      std::cout << "Done with kc3=" << kc3 << "k3=" << k3 <<  "\n";
      k3 += pointDistance;
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
    std::cout << "First: " << first << " Last: " << last << " Total: " << total << " Nan: " << nan << " Duplicates: " << dupes << "\n";
  };

  // Interpolation functions

  std::pair<std::array<std::array<floatType, 3>, 3>, floatType> getNextNeighbor(floatType i1, floatType i2, floatType i3, floatType j1, floatType j2, floatType j3, floatType k1, floatType k2, floatType k3) {
    auto accurate = ATFunctor(i1, i2, i3, j1, j2, j3, k1, k2, k3);
    if constexpr (intervalType == evenSpacing) {
      i1 += cutoffSquared;
      i2 += cutoffSquared;
      i3 += cutoffSquared;
      j1 += cutoffSquared;
      j2 += cutoffSquared;
      j3 += cutoffSquared;
      k1 += cutoffSquared;
      k2 += cutoffSquared;
      k3 += cutoffSquared;
      size_t ii1 = std::floor(i1 / pointDistance);
      size_t ii2 = std::floor(i2 / pointDistance);
      size_t ii3 = std::floor(i3 / pointDistance);
      size_t ij1 = std::floor(j1 / pointDistance);
      size_t ij2 = std::floor(j2 / pointDistance);
      size_t ij3 = std::floor(j3 / pointDistance);
      size_t ik1 = std::floor(k1 / pointDistance);
      size_t ik2 = std::floor(k2 / pointDistance);
      size_t ik3 = std::floor(k3 / pointDistance);
      if (ii1 == numberOfPoints)
        ii1--;
      if (ii2 == numberOfPoints)
        ii2--;
      if (ii3 == numberOfPoints)
        ii3--;
      if (ij1 == numberOfPoints)
        ij1--;
      if (ij2 == numberOfPoints)
        ij2--;
      if (ij3 == numberOfPoints)
        ij3--;
      if (ik1 == numberOfPoints)
        ik1--;
      if (ik2 == numberOfPoints)
        ik2--;
      if (ik3 == numberOfPoints)
        ik3--;

      auto ret = lut.at(getIndexNoP(ii1, ii2, ii3, ij1, ij2, ij3, ik1, ik2, ik3)); // How slow is std::floor?
      //AutoPasLog(DEBUG, "Return {} instead of {}", ret, accurate);
      return ret;
    }
  }

  // Functor stub

  std::pair<std::array<std::array<floatType, 3>, 3>, floatType> ATFunctor(floatType i1, floatType i2, floatType i3, floatType j1, floatType j2, floatType j3, floatType k1, floatType k2, floatType k3) {
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


};

}


#endif  // AUTOPAS_ATLOOKUPTABLE_H
