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

template <PositionType positionType, IntervalType intervalType, InterpolationType interpolationType, typename floatType = double, typename intType = unsigned long>
class ATLookUpTable{};

template <IntervalType intervalType, InterpolationType interpolationType, typename floatType, typename intType>
class ATLookUpTable<absolute, intervalType, interpolationType, floatType, intType> {

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
      pointDistance = (cutoffSquared * 4) / numberOfPoints;
      two = numberOfPoints;
      three = two * numberOfPoints;
      four = three * numberOfPoints;
      five = four * numberOfPoints;
      six = five * numberOfPoints;
      fillTableEvenSpacing();
    }
  }

  Entry retrieveValue(const std::array<double, 3>& displacementIJ, const std::array<double, 3>& displacementJK, const std::array<double, 3>& displacementKI) {
    AutoPasLog(DEBUG, "Retrieved value from AT-LUT");
    if constexpr (interpolationType == nextNeighbor) {
      return getNextNeighbor(displacementIJ[0], displacementIJ[1], displacementIJ[2], displacementJK[0], displacementJK[1], displacementJK[2], displacementKI[0], displacementKI[1], displacementKI[2]);
    }
  }

 private:
  std::vector<Entry> lut; // Pair of potential energy and Triplets of forces for newton3, each of which consist of 3 values fpr the 3 dimensions
  intType numberOfPoints; // For even spacing
  floatType pointDistance; // For even spacing
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

  void fillTableEvenSpacing () {
    std::cout << "Building Table\n";
    std::cout << "Number of points: " << numberOfPoints << " Point distance: " << pointDistance << "\n";
    floatType i1, i2, i3, j1, j2, j3, k1, k2, k3;
    j1 = j2 = j3 = k1 = k2 = k3 = (pointDistance / 2) - 2 * cutoffSquared; // Why 2 * cutoffSquared?
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
    std::cout << "First: " << first << " Last: " << last << " Total: " << total << " Nan: " << nan << " Duplicates: " << dupes << " Size: "<< lut.size() << "\n";
  };

  // Interpolation functions

  Entry getNextNeighbor(floatType i1, floatType i2, floatType i3, floatType j1, floatType j2, floatType j3, floatType k1, floatType k2, floatType k3) {
    //using namespace autopas::utils::ArrayMath::literals;
    //auto accurate = ATFunctor(i1, i2, i3, j1, j2, j3, k1, k2, k3);
    AutoPasLog(DEBUG, "Input was {} {} {} | {} {} {} | {} {} {}", i1, i2, i3, j1, j2, j3, k1, k2, k3);
    //static Entry totalRelError;
    if constexpr (intervalType == evenSpacing) {
      j1 -= i1;
      j2 -= i2;
      j3 -= i3;
      k1 -= i1;
      k2 -= i2;
      k3 -= i3;
      //AutoPasLog(DEBUG, "Subtracted: {} {} {} | {} {} {}", j1, j2, j3, k1, k2, k3);
      // 2 * cutoffSquared als member
      j1 += 2 * cutoffSquared;
      j2 += 2 * cutoffSquared;
      j3 += 2 * cutoffSquared;
      k1 += 2 * cutoffSquared;
      k2 += 2 * cutoffSquared;
      k3 += 2 * cutoffSquared;
      // / pointDistance als member
      size_t ij1 = std::floor(j1 / pointDistance);
      size_t ij2 = std::floor(j2 / pointDistance);
      size_t ij3 = std::floor(j3 / pointDistance);
      size_t ik1 = std::floor(k1 / pointDistance);
      size_t ik2 = std::floor(k2 / pointDistance);
      size_t ik3 = std::floor(k3 / pointDistance);
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

      auto ret = lut[(getIndexNoP(ij1, ij2, ij3, ik1, ik2, ik3))]; // How slow is std::floor?
//      auto relErr = relError(ret, accurate);
//      totalRelError.first[0] += relErr.first[0];
//      totalRelError.first[1] += relErr.first[1];
//      totalRelError.first[2] += relErr.first[2];
//      totalRelError.second += relErr.second;
//      AutoPasLog(DEBUG, "Total rel Error: {}", totalRelError);
      //AutoPasLog(DEBUG, "Return {} instead of {}\nAbs: {} Rel: {}", ret, accurate, absError(ret, accurate), relError(ret, accurate));
      //AutoPasLog(DEBUG, "Used {} {} {} | {} {} {}", (pointDistance / 2) - 2*cutoffSquared + ij1 * pointDistance, (pointDistance / 2) - 2*cutoffSquared + ij2 * pointDistance, (pointDistance / 2) - 2*cutoffSquared + ij3 * pointDistance, (pointDistance / 2) - 2*cutoffSquared + ik1 * pointDistance, (pointDistance / 2) - 2*cutoffSquared + ik2 * pointDistance, (pointDistance / 2) - 2*cutoffSquared + ik3 * pointDistance);
      return ret;
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


};

}


#endif  // AUTOPAS_ATLOOKUPTABLE_H
