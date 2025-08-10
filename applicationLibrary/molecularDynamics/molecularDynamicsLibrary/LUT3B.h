//
// Created by sliep on 05.05.2025.
//
#pragma once

#include "../examples/md-flexible/src/TypeDefinitions.h"
#include "ParentLUT.h"
#include "autopas/utils/Math.h"
#include "autopas/utils/Timer.h"

namespace mdLib {
class LUT3B : public ParentLUT {
 public:
  //  LUT3B(int resolution, double cutoffSquared, double Nu) : ParentLUT(resolution, cutoffSquared) {
  //    nu = Nu;
  //
  //    fill(cutoffSquared, false);
  //
  //  }
  //
  //  LUT3B(int resolution, double cutoffSquared) : ParentLUT(resolution, cutoffSquared) {
  //    // purposely not fill the lut so it can exist in the PPL
  //    if (resolution == 0) {
  //      return;
  //    }
  //
  //    fill(cutoffSquared, false);
  //
  //  }

  LUT3B(int resolution, double cutoffSquared, bool global = false, mdLib::InterpolationMethod method = NN,
        bool krypton = false, double _nu = 0)
      : ParentLUT(resolution, cutoffSquared) {
    KRYPTON = krypton;
    nu = _nu;
    useGlobal = global;
    interpolationMethod = method;
    _numberOfPoints = calculateNumberOfPoints3B(resolution);

    if (resolution == 0) {
      return;
    }

    fill();
  }

  void fill() {
    if (!useGlobal) {
      fill_plain();
    }
    if (useGlobal) {
      if (KRYPTON) {
        fill_global_krypton();
      } else {
        fill_global();
      }
    }
  }

  void fill_plain() {
    _lut3B.resize(_resolution + 1);
    for (auto a = 0; a <= _resolution; a++) {
      double distA = _lutCutoff + a * _pointDistance;
      _lut3B[a].resize(a + 1);
      for (auto b = 0; b <= a; b++) {
        double distB = _lutCutoff + b * _pointDistance;
        _lut3B[a][b].reserve(b + 1);
        for (auto c = 0; c <= b; c++) {
          double distC = _lutCutoff + c * _pointDistance;

          _lut3B[a][b].push_back(getLUTValuesAT(distA, distB, distC));
        }
      }
    }
  }

  void fill_global_krypton() {
    _lut3B_global.resize(_resolution + 1);
    for (auto a = 0; a <= _resolution; a++) {
      double distA = _lutCutoff + a * _pointDistance;
      _lut3B_global[a].resize(a + 1);
      for (auto b = 0; b <= a; b++) {
        double distB = _lutCutoff + b * _pointDistance;
        _lut3B_global[a][b].reserve(b + 1);
        for (auto c = 0; c <= b; c++) {
          double distC = _lutCutoff + c * _pointDistance;
          auto x = getLUTValuesKrypton(distA, distB, distC);
          _lut3B_global[a][b].push_back(getLUTValuesKrypton(distA, distB, distC));
          //          std::cout << "index  :" << a << ", " << b << ", " << c << "  "<< distA << ", " << distB << ", " <<
          //          distC << std::endl;
        }
      }
    }
  }

  void fill_global() {
    _lut3B_global.resize(_resolution + 1);
    for (auto a = 0; a <= _resolution; a++) {
      double distA = _lutCutoff + a * _pointDistance;
      _lut3B_global[a].resize(a + 1);
      for (auto b = 0; b <= a; b++) {
        double distB = _lutCutoff + b * _pointDistance;
        _lut3B_global[a][b].reserve(b + 1);
        for (auto c = 0; c <= b; c++) {
          double distC = _lutCutoff + c * _pointDistance;

          auto global = calculatePotEnergy(distA, distB, distC);

          auto lutVal = getLUTValuesAT(distA, distB, distC);
          std::array<double, 4> res = {lutVal[0], lutVal[1], lutVal[2], global};
          _lut3B_global[a][b].push_back(res);
        }
      }
    }
  }

  std::array<double, 3> getNextNeighbour(size_t index1, size_t index2, size_t index3) const {
    // just a function for consistency, floored values are basically next neighbours

    return _lut3B[index1][index2][index3];
  }

  std::array<double, 4> getNextNeighbour_global(size_t index1, size_t index2, size_t index3) const {
    // just a function for consistency, floored values are basically next neighbours
    return _lut3B_global[index1][index2][index3];
  }

  std::array<double, 3> getMidpoint(size_t index1, size_t index2, size_t index3) const {
    using namespace autopas::utils::ArrayMath::literals;
    // Retrieve values from the LUT by linear interpolation

    auto resultLeft = _lut3B[index1][index2][index3];

    std::array<double, 3> resultRight{};

    auto z_max = _lut3B[index1][index2].size();
    auto y_max = _lut3B[index1].size();
    auto x_max = _lut3B.size();

    if (index3 < z_max - 1) {
      resultRight = _lut3B[index1][index2][index3 + 1];
      return (resultRight + resultLeft) * 0.5;
    } else if (index2 < y_max - 1) {
      resultRight = _lut3B[index1][index2 + 1][index3];
      return (resultRight + resultLeft) * 0.5;
    }

    // handle case that we are on the very last entry
    if (index1 > x_max) {
      return _lut3B[index1][index2][index3];
    }
    resultRight = _lut3B[index1 + 1][index2][index3];

    return (resultRight + resultLeft) * 0.5;
  }

  std::array<double, 4> getMidpoint_global(size_t index1, size_t index2, size_t index3) const {
    using namespace autopas::utils::ArrayMath::literals;
    // Retrieve values from the LUT by linear interpolation

    auto resultLeft = _lut3B_global[index1][index2][index3];

    std::array<double, 4> resultRight{};

    auto z_max = _lut3B_global[index1][index2].size();
    auto y_max = _lut3B_global[index1].size();
    auto x_max = _lut3B_global.size();

    if (index3 < z_max - 1) {
      resultRight = _lut3B_global[index1][index2][index3 + 1];
      return (resultRight + resultLeft) * 0.5;
    } else if (index2 < y_max - 1) {
      resultRight = _lut3B_global[index1][index2 + 1][index3];
      return (resultRight + resultLeft) * 0.5;
    }
    // handle case that we are on the very last entry
    if (index1 > x_max) {
      return _lut3B_global[index1][index2][index3];
    }
    resultRight = _lut3B_global[index1 + 1][index2][index3];

    return (resultRight + resultLeft) * 0.5;
  }

  // potential energies

  double calculatePotEnergy(double distSquaredIJ, double distSquaredJK, double distSquaredKI) const {
    double distIJ = std::sqrt(distSquaredIJ);
    double distJK = std::sqrt(distSquaredJK);
    double distKI = std::sqrt(distSquaredKI);

    double distIJ3 = distSquaredIJ * distIJ;
    double distJK3 = distSquaredJK * distJK;
    double distKI3 = distSquaredKI * distKI;
    double KIcosIJ = (distSquaredIJ + distSquaredKI - distSquaredJK) / (2 * distIJ * distKI);
    double IJcosJK = (distSquaredIJ + distSquaredJK - distSquaredKI) / (2 * distIJ * distJK);
    double JKcosKI = (distSquaredJK + distSquaredKI - distSquaredIJ) / (2 * distJK * distKI);
    double res = 3 * (nu * (1 + (3 * KIcosIJ * IJcosJK * JKcosKI))) / (distIJ3 * distKI3 * distJK3);

    return res;
  }

  double calculatePotEnergy_Krypton(double cosines, double _nu, double allDistsTripled, double sum,
                                    double expTerm) const {
    // Add 3 * potential energy to every owned particle of the interaction.
    // Division to the correct value is handled in endTraversal().
    const double potentialEnergy = (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);
    //            return {potentialEnergy, 0};
    return potentialEnergy;
  }

  // retrieve values

  template <class Functor>
  [[nodiscard]] std::pair<const std::array<double, 3>, std::array<u_int8_t, 3>> retrieveValues(const Functor &functor,
                                                                                               double dist1,
                                                                                               double dist2,
                                                                                               double dist3) const {
    if (dist1 < _lutCutoff or dist2 < _lutCutoff or dist3 < _lutCutoff) {
      //      std::cout << "LUT3B::retrieveValues(): Particle distance smaller than lower lutCutoff(" <<
      //      std::sqrt(_lutCutoff)
      //                << ") - Calculating force factors directly!" << std::endl;
      return std::make_pair(functor.getLUTValues(dist1, dist2, dist3), std::array<u_int8_t, 3>({0, 1, 2}));
    }

    size_t index1, index2, index3;
    std::array<u_int8_t, 3> order{};  // Array to track the original indices of the distances

    if (dist1 >= dist2 && dist2 >= dist3) {
      index1 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor);
      index2 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor);
      index3 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor);
      order[0] = 0;
      order[1] = 1;
      order[2] = 2;
    } else if (dist1 >= dist3 && dist3 >= dist2) {
      index1 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor);
      index2 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor);
      index3 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor);
      order[0] = 0;
      order[1] = 2;
      order[2] = 1;
    } else if (dist2 >= dist1 && dist1 >= dist3) {
      index1 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor);
      index2 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor);
      index3 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor);
      order[0] = 1;
      order[1] = 0;
      order[2] = 2;
    } else if (dist2 >= dist3 && dist3 >= dist1) {
      index1 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor);
      index2 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor);
      index3 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor);
      order[0] = 2;
      order[1] = 0;
      order[2] = 1;
    } else if (dist3 >= dist1 && dist1 >= dist2) {
      index1 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor);
      index2 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor);
      index3 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor);
      order[0] = 1;
      order[1] = 2;
      order[2] = 0;
    } else {
      index1 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor);
      index2 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor);
      index3 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor);
      order[0] = 2;
      order[1] = 1;
      order[2] = 0;
    }

    std::array<double, 3> result;
    switch (interpolationMethod) {
      case NN:
        result = getNextNeighbour(index1, index2, index3);
        return std::make_pair(result, order);
        break;
      case MP:
        result = getMidpoint(index1, index2, index3);
        return std::make_pair(result, order);
      default:
        throw autopas::utils::ExceptionHandler::AutoPasException(
            "Choose one of the available interpolation methods for 3B LUTS");

    }
  }

  template <class Functor>
  [[nodiscard]] std::pair<const std::array<double, 4>, std::array<u_int8_t, 3>> retrieveValues_global(
      const Functor &functor, double dist1, double dist2, double dist3) const {
    if (dist1 < _lutCutoff or dist2 < _lutCutoff or dist3 < _lutCutoff) {
      //        std::cout << "LUT3B::retrieveValues(): Particle distance smaller than lower lutCutoff(" <<
      //        std::sqrt(_lutCutoff)
      //                  << ") - Calculating force factors directly!" << std::endl;

#if defined(MD_FLEXIBLE_FUNCTOR_AT_AUTOVEC) && defined(MD_FLEXIBLE_FUNCTOR_KRYPTON)
      throw std::runtime_error(
          "LUT3B can either run with MD_FLEXIBLE_FUNCTOR_KRYPTON or MD_FLEXIBLE_FUNCTOR_AT_AUTOVEC but not with both "
          "at the same time. Please turn one off");
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_AT_AUTOVEC)
      auto globs = calculatePotEnergy(dist1, dist2, dist3);
      auto vals = functor.getLUTValues(dist1, dist2, dist3);

      std::array<double, 4> res = {vals[0], vals[1], vals[2], globs};
      return std::make_pair(res, std::array<u_int8_t, 3>({0, 1, 2}));
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_KRYPTON)
      auto vals = functor.getLUTValues(dist1, dist2, dist3);
      return std::make_pair(vals, std::array<u_int8_t, 3>({0, 1, 2}));
#endif
    }

    size_t index1, index2, index3;
    std::array<u_int8_t, 3> order{};  // Array to track the original indices of the distances

    if (dist1 >= dist2 && dist2 >= dist3) {
      index1 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor);
      index2 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor);
      index3 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor);
      order[0] = 0;
      order[1] = 1;
      order[2] = 2;
    } else if (dist1 >= dist3 && dist3 >= dist2) {
      index1 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor);
      index2 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor);
      index3 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor);
      order[0] = 0;
      order[1] = 2;
      order[2] = 1;
    } else if (dist2 >= dist1 && dist1 >= dist3) {
      index1 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor);
      index2 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor);
      index3 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor);
      order[0] = 1;
      order[1] = 0;
      order[2] = 2;
    } else if (dist2 >= dist3 && dist3 >= dist1) {
      index1 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor);
      index2 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor);
      index3 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor);
      order[0] = 2;
      order[1] = 0;
      order[2] = 1;
    } else if (dist3 >= dist1 && dist1 >= dist2) {
      index1 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor);
      index2 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor);
      index3 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor);
      order[0] = 1;
      order[1] = 2;
      order[2] = 0;
    } else {
      index1 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor);
      index2 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor);
      index3 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor);
      order[0] = 2;
      order[1] = 1;
      order[2] = 0;
    }

    std::array<double, 4> result;
    switch (interpolationMethod) {
      case NN:
        result = getNextNeighbour_global(index1, index2, index3);
        return std::make_pair(result, order);
        break;
      case MP:
        result = getMidpoint_global(index1, index2, index3);
        return std::make_pair(result, order);
      default:
        throw autopas::utils::ExceptionHandler::AutoPasException(
            "Choose one of the available interpolation methods for 3B LUTS");
    }
  }

 public:
  [[nodiscard]] std::array<double, 3> getLUTValuesAT(double dist1Squared, double dist2Squared,
                                                     double dist3Squared) const {
    // Calculate prefactor
    const double allDistsSquared = dist1Squared * dist2Squared * dist3Squared;
    const double allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
    const double factor = 3.0 * nu / allDistsTo5;

    // Dot products of both distance vectors going from one particle
    const double IJDotKI = -0.5 * (dist1Squared + dist3Squared - dist2Squared);
    const double IJDotJK = -0.5 * (dist1Squared + dist2Squared - dist3Squared);
    const double JKDotKI = -0.5 * (dist2Squared + dist3Squared - dist1Squared);

    const double allDotProducts = IJDotKI * IJDotJK * JKDotKI;

    const auto forceIDirIJ = IJDotJK * JKDotKI - dist2Squared * dist3Squared + 5.0 * allDotProducts / dist1Squared;
    const auto forceJDirJK = IJDotKI * JKDotKI - dist1Squared * dist3Squared + 5.0 * allDotProducts / dist2Squared;
    const auto forceKDirKI = IJDotJK * JKDotKI - dist1Squared * dist2Squared + 5.0 * allDotProducts / dist3Squared;
    const auto forceIDirJK = IJDotKI * (IJDotJK - JKDotKI);
    const auto forceJDirKI = IJDotJK * (JKDotKI - IJDotKI);

    const auto forceIJ = (forceIDirIJ - forceIDirJK) * factor;
    const auto forceJK = (forceJDirJK - forceJDirKI) * factor;
    const auto forceKI = (forceKDirKI + forceIDirJK) * factor;

    return {forceIJ, forceJK, forceKI};
  }

  [[nodiscard]] std::array<double, 4> getLUTValuesKrypton(double distSquaredIJ, double distSquaredJK,
                                                          double distSquaredKI) const {
    // Actual distances
    const double distIJ = std::sqrt(distSquaredIJ);
    const double distJK = std::sqrt(distSquaredJK);
    const double distKI = std::sqrt(distSquaredKI);

    // Numerators of cosine representation (cos_i = (r_ij^2 + r_ik^2 - r_jk^2) / (2 * r_ij * r_ik)
    const double numKI = distSquaredIJ + distSquaredJK - distSquaredKI;
    const double numJK = distSquaredIJ + distSquaredKI - distSquaredJK;
    const double numIJ = distSquaredJK + distSquaredKI - distSquaredIJ;

    const double numerator = numKI * numJK * numIJ;

    const double allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
    const double allDists = distIJ * distJK * distKI;
    const auto allDistsTripled = allDistsSquared * allDists;

    // Gradient factors of 1. / (rrr)^3
    const double allDistsTriplesGradientIJ = 3. / (allDistsTripled * distSquaredIJ);
    const double allDistsTriplesGradientKI = -3. / (allDistsTripled * distSquaredKI);

    // Product of all cosines multiplied with 3: 3 * cos(a)cos(b)cos(c)
    const auto cosines = (3. / 8.) * numerator / allDistsSquared;
    const double cosinesGradientIJ =
        (3. / 4.) * ((numerator / distSquaredIJ - numKI * numIJ - numJK * numIJ + numJK * numKI) / allDistsSquared);
    const double cosinesGradientKI =
        (3. / 4.) * ((-numerator / distSquaredKI + numKI * numIJ - numJK * numIJ + numJK * numKI) / allDistsSquared);

    // Gradient factors corresponding to the normal ATM term
    const auto fullATMGradientIJ =
        _nu_krypton * ((1. + cosines) * allDistsTriplesGradientIJ + cosinesGradientIJ / allDistsTripled);
    const auto fullATMGradientKI =
        _nu_krypton * ((1. + cosines) * allDistsTriplesGradientKI + cosinesGradientKI / allDistsTripled);

    const auto expTerm = std::exp(-_alpha_krypton * (distIJ + distJK + distKI));

    // Calculate factors and sum for: \sum_{n=0}^5 A_{2n}(r_ij*r_jk*r_ki)^(2n/3)
    std::array<double, 6> sumFactors{};
    double sum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      sumFactors[n] = _constantsA_krypton[n] * std::pow((distIJ * distJK * distKI), 2. * n / 3.);
      sum += sumFactors[n];
    }

    // Gradient factor of the sum in ij-direction
    double ijSum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      ijSum += sumFactors[n] * (2. * n / (3. * distIJ) - _alpha_krypton);
    }

    // Gradient factor of the sum in ki-direction
    double kiSum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      kiSum += sumFactors[n] * (2. * n / (3. * distKI) - _alpha_krypton);
    }

    // Total gradient factors for the exponential term times the cosines term
    const double fullExpGradientIJ = expTerm * (-(1. + cosines) * ijSum / distIJ + cosinesGradientIJ * sum);
    const double fullExpGradientKI = expTerm * ((1. + cosines) * kiSum / distKI + cosinesGradientKI * sum);

    const auto factorIJ = -fullATMGradientIJ - fullExpGradientIJ;
    const auto factorKI = fullATMGradientKI + fullExpGradientKI;

    // for the newton part

    // Calculate all components for jk-direction
    const double allDistsTriplesGradientJK = 3. / (allDistsTripled * distSquaredJK);
    const double cosinesGradientJK =
        (3. / 4.) * ((numerator / distSquaredJK + numKI * numIJ - numJK * numIJ - numJK * numKI) / allDistsSquared);
    const auto fullATMGradientJK =
        _nu_krypton * ((1. + cosines) * allDistsTriplesGradientJK + cosinesGradientJK / allDistsTripled);

    double jkSum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      jkSum += sumFactors[n] * (2. * n / (3. * distJK) - _alpha_krypton);
    }
    double fullExpGradientJK = expTerm * (-(1. + cosines) * jkSum / distJK + cosinesGradientJK * sum);

    const auto factorJK = -fullATMGradientJK - fullExpGradientJK;

    double epot = calculatePotEnergy_Krypton(cosines, _nu_krypton, allDistsTripled, sum, expTerm);

    return {factorIJ, factorJK, factorKI, epot};
  }

  std::array<double, 3> getWeightedAverage(size_t index1, size_t index2, size_t index3, size_t actual_dist0,
                                           size_t actual_dist1, size_t actual_dist2) const {
    auto V_Middle = _lut3B[index1][index2][index3];
    auto distances_middle = getDistances(index1, index2, index3);

    std::array<double, 3> V_Above{};
    std::array<double, 3> distances_above;
    auto z_max = _lut3B[index1][index2].size();
    auto y_max = _lut3B[index1].size();
    auto x_max = _lut3B.size();

    // index is enough doesnt need specific distance since the can be directly mapped to each other and bigger index
    // always == bigger distance no use the actual distances cause need to comapre to transported distances
    if (index3 < z_max - 1) {
      V_Above = _lut3B[index1][index2][index3 + 1];
      distances_above = getDistances(index1, index2, index3 + 1);

    } else if (index2 < y_max - 1) {
      V_Above = _lut3B[index1][index2 + 1][index3];
      distances_above = getDistances(index1, index2 + 1, index3);

    } else {
      V_Above = _lut3B[index1 + 1][index2][index3];
      distances_above = getDistances(index1 + 1, index2, index3);
    }

    std::array<double, 3> V_Below{};
    std::array<double, 3> distances_below;

    if (index3 > 0) {
      V_Below = _lut3B[index1][index2][index3 - 1];
      distances_below = getDistances(index1, index2, index3 - 1);

    } else if (index2 > 0) {
      V_Below = _lut3B[index1][index2 - 1][index3];
      distances_below = getDistances(index1, index2 - 1, index3);

    } else {
      if (index1 > 0) {
        V_Below = _lut3B[index1 - 1][index2][index3];
        distances_below = getDistances(index1 - 1, index2, index3);

      } else {
        V_Below = _lut3B[0][0][0];
        distances_below = getDistances(0, 0, 0);
      }
    }

    // weights
    auto a = pow(actual_dist0 - distances_middle[0], 2) + pow(actual_dist1 - distances_middle[1], 2) +
             pow(actual_dist2 - distances_middle[2], 2);
    auto b = pow(actual_dist0 - distances_below[0], 2) + pow(actual_dist1 - distances_below[1], 2) +
             pow(actual_dist2 - distances_below[2], 2);
    auto c = pow(actual_dist0 - distances_above[0], 2) + pow(actual_dist1 - distances_above[1], 2) +
             pow(actual_dist2 - distances_above[2], 2);

    auto w1 = 1 / a;
    auto w2 = 1 / b;
    auto w3 = 1 / c;

    auto w_total = w1 + w2 + w3;

    w1 = w1 / w_total;
    w2 = w2 / w_total;
    w3 = w3 / w_total;

    auto v0 = (V_Middle[0] * w1) + (V_Below[0] * w2) + (V_Above[0] * w3);
    auto v1 = (V_Middle[1] * w1) + (V_Below[1] * w2) + (V_Above[1] * w3);
    auto v2 = (V_Middle[2] * w1) + (V_Below[2] * w2) + (V_Above[2] * w3);

    return std::array<double, 3>{v0, v1, v2};
  }

  std::array<double, 4> getWeightedAverage_global(size_t index1, size_t index2, size_t index3, size_t actual_dist0,
                                                  size_t actual_dist1, size_t actual_dist2) const {
    auto V_Middle = _lut3B_global[index1][index2][index3];
    auto distances_middle = getDistances(index1, index2, index3);

    std::array<double, 4> V_Above{};
    std::array<double, 3> distances_above;
    auto z_max = _lut3B_global[index1][index2].size();
    auto y_max = _lut3B_global[index1].size();
    auto x_max = _lut3B_global.size();

    // index is enough doesnt need specific distance since the can be directly mapped to each other and bigger index
    // always == bigger distance no use the actual distances cause need to comapre to transported distances
    if (index3 < z_max - 1) {
      V_Above = _lut3B_global[index1][index2][index3 + 1];
      distances_above = getDistances(index1, index2, index3 + 1);

    } else if (index2 < y_max - 1) {
      V_Above = _lut3B_global[index1][index2 + 1][index3];
      distances_above = getDistances(index1, index2 + 1, index3);

    } else {
      V_Above = _lut3B_global[index1 + 1][index2][index3];
      distances_above = getDistances(index1 + 1, index2, index3);
    }

    std::array<double, 4> V_Below{};
    std::array<double, 3> distances_below;

    if (index3 > 0) {
      V_Below = _lut3B_global[index1][index2][index3 - 1];
      distances_below = getDistances(index1, index2, index3 - 1);

    } else if (index2 > 0) {
      V_Below = _lut3B_global[index1][index2 - 1][index3];
      distances_below = getDistances(index1, index2 - 1, index3);

    } else {
      if (index1 > 0) {
        V_Below = _lut3B_global[index1 - 1][index2][index3];
        distances_below = getDistances(index1 - 1, index2, index3);

      } else {
        V_Below = _lut3B_global[0][0][0];
        distances_below = getDistances(0, 0, 0);
      }
    }

    // weights
    auto a = pow(actual_dist0 - distances_middle[0], 2) + pow(actual_dist1 - distances_middle[1], 2) +
             pow(actual_dist2 - distances_middle[2], 2);
    auto b = pow(actual_dist0 - distances_below[0], 2) + pow(actual_dist1 - distances_below[1], 2) +
             pow(actual_dist2 - distances_below[2], 2);
    auto c = pow(actual_dist0 - distances_above[0], 2) + pow(actual_dist1 - distances_above[1], 2) +
             pow(actual_dist2 - distances_above[2], 2);

    auto w1 = 1 / a;
    auto w2 = 1 / b;
    auto w3 = 1 / c;

    auto w_total = w1 + w2 + w3;

    w1 = w1 / w_total;
    w2 = w2 / w_total;
    w3 = w3 / w_total;

    auto v0 = (V_Middle[0] * w1) + (V_Below[0] * w2) + (V_Above[0] * w3);
    auto v1 = (V_Middle[1] * w1) + (V_Below[1] * w2) + (V_Above[1] * w3);
    auto v2 = (V_Middle[2] * w1) + (V_Below[2] * w2) + (V_Above[2] * w3);
    auto glob = (V_Middle[3] * w1) + (V_Below[3] * w2) + (V_Above[3] * w3);

    return std::array<double, 4>{v0, v1, v2, glob};
  }

 private:
  bool KRYPTON;
  bool useGlobal;
  mdLib::InterpolationMethod interpolationMethod;

  std::vector<std::vector<std::vector<std::array<double, 3>>>> _lut3B;

  std::vector<std::vector<std::vector<std::array<double, 4>>>> _lut3B_global;
  double nu;

  double _nu_krypton = 1.61525e6;

  double _alpha_krypton = 1.378382;
  std::vector<autopas::utils::Timer> lutTimer;

  // Units: {K, K*A^-2, K*A^-4, K*A^-6, K*A^-8, K*A^-10}
  std::array<double, 6> _constantsA_krypton = {-0.3081304e8, -0.3519442e8, 0.4928052e7, -0.2182411e6, 0.343088e4, 0.0};

  // helper functions & getters & setters

  std::array<double, 3> getDistances(size_t i1, size_t i2, size_t i3) const {
    double d1 = (i1 / _lutFactor) + _lutCutoff;
    double d2 = (i2 / _lutFactor) + _lutCutoff;
    double d3 = (i3 / _lutFactor) + _lutCutoff;

    return std::array<double, 3>{d1, d2, d3};
  }

 public:
  void setNu(double _nu) { nu = _nu; }
  void setLutTimer(const std::vector<autopas::utils::Timer> &lutTimer) { this->lutTimer = lutTimer; }
  void setAlphaKrypton(double alphaKrypton) { _alpha_krypton = alphaKrypton; }
  void setNuKrypton(double nuKrypton) { _nu_krypton = nuKrypton; }
  void setConstantsAKrypton(const std::array<double, 6> &constantsAKrypton) { _constantsA_krypton = constantsAKrypton; }
};

}  // namespace mdLib