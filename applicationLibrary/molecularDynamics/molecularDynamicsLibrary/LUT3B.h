//
// Created by sliep on 05.05.2025.
//
#pragma once

#include "ParentLUT.h"
#include "autopas/utils/Math.h"
#include "autopas/utils/Timer.h"

namespace mdLib {
class LUT3B : public ParentLUT {
 public:
  LUT3B(int resolution, double cutoffSquared, double Nu) : ParentLUT(resolution, cutoffSquared) {
    nu = Nu;

    fill(cutoffSquared, false);
    std::cout << "FILLING LUT 3 args " << std::endl;
  }

  LUT3B(int resolution, double cutoffSquared) : ParentLUT(resolution, cutoffSquared) {
    if (resolution == 0) {
      return;
    }

    fill(cutoffSquared, false);
    std::cout << "FILLING LUT 2 rgs" << std::endl;
  }

  LUT3B(int resolution, double cutoffSquared, bool global, double _nu = 0) : ParentLUT(resolution, cutoffSquared) {
    nu = _nu;
    if (resolution == 0) {
      return;
    }

    fill(cutoffSquared, global);
    std::cout << "FILLING LUT 2 rgs" << std::endl;
  }

  template <class Functor>
  [[nodiscard]] std::pair<const std::array<double, 3>, std::array<u_int8_t, 3>>
  retrieveValues(const Functor &functor, double dist1, double dist2, double dist3) const {

    if (dist1 < _lutCutoff or dist2 < _lutCutoff or dist3 < _lutCutoff) {
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

    // if Linear Interpolation
    // auto result = getLinear(index1, index2, index3);
    auto result = getWeightedAverage(index1, index2, index3, dist1, dist2, dist3);

    // space for other interpolation
    // TODO NN does not work for Krypton yet
    //            auto result = getNextNeighbour(index1, index2, index3);
    //            auto result = getNextNeighbour_krypton(index1, index2, index3);

    auto resultPair = std::make_pair(result, order);
    return resultPair;
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

    std::array<double, 3> resultBelow{};
    if (index3 > 0) {
      resultBelow = _lut3B[index1][index2][index3 - 1];
      V_Below = _lut3B[index1][index2][index3 - 1];
      distances_below = getDistances(index1, index2, index3 - 1);

    } else if (index2 > 0) {
      resultBelow = _lut3B[index1][index2 - 1][index3];
      V_Below = _lut3B[index1][index2 - 1][index3];
      distances_below = getDistances(index1, index2 - 1, index3);

    } else {
      if (index1 > 0) {
        resultBelow = _lut3B[index1 - 1][index2][index3];
        V_Below = _lut3B[index1 - 1][index2][index3];
        distances_below = getDistances(index1 - 1, index2, index3);

      } else {
        resultBelow = _lut3B[0][0][0];
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

  std::array<double, 3> getLinear(size_t index1, size_t index2, size_t index3) const {
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
    resultRight = _lut3B[index1 + 1][index2][index3];

    std::cout << "in LINEAR ";

    return (resultRight + resultLeft) * 0.5;
  }

  //        template<class Functor>
  //        void fill(const Functor &functor, double cutoffSquared, bool withGlobals) {
  void fill(double cutoffSquared, bool withGlobals) {
    //            bool withGlobals = false;

    //        lutTimer[0].start(); //0th being for fill time
    if (!withGlobals) {
      //                fill_plain(functor, cutoffSquared);
      fill_plain(cutoffSquared);

      // fill_plain_krypton( cutoffSquared);
    }
    if (withGlobals) {
      //  fill_gobal(functor, cutoffSquared);
      //                fill_global_all(functor, cutoffSquared);
      fill_global_all(cutoffSquared);
    }
    //      lutTimer[0].stop();
    std::cout << "filled ";
  }

  std::array<double, 3> getNextNeighbour(size_t index1, size_t index2, size_t index3) const {
    // just a function for consistency, floored values are basically next neighbours
    return _lut3B[index1][index2][index3];
  }

  std::array<double, 4> getNextNeighbour_krypton(size_t index1, size_t index2, size_t index3) const {
    // just a function for consistency, floored values are basically next neighbours
    return _lut3B_krypton[index1][index2][index3];
  }

  // alternative global LUT

  std::array<double, 2> calculateGlobalFac(double distSquaredIJ, double distSquaredJK, double distSquaredKI) const {
    const double allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
    const double allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
    //      const double factor = 3.0 * nu / allDistsTo5;
    const double factor = 3.0 * 1.6152500E-3 / allDistsTo5;

    return {allDistsSquared, factor};
  }

  std::array<double, 2> calculateFullGlobal(double distSquaredIJ, double distSquaredJK, double distSquaredKI) const {
    const double allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
    const double allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
    //      const double factor = 3.0 * nu / allDistsTo5;
    const double factor = 3.0 * 1.6152500E-3 / allDistsTo5;

    //        // Dot products of both distance vectors going from one particle  -> over law of cosine
    const double IJDotKI = 0.5 * (distSquaredIJ + distSquaredKI - distSquaredJK);
    const double IJDotJK = 0.5 * (distSquaredIJ + distSquaredJK - distSquaredKI);
    const double JKDotKI = 0.5 * (distSquaredJK + distSquaredKI - distSquaredIJ);
    const double allDotProducts = IJDotKI * IJDotJK * JKDotKI;
    const double potentialEnergy3 = factor * (allDistsSquared - 3.0 * allDotProducts);

    return {potentialEnergy3, 0};
  }

  std::array<double, 2> calculatePotEnergyTest(double distSquaredIJ, double distSquaredJK, double distSquaredKI) const {
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

    return {res, 0};
  }

  double
  //        calculatePotEnergyTest_krypton(double distSquaredIJ, double distSquaredJK, double distSquaredKI) const {
  calculatePotEnergyTest_krypton(double cosines, double _nu, double allDistsTripled, double sum, double expTerm) const {
    // Add 3 * potential energy to every owned particle of the interaction.
    // Division to the correct value is handled in endTraversal().
    const double potentialEnergy = (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);
    //            return {potentialEnergy, 0};
    return potentialEnergy;
  }

  //        template<class Functor>
  //        void fill_plain(const Functor &functor, double cutoffSquared) {
  void fill_plain(double cutoffSquared) {
    _lut3B.resize(_resolution + 1);
    for (auto a = 0; a <= _resolution; a++) {
      double distA = _lutCutoff + a * _pointDistance;
      _lut3B[a].resize(a + 1);
      for (auto b = 0; b <= a; b++) {
        double distB = _lutCutoff + b * _pointDistance;
        _lut3B[a][b].reserve(b + 1);
        for (auto c = 0; c <= b; c++) {
          double distC = _lutCutoff + c * _pointDistance;
          //                        _lut3B[a][b].push_back(functor.getLUTValues(distA, distB, distC));
          _lut3B[a][b].push_back(getLUTValuesAT(distA, distB, distC));
        }
      }
    }
  }

  //        template<class Functor>
  //        void fill_plain_krypton(const Functor &functor, double cutoffSquared) {
  void fill_plain_krypton(double cutoffSquared) {
    _lut3B.resize(_resolution + 1);
    for (auto a = 0; a <= _resolution; a++) {
      double distA = _lutCutoff + a * _pointDistance;
      _lut3B[a].resize(a + 1);
      for (auto b = 0; b <= a; b++) {
        double distB = _lutCutoff + b * _pointDistance;
        _lut3B[a][b].reserve(b + 1);
        for (auto c = 0; c <= b; c++) {
          double distC = _lutCutoff + c * _pointDistance;
          auto x = getLUTValuesKrypton(distA, distB, distC);
          _lut3B[a][b].push_back(getLUTValuesKrypton(distA, distB, distC));
          //                _lut3B_krypton[a][b].push_back(getLUTValuesKrypton(distA, distB, distC));
        }
      }
    }
  }

  //        template<class Functor>
  //        void fill_gobal(const Functor &functor, double cutoffSquared) {
  void fill_gobal(double cutoffSquared) {
    _lut3B_global.resize(_resolution + 1);
    for (auto a = 0; a <= _resolution; a++) {
      double distA = _lutCutoff + a * _pointDistance;
      _lut3B_global[a].resize(a + 1);
      for (auto b = 0; b <= a; b++) {
        double distB = _lutCutoff + b * _pointDistance;
        _lut3B_global[a][b].reserve(b + 1);
        for (auto c = 0; c <= b; c++) {
          double distC = _lutCutoff + c * _pointDistance;
          auto global = calculateGlobalFac(distA, distB, distC);
          //                        auto lutVal = functor.getLUTValues(distA, distB, distC);
          auto lutVal = getLUTValuesAT(distA, distB, distC);
          std::array<double, 5> res = {lutVal[0], lutVal[1], lutVal[2], global[0], global[1]};
          _lut3B_global[a][b].push_back(res);
        }
      }
    }

    std::cout << "filled ";
  }

  //        template<class Functor>
  //        void fill_global_all(const Functor &functor, double cutoffSquared) {
  void fill_global_all(double cutoffSquared) {
    _lut3B_global.resize(_resolution + 1);
    for (auto a = 0; a <= _resolution; a++) {
      double distA = _lutCutoff + a * _pointDistance;
      _lut3B_global[a].resize(a + 1);
      for (auto b = 0; b <= a; b++) {
        double distB = _lutCutoff + b * _pointDistance;
        _lut3B_global[a][b].reserve(b + 1);
        for (auto c = 0; c <= b; c++) {
          double distC = _lutCutoff + c * _pointDistance;
          //            auto global = calculateFullGlobal(distA, distB, distC);
          // testing
          auto global = calculatePotEnergyTest(distA, distB, distC);
          //                        auto lutVal = functor.getLUTValues(distA, distB, distC);
          auto lutVal = getLUTValuesAT(distA, distB, distC);
          std::array<double, 5> res = {lutVal[0], lutVal[1], lutVal[2], global[0], 0};
          _lut3B_global[a][b].push_back(res);
        }
      }
    }

    std::cout << "filled ";
  }

  template <class Functor>
  //    [[nodiscard]] std::pair<const std::array<double, 3>, std::pair<std::array<u_int8_t, 3>, std::array<double,2>>
  //    retrieveValues_global(const Functor &functor, double dist1, double dist2, double dist3) const {
  [[nodiscard]] std::pair<const std::array<double, 5>, std::array<u_int8_t, 3>> retrieveValues_global(
      const Functor &functor, double dist1, double dist2, double dist3) const {
    //    using namespace autopas::utils::ArrayMath::literals;
    if (dist1 < _lutCutoff or dist2 < _lutCutoff or dist3 < _lutCutoff) {
      auto globs = calculateGlobalFac(dist1, dist2, dist3);
      //                auto vals = functor.getLUTValues(dist1, dist2, dist3);
      auto vals = getLUTValuesAT(dist1, dist2, dist3);
    std:
      std::array<double, 5> res = {vals[0], vals[1], vals[2], globs[0], globs[1]};
      return std::make_pair(res, std::array<u_int8_t, 3>({0, 1, 2}));
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

    // if Linear Interpolation
    //     auto result = getLinear(index1, index2, index3);

    // space for other interpolation
    auto result = getNextNeighbour_global(index1, index2, index3);
    //      auto val ={result[0], result[1], result[2]};
    //      auto fac = {result[3], result[4]};
    //
    //      auto pait = std::make_pair(order, fac);
    //      auto resultPair = std::make_pair(val, pait);
    auto resultPair = std::make_pair(result, order);

    //      std::cout << "RETRIEVE VALUE ";

    return resultPair;
  }

  std::array<double, 5> getNextNeighbour_global(size_t index1, size_t index2, size_t index3) const {
    // just a function for consistency, floored values are basically next neighbours
    return _lut3B_global[index1][index2][index3];
  }

  //        template<class Functor>
  void setNu(double _nu) { nu = _nu; }

  void setLutTimer(const std::vector<autopas::utils::Timer> &lutTimer) { this->lutTimer = lutTimer; }
  void setAlphaKrypton(double alphaKrypton) { _alpha_krypton = alphaKrypton; }
  void setNuKrypton(double nuKrypton) { _nu_krypton = nuKrypton; }
  void setConstantsAKrypton(const std::array<double, 6> &constantsAKrypton) { _constantsA_krypton = constantsAKrypton; }

 private:
  std::vector<std::vector<std::vector<std::array<double, 3>>>> _lut3B;
  std::vector<std::vector<std::vector<std::array<double, 4>>>> _lut3B_krypton;
  std::vector<std::vector<std::vector<std::array<double, 5>>>> _lut3B_global;
  double nu;
  //         double _nu_krypton = 1.61525e-3;   // K·nm^9  // K*A^9
  double _nu_krypton = 1.61525e6;
  ;
  //        const double _alpha = 1.378382;  // A^-1
  //         double _alpha_krypton = 13.78382;  // A^-1  changed to nm
  double _alpha_krypton = 1.378382;
  std::vector<autopas::utils::Timer> lutTimer;

  // Units: {K, K*A^-2, K*A^-4, K*A^-6, K*A^-8, K*A^-10}
  std::array<double, 6> _constantsA_krypton = {-0.3081304e8, -0.3519442e8, 0.4928052e7, -0.2182411e6, 0.343088e4, 0.0};

  // helper functions

  std::array<double, 3> getDistances(size_t i1, size_t i2, size_t i3) const {
    double d1 = (i1 / _lutFactor) + _lutCutoff;
    double d2 = (i2 / _lutFactor) + _lutCutoff;
    double d3 = (i3 / _lutFactor) + _lutCutoff;

    return std::array<double, 3>{d1, d2, d3};
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

  [[nodiscard]] std::array<double, 3> getLUTValuesKrypton2(double dist1Squared, double dist2Squared,
                                                           double dist3Squared) const {
    double distIJ = std::sqrt(dist1Squared);
    double distJK = std::sqrt(dist2Squared);
    double distKI = std::sqrt(dist3Squared);

    // Calculate prefactor
    const auto allDist = distIJ * distJK * distKI;
    const double allDistsSquared = dist1Squared * dist2Squared * dist3Squared;
    const double allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
    //            const auto allDistsTripled = allDistsSquared * allDists;

    const auto allDistCubed = (allDist) * (allDist) * (allDist);

    const double factor = _nu_krypton / allDistCubed;  // C_ATM / (R1R2R3)^3

    // Gradient factors of 1. / (rrr)^3
    const double allDistsTriplesGradientIJ = -3. * (distJK * distKI) / (allDist);
    const double allDistsTriplesGradientKI = -3. * (distJK * distIJ) / (allDist);
    const double allDistsTriplesGradientJJ = -3. * (distIJ * distIJ) / (allDist);

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

    double dist_sum = distIJ + distKI + distJK;
    double dist_prod = distIJ * distKI * distJK;

    // Calculate factors and sum for: \sum_{n=0}^5 A_{2n}(r_ij*r_jk*r_ki)^(2n/3)
    std::array<double, 6> sumFactors{};
    double sum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      sumFactors[n] = _constantsA_krypton[n] * std::pow((distIJ * distJK * distKI), 2. * n / 3.);
      sum += sumFactors[n];
    }

    const auto expTerm = std::exp(-1 * _alpha_krypton * (distIJ + distJK + distKI));

    const auto expandedTerm = expTerm * sum;
    const auto A = factor + expandedTerm;

    auto devAatm_IJ = -3 * (_nu_krypton / (allDistCubed * distIJ));
    auto devAatm_KI = -3 * (_nu_krypton / (allDistCubed * distKI));
    auto devAatm_JK = -3 * (_nu_krypton / (allDistCubed * distJK));

    const double numKI = dist1Squared + dist2Squared - dist3Squared;
    const double numJK = dist1Squared + dist2Squared - dist3Squared;
    const double numIJ = dist3Squared + dist2Squared - dist1Squared;

    const double numerator = numKI * numJK * numIJ;
    auto cosines = (3. / 8.) * numerator / allDistsSquared;

    // Derivatives of expandedTerm
    //  Gradient factor of the sum in ij-direction
    double ijSum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      ijSum += sumFactors[n] * (2. * n / (3. * distIJ) - _alpha_krypton);
    }
    // Gradient factor of the sum in ki-direction
    double kiSum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      kiSum += sumFactors[n] * (2. * n / (3. * distKI) - _alpha_krypton);
    }

    // Gradient factor of the sum in jk-direction
    double jkSum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      jkSum += sumFactors[n] * (2. * n / (3. * distJK) - _alpha_krypton);
    }

    double ijSum_dev = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      ijSum_dev += sumFactors[n] * (-1 * _alpha_krypton - (2. * n / (3. * distIJ)));
    }
    double kiSum_dev = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      kiSum_dev += sumFactors[n] * (-1 * _alpha_krypton - (2. * n / (3. * distKI)));
    }

    double jkSum_dev = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      jkSum_dev += sumFactors[n] * (-1 * _alpha_krypton - (2. * n / (3. * distJK)));
    }

    //            auto devCorr_IJ =  -1  * _nu * expTerm* sum + expTerm *

    //            auto devIJ = -1*  _alpha*expTerm * sum + expTerm *ijSum;
    auto devIJ = devAatm_IJ + expTerm * ijSum_dev;
    //            auto devIJ =  (-(1. + cosines) * ijSum / distIJ);
    //            auto devKI = -_alpha*expTerm * sum + expTerm *kiSum;
    auto devKI = devAatm_KI + expTerm * kiSum_dev;
    //            auto devKI =  ((1. + cosines) * kiSum / distKI);
    //            auto devJK = -_alpha*expTerm * sum + expTerm *jkSum;
    auto devJK = devAatm_JK + expTerm * jkSum_dev;
    //            auto devJK =  (-(1. + cosines) * jkSum / distJK);
    return {devIJ, devKI, devJK};
  }

  [[nodiscard]] std::array<double, 3> getLUTValuesKrypton(double distSquaredIJ, double distSquaredKI,
                                                          double distSquaredJK) const {
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

    const auto factorForceJDirectionIJ = (fullATMGradientIJ + fullExpGradientIJ);
    const auto factorForceJDirectionKI = (fullATMGradientKI + fullExpGradientKI);

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

    auto factorForceJDirectionJK = (fullATMGradientJK + fullExpGradientJK);

    //             double epot =     calculatePotEnergyTest_krypton( cosines,  _nu_krypton,  allDistsTripled,  sum,
    //             expTerm);

    return {-factorForceJDirectionIJ, factorForceJDirectionKI, factorForceJDirectionJK};
    //          return { epot, 0, 0};
  }
};

}  // namespace mdLib