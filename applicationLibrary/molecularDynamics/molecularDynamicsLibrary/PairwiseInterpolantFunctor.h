/**
 * @file PairwiseInterpolantFunctor
 *
 * @date 25 April 2025
 * @author Luis Gall
 */

#pragma once

#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace mdLib {

template <class PairwiseKernel, class Particle_T,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>
class PairwiseInterpolantFunctor
    : public autopas::PairwiseFunctor<
          Particle_T, PairwiseInterpolantFunctor<PairwiseKernel, Particle_T, useNewton3,
                                                 calculateGlobals, countFLOPs, relevantForTuning>> {
  using SoAArraysType = typename Particle_T::SoAArraysType;

  using SoAFloatPrecision = typename Particle_T::ParticleSoAFloatPrecision;

 public:
  PairwiseInterpolantFunctor() = delete;

 private:
  explicit PairwiseInterpolantFunctor(PairwiseKernel f, double cutoff, double a, std::vector<size_t> numNodes,
                                      std::vector<double> intervalSplits, void * /*dummy*/)
      : autopas::PairwiseFunctor<
            Particle_T, PairwiseInterpolantFunctor<PairwiseKernel, Particle_T, useNewton3,
                                                   calculateGlobals, countFLOPs, relevantForTuning>>(cutoff),
        _cutoffSquared{cutoff * cutoff},
        _numNodes{numNodes},
        _intervalSplits{intervalSplits},
        _a{a},
        _b{cutoff},
        _potentialEnergySum{0.},
        _absErrorSum{0.},
        _relErrorSum{0.},
        _virialSum{0., 0., 0.},
        _postProcessed{false},
        _kernel{f} {
    if (_numNodes.size() != _intervalSplits.size() + 1) {
      throw autopas::utils::ExceptionHandler::AutoPasException("Size of node list must be = size of intervals + 1");
    }

    if constexpr (calculateGlobals) {
      _aosThreadDataGlobals.resize(autopas::autopas_get_max_threads());
    }

    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs.resize(autopas::autopas_get_max_threads());
    }

    _distanceStatistics.resize(autopas::autopas_get_max_threads());

    int interval = 0;
    for (size_t node : _numNodes) {
      _coefficients.insert(std::pair<size_t, std::vector<double>>{interval++, std::vector<double>(node)});
    }
    constructPolynomial();
  }

  double evaluateChebyshev(int n, double x) {
    if (n == 0) return 1.;
    if (n == 1) return x;

    return 2 * x * evaluateChebyshev(n - 1, x) - evaluateChebyshev(n - 2, x);
  }

  double evaluateChebyshev(int n, double x, const std::vector<double> &coeffs) {
    if (n == 0) return coeffs[0];
    if (n == 1) return x * coeffs[1];

    return 2 * x * evaluateChebyshev(n - 1, x) * coeffs[n - 1] - evaluateChebyshev(n - 2, x) * coeffs[n - 2];
  }

  float evaluateChebyshevFast(double x, int n, const std::vector<double> &coeff) {
    /*Following the Clenshaw Algorithm*/
    double b_n = 0;
    double b_n1 = 0;
    double b_n2 = 0;
    for (int k = n - 1; k > 0; --k) {
      b_n = coeff[k] + 2. * x * b_n1 - b_n2;
      b_n2 = b_n1;
      b_n1 = b_n;
    }

    /*Special Treatment of k=0 because of c0 (has 1/2 in the formula)*/
    b_n = 2. * coeff[0] + 2. * x * b_n1 - b_n2;
    return (b_n - b_n2) / 2.;
  }

  double mapToInterval(double x, double a, double b) {
    double intermediate = x + 1;
    double newX = a + ((b - a) * intermediate) / 2.;
    return newX;
  }

  double mapToCheb(double x, double a, double b) {
    double intermediate = x - a;
    double newX = 2. * intermediate / (b - a) - 1.;
    return newX;
  }

  void dct(const std::vector<double> &values, int interval) {
    for (int i = 0; i < _numNodes[interval]; ++i) {
      float coefficient = 0.;
      for (int k = 0; k < _numNodes[interval]; ++k) {
        coefficient += values[k] * std::cos(PI * i * (2 * k + 1) / (2 * _numNodes[interval]));
      }
      _coefficients.at(interval)[i] = coefficient * 2. / _numNodes[interval];
    }
    _coefficients.at(interval)[0] = _coefficients.at(interval)[0] / 2.;
  }
 
  void constructPolynomial() {

    double a = _a;
    for (int interval = 0; interval < _numNodes.size(); ++interval) {
      std::vector<double> values{};
      values.resize(_numNodes[interval]);

      double b = _intervalSplits.size() > interval ? _intervalSplits[interval] : _b;

      // evaluate Functor at optimal Chebyshev nodes
      for (int i = 0; i < _numNodes[interval]; ++i) {
        double x = ((2 * i + 1.) / (2. * (_numNodes[interval]))) * PI;
        double node = std::cos(x);
        double d = mapToInterval(node, a, b);

        double value = _kernel.calculatePairDerivative(d);
        values[i] = value;
      }

      dct(values, interval);
      a = b;
    }
  }

 public:
  explicit PairwiseInterpolantFunctor(PairwiseKernel f, double cutoff, double a, std::vector<size_t> nNodes,
                                      std::vector<double> intervalSplits)
      : PairwiseInterpolantFunctor(f, cutoff, a, nNodes, intervalSplits, nullptr) {
  }

  std::string getName() final { return "PairwiseInterpolantFunctor"; }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;

    /* Filter unnecessary force computations */
    if (i.isDummy() or j.isDummy()) {
      return;
    }

    const auto threadnum = autopas::autopas_get_thread_num();

    if constexpr (countFLOPs) {
      ++_aosThreadDataFLOPs[threadnum].numDistCalls;
    }

    auto dr = i.getR() - j.getR();
    double dr2 = autopas::utils::ArrayMath::dot(dr, dr);

    if (dr2 > _cutoffSquared) {
      return;
    }

    /* Evaluate Polynomial */
    double d = std::sqrt(dr2);

#if defined(MD_FLEXIBLE_BENCHMARK_INTERPOLANT_ACCURACY)
    _distanceStatistics[threadnum].push_back(d);
#endif

    int interval = 0;
    double a = _a;
    double b = _b;
    /* determine interval*/
    for (; interval < _intervalSplits.size(); ++interval) {
      if (d < _intervalSplits[interval]) {
        b = _intervalSplits[interval];
        break;
      } else {
        a = _intervalSplits[interval];
      }
    }

    double fac = evaluateChebyshevFast(mapToCheb(d, a, b), _numNodes.at(interval), _coefficients.at(interval));
    double interpolationError = 0.;
    double relInterpolationError = 0.;
    auto f = dr * fac;

#if defined(MD_FLEXIBLE_BENCHMARK_INTERPOLANT_ACCURACY)
    double real_fac = _kernel.calculatePairDerivative(d);
    interpolationError = std::abs(real_fac - fac);
    if (real_fac != 0.) {
      relInterpolationError = std::abs(interpolationError / real_fac);
    }
#endif

    i.addF(f);
    if (newton3) {
      j.subF(f);
    }

    if constexpr (countFLOPs) {
      if (newton3) {
        ++_aosThreadDataFLOPs[threadnum].numKernelCallsN3;
      } else {
        ++_aosThreadDataFLOPs[threadnum].numKernelCallsNoN3;
      }
    }

    if constexpr (calculateGlobals) {
      // TODO: delegate potential energy to kernel?
      double potentialEnergy = 0.;
      auto virial = dr * f;
      if (i.isOwned()) {
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadDataGlobals[threadnum].virialSum += virial;
        _aosThreadDataGlobals[threadnum].absErrorSum += interpolationError;
        _aosThreadDataGlobals[threadnum].relErrorSum += relInterpolationError;
      }
      if (newton3 and j.isOwned()) {
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadDataGlobals[threadnum].virialSum += virial;
        _aosThreadDataGlobals[threadnum].absErrorSum += interpolationError;
        _aosThreadDataGlobals[threadnum].relErrorSum += relInterpolationError;
      }
    }
  }

  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {}

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {}

  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        const bool newton3) final {}

  void setParticleProperties(SoAFloatPrecision epsilon24, SoAFloatPrecision sigmaSquared) {}

  constexpr static auto getNeededAttr() {
    return std::array<typename Particle_T::AttributeNames, 9>{Particle_T::AttributeNames::id,
                                                              Particle_T::AttributeNames::posX,
                                                              Particle_T::AttributeNames::posY,
                                                              Particle_T::AttributeNames::posZ,
                                                              Particle_T::AttributeNames::forceX,
                                                              Particle_T::AttributeNames::forceY,
                                                              Particle_T::AttributeNames::forceZ,
                                                              Particle_T::AttributeNames::typeId,
                                                              Particle_T::AttributeNames::ownershipState};
  }

  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle_T::AttributeNames, 6>{
        Particle_T::AttributeNames::id,     Particle_T::AttributeNames::posX,
        Particle_T::AttributeNames::posY,   Particle_T::AttributeNames::posZ,
        Particle_T::AttributeNames::typeId, Particle_T::AttributeNames::ownershipState};
  }

  constexpr static auto getComputedAttr() {
    return std::array<typename Particle_T::AttributeNames, 3>{
        Particle_T::AttributeNames::forceX, Particle_T::AttributeNames::forceY, Particle_T::AttributeNames::forceZ};
  }

  constexpr static bool getMixing() { return false; }

  void initTraversal() final {
    _potentialEnergySum = 0.;
    _absErrorSum = 0.;
    _relErrorSum = 0.;
    _virialSum = {0., 0., 0.};
    _postProcessed = false;
    if constexpr (calculateGlobals) {
      for (auto &data : _aosThreadDataGlobals) {
        data.setZero();
      }
    }
    if constexpr (countFLOPs) {
      for (auto &data : _aosThreadDataFLOPs) {
        data.setZero();
      }
    }
    for (auto &stat : _distanceStatistics) {
      stat.clear();
    }
  }

  void endTraversal(bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;

    if (_postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
    }
    if (calculateGlobals) {
      for (const auto &data : _aosThreadDataGlobals) {
        _potentialEnergySum += data.potentialEnergySum;
        _virialSum += data.virialSum;
        _absErrorSum += data.absErrorSum;
        _relErrorSum += data.relErrorSum;
      }
      // For each interaction, we added the full contribution for both particles. Divide by 2 here, so that each
      // contribution is only counted once per pair.
      _potentialEnergySum *= 0.5;
      _virialSum *= 0.5;

      size_t numKernelCalls = 0;

      for (const auto &data : _aosThreadDataFLOPs) {
        numKernelCalls += data.numKernelCallsN3;
        numKernelCalls += data.numKernelCallsNoN3;
      }

      // We have always calculated 6*potentialEnergy, so we divide by 6 here!
      _potentialEnergySum /= 6.;
      _postProcessed = true;

      size_t number = 0;
      std::vector<double> cuts{};
      for (double cut = _a; cut <= _b; cut+= 0.25) {
        cuts.push_back(cut);
      }
      double prev_cut = 0.;
      for (auto cut : cuts) {
        number = 0;
        for (const auto &stat : _distanceStatistics) {
          number += std::count_if(stat.begin(), stat.end(),
                                  [cut, prev_cut](double value) { return value <= cut && value > prev_cut; });
        }
        AutoPasLog(DEBUG, "Distances below {} : {}", cut, number);
        prev_cut = cut;
      }
      size_t total = 0;
      for (auto &stat : _distanceStatistics) {
        total += stat.size();
      }
      AutoPasLog(DEBUG, "Total number distances: {}", total);

      AutoPasLog(DEBUG, "Final potential energy {}", _potentialEnergySum);
      AutoPasLog(DEBUG, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
      AutoPasLog(DEBUG, "Final absolute error   {}", _absErrorSum);
      AutoPasLog(DEBUG, "Final relative error   {}", _relErrorSum);
      AutoPasLog(DEBUG, "Number of kernel calls {}", numKernelCalls);
    }
  }

  double getPotentialEnergy() {
    if (not calculateGlobals) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Trying to get potential energy even though calculateGlobals is false. If you want this functor to calculate "
          "global "
          "values, please specify calculateGlobals to be true.");
    }
    if (not _postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Cannot get potential energy, because endTraversal was not called.");
    }
    return _potentialEnergySum;
  }

  double getVirial() {
    if (not calculateGlobals) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Trying to get virial even though calculateGlobals is false. If you want this functor to calculate global "
          "values, please specify calculateGlobals to be true.");
    }
    if (not _postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Cannot get virial, because endTraversal was not called.");
    }
    return _virialSum[0] + _virialSum[1] + _virialSum[2];
  }

  [[nodiscard]] size_t getNumFLOPs() const override {
    if constexpr (countFLOPs) {
      // TODO: implement
      return 0;
    } else {
      // This is needed because this function still gets called with FLOP logging disabled, just nothing is done with it
      return std::numeric_limits<size_t>::max();
    }
  }

  [[nodiscard]] double getHitRate() const override {
    if constexpr (countFLOPs) {
      const size_t numDistCallsAcc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numDistCalls; });
      const size_t numKernelCallsN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsN3; });
      const size_t numKernelCallsNoN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsNoN3; });

      return (static_cast<double>(numKernelCallsNoN3Acc) + static_cast<double>(numKernelCallsN3Acc)) /
             (static_cast<double>(numDistCallsAcc));
    } else {
      // This is needed because this function still gets called with FLOP logging disabled, just nothing is done with it
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

  class AoSThreadDataGlobals {
   public:
    AoSThreadDataGlobals() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, absErrorSum{0.}, relErrorSum{0.}, __remainingTo64{} {}

    void setZero() {
      virialSum = {0., 0., 0.};
      potentialEnergySum = 0.;
      absErrorSum = 0.;
      relErrorSum = 0.;
    }

    // variables
    std::array<double, 3> virialSum;
    double potentialEnergySum;
    double absErrorSum;
    double relErrorSum;

   private:
    // dummy parameter to get the right size (64 bytes)
    double __remainingTo64[(64 - 6 * sizeof(double)) / sizeof(double)];
  };

  class AoSThreadDataFLOPs {
   public:
    AoSThreadDataFLOPs() : __remainingTo64{} {}

    void setZero() {
      numKernelCallsNoN3 = 0;
      numKernelCallsN3 = 0;
      numDistCalls = 0;
      numGlobalCalcsNoN3 = 0;
      numGlobalCalcsN3 = 0;
    }

    size_t numKernelCallsNoN3 = 0;

    size_t numKernelCallsN3 = 0;

    size_t numDistCalls = 0;

    size_t numGlobalCalcsN3 = 0;

    size_t numGlobalCalcsNoN3 = 0;

   private:
    double __remainingTo64[(64 - 5 * sizeof(size_t)) / sizeof(size_t)];
  };

  /* Interpolation Parameters */
  const float PI = 2 * std::acos(0.0);
  const double _a;
  const double _b;
  const std::vector<size_t> _numNodes;
  const std::vector<double> _intervalSplits;
  std::map<size_t, std::vector<double>> _coefficients{};
  PairwiseKernel _kernel;

  const double _cutoffSquared;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _potentialEnergySum;
  double _absErrorSum;
  double _relErrorSum;

  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadDataGlobals> _aosThreadDataGlobals{};
  std::vector<AoSThreadDataFLOPs> _aosThreadDataFLOPs{};

  std::vector<std::vector<double>> _distanceStatistics{};

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;
};
}  // namespace mdLib