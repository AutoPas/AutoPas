/**
 * @file LuTFunctor.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "LuT.h"
#include "ParticlePropertiesLibrary.h"
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

/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * @tparam Particle The type of particle.
 * @tparam applyShift Switch for the lj potential to be truncated shifted.
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam countFLOPs counts FLOPs and hitrate
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>
class LuTFunctor
    : public autopas::PairwiseFunctor<Particle, LuTFunctor<Particle, applyShift, useMixing, useNewton3,
                                                           calculateGlobals, countFLOPs, relevantForTuning>> {
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle::SoAArraysType;

  /**
   * Precision of SoA entries.
   */
  using SoAFloatPrecision = typename Particle::ParticleSoAFloatPrecision;

 public:
  /**
   * Deleted default constructor
   */
  LuTFunctor() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit LuTFunctor(double cutoff, void * /*dummy*/)
      : autopas::PairwiseFunctor<Particle, LuTFunctor<Particle, applyShift, useMixing, useNewton3, calculateGlobals,
                                                      countFLOPs, relevantForTuning>>(cutoff),
        _cutoffSquared{cutoff * cutoff},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _postProcessed{false} {
    if constexpr (calculateGlobals) {
      _aosThreadDataGlobals.resize(autopas::autopas_get_max_threads());
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs.resize(autopas::autopas_get_max_threads());
    }
  }

 public:
  /**
   * Constructor for Functor with mixing disabled. When using this functor it is necessary to call
   * setParticleProperties() to set internal constants because it does not use a particle properties library.
   *
   * @note Only to be used with mixing == false.
   *
   * @param cutoff
   */
  explicit LuTFunctor(double cutoff) : LuTFunctor(cutoff, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  /**
   * Constructor for Functor with mixing active. This functor takes a ParticlePropertiesLibrary to look up (mixed)
   * properties like sigma, epsilon and shift.
   * @param cutoff
   * @param particlePropertiesLibrary
   */
  explicit LuTFunctor(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LuTFunctor(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  void setLuT(std::shared_ptr<LookupTable> lookupTable) { _lookupTable = lookupTable; }

  std::string getName() final { return "LuTFunctorAutoVec"; }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  void AoSFunctor(Particle &i, Particle &j, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;

    if (i.isDummy() or j.isDummy()) {
      return;
    }

    auto dr = j.getR() - i.getR();
    double dr2 = autopas::utils::ArrayMath::dot(dr, dr);

    if (dr2 > _cutoffSquared) {
      return;
    }

    double r = std::sqrt(dr2);
    double potential = _lookupTable->interpolate(r);

    // Calculation of the force: F = -dU/dr, here approximation with difference quotient
    double h = 1e-5;
    double potentialPlus = _lookupTable->interpolate(r + h);
    double potentialMinus = _lookupTable->interpolate(r - h);
    double forceMagnitude = -(potentialPlus - potentialMinus) / (2 * h);

    auto forceVec = (forceMagnitude / r) * dr;
    i.addF(forceVec);
    if (newton3) {
      j.subF(forceVec);
    }
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorSingle()
   * This functor will always use a newton3 like traversal of the soa.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {}

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorPair()
   */
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {}

 private:
  /**
   * Implementation function of SoAFunctorPair(soa1, soa2, newton3)
   *
   * @tparam newton3
   * @param soa1
   * @param soa2
   */
  template <bool newton3>
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {}

 public:
  // clang-format off
  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorVerlet()
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors!
   */
  // clang-format on
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {}

  /**
   * Sets the particle properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   *
   * @param epsilon24
   * @param sigmaSquared
   */
  void setParticleProperties(SoAFloatPrecision epsilon24, SoAFloatPrecision sigmaSquared) {}

  /**
   * @copydoc autopas::Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 9>{
        Particle::AttributeNames::id,     Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,   Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ, Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 6>{
        Particle::AttributeNames::id,   Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ, Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 3>{
        Particle::AttributeNames::forceX, Particle::AttributeNames::forceY, Particle::AttributeNames::forceZ};
  }

  /**
   *
   * @return useMixing
   */
  constexpr static bool getMixing() { return useMixing; }

  /**
   * Reset the global values.
   * Will set the global values to zero to prepare for the next iteration.
   */
  void initTraversal() final {
    _potentialEnergySum = 0.;
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
  }

  /**
   * Accumulates global values, e.g. potential energy and virial.
   * @param newton3
   */
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
      }
      // For each interaction, we added the full contribution for both particles. Divide by 2 here, so that each
      // contribution is only counted once per pair.
      _potentialEnergySum *= 0.5;
      _virialSum *= 0.5;

      // We have always calculated 6*potentialEnergy, so we divide by 6 here!
      _potentialEnergySum /= 6.;
      _postProcessed = true;

      AutoPasLog(DEBUG, "Final potential energy {}", _potentialEnergySum);
      AutoPasLog(DEBUG, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
    }
  }

  /**
   * Get the potential Energy.
   * @return the potential Energy
   */
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

  /**
   * Get the virial.
   * @return
   */
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

  /**
   * Gets the number of useful FLOPs.
   *
   * For the distance squared calculation, this is:
   * - Displacement: 3
   * - DistanceSquared in each dimension: 3
   * - DistanceSquared Total: 2
   * - Total: 8
   *
   * For the force kernel, this is:
   * - inverse distance squared: 1 (assume division is 1 FLOP)
   * - lj2: 1
   * - lj6: 2
   * - lj12: 1
   * - lj12m6: 1
   * - scalar factor: 3
   * - force: 3
   * - accumulate force on mol i: 3
   * - accumulate force on mol j if n3: 3
   * - Total: 15 without n3, 18 with n3
   *
   * For the globals calculation, this is:
   * - virial: 3
   * - potential: 1, or 2 with shift
   * - accumulation: 4 without n3, 8 with n3
   * - Total: 8 or 9 without n3, 12 or 13 with n3
   *
   * Caveats:
   *
   * This function is supposed to return useful FLOPs, e.g. without counting masked vector instructions.
   * You could also argue that, strictly speaking, we redundantly calculate forces and globals twice in the newton3 case
   * on a owned/halo boundary. This function does not treat such "redundant" calculations as useless. Similarly, this
   * function does not treat halo-halo interactions as redundant useless calculations.
   *
   * @return number of FLOPs since initTraversal() is called.
   */
  [[nodiscard]] size_t getNumFLOPs() const override {
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
      const size_t numGlobalCalcsN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numGlobalCalcsN3; });
      const size_t numGlobalCalcsNoN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numGlobalCalcsNoN3; });

      constexpr size_t numFLOPsPerDistanceCall = 8;
      constexpr size_t numFLOPsPerN3KernelCall = 18;
      constexpr size_t numFLOPsPerNoN3KernelCall = 15;
      constexpr size_t numFLOPsPerN3GlobalCalc = applyShift ? 13 : 12;
      constexpr size_t numFLOPsPerNoN3GlobalCalc = applyShift ? 9 : 8;

      return numDistCallsAcc * numFLOPsPerDistanceCall + numKernelCallsN3Acc * numFLOPsPerN3KernelCall +
             numKernelCallsNoN3Acc * numFLOPsPerNoN3KernelCall + numGlobalCalcsN3Acc * numFLOPsPerN3GlobalCalc +
             numGlobalCalcsNoN3Acc * numFLOPsPerNoN3GlobalCalc;
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

 private:
  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {}

  /**
   * This class stores internal data for global calculations for each thread. Make sure that this data has proper size,
   * i.e. k*64 Bytes!
   */
  class AoSThreadDataGlobals {
   public:
    AoSThreadDataGlobals() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, __remainingTo64{} {}
    void setZero() {
      virialSum = {0., 0., 0.};
      potentialEnergySum = 0.;
    }

    // variables
    std::array<double, 3> virialSum;
    double potentialEnergySum;

   private:
    // dummy parameter to get the right size (64 bytes)
    double __remainingTo64[(64 - 4 * sizeof(double)) / sizeof(double)];
  };

  /**
   * This class stores internal data for FLOP counters for each thread. Make sure that this data has proper size, i.e.
   * k*64 Bytes!
   * The FLOP count and HitRate are not counted/calculated directly, but through helper counters (numKernelCallsNoN3,
   * numKernelCallsN3, numDistCalls, numGlobalCalcs) to reduce computational cost in the functors themselves and to
   * improve maintainability (e.g. if the cost of a kernel call changes).
   */
  class AoSThreadDataFLOPs {
   public:
    AoSThreadDataFLOPs() : __remainingTo64{} {}

    /**
     * Set all counters to zero.
     */
    void setZero() {
      numKernelCallsNoN3 = 0;
      numKernelCallsN3 = 0;
      numDistCalls = 0;
      numGlobalCalcsNoN3 = 0;
      numGlobalCalcsN3 = 0;
    }

    /**
     * Number of calls to Lennard-Jones Kernel with newton3 disabled.
     * Used for calculating number of FLOPs and hit rate.
     */
    size_t numKernelCallsNoN3 = 0;

    /**
     * Number of calls to Lennard-Jones Kernel with newton3 enabled.
     * Used for calculating number of FLOPs and hit rate.
     */
    size_t numKernelCallsN3 = 0;

    /**
     * Number of distance calculations.
     * Used for calculating number of FLOPs and hit rate.
     */
    size_t numDistCalls = 0;

    /**
     * Counter for the number of times the globals have been calculated. Excludes the special case that N3 is enabled
     * and we calculate globals for an owned-halo particle pair.
     */
    size_t numGlobalCalcsN3 = 0;

    /**
     * Counter for the number of times the globals have been calculated. Excludes the special case that N3 is enabled
     * and we calculate globals for an owned-halo particle pair.
     */
    size_t numGlobalCalcsNoN3 = 0;

   private:
    /**
     * dummy parameter to get the right size (64 bytes)
     */
    double __remainingTo64[(64 - 5 * sizeof(size_t)) / sizeof(size_t)];
  };

  // make sure of the size of AoSThreadDataGlobals and AoSThreadDataFLOPs
  static_assert(sizeof(AoSThreadDataGlobals) % 64 == 0, "AoSThreadDataGlobals has wrong size");
  static_assert(sizeof(AoSThreadDataFLOPs) % 64 == 0, "AoSThreadDataFLOPs has wrong size");

  const double _cutoffSquared;
  // not const because they might be reset through PPL
  double _epsilon24, _sigmaSquared, _shift6 = 0;

  ParticlePropertiesLibrary<SoAFloatPrecision, size_t> *_PPLibrary = nullptr;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _potentialEnergySum;

  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadDataGlobals> _aosThreadDataGlobals{};
  std::vector<AoSThreadDataFLOPs> _aosThreadDataFLOPs{};

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

  std::shared_ptr<LookupTable> _lookupTable;
};
}  // namespace mdLib
