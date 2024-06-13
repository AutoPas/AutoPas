/**
 * @file KryptonFunctor.h
 *
 * @date 27.04.2024
 * @author muehlhaeusser
 */

#pragma once

#include "autopas/pairwiseFunctors/Functor.h"
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
 * A functor to handle interactions between two krypton atoms.
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * @tparam Particle The type of particle.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool relevantForTuning = true>
class KryptonPairFunctor
    : public autopas::Functor<
          Particle, KryptonPairFunctor<Particle, useNewton3, calculateGlobals, relevantForTuning>> {
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
  KryptonPairFunctor() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit KryptonPairFunctor(double cutoff, void * /*dummy*/)
      : autopas::Functor<Particle,
                         KryptonPairFunctor<Particle, useNewton3, calculateGlobals, relevantForTuning>>(
            cutoff),
        _cutoffSquared{cutoff * cutoff},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false} {
    if constexpr (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
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
  explicit KryptonPairFunctor(double cutoff) : KryptonPairFunctor(cutoff, nullptr) {}

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

    const auto displacement = i.getR() - j.getR();
    const double dist2 = autopas::utils::ArrayMath::dot(displacement, displacement);

    if (dist2 > _cutoffSquared) {
      return;
    }

    // Pre calculations
    const double dist = std::sqrt(dist2);
    const double distNeg = 1. / dist;
    const double distNeg2 = 1. / dist2;
    const double distNeg6 = distNeg2 * distNeg2 * distNeg2;
    const double distNeg8 = distNeg6 * distNeg2;
    const double distNeg10 = distNeg8 * distNeg2;
    const double distNeg12 = distNeg10 * distNeg2;
    const double distNeg14 = distNeg12 * distNeg2;
    const double distNeg16 = distNeg14 * distNeg2;

    const double expAlphaTerm = std::exp(_alpha1 * dist + _alpha2 * dist2 + _alphaneg1 * distNeg);
    const double alphaTerm = _alpha1 + 2 * _alpha2 * dist - _alphaneg1 * distNeg2;

    const double firstTerm = (- _A * alphaTerm * expAlphaTerm) * distNeg;

    const double bdist = _b * dist;
    const double bdist2 = bdist * bdist;
    const double bdist3 = bdist2 * bdist;
    const double bdist4 = bdist3 * bdist;
    const double bdist5 = bdist4 * bdist;
    const double bdist6 = bdist5 * bdist;
    const double bdist7 = bdist6 * bdist;
    const double bdist8 = bdist7 * bdist;
    const double bdist9 = bdist8 * bdist;
    const double bdist10 = bdist9 * bdist;
    const double bdist11 = bdist10 * bdist;
    const double bdist12 = bdist11 * bdist;
    const double bdist13 = bdist12 * bdist;
    const double bdist14 = bdist13 * bdist;
    const double bdist15 = bdist14 * bdist;
    const double bdist16 = bdist15 * bdist;

    const double ksum0 = 1.0;
    const double ksum1 = bdist;
    const double ksum2 = bdist2 * 0.5;
    const double ksum3 = bdist3 * _fac3Inv;
    const double ksum4 = bdist4 * _fac4Inv;
    const double ksum5 = bdist5 * _fac5Inv;
    const double ksum6 = bdist6 * _fac6Inv;
    const double ksum7 = bdist7 * _fac7Inv;
    const double ksum8 = bdist8 * _fac8Inv;
    const double ksum9 = bdist9 * _fac9Inv;
    const double ksum10 = bdist10 * _fac10Inv;
    const double ksum11 = bdist11 * _fac11Inv;
    const double ksum12 = bdist12 * _fac12Inv;
    const double ksum13 = bdist13 * _fac13Inv;
    const double ksum14 = bdist14 * _fac14Inv;
    const double ksum15 = bdist15 * _fac15Inv;
    const double ksum16 = bdist16 * _fac16Inv;

    const double ksumacc6 = ksum0 + ksum1 + ksum2 + ksum3 + ksum4 + ksum5 + ksum6;
    const double ksumacc8 = ksumacc6 + ksum7 + ksum8;
    const double ksumacc10 = ksumacc8 + ksum9 + ksum10;
    const double ksumacc12 = ksumacc10 + ksum11 + ksum12;
    const double ksumacc14 = ksumacc12 + ksum13 + ksum14;
    const double ksumacc16 = ksumacc14 + ksum15 + ksum16;

    const double kksumacc6 = ksum1 + ksum2 * 2. + ksum3 * 3. + ksum4 * 4. + ksum5 * 5. + ksum6 * 6.;
    const double kksumacc8 = kksumacc6 + ksum7 * 7. + ksum8 * 8.;
    const double kksumacc10 = kksumacc8 + ksum9 * 9. + ksum10 * 10.;
    const double kksumacc12 = kksumacc10 + ksum11 * 11. + ksum12 * 12.;
    const double kksumacc14 = kksumacc12 + ksum13 * 13. + ksum14 * 14.;
    const double kksumacc16 = kksumacc14 + ksum15 * 15. + ksum16 * 16.;

    const double expbr = std::exp(-_b * dist);

    const double term6 = _C6 * distNeg6 * (-6.0 + expbr * ((6.0 + bdist) * ksumacc6 - kksumacc6));
    const double term8 = _C8 * distNeg8 * (-8.0 + expbr * ((8.0 + bdist) * ksumacc8 - kksumacc8));
    const double term10 = _C10 * distNeg10 * (-10.0 + expbr * ((10.0 + bdist) * ksumacc10 - kksumacc10));
    const double term12 = _C12 * distNeg12 * (-12.0 + expbr * ((12.0 + bdist) * ksumacc12 - kksumacc12));
    const double term14 = _C14 * distNeg14 * (-14.0 + expbr * ((14.0 + bdist) * ksumacc14 - kksumacc14));
    const double term16 = _C16 * distNeg16 * (-16.0 + expbr * ((16.0 + bdist) * ksumacc16 - kksumacc16));

    const double secondTerm = (term6 + term8 + term10 + term12 + term14 + term16) * distNeg2;
    const auto f = displacement * (firstTerm + secondTerm);

    i.addF(f);
    if (newton3) {
      // only if we use newton 3 here, we want to
      j.subF(f);
    }

    if (calculateGlobals) {
      auto virial = displacement * f;
      double potentialEnergy = _A * expAlphaTerm
                               - _C6 * distNeg6 * (1 - expbr * ksumacc6)
                               - _C8 * distNeg8 * (1 - expbr * ksumacc8)
                               - _C10 * distNeg10 * (1 - expbr * ksumacc10)
                               - _C12 * distNeg12 * (1 - expbr * ksumacc12)
                               - _C14 * distNeg14 * (1 - expbr * ksumacc14)
                               - _C16 * distNeg16 * (1 - expbr * ksumacc16);

      const int threadnum = autopas::autopas_get_thread_num();
      if (i.isOwned()) {
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadData[threadnum].virialSum += virial;
      }
      // for non-newton3 the second particle will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadData[threadnum].virialSum += virial;
      }
    }
  }

  /**
   * @copydoc autopas::Functor::SoAFunctorSingle()
   * This functor will always use a newton3 like traversal of the soa.
   * However, it still needs to know about newton3 to correctly add up the global values.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {}

  /**
   * @copydoc autopas::Functor::SoAFunctorPair()
   */
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {}

 public:
  // clang-format off
  /**
   * @copydoc autopas::Functor::SoAFunctorVerlet()
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors!
   */
  // clang-format on
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {
    if (soa.size() == 0 or neighborList.empty()) return;
  }

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
   * Get the number of flops used per kernel call for a given particle pair. This should count the
   * floating point operations needed for two particles that lie within a cutoff radius, having already calculated the
   * distance.
   * @param molAType molecule A's type id
   * @param molBType molecule B's type id
   * @param newton3 is newton3 applied.
   * @note molAType and molBType should always be the same for KryptonPairFunctor (i.e. Krypton atoms), but the method is kept to have a consistent interface with other
   * functors.
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall(size_t molAType, size_t molBType, bool newton3) {
    // Kernel:
    // Total =
    return newton3 ? 0ul : 0ul;
  }

  /**
   * Reset the global values.
   * Will set the global values to zero to prepare for the next iteration.
   */
  void initTraversal() final {
    _potentialEnergySum = 0.;
    _virialSum = {0., 0., 0.};
    _postProcessed = false;
    for (size_t i = 0; i < _aosThreadData.size(); ++i) {
      _aosThreadData[i].setZero();
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
      for (size_t i = 0; i < _aosThreadData.size(); ++i) {
        _virialSum += _aosThreadData[i].potentialEnergySum;
        _potentialEnergySum += _aosThreadData[i].potentialEnergySum;
      }

      _postProcessed = true;

      AutoPasLog(TRACE, "Final potential energy {}", _potentialEnergySum);
      AutoPasLog(TRACE, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
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

 private:
  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
   */
  class AoSThreadData {
   public:
    AoSThreadData()
        : virialSum{0., 0., 0.},
          potentialEnergySum{0.},
          __remainingTo64{} {}
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
  // make sure of the size of AoSThreadData
  static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size");

  const double _cutoffSquared;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _potentialEnergySum;

  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadData> _aosThreadData;

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

  // helper constants

  const double _A = 0.3200711798e8;
  const double _alpha1 = -0.2430565544e2;
  const double _alpha2 = -0.1435536209e2;
  const double _alphaneg1 = -0.4532273868e-1;
  const double _b = 0.2786344368e2;
  const double _C6 = 0.8992209265;
  const double _C8 = 0.7316713603e-1;
  const double _C10 = 0.7835488511e-2;
  const double _C12 = 1.1043747590e-3;
  const double _C14 = 2.0486474980e-4;
  const double _C16 = 5.0017084700e-5;

  const double _fac0Inv = 1.;
  const double _fac1Inv = 1.;
  const double _fac2Inv = 0.5;
  const double _fac3Inv = 1. / 6.;
  const double _fac4Inv = 1. / 24.;
  const double _fac5Inv = 1. / 120.;
  const double _fac6Inv = 1. / 720.;
  const double _fac7Inv = 1. / 5040.;
  const double _fac8Inv = 1. / 40320.;
  const double _fac9Inv = 1. / 362880.;
  const double _fac10Inv = 1. / 3628800.;
  const double _fac11Inv = 1. / 39916800.;
  const double _fac12Inv = 1. / 479001600.;
  const double _fac13Inv = 1. / 6227020800.;
  const double _fac14Inv = 1. / 87178291200.;
  const double _fac15Inv = 1. / 1307674368000.;
  const double _fac16Inv = 1. / 20922789888000.;
};
}  // namespace mdLib
