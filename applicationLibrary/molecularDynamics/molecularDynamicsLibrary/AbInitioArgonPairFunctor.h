/**
 * @file AbInitioArgonPairFunctor.h
 *
 * @date 22.08.2024
 * @author muehlhaeusser
 */
#pragma once

#include <array>

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace mdLib {

/**
 * A functor to handle argon pair interactions based on the ab-initio potential from Jäger et al. (https://doi.org/10.1080/00268976.2010.507557).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * @tparam Particle The type of particle.
 * @tparam ParticleCell The type of particle cell.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam countFLOPs counts FLOPs and hitrate. Not implemented for this functor.
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>
class AbInitioArgonPairFunctor : public autopas::PairwiseFunctor<Particle, AbInitioArgonPairFunctor<Particle, useNewton3,
                                                                    calculateGlobals, countFLOPs, relevantForTuning>> {
  using SoAArraysType = typename Particle::SoAArraysType;

 public:
  /**
   * Deleted default constructor
   */
  AbInitioArgonPairFunctor() = delete;

  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy unused, only there to make the signature different from the public constructor.
   */
  explicit AbInitioArgonPairFunctor(double cutoff)
      : autopas::PairwiseFunctor<Particle, AbInitioArgonPairFunctor<Particle, useNewton3, calculateGlobals,
                                                countFLOPs, relevantForTuning>>(cutoff),
        _cutoffSquared{cutoff * cutoff},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false} {
    if (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
    if constexpr (countFLOPs) {
      AutoPasLog(DEBUG, "Using AbInitioArgonPairFunctor with countFLOPs but FLOP counting is not implemented.");
    }
  }

  std::string getName() final { return "ArgonFunctorPairwise"; }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  inline void AoSFunctor(Particle &i, Particle &j, bool newton3) {
    using namespace autopas::utils::ArrayMath::literals;
    if (i.isDummy() or j.isDummy()) {
      return;
    }
    auto displacement = i.getR() - j.getR();
    double dist2 = autopas::utils::ArrayMath::dot(displacement, displacement);

    if (dist2 > _cutoffSquared) {
      return;
    }

    // Pre calculations
    const double dist = std::sqrt(dist2);
    const double distInv = 1. / dist;
    const double distInv2 = distInv * distInv;
    const double distInv6 = distInv2 * distInv2 * distInv2;
    const double distNeg8 = distInv6 * distInv2;
    const double distNeg10 = distNeg8 * distInv2;
    const double distNeg12 = distNeg10 * distInv2;
    const double distNeg14 = distNeg12 * distInv2;
    const double distNeg16 = distNeg14 * distInv2;

    const double expAlphaTerm = std::exp(_alpha1 * dist + _alpha2 * dist2 + _alphaneg1 * distInv + _alphaneg2 * distInv2);
    const double alphaTerm = _alpha1 + 2 * _alpha2 * dist - _alphaneg1 * distInv2 - 2 * _alphaneg2 * distInv2 * distInv;

    const double firstTerm = (-_A * alphaTerm * expAlphaTerm) * distInv;

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
    const double ksum2 = bdist2 * _invFactorials[2];
    const double ksum3 = bdist3 * _invFactorials[3];
    const double ksum4 = bdist4 * _invFactorials[4];
    const double ksum5 = bdist5 * _invFactorials[5];
    const double ksum6 = bdist6 * _invFactorials[6];
    const double ksum7 = bdist7 * _invFactorials[7];
    const double ksum8 = bdist8 * _invFactorials[8];
    const double ksum9 = bdist9 * _invFactorials[9];
    const double ksum10 = bdist10 * _invFactorials[10];
    const double ksum11 = bdist11 * _invFactorials[11];
    const double ksum12 = bdist12 * _invFactorials[12];
    const double ksum13 = bdist13 * _invFactorials[13];
    const double ksum14 = bdist14 * _invFactorials[14];
    const double ksum15 = bdist15 * _invFactorials[15];
    const double ksum16 = bdist16 * _invFactorials[16];

    const double ksumacc6 = ksum0 + ksum1 + ksum2 + ksum3 + ksum4 + ksum5 + ksum6;
    const double ksumacc8 = ksumacc6 + ksum7 + ksum8;
    const double ksumacc10 = ksumacc8 + ksum9 + ksum10;
    const double ksumacc12 = ksumacc10 + ksum11 + ksum12;
    const double ksumacc14 = ksumacc12 + ksum13 + ksum14;
    const double ksumacc16 = ksumacc14 + ksum15 + ksum16;

    const double expbr = std::exp(-bdist);

    const double term6 = _cConstants[0] * distInv6 * (-6.0 + expbr * (6.0 * ksumacc6 + bdist * ksum6));
    const double term8 = _cConstants[1] * distNeg8 * (-8.0 + expbr * (8.0 * ksumacc8 + bdist * ksum8));
    const double term10 = _cConstants[2] * distNeg10 * (-10.0 + expbr * (10.0 * ksumacc10 + bdist * ksum10));
    const double term12 = _cConstants[3] * distNeg12 * (-12.0 + expbr * (12.0 * ksumacc12 + bdist * ksum12));
    const double term14 = _cConstants[4] * distNeg14 * (-14.0 + expbr * (14.0 * ksumacc14 + bdist * ksum14));
    const double term16 = _cConstants[5] * distNeg16 * (-16.0 + expbr * (16.0 * ksumacc16 + bdist * ksum16));

    const double secondTerm = (term6 + term8 + term10 + term12 + term14 + term16) * distInv2;
    auto f = displacement * (firstTerm + secondTerm);

    i.addF(f);
    if (newton3) {
      // only if we use newton 3 here, we want to
      j.subF(f);
    }

    if (calculateGlobals) {
      double potentialEnergy =
          _A * expAlphaTerm - _cConstants[0] * distInv6 * (1. - expbr * ksumacc6) - _cConstants[1] * distNeg8 * (1. - expbr * ksumacc8) -
          _cConstants[2] * distNeg10 * (1. - expbr * ksumacc10) - _cConstants[3] * distNeg12 * (1. - expbr * ksumacc12) -
          _cConstants[4] * distNeg14 * (1. - expbr * ksumacc14) - _cConstants[5] * distNeg16 * (1. - expbr * ksumacc16);
      auto virial = displacement * f;

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
   * This functor will always do a newton3 like traversal of the soa.
   * However, it still needs to know about newton3 to correctly add up the global values.
   */
  inline void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) {}

  // clang-format off
  /**
   * @copydoc autopas::Functor::SoAFunctorPair()
   */
  // clang-format on
  inline void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                             const bool newton3) { }

  // clang-format off
  /**
   * @copydoc autopas::Functor::SoAFunctorVerlet()
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors and specify iFrom and iTo accordingly.
   */
  // clang-format on
  inline void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                               const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                               bool newton3) {}

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
        _potentialEnergySum += _aosThreadData[i].potentialEnergySum;
        _virialSum += _aosThreadData[i].virialSum;
      }
      // For each interaction, we added the full contribution for both particles. Divide by 2 here, so that each
      // contribution is only counted once per pair.
      _potentialEnergySum *= 0.5;
      _virialSum *= 0.5;

      _postProcessed = true;

      AutoPasLog(TRACE, "Final potential energy {}", _potentialEnergySum);
      AutoPasLog(TRACE, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
    }
  }

  /**
   * Get the potential Energy
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
   * Get the virial
   * @return the virial
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
    AoSThreadData() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, __remainingTo64{} {}
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

  const double _cutoffSquared = 0;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _potentialEnergySum;

  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadData> _aosThreadData;

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

  // constants
  const double _A = 4.61330146e7;
  const double _alpha1 = -2.98337630e1;
  const double _alpha2 = -9.71208881;
  const double _alphaneg1 = 2.75206827e-2;
  const double _alphaneg2 = -1.01489050e-2;
  const double _b = 4.02517211e1;

  // Constants C6, C8, C10, C12, C14, C16
  const std::array<double, 6> _cConstants = {4.42812017e-1, 3.26707684e-2, 2.45656537e-3,
                                             1.88246247e-4, 1.47012192e-5, 1.17006343e-6};

  const std::array<double, 17> _bPots = [&] {
    std::array<double, 17> pots{1.0};
    for (size_t k = 1; k < 17; k++) {
      pots[k] = pots[k - 1] * _b;
    }
    return pots;
  }();

  const std::array<double, 17> _invFactorials = {1.,
                                                 1.,
                                                 0.5,
                                                 1. / 6.,
                                                 1. / 24.,
                                                 1. / 120.,
                                                 1. / 720.,
                                                 1. / 5040.,
                                                 1. / 40320.,
                                                 1. / 362880.,
                                                 1. / 3628800.,
                                                 1. / 39916800.,
                                                 1. / 479001600.,
                                                 1. / 6227020800.,
                                                 1. / 87178291200.,
                                                 1. / 1307674368000.,
                                                 1. / 20922789888000.};

};
}  // namespace mdLib
