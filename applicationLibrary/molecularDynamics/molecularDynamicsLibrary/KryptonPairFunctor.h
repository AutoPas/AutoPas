/**
 * @file KryptonPairFunctor.h
 *
 * @date 27.04.2024
 * @author muehlhaeusser
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
/**
 * A functor to handle interactions between two krypton atoms based on the ab-initio potential given by Jäger et al.
 * (2016) https://doi.org/10.1063/1.4943959. The parameters are given in Kelvin (K) and Angström (A).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * @tparam Particle_T The type of particle.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle_T, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGlobals = false, bool relevantForTuning = true>
class KryptonPairFunctor
    : public autopas::PairwiseFunctor<Particle_T,
                                      KryptonPairFunctor<Particle_T, useNewton3, calculateGlobals, relevantForTuning>> {
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle_T::SoAArraysType;

  /**
   * Precision of SoA entries.
   */
  using SoAFloatPrecision = typename Particle_T::ParticleSoAFloatPrecision;

 public:
  /**
   * Deleted default constructor
   */
  KryptonPairFunctor() = delete;

  /**
   * Constructor.
   * @param cutoff
   */
  explicit KryptonPairFunctor(double cutoff)
      : autopas::PairwiseFunctor<Particle_T,
                                 KryptonPairFunctor<Particle_T, useNewton3, calculateGlobals, relevantForTuning>>(
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

  std::string getName() final { return "KryptonPairFunctor"; }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) final {
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
    const double distInv = 1. / dist;
    const double distInv2 = distInv * distInv;
    const double distInv6 = distInv2 * distInv2 * distInv2;
    const double distNeg8 = distInv6 * distInv2;
    const double distNeg10 = distNeg8 * distInv2;
    const double distNeg12 = distNeg10 * distInv2;
    const double distNeg14 = distNeg12 * distInv2;
    const double distNeg16 = distNeg14 * distInv2;

    const double expAlphaTerm = std::exp(_alpha1 * dist + _alpha2 * dist2 + _alphaInv1 * distInv);
    const double alphaTerm = _alpha1 + 2 * _alpha2 * dist - _alphaInv1 * distInv2;

    const double firstTerm = (-_constA * alphaTerm * expAlphaTerm) * distInv;

    const double bdist = _constb * dist;
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

    const double term6 = _constC6 * distInv6 * (-6.0 + expbr * (6.0 * ksumacc6 + bdist * ksum6));
    const double term8 = _constC8 * distNeg8 * (-8.0 + expbr * (8.0 * ksumacc8 + bdist * ksum8));
    const double term10 = _constC10 * distNeg10 * (-10.0 + expbr * (10.0 * ksumacc10 + bdist * ksum10));
    const double term12 = _constC12 * distNeg12 * (-12.0 + expbr * (12.0 * ksumacc12 + bdist * ksum12));
    const double term14 = _constC14 * distNeg14 * (-14.0 + expbr * (14.0 * ksumacc14 + bdist * ksum14));
    const double term16 = _constC16 * distNeg16 * (-16.0 + expbr * (16.0 * ksumacc16 + bdist * ksum16));

    const double secondTerm = (term6 + term8 + term10 + term12 + term14 + term16) * distInv2;
    auto f = displacement * (dist >= _minDistance
                                 ? firstTerm + secondTerm
                                 : _constATilde * distInv2 * std::exp(-_alphaTilde * dist) * (distInv + _alphaTilde));

    i.addF(f);
    if (newton3) {
      // only if we use newton 3 here, we want to
      j.subF(f);
    }

    if (calculateGlobals) {
      double potentialEnergy =
          dist >= _minDistance
              ? _constA * expAlphaTerm - _constC6 * distInv6 * (1. - expbr * ksumacc6) -
                    _constC8 * distNeg8 * (1. - expbr * ksumacc8) - _constC10 * distNeg10 * (1. - expbr * ksumacc10) -
                    _constC12 * distNeg12 * (1. - expbr * ksumacc12) -
                    _constC14 * distNeg14 * (1. - expbr * ksumacc14) - _constC16 * distNeg16 * (1. - expbr * ksumacc16)
              : (_constATilde * distInv) * std::exp(-_alphaTilde * dist);

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
   * @copydoc autopas::PairwiseFunctor::SoAFunctorSingle()
   * This functor will always use a newton3 like traversal of the soa.
   * However, it still needs to know about newton3 to correctly add up the global values.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {}

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorPair()
   */
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {}

 public:
  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorVerlet()
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors!
   */
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {
    if (soa.size() == 0 or neighborList.empty()) return;
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr()
   */
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

  /**
   * @copydoc autopas::Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle_T::AttributeNames, 6>{
        Particle_T::AttributeNames::id,     Particle_T::AttributeNames::posX,
        Particle_T::AttributeNames::posY,   Particle_T::AttributeNames::posZ,
        Particle_T::AttributeNames::typeId, Particle_T::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {
    return std::array<typename Particle_T::AttributeNames, 3>{
        Particle_T::AttributeNames::forceX, Particle_T::AttributeNames::forceY, Particle_T::AttributeNames::forceZ};
  }

  /**
   * Get the number of flops used per kernel call for a given particle pair. This should count the
   * floating point operations needed for two particles that lie within a cutoff radius, having already calculated the
   * distance.
   * @param molAType molecule A's type id
   * @param molBType molecule B's type id
   * @param newton3 is newton3 applied.
   * @note molAType and molBType should always be the same for KryptonPairFunctor (i.e. Krypton atoms), but the method
   * is kept to have a consistent interface with other functors.
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
        _virialSum += _aosThreadData[i].virialSum;
        _potentialEnergySum += _aosThreadData[i].potentialEnergySum;
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

  const double _constA = 0.3200711798e8;
  const double _alpha1 = -0.2430565544e1;
  const double _alpha2 = -0.1435536209;
  const double _alphaInv1 = -0.4532273868;
  const double _constb = 0.2786344368e1;
  const double _constC6 = 0.8992209265e6;
  const double _constC8 = 0.7316713603e7;
  const double _constC10 = 0.7835488511e8;
  const double _constC12 = 1.1043747590e9;
  const double _constC14 = 2.0486474980e10;
  const double _constC16 = 5.0017084700e11;

  const double _constATilde = 0.8268005465e7;
  const double _alphaTilde = 0.1682493666e1;
  // Use a different formula in the extreme repulsive region (r < _minDistance)
  const double _minDistance = 1.2047406;

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
