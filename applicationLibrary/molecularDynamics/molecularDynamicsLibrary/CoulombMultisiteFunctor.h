/**
 * @file CoulombMultisiteFunctor.h
 * @date 21/02/2022
 * @author S. Newcome
 */

#pragma once

#include "MoleculeLJ.h"
#include "MultisiteMoleculeLJ.h"
#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Math.h"
#include "autopas/utils/Quaternion.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace mdLib {

/**
 * A functor to handle Lennard-Jones interactions between two Multisite Molecules.
 *
 * @tparam Particle The type of particle.
 * @tparam applyShift Flag for the LJ potential to have a truncated shift.
 * @tparam useMixing Flag for if the functor is to be used with multiple particle types. If set to false, _epsilon and
 * _sigma need to be set and the constructor with PPL can be omitted.
 * @warning: Whilst this class allows for mixing to be disabled, this feature is not of much value in real applications
 * and as such is experimental only and untested!
 * @tparam useNewton3 Switch for the functor to support newton3 on, off, or both. See FunctorN3Nodes for possible
 * values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 * @tparam countFLOPs counts FLOPs and hitrate. Currently not implemented as this functor is problematically bad and
 * will be replaced.
 */
template <class Particle, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>
class CoulombMultisiteFunctor
    : public autopas::PairwiseFunctor<Particle,
                                      CoulombMultisiteFunctor<Particle, applyShift, useMixing, useNewton3,
                                                              calculateGlobals, relevantForTuning, countFLOPs>> {
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle::SoAArraysType;

  /**
   * Precision of SoA entries
   */
  using SoAFloatPrecision = typename Particle::ParticleSoAFloatPrecision;

  /**
   * cutoff^2
   */
  const double _cutoffSquared;

  /**
   * epsilon x 24. Not constant as may be reset through PPL.
   */
  double _epsilon24;

  /**
   * sigma^2. Not constant as may be reset through PPL.
   */
  double _sigmaSquared;

  /**
   * Not constant as may be reset through PPL.
   */
  double _shift6 = 0;

  /**
   * List of relative unrotated LJ Site Positions. This is to be used when there is no mixing of molecules.
   */
  const std::vector<std::array<double, 3>> _sitePositionsLJ{};

  /**
   * Particle property library. Not used if all sites are of the same species.
   */
  ParticlePropertiesLibrary<SoAFloatPrecision, size_t> *_PPLibrary = nullptr;

  /**
   * Sum of potential energy. Only calculated if calculateGlobals is true.
   */
  double _potentialEnergySum;

  /**
   * Sum of the virial. Only calculated if calculateGlobals is true.
   */
  std::array<double, 3> _virialSum;

  /**
   * Defines whether or whether not the global values are already processed
   */
  bool _postProcessed;

 public:
  /**
   * Delete Default constructor
   */
  CoulombMultisiteFunctor() = delete;

 private:
  /**
   * Internal (actually used) constructor
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit CoulombMultisiteFunctor(SoAFloatPrecision cutoff, void * /*dummy*/)
      : autopas::PairwiseFunctor<Particle, CoulombMultisiteFunctor<Particle, applyShift, useMixing, useNewton3,
                                                                   calculateGlobals, relevantForTuning, countFLOPs>>(
            cutoff),
        _cutoffSquared{cutoff * cutoff},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false} {
    if constexpr (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
    if constexpr (countFLOPs) {
      AutoPasLog(DEBUG, "FLOP counting is enabled but is not supported for multi-site functors yet!");
    }
  }

 public:
  /**
   * Constructor for Functor with particle mixing disabled. setParticleProperties() must be called.
   * @note Only to be used with mixing == false
   * @param cutoff
   */
  explicit CoulombMultisiteFunctor(double cutoff) : CoulombMultisiteFunctor(cutoff, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
    AutoPasLog(WARN, "Using CoulombMultisiteFunctor with mixing disabled is untested!");
  }

  /**
   * Constructor for Functor with particle mixing enabled.
   * Calculating global attributes is done with the center of mass and overall forces applied
   * @param cutoff
   * @param particlePropertiesLibrary Library used to look up the properties of each type of particle e.g. sigma,
   * epsilon, shift.
   */
  explicit CoulombMultisiteFunctor(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : CoulombMultisiteFunctor(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  std::string getName() final { return "CoulombMultisiteFunctor"; }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  /**
   * Functor for arrays of structures (AoS).
   *
   * @param particleA Particle A
   * @param particleB Particle B
   * @param newton3 Flag for if newton3 is used.
   */
  void AoSFunctor(Particle &particleA, Particle &particleB, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;
    if (particleA.isDummy() or particleB.isDummy()) {
      return;
    }

    // Don't calculate force if particleB outside cutoff of particleA
    const auto displacementCoM = autopas::utils::ArrayMath::sub(particleA.getR(), particleB.getR());
    const auto distanceSquaredCoM = autopas::utils::ArrayMath::dot(displacementCoM, displacementCoM);

    if (distanceSquaredCoM > _cutoffSquared) {
      return;
    }

    // get number of sites
    const size_t numSitesA = useMixing ? _PPLibrary->getNumSites(particleA.getTypeId()) : _sitePositionsLJ.size();
    const size_t numSitesB = useMixing ? _PPLibrary->getNumSites(particleB.getTypeId()) : _sitePositionsLJ.size();

    // get siteIds
    const std::vector<size_t> siteIdsA =
        useMixing ? _PPLibrary->getSiteTypes(particleA.getTypeId()) : std::vector<unsigned long>();
    const std::vector<size_t> siteIdsB =
        useMixing ? _PPLibrary->getSiteTypes(particleB.getTypeId()) : std::vector<unsigned long>();

    // get unrotated relative site positions
    const std::vector<std::array<double, 3>> unrotatedSitePositionsA =
        useMixing ? _PPLibrary->getSitePositions(particleA.getTypeId()) : _sitePositionsLJ;
    const std::vector<std::array<double, 3>> unrotatedSitePositionsB =
        useMixing ? _PPLibrary->getSitePositions(particleB.getTypeId()) : _sitePositionsLJ;

    // calculate correctly rotated relative site positions
    const auto rotatedSitePositionsA =
        autopas::utils::quaternion::rotateVectorOfPositions(particleA.getQuaternion(), unrotatedSitePositionsA);
    const auto rotatedSitePositionsB =
        autopas::utils::quaternion::rotateVectorOfPositions(particleB.getQuaternion(), unrotatedSitePositionsB);

    for (int i = 0; i < numSitesA; i++) {
      for (int j = 0; j < numSitesB; j++) {
        const auto displacement = autopas::utils::ArrayMath::add(
            autopas::utils::ArrayMath::sub(displacementCoM, rotatedSitePositionsB[j]), rotatedSitePositionsA[i]);
        const auto distanceSquared = autopas::utils::ArrayMath::dot(displacement, displacement);

        if (distanceSquared > _cutoffSquared) {
          continue;
        }

        const auto coulombEpsilon = useMixing ? _PPLibrary->getCoulombEpsilon(siteIdsA[i]) : 0;

        const double invDist = 1.0 / std::sqrt(distanceSquared);
        const double coulombFactor = 1.0 / (4.0 * M_PI * coulombEpsilon);
        const double chargeA = _PPLibrary->getCharge(siteIdsA[i]);
        const double chargeB = _PPLibrary->getCharge(siteIdsB[j]);

        const double forceMagnitude = coulombFactor * chargeA * chargeB * invDist * invDist;
        const auto force = autopas::utils::ArrayMath::mulScalar(displacement, forceMagnitude);

        particleA.addF(force);
        if (newton3) {
          particleB.subF(force);
        }

        particleA.addTorque(autopas::utils::ArrayMath::cross(rotatedSitePositionsA[i], force));
        if (newton3) {
          particleB.subTorque(autopas::utils::ArrayMath::cross(rotatedSitePositionsB[j], force));
        }

        if (calculateGlobals) {
          const double potentialEnergy = coulombFactor * chargeA * chargeB * invDist;
          const auto virial = displacement * force;
          const auto threadNum = autopas::autopas_get_thread_num();

          if (particleA.isOwned()) {
            _aosThreadData[threadNum].potentialEnergySum += potentialEnergy;
            _aosThreadData[threadNum].virialSum += virial;
          }
          if (newton3 and particleB.isOwned()) {
            _aosThreadData[threadNum].potentialEnergySum += potentialEnergy;
            _aosThreadData[threadNum].virialSum += virial;
          }
        }
      }
    }
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorSingle()
   * This functor will always use a newton3 like traversing of the soa, however, it still needs to know about newton3
   * to use it correctly for the global values.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {}
  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorPair()
   */
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {}

  // clang-format off
  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorVerlet()
   */
  // clang-format on
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {}

  /**
   * Sets the molecule properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   * @param epsilon24 epsilon * 24
   * @param sigmaSquared sigma^2
   * @param sitePositionsLJ vector of 3D relative unrotated untranslated site positions
   */
  void setParticleProperties(SoAFloatPrecision epsilon24, SoAFloatPrecision sigmaSquared,
                             std::vector<std::array<SoAFloatPrecision, 3>> sitePositionsLJ) {}

  /**
   * @copydoc autopas::Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 16>{
        Particle::AttributeNames::id,          Particle::AttributeNames::posX,
        Particle::AttributeNames::posY,        Particle::AttributeNames::posZ,
        Particle::AttributeNames::forceX,      Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ,      Particle::AttributeNames::quaternion0,
        Particle::AttributeNames::quaternion1, Particle::AttributeNames::quaternion2,
        Particle::AttributeNames::quaternion3, Particle::AttributeNames::torqueX,
        Particle::AttributeNames::torqueY,     Particle::AttributeNames::torqueZ,
        Particle::AttributeNames::typeId,      Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 16>{
        Particle::AttributeNames::id,          Particle::AttributeNames::posX,
        Particle::AttributeNames::posY,        Particle::AttributeNames::posZ,
        Particle::AttributeNames::forceX,      Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ,      Particle::AttributeNames::quaternion0,
        Particle::AttributeNames::quaternion1, Particle::AttributeNames::quaternion2,
        Particle::AttributeNames::quaternion3, Particle::AttributeNames::torqueX,
        Particle::AttributeNames::torqueY,     Particle::AttributeNames::torqueZ,
        Particle::AttributeNames::typeId,      Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 6>{
        Particle::AttributeNames::forceX,  Particle::AttributeNames::forceY,  Particle::AttributeNames::forceZ,
        Particle::AttributeNames::torqueX, Particle::AttributeNames::torqueY, Particle::AttributeNames::torqueZ};
  }

  /**
   * @return useMixing
   */
  constexpr static bool getMixing() { return useMixing; }

  /**
   * Get the number of flops used per kernel call - i.e. number of flops to calculate kernel *given* the two particles
   * lie within the cutoff (i.e. distance^2 / cutoff has been already been calculated).
   * Note: there is currently a large difference between AoS & SoA number of flops. This function returns the AoS
   * number of flops.
   * @param molAType molecule A's type id
   * @param molBType molecule B's type id
   * @param newton3 true if newton3 optimizations enabled
   * @return Number of FLOPs
   */
  unsigned long getNumFlopsPerKernelCall(size_t molAType, size_t molBType, bool newton3) {
    // Site-to-site displacement: 6 (3 in the SoA case, but this requires O(N) precomputing site positions)
    // Site-to-site distance squared: 4
    // Compute scale: 9
    // Apply scale to force: With newton3: 6, Without: 3
    // Apply scale to torque: With newton3 18, Without: 9 (0 in SoA case, with O(N) post computing)
    // Site-to-site total: With newton3: 33, Without: 26
    // (SoA total: With N3L: 22, Without N3L: 19)
    // Above multiplied by number sites of i * number sites of j
    const unsigned long siteToSiteFlops = newton3 ? 33ul : 26ul;
    return _PPLibrary->getNumSites(molAType) * _PPLibrary->getNumSites(molBType) * siteToSiteFlops;
  }

  /**
   * Reset the global values.
   * Will set the global values to zero to prepare for the next iteration.
   */
  void initTraversal() final {
    _potentialEnergySum = 0;
    _virialSum = {0., 0., 0.};
    _postProcessed = false;
    for (size_t i = 0; i < _aosThreadData.size(); i++) {
      _aosThreadData[i].setZero();
    }
  }

  /**
   * Postprocesses global values, e.g. potential energy & virial
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

      // We have always calculated 6*potentialEnergy, so we divide by 6 here!
      _potentialEnergySum /= 6.;
      _postProcessed = true;

      AutoPasLog(DEBUG, "Final potential energy {}", _potentialEnergySum);
      AutoPasLog(DEBUG, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
    }
  }

  /**
   * Get the potential energy.
   *
   * @return the potential energy
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
   * Implementation function of SoAFunctorPair(soa1, soa2, newton3)
   * @tparam newton3 flag for if newton's third law is used
   * @param soaA structure of arrays A
   * @param soaB structure of arrays B
   */
  template <bool newton3>
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soaA, autopas::SoAView<SoAArraysType> soaB) {}

  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexPrime,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {}

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

  static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size (should be multiple of 64)");

  /**
   * Thread buffer for AoS
   */
  std::vector<AoSThreadData> _aosThreadData;
};
}  // namespace mdLib