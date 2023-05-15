/**
 * @file LJMultisiteFunctorCTS.h
 * @date 06/05/2022
 * @author Q. Behrami
 */
#pragma once

#ifndef __AVX__
#pragma message "LJMultisiteFunctorCTS.h included, but AVX is not supported by the compiler."
#else
#include "immintrin.h"
#endif

#include <array>

#include "../pairwiseFunctors/Functor.h"
#include "../utils/ExceptionHandler.h"
#include "MultisiteMoleculeLJ.h"
#include "ParticlePropertiesLibrary.h"
#include "autopas/utils/Quaternion.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * A functor to handle Lennard-Jones interactions between two Multisite Molecules.
 * This functor utilizes AVX instructions to speed up the computation.
 *
 * @tparam Particle The type of particle.
 * @tparam applyShift Flag for the LJ potential to have a truncated shift.
 * @tparam useMixing Flag for if the functor is to be used with multiple particle types. If set to false, _epsilon and
 * _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off, or both. See FunctorN3Nodes for possible
 * values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool relevantForTuning = true>
class LJMultisiteFunctorCTS
    : public autopas::Functor<Particle, LJMultisiteFunctorCTS<Particle, applyShift, useMixing, useNewton3,
                                                              calculateGlobals, relevantForTuning>> {
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle::SoAArraysType;

  /**
   * Precision of SoA entries
   */
  using SoAFloatPrecision = typename Particle::ParticleSoAFloatPrecision;

  /**
   * cutoff^2 for AoS Functor
   */
  const double _cutoffSquaredAoS = 0;

  /**
   * epsilon x 24 for AoS Functor. Not constant as may be reset through PPL.
   */
  double _epsilon24AoS = 0;

  /**
   * sigma^2 for AoS Functor. Not constant as may be reset through PPL.
   */
  double _sigmaSquaredAoS = 0;

  /**
   * Shift for AoS Functor. Not constant as may be reset through PPL.
   */
  double _shift6AoS = 0;

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

  /**
   * List of relative unrotated LJ Site Positions. This is to be used when there is no mixing of molecules.
   */
  const std::vector<std::array<double, 3>> _sitePositionsLJ{};

  // number of double values that fit into a vector register.
  // MUST be power of 2 because some optimizations make this assumption
  // @Note: The functor will probably not work when vecLength is not 4.
  constexpr static size_t vecLength = 4;

#ifdef __AVX__
  const __m256d _cutoffSquared{};
  const __m256d _zero{_mm256_set1_pd(0.)};
  const __m256d _one{_mm256_set1_pd(1.)};
  const __m256i _masks[3]{
      _mm256_set_epi64x(0, 0, 0, -1),
      _mm256_set_epi64x(0, 0, -1, -1),
      _mm256_set_epi64x(0, -1, -1, -1),
  };
  const __m256i _ownedStateOwnedMM256i{_mm256_set1_epi64x(static_cast<int64_t>(OwnershipState::owned))};

#endif

 public:
  /**
   * Deleted default constructor.
   */
  LJMultisiteFunctorCTS() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy unused, only there to make the signature different from the public constructor.
   */
  explicit LJMultisiteFunctorCTS(double cutoff, void * /*dummy*/)
#ifdef __AVX__
      : Functor<Particle, LJMultisiteFunctorCTS<Particle, applyShift, useMixing, useNewton3, calculateGlobals,
                                                relevantForTuning>>(cutoff),
        _cutoffSquared{_mm256_set1_pd(cutoff * cutoff)},
        _cutoffSquaredAoS(cutoff * cutoff),
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false} {
    if (calculateGlobals) {
      _aosThreadData.resize(autopas_get_max_threads());
    }
  }
#else
      : Functor<Particle, LJMultisiteFunctorCTS<Particle, applyShift, useMixing, useNewton3, calculateGlobals,
                                                relevantForTuning>>(cutoff) {
    utils::ExceptionHandler::exception("AutoPas was compiled without AVX support!");
  }
#endif

 public:
  /**
   * Constructor for Functor with particle mixing disabled. setParticleProperties() must be called.
   * @note Only to be used with mixing == false
   * @param cutoff
   */
  explicit LJMultisiteFunctorCTS(double cutoff) : LJMultisiteFunctorCTS(cutoff, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  /**
   * Constructor for Functor with particle mixing enabled.
   * Calculating global attributes is done with CoM and overall forces applied
   * @param cutoff
   * @param particlePropertiesLibrary Library used to look up the properties of each type of particle e.g. sigma,
   * epsilon, shift.
   */
  explicit LJMultisiteFunctorCTS(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJMultisiteFunctorCTS(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

 public:
  /**
   * Sets the particle properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   * @param epsilon24 epsilon * 24
   * @param sigmaSquared sigma^2
   */
  void setParticleProperties(SoAFloatPrecision epsilon24, SoAFloatPrecision sigmaSquared) {
    _epsilon24AoS = epsilon24;
    _sigmaSquaredAoS = sigmaSquared;
    if (applyShift) {
      _shift6AoS =
          ParticlePropertiesLibrary<double, size_t>::calcShift6(_epsilon24AoS, _sigmaSquaredAoS, _cutoffSquaredAoS);
    } else {
      _shift6AoS = 0;
    }
  }

  /**
   * @copydoc Functor::getNeededAttr()
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
   * @copydoc Functor::getNeededAttr(std::false_type)
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
   * @copydoc Functor::getComputedAttr()
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
   * @param numB number of sites in molecule B
   * @return #FLOPs
   * @TODO this will probably need a change for vectorized functors
   */
  unsigned long getNumFlopsPerKernelCall(size_t molAType, size_t molBType, bool newton3) {
    // Site-to-site displacement: 6 (3 in the SoA case, but this requires O(N) precomputing site positions)
    // Site-to-site distance squared: 4
    // Compute scale: 9
    // Apply scale to force: With newton3: 6, Without: 3
    // Apply scale to torque: With newton3 18, Without: 9 (0 in SoA case, with O(N) post computing)
    // Site-to-site total: With newton3: 33, Without: 26
    // (SoA total: With N3L: 19)
    // Above multiplied by number sites of i * number sites of j
    const unsigned long siteToSiteFlops = newton3 ? 33ul : 26ul;
    return _PPLibrary->getNumSites(molAType) * _PPLibrary->getNumSites(molBType) * siteToSiteFlops;
  }

  /**
   * Reset the global values.
   * Will set the global values to zero to prepare for the next iteration.
   * @TODO maybe needs refactoring for vectorization
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
   * @TODO maybe needs refactoring for vectorization
   * @param newton3
   */
  void endTraversal(bool newton3) final {
    if (_postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
    }
    if (calculateGlobals) {
      for (size_t i = 0; i < _aosThreadData.size(); ++i) {
        _potentialEnergySum += _aosThreadData[i].potentialEnergySum;
        _virialSum = autopas::utils::ArrayMath::add(_virialSum, _aosThreadData[i].virialSum);
      }
      if (not newton3) {
        // if the newton3 optimization is disabled we have added every energy contribution twice, so we divide by 2
        // here.
        _potentialEnergySum *= 0.5;
        _virialSum = autopas::utils::ArrayMath::mulScalar(_virialSum, 0.5);
      }
      // we have always calculated 6*potential, so we divide by 6 here!
      _potentialEnergySum /= 6.;
      _postProcessed = true;
    }
  }

  /**
   * Get the potential energy.
   * @return the potential energy
   */
  double getPotentialEnergy() {
    if (not calculateGlobals) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Trying to get potential energy even though calculateGlobals is false. If you want this functor to "
          "calculate "
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

  /**
   * Functors
   */

  /**
   * @copydoc Functor::AoSFunctor(Particle&, Particle&, bool)
   */
  void AoSFunctor(Particle &particleA, Particle &particleB, bool newton3) final {
    if (particleA.isDummy() or particleB.isDummy()) {
      return;
    }

    // Don't calculate force if particleB outside cutoff of particleA
    const auto displacementCoM = autopas::utils::ArrayMath::sub(particleA.getR(), particleB.getR());
    const auto distanceSquaredCoM = autopas::utils::ArrayMath::dot(displacementCoM, displacementCoM);

    if (distanceSquaredCoM > _cutoffSquaredAoS) {
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

    // calculate correctly rotated relative site positions (rotated correctly)
    const auto rotatedSitePositionsA =
        autopas::utils::quaternion::rotateVectorOfPositions(particleA.getQ(), unrotatedSitePositionsA);
    const auto rotatedSitePositionsB =
        autopas::utils::quaternion::rotateVectorOfPositions(particleB.getQ(), unrotatedSitePositionsB);

    for (int i = 0; i < numSitesA; i++) {
      for (int j = 0; j < numSitesB; j++) {
        const auto displacement = autopas::utils::ArrayMath::add(
            autopas::utils::ArrayMath::sub(displacementCoM, rotatedSitePositionsB[j]), rotatedSitePositionsA[i]);
        const auto distanceSquared = autopas::utils::ArrayMath::dot(displacement, displacement);

        const auto sigmaSquared =
            useMixing ? _PPLibrary->getMixingSigmaSquared(siteIdsA[i], siteIdsB[j]) : _sigmaSquaredAoS;
        const auto epsilon24 = useMixing ? _PPLibrary->getMixing24Epsilon(siteIdsA[i], siteIdsB[j]) : _epsilon24AoS;
        const auto shift6 = applyShift ? _PPLibrary->getMixingShift6(siteIdsA[i], siteIdsB[j]) : _shift6AoS;

        // clang-format off
        // Calculate potential between sites and thus force
        // Force = 24 * epsilon * (2*(sigma/distance)^12 - (sigma/distance)^6) * (1/distance)^2 * [x_displacement, y_displacement, z_displacement]
        //         {                         scalarMultiple                                   } * {                     displacement             }
        // clang-format on
        const auto invDistSquared = 1. / distanceSquared;
        const auto lj2 = sigmaSquared * invDistSquared;
        const auto lj6 = lj2 * lj2 * lj2;
        const auto lj12 = lj6 * lj6;
        const auto lj12m6 = lj12 - lj6;  // = LJ potential / (4x epsilon)
        const auto scalarMultiple = epsilon24 * (lj12 + lj12m6) * invDistSquared;
        const auto force = autopas::utils::ArrayMath::mulScalar(displacement, scalarMultiple);

        // Add force on site to net force
        particleA.addF(force);
        if (newton3) {
          particleB.subF(force);
        }

        // Add torque applied by force
        particleA.addTorque(autopas::utils::ArrayMath::cross(rotatedSitePositionsA[i], force));
        if (newton3) {
          particleB.subTorque(autopas::utils::ArrayMath::cross(rotatedSitePositionsB[j], force));
        }

        if (calculateGlobals) {
          // Here we calculate either the potential energy * 6 or the potential energy * 12.
          // For newton3, this potential energy contribution is distributed evenly to the two molecules.
          // For non-newton3, the full potential energy is added to the one molecule.
          // I.e. double the potential energy will be added in this case.
          // The division by 6 is handled in endTraversal, as well as the division by two needed if newton3 is not used.
          // There is a similar handling of the virial, but without the mutliplication/division by 6.
          const auto potentialEnergy6 = newton3 ? 0.5 * (epsilon24 * lj12m6 + shift6) : (epsilon24 * lj12m6 + shift6);
          const auto virial = newton3 ? utils::ArrayMath::mulScalar(utils::ArrayMath::mul(displacement, force), 0.5)
                                      : utils::ArrayMath::mul(displacement, force);

          const auto threadNum = autopas_get_thread_num();

          if (particleA.isOwned()) {
            _aosThreadData[threadNum].potentialEnergySum += potentialEnergy6;
            _aosThreadData[threadNum].virialSum = utils::ArrayMath::add(_aosThreadData[threadNum].virialSum, virial);
          }
          if (newton3 and particleB.isOwned()) {
            _aosThreadData[threadNum].potentialEnergySum += potentialEnergy6;
            _aosThreadData[threadNum].virialSum = utils::ArrayMath::add(_aosThreadData[threadNum].virialSum, virial);
          }
        }
      }
    }
  }

  /**
   * @copydoc Functor::SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3)
   * @note This functor ignores the newton3 value, as we do not expect any benefit from disabling newton3.
   */
  void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) final {
    if (newton3) {
      SoAFunctorSingleImpl<true>(soa);
    } else {
      SoAFunctorSingleImpl<false>(soa);
    }
  }

  /**
   * @copydoc Functor::SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3)
   */
  void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, const bool newton3) final {
    if (newton3) {
      SoAFunctorPairImpl<true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1, soa2);
    }
  }

  /**
   * @copydoc Functor::SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst, const std::vector<size_t,
   * autopas::AlignedAllocator<size_t>> &neighborList, bool newton3)
   */
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {
    if (soa.getNumberOfParticles() == 0 or neighborList.empty()) return;
    if (newton3) {
      SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
    }
  }

  /**
   * @TODO Documentation
   */
 private:
  template <bool newton3>
  void SoAFunctorSingleImpl(SoAView<SoAArraysType> soa) {
#ifndef __AVX__
#pragma message "SoAFunctorCTS called without AVX support!"
#endif
    if (soa.getNumberOfParticles() == 0) return;

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    const auto *const __restrict q0ptr = soa.template begin<Particle::AttributeNames::quaternion0>();
    const auto *const __restrict q1ptr = soa.template begin<Particle::AttributeNames::quaternion1>();
    const auto *const __restrict q2ptr = soa.template begin<Particle::AttributeNames::quaternion2>();
    const auto *const __restrict q3ptr = soa.template begin<Particle::AttributeNames::quaternion3>();

    SoAFloatPrecision *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    SoAFloatPrecision *const __restrict txptr = soa.template begin<Particle::AttributeNames::torqueX>();
    SoAFloatPrecision *const __restrict typtr = soa.template begin<Particle::AttributeNames::torqueY>();
    SoAFloatPrecision *const __restrict tzptr = soa.template begin<Particle::AttributeNames::torqueZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();

    __m256d potentialEnergySum = _zero;
    __m256d virialSumX = _zero;
    __m256d virialSumY = _zero;
    __m256d virialSumZ = _zero;

    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const SoAFloatPrecision cutoffSquared = _cutoffSquaredAoS;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionZ;

    // we require arrays for forces for sites to maintain SIMD in site-site calculations
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceZ;

    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteTypes;
    std::vector<size_t, autopas::AlignedAllocator<size_t>> isSiteOwned;

    const SoAFloatPrecision const_sigmaSquared = _sigmaSquaredAoS;
    const SoAFloatPrecision const_epsilon24 = _epsilon24AoS;
    const SoAFloatPrecision const_shift6 = _shift6AoS;

    const auto const_unrotatedSitePositions = _sitePositionsLJ;

    // count number of sites in SoA
    size_t siteCount = 0;
    if constexpr (useMixing) {
      for (size_t mol = 0; mol < soa.getNumberOfParticles(); ++mol) {
        siteCount += _PPLibrary->getNumSites(typeptr[mol]);
      }
    } else {
      siteCount = const_unrotatedSitePositions.size() * soa.getNumberOfParticles();
    }

    // pre-reserve site std::vectors
    exactSitePositionX.reserve(siteCount);
    exactSitePositionY.reserve(siteCount);
    exactSitePositionZ.reserve(siteCount);

    if constexpr (useMixing) {
      siteTypes.reserve(siteCount);
    }

    siteForceX.reserve((siteCount));
    siteForceY.reserve((siteCount));
    siteForceZ.reserve((siteCount));

    if constexpr (calculateGlobals) {
      // this is only needed for vectorization when calculating globals
      isSiteOwned.reserve(siteCount);
    }

    // Fill site-wise std::vectors for SIMD
    std::vector<std::array<double, 3>> rotatedSitePositions;
    std::fill_n(siteForceX.begin(), siteCount, 0.);
    std::fill_n(siteForceY.begin(), siteCount, 0.);
    std::fill_n(siteForceZ.begin(), siteCount, 0.);
    for (size_t mol = 0; mol < soa.getNumberOfParticles(); ++mol) {
      if constexpr (useMixing) {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, _PPLibrary->getSitePositions(typeptr[mol]));
      } else {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, const_unrotatedSitePositions);
      }
      size_t localSiteCount = useMixing ? _PPLibrary->getNumSites(typeptr[mol]) : const_unrotatedSitePositions.size();

      for (size_t site = 0; site < localSiteCount; ++site) {
        exactSitePositionX.push_back(rotatedSitePositions[site][0] + xptr[mol]);
        exactSitePositionY.push_back(rotatedSitePositions[site][1] + yptr[mol]);
        exactSitePositionZ.push_back(rotatedSitePositions[site][2] + zptr[mol]);
        if constexpr (calculateGlobals) {
          isSiteOwned.push_back(ownedStatePtr[mol] == OwnershipState::owned);
        }
        if constexpr (useMixing) {
          siteTypes.push_back(_PPLibrary->getSiteTypes(typeptr[mol])[site]);
        }
      }
    }

    // main force calculation loop
    size_t siteIndexMolA = 0;  // index of first site in molA
    for (size_t molA = 0; molA < soa.getNumberOfParticles(); ++molA) {
      const auto ownedStateA = ownedStatePtr[molA];
      if (ownedStateA == autopas::OwnershipState::dummy) {
        continue;
      }

      const size_t noSitesInMolA = useMixing ? _PPLibrary->getNumSites(typeptr[molA])
                                             : const_unrotatedSitePositions.size();  // Number of sites in molecule A
      const size_t siteIndexMolB = siteIndexMolA + noSitesInMolA;                    // index of first site in molB
      const size_t noSitesB = (siteCount - siteIndexMolB);  // Number of sites in molecules that A interacts with

      // create mask over every mol 'above' molA  (doubles because they are easy to vectorize)
      std::vector<double, autopas::AlignedAllocator<double>> molMask;
      molMask.reserve(soa.getNumberOfParticles() - (molA + 1));

      const __m256d xposA = _mm256_broadcast_sd(&xptr[molA]);
      const __m256d yposA = _mm256_broadcast_sd(&yptr[molA]);
      const __m256d zposA = _mm256_broadcast_sd(&zptr[molA]);

      // Build the molMask
      for (size_t molB = molA + 1; molB < soa.getNumberOfParticles(); molB += vecLength) {
        const size_t rest = soa.getNumberOfParticles() - molB;
        const bool remainderCase = rest < vecLength;
        __m256i remainderMask = remainderCase ? _masks[rest - 1] : _mm256_set1_epi64x(-1);

        const __m256d xposB = load_pd(remainderCase, &xptr[molB], remainderMask);
        const __m256d yposB = load_pd(remainderCase, &yptr[molB], remainderMask);
        const __m256d zposB = load_pd(remainderCase, &zptr[molB], remainderMask);

        const __m256d displacementCoMX = _mm256_sub_pd(xposA, xposB);
        const __m256d displacementCoMY = _mm256_sub_pd(yposA, yposB);
        const __m256d displacementCoMZ = _mm256_sub_pd(zposA, zposB);

        const __m256d distanceSquaredCoMX = _mm256_mul_pd(displacementCoMX, displacementCoMX);
        const __m256d distanceSquaredCoMY = _mm256_mul_pd(displacementCoMY, displacementCoMY);
        const __m256d distanceSquaredCoMZ = _mm256_mul_pd(displacementCoMZ, displacementCoMZ);

        const __m256d distanceSquaredCoM =
            _mm256_add_pd(_mm256_add_pd(distanceSquaredCoMX, distanceSquaredCoMY), distanceSquaredCoMZ);

        const __m256d cutoffMask = _mm256_cmp_pd(distanceSquaredCoM, _cutoffSquared, _CMP_LE_OS);
        const __m256i ownedStateB =
            load_epi64(remainderCase, reinterpret_cast<const size_t *>(&ownedStatePtr[molB]), remainderMask);
        const __m256d dummyMask =
            _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateB), _zero, _CMP_NEQ_OS);  // Assuming that dummy = 0
        const __m256d totalMask = _mm256_and_pd(cutoffMask, dummyMask);

        if (remainderCase) {
          _mm256_maskstore_pd((&molMask[molB - (molA + 1)]), remainderMask, totalMask);
        } else {
          _mm256_storeu_pd((&molMask[molB - (molA + 1)]), totalMask);
        }
      }

      // generate mask for each site in the mols 'above' molA from molecular mask
      std::vector<size_t, autopas::AlignedAllocator<size_t>> siteVector;
      siteVector.reserve(noSitesB);

      for (size_t molB = molA + 1, siteIndex = 0; molB < soa.getNumberOfParticles(); ++molB) {
        bool condition = static_cast<bool>(molMask[molB - (molA + 1)]);
        size_t siteCountB = _PPLibrary->getNumSites(typeptr[molB]);  // Maybe needs to change if no mixing
        if (condition) {
          for (size_t siteB = 0; siteB < siteCountB; ++siteB) {
            siteVector.emplace_back(siteIndex++);
          }
        } else {
          siteIndex += siteCountB;
        }
      }

      // ------- Main force calculation loop -------
      for (size_t siteA = siteIndexMolA; siteA < siteIndexMolB; ++siteA) {
        const double *mixingPtr = useMixing ? _PPLibrary->getMixingDataPtr(siteTypes[siteA], 0) : nullptr;

        // sums used for siteA
        __m256d forceSumX = _zero;
        __m256d forceSumY = _zero;
        __m256d forceSumZ = _zero;

        const __m256d exactSitePositionsAX = _mm256_broadcast_sd(&exactSitePositionX[siteA]);
        const __m256d exactSitePositionsAY = _mm256_broadcast_sd(&exactSitePositionY[siteA]);
        const __m256d exactSitePositionsAZ = _mm256_broadcast_sd(&exactSitePositionZ[siteA]);

        __m256d epsilon24 = _mm256_set1_pd(const_epsilon24);
        __m256d sigmaSquared = _mm256_set1_pd(const_sigmaSquared);
        __m256d shift6 = applyShift ? _mm256_set1_pd(const_shift6) : _zero;

        for (size_t siteVectorIndex = 0; siteVectorIndex < siteVector.size(); siteVectorIndex += vecLength) {
          const size_t rest = siteVector.size() - siteVectorIndex;
          const bool remainderCase = rest < vecLength;
          const __m256i remainderMask = remainderCase ? _masks[rest - 1] : _mm256_set1_epi64x(-1);

          const __m256i siteIndices =
              load_epi64(remainderCase, &siteVector[siteVectorIndex], remainderMask);

          if constexpr (useMixing) {
            __m256i siteTypeValues = gather_epi64(remainderCase, &siteTypes[siteIndexMolB], siteIndices, remainderMask);
            siteTypeValues =
                _mm256_mul_epu32(siteTypeValues, _mm256_set1_epi64x(3));  // This could maybe be problematic!

            epsilon24 = _mm256_i64gather_pd(mixingPtr, siteTypeValues, 8);
            sigmaSquared = _mm256_i64gather_pd(mixingPtr + 1, siteTypeValues, 8);
            if constexpr (applyShift) {
              shift6 = _mm256_i64gather_pd(mixingPtr + 2, siteTypeValues, 8);
            }
          }

          const __m256d exactSitePositionsBX =
              gather_pd(remainderCase, &exactSitePositionX[siteIndexMolB], siteIndices, remainderMask);
          const __m256d exactSitePositionsBY =
              gather_pd(remainderCase, &exactSitePositionY[siteIndexMolB], siteIndices, remainderMask);
          const __m256d exactSitePositionsBZ =
              gather_pd(remainderCase, &exactSitePositionZ[siteIndexMolB], siteIndices, remainderMask);
          const __m256d displacementX = _mm256_sub_pd(exactSitePositionsAX, exactSitePositionsBX);
          const __m256d displacementY = _mm256_sub_pd(exactSitePositionsAY, exactSitePositionsBY);
          const __m256d displacementZ = _mm256_sub_pd(exactSitePositionsAZ, exactSitePositionsBZ);

          const __m256d distanceSquaredX = _mm256_mul_pd(displacementX, displacementX);
          const __m256d distanceSquaredY = _mm256_mul_pd(displacementY, displacementY);
          const __m256d distanceSquaredZ = _mm256_mul_pd(displacementZ, displacementZ);

          const __m256d distanceSquared =
              _mm256_add_pd(distanceSquaredX, _mm256_add_pd(distanceSquaredY, distanceSquaredZ));

          const __m256d invDistSquared = _mm256_div_pd(_one, distanceSquared);
          const __m256d lj2 = _mm256_mul_pd(sigmaSquared, invDistSquared);
          const __m256d lj6 = _mm256_mul_pd(_mm256_mul_pd(lj2, lj2), lj2);
          const __m256d lj12 = _mm256_mul_pd(lj6, lj6);
          const __m256d lj12m6 = _mm256_sub_pd(lj12, lj6);
          const __m256d scalar = _mm256_mul_pd(epsilon24, _mm256_mul_pd(_mm256_add_pd(lj12, lj12m6), invDistSquared));
          const __m256d scalarMultiple =
              _mm256_and_pd(_mm256_castsi256_pd(remainderMask), scalar);  // Idk why this is needed

          const __m256d forceX = _mm256_mul_pd(scalarMultiple, displacementX);
          const __m256d forceY = _mm256_mul_pd(scalarMultiple, displacementY);
          const __m256d forceZ = _mm256_mul_pd(scalarMultiple, displacementZ);

          forceSumX = _mm256_add_pd(forceSumX, forceX);
          forceSumY = _mm256_add_pd(forceSumY, forceY);
          forceSumZ = _mm256_add_pd(forceSumZ, forceZ);

          // N3

          __m256d forceSumBX = gather_pd(remainderCase, &siteForceX[siteIndexMolB], siteIndices, remainderMask);
          __m256d forceSumBY = gather_pd(remainderCase, &siteForceY[siteIndexMolB], siteIndices, remainderMask);
          __m256d forceSumBZ = gather_pd(remainderCase, &siteForceZ[siteIndexMolB], siteIndices, remainderMask);

          forceSumBX = _mm256_sub_pd(forceSumBX, forceX);
          forceSumBY = _mm256_sub_pd(forceSumBY, forceY);
          forceSumBZ = _mm256_sub_pd(forceSumBZ, forceZ);

          scatter_pd(remainderCase, &siteForceX[siteIndexMolB], siteIndices, remainderMask, forceSumBX);
          scatter_pd(remainderCase, &siteForceY[siteIndexMolB], siteIndices, remainderMask, forceSumBY);
          scatter_pd(remainderCase, &siteForceZ[siteIndexMolB], siteIndices, remainderMask, forceSumBZ);

          if constexpr (calculateGlobals) {
            const __m256d virialX = _mm256_mul_pd(displacementX, forceX);
            const __m256d virialY = _mm256_mul_pd(displacementY, forceY);
            const __m256d virialZ = _mm256_mul_pd(displacementZ, forceZ);

            const __m256d potentialEnergy6 =
                _mm256_fmadd_pd(epsilon24, lj12m6, shift6);  // FMA may not be supported on all CPUs
            const __m256d potentialEnergy6Masked = _mm256_and_pd(remainderMask, potentialEnergy6);

            __m256i ownedStateA4 = _mm256_set1_epi64x(static_cast<int64_t>(ownedStateA));
            __m256d ownedMaskA = _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateA4),
                                               _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
            __m256d energyFactor = _mm256_blendv_pd(_zero, _one, ownedMaskA);
            if constexpr (newton3) {
              const __m256i ownedStateB =
                  gather_epi64(remainderCase, &isSiteOwned[siteIndexMolB], siteIndices, remainderMask);
              __m256d ownedMaskB = _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateB),
                                                 _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
              energyFactor = _mm256_add_pd(energyFactor, _mm256_blendv_pd(_zero, _one, ownedMaskB));
            }
            potentialEnergySum = _mm256_fmadd_pd(energyFactor, potentialEnergy6Masked, potentialEnergySum);
            virialSumX = _mm256_fmadd_pd(energyFactor, virialX, virialSumX);
            virialSumY = _mm256_fmadd_pd(energyFactor, virialY, virialSumY);
            virialSumZ = _mm256_fmadd_pd(energyFactor, virialZ, virialSumZ);
          }
        }
        // sum forces on single site in mol A
        siteForceX[siteA] += horizontalSum(forceSumX);
        siteForceY[siteA] += horizontalSum(forceSumY);
        siteForceZ[siteA] += horizontalSum(forceSumZ);
      }
      // ------------------- end of vectorized part -------------------
      siteIndexMolA += noSitesInMolA;
    }

    // reduce the forces on individual sites to forces & torques on whole molecules.
    if constexpr (useMixing) {
      size_t siteIndex = 0;
      for (size_t mol = 0; mol < soa.getNumberOfParticles(); ++mol) {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, _PPLibrary->getSitePositions(typeptr[mol]));
        for (size_t site = 0; site < _PPLibrary->getNumSites(typeptr[mol]); ++site) {
          fxptr[mol] += siteForceX[siteIndex];
          fyptr[mol] += siteForceY[siteIndex];
          fzptr[mol] += siteForceZ[siteIndex];
          txptr[mol] += rotatedSitePositions[site][1] * siteForceZ[siteIndex] -
                        rotatedSitePositions[site][2] * siteForceY[siteIndex];
          typtr[mol] += rotatedSitePositions[site][2] * siteForceX[siteIndex] -
                        rotatedSitePositions[site][0] * siteForceZ[siteIndex];
          tzptr[mol] += rotatedSitePositions[site][0] * siteForceY[siteIndex] -
                        rotatedSitePositions[site][1] * siteForceX[siteIndex];
          ++siteIndex;
        }
      }
    } else {
      size_t siteIndex = 0;
      for (size_t mol = 0; mol < soa.getNumberOfParticles(); mol++) {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, const_unrotatedSitePositions);
        for (size_t site = 0; site < const_unrotatedSitePositions.size(); ++site) {
          fxptr[mol] += siteForceX[siteIndex];
          fyptr[mol] += siteForceY[siteIndex];
          fzptr[mol] += siteForceZ[siteIndex];
          txptr[mol] += rotatedSitePositions[site][1] * siteForceZ[siteIndex] -
                        rotatedSitePositions[site][2] * siteForceY[siteIndex];
          typtr[mol] += rotatedSitePositions[site][2] * siteForceX[siteIndex] -
                        rotatedSitePositions[site][0] * siteForceZ[siteIndex];
          tzptr[mol] += rotatedSitePositions[site][0] * siteForceY[siteIndex] -
                        rotatedSitePositions[site][1] * siteForceX[siteIndex];
          ++siteIndex;
        }
      }
    }

    if constexpr (calculateGlobals) {
      const auto threadNum = autopas_get_thread_num();
      // SoAFunctorSingle obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
      // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.
      const auto newton3Factor = newton3 ? .5 : 1.;

      _aosThreadData[threadNum].potentialEnergySum += horizontalSum(potentialEnergySum) * newton3Factor;
      _aosThreadData[threadNum].virialSum[0] += horizontalSum(virialSumX) * newton3Factor;
      _aosThreadData[threadNum].virialSum[1] += horizontalSum(virialSumY) * newton3Factor;
      _aosThreadData[threadNum].virialSum[2] += horizontalSum(virialSumZ) * newton3Factor;
    }
  }

  /**
   * Implementation function of SoAFunctorPair(soa1, soa2, newton3)
   * @tparam newton3 flag for if newton's third law is used
   * @param soaA structure of arrays A
   * @param soaB structure of arrays B
   */
  template <bool newton3>
  void SoAFunctorPairImpl(SoAView<SoAArraysType> soaA, SoAView<SoAArraysType> soaB) {
#ifndef __AVX__
#pragma message "SoAFunctorPair functor called without AVX support!"
#endif
    if (soaA.getNumberOfParticles() == 0 || soaB.getNumberOfParticles() == 0) return;

    const auto *const __restrict xAptr = soaA.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yAptr = soaA.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zAptr = soaA.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict xBptr = soaB.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yBptr = soaB.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zBptr = soaB.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtrA = soaA.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtrB = soaB.template begin<Particle::AttributeNames::ownershipState>();

    const auto *const __restrict q0Aptr = soaA.template begin<Particle::AttributeNames::quaternion0>();
    const auto *const __restrict q1Aptr = soaA.template begin<Particle::AttributeNames::quaternion1>();
    const auto *const __restrict q2Aptr = soaA.template begin<Particle::AttributeNames::quaternion2>();
    const auto *const __restrict q3Aptr = soaA.template begin<Particle::AttributeNames::quaternion3>();
    const auto *const __restrict q0Bptr = soaB.template begin<Particle::AttributeNames::quaternion0>();
    const auto *const __restrict q1Bptr = soaB.template begin<Particle::AttributeNames::quaternion1>();
    const auto *const __restrict q2Bptr = soaB.template begin<Particle::AttributeNames::quaternion2>();
    const auto *const __restrict q3Bptr = soaB.template begin<Particle::AttributeNames::quaternion3>();

    SoAFloatPrecision *const __restrict fxAptr = soaA.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyAptr = soaA.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzAptr = soaA.template begin<Particle::AttributeNames::forceZ>();
    SoAFloatPrecision *const __restrict fxBptr = soaB.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyBptr = soaB.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzBptr = soaB.template begin<Particle::AttributeNames::forceZ>();

    SoAFloatPrecision *const __restrict txAptr = soaA.template begin<Particle::AttributeNames::torqueX>();
    SoAFloatPrecision *const __restrict tyAptr = soaA.template begin<Particle::AttributeNames::torqueY>();
    SoAFloatPrecision *const __restrict tzAptr = soaA.template begin<Particle::AttributeNames::torqueZ>();
    SoAFloatPrecision *const __restrict txBptr = soaB.template begin<Particle::AttributeNames::torqueX>();
    SoAFloatPrecision *const __restrict tyBptr = soaB.template begin<Particle::AttributeNames::torqueY>();
    SoAFloatPrecision *const __restrict tzBptr = soaB.template begin<Particle::AttributeNames::torqueZ>();

    [[maybe_unused]] auto *const __restrict typeptrA = soaA.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptrB = soaB.template begin<Particle::AttributeNames::typeId>();

    __m256d virialSumX = _mm256_setzero_pd();
    __m256d virialSumY = _mm256_setzero_pd();
    __m256d virialSumZ = _mm256_setzero_pd();
    __m256d potentialEnergySum = _mm256_setzero_pd();

    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const SoAFloatPrecision cutoffSquared = _cutoffSquaredAoS;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionBx;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionBy;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionBz;

    // we require arrays for forces for sites to maintain SIMD in site-site calculations
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceBx;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceBy;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceBz;

    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteTypesB;
    std::vector<size_t, autopas::AlignedAllocator<size_t>> isSiteOwnedBArr;

    const SoAFloatPrecision const_sigmaSquared = _sigmaSquaredAoS;
    const SoAFloatPrecision const_epsilon24 = _epsilon24AoS;
    const SoAFloatPrecision const_shift6 = _shift6AoS;

    const auto const_unrotatedSitePositions = _sitePositionsLJ;

    // count number of sites in both SoAs
    size_t siteCountB = 0;
    if constexpr (useMixing) {
      for (size_t mol = 0; mol < soaB.getNumberOfParticles(); ++mol) {
        siteCountB += _PPLibrary->getNumSites(typeptrB[mol]);
      }
    } else {
      siteCountB = const_unrotatedSitePositions.size() * soaB.getNumberOfParticles();
    }

    // pre-reserve std::vectors
    exactSitePositionBx.reserve(siteCountB);
    exactSitePositionBy.reserve(siteCountB);
    exactSitePositionBz.reserve(siteCountB);

    if constexpr (useMixing) {
      siteTypesB.reserve(siteCountB);
    }

    siteForceBx.reserve(siteCountB);
    siteForceBy.reserve(siteCountB);
    siteForceBz.reserve(siteCountB);

    if constexpr (calculateGlobals) {
      // this is only needed for vectorization when calculating globals
      isSiteOwnedBArr.reserve(siteCountB);
    }

    // Fill site-wise std::vectors for SIMD
    std::vector<std::array<double, 3>> rotatedSitePositions;
    std::fill_n(siteForceBx.begin(), siteCountB, 0.);
    std::fill_n(siteForceBy.begin(), siteCountB, 0.);
    std::fill_n(siteForceBz.begin(), siteCountB, 0.);
    for (size_t mol = 0; mol < soaB.getNumberOfParticles(); ++mol) {
      if constexpr (useMixing) {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0Bptr[mol], q1Bptr[mol], q2Bptr[mol], q3Bptr[mol]}, _PPLibrary->getSitePositions(typeptrB[mol]));
      } else {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0Bptr[mol], q1Bptr[mol], q2Bptr[mol], q3Bptr[mol]}, const_unrotatedSitePositions);
      }

      size_t localSiteCount = useMixing ? _PPLibrary->getNumSites(typeptrB[mol]) : const_unrotatedSitePositions.size();
      for (size_t site = 0; site < localSiteCount; ++site) {
        exactSitePositionBx.push_back(rotatedSitePositions[site][0] + xBptr[mol]);
        exactSitePositionBy.push_back(rotatedSitePositions[site][1] + yBptr[mol]);
        exactSitePositionBz.push_back(rotatedSitePositions[site][2] + zBptr[mol]);
        if constexpr (calculateGlobals) {
          isSiteOwnedBArr.push_back(ownedStatePtrB[mol] == OwnershipState::owned);
        }
        if constexpr (useMixing) {
          siteTypesB.push_back(_PPLibrary->getSiteTypes(typeptrB[mol])[site]);
        }
      }
    }

    // main force calculation loop
    for (size_t molA = 0; molA < soaA.getNumberOfParticles(); ++molA) {
      const auto ownedStateA = ownedStatePtrA[molA];
      if (ownedStateA == autopas::OwnershipState::dummy) {
        continue;
      }

      const auto noSitesInMolA =
          useMixing ? _PPLibrary->getNumSites(typeptrA[molA]) : const_unrotatedSitePositions.size();
      const auto unrotatedSitePositionsA =
          useMixing ? _PPLibrary->getSitePositions(typeptrA[molA]) : const_unrotatedSitePositions;

      const auto rotatedSitePositionsA = autopas::utils::quaternion::rotateVectorOfPositions(
          {q0Aptr[molA], q1Aptr[molA], q2Aptr[molA], q3Aptr[molA]}, unrotatedSitePositionsA);

      // create mask over every mol in cell B (char to keep arrays aligned)
      std::vector<double, autopas::AlignedAllocator<double>> molMask;
      molMask.reserve(soaB.getNumberOfParticles());

      const __m256d xposA = _mm256_broadcast_sd(&xAptr[molA]);
      const __m256d yposA = _mm256_broadcast_sd(&yAptr[molA]);
      const __m256d zposA = _mm256_broadcast_sd(&zAptr[molA]);

      for (size_t molB = 0; molB < soaB.getNumberOfParticles(); molB += vecLength) {
        const size_t rest = soaB.getNumberOfParticles() - molB;
        const bool remainderCase = rest < vecLength;
        __m256i remainderMask = remainderCase ? _masks[rest - 1] : _mm256_set1_epi64x(-1);

        const __m256d xposB = load_pd(remainderCase, &xBptr[molB], remainderMask);
        const __m256d yposB = load_pd(remainderCase, &yBptr[molB], remainderMask);
        const __m256d zposB = load_pd(remainderCase, &zBptr[molB], remainderMask);

        const __m256d displacementCoMX = _mm256_sub_pd(xposA, xposB);
        const __m256d displacementCoMY = _mm256_sub_pd(yposA, yposB);
        const __m256d displacementCoMZ = _mm256_sub_pd(zposA, zposB);

        const __m256d distanceSquaredCoMX = _mm256_mul_pd(displacementCoMX, displacementCoMX);
        const __m256d distanceSquaredCoMY = _mm256_mul_pd(displacementCoMY, displacementCoMY);
        const __m256d distanceSquaredCoMZ = _mm256_mul_pd(displacementCoMZ, displacementCoMZ);

        const __m256d distanceSquaredCoM =
            _mm256_add_pd(_mm256_add_pd(distanceSquaredCoMX, distanceSquaredCoMY), distanceSquaredCoMZ);

        const __m256d cutoffMask = _mm256_cmp_pd(distanceSquaredCoM, _cutoffSquared, _CMP_LE_OS);
        const __m256i ownedStateB =
            remainderCase
                ? _mm256_maskload_epi64(reinterpret_cast<const long long *>(&ownedStatePtrB[molB]), remainderMask)
                : _mm256_loadu_si256((__m256i *)&ownedStatePtrB[molB]);
        const __m256d dummyMask =
            _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateB), _zero, _CMP_NEQ_OS);  // Assuming that dummy = 0
        const __m256d totalMask = _mm256_and_pd(cutoffMask, dummyMask);

        if (remainderCase) {
          _mm256_maskstore_pd((&molMask[molB]), remainderMask, totalMask);
        } else {
          _mm256_storeu_pd((&molMask[molB]), totalMask);
        }
      }

      // generate mask for each site in cell B from molecular mask
      std::vector<size_t, autopas::AlignedAllocator<size_t>> siteVector;
      siteVector.reserve(siteCountB);

      for (size_t molB = 0, siteIndex = 0; molB < soaB.getNumberOfParticles(); ++molB) {
        bool condition = static_cast<bool>(molMask[molB]);
        size_t siteCount = _PPLibrary->getNumSites(typeptrB[molB]);  // TODO: change if no mixing is used
        if (condition) {
          for (size_t siteB = 0; siteB < siteCount; ++siteB) {
            siteVector.emplace_back(siteIndex++);
          }
        } else {
          siteIndex += siteCount;
        }
      }

      // sums used for molA
      __m256d forceSumX = _zero;
      __m256d forceSumY = _zero;
      __m256d forceSumZ = _zero;
      __m256d torqueSumX = _zero;
      __m256d torqueSumY = _zero;
      __m256d torqueSumZ = _zero;

      for (size_t siteA = 0; siteA < noSitesInMolA; ++siteA) {
        const double *mixingPtr =
            useMixing ? _PPLibrary->getMixingDataPtr(_PPLibrary->getSiteTypes(typeptrA[molA])[siteA], 0) : nullptr;

        const __m256d rotatedSitePositionsAX = _mm256_broadcast_sd(&rotatedSitePositionsA[siteA][0]);
        const __m256d rotatedSitePositionsAY = _mm256_broadcast_sd(&rotatedSitePositionsA[siteA][1]);
        const __m256d rotatedSitePositionsAZ = _mm256_broadcast_sd(&rotatedSitePositionsA[siteA][2]);

        const __m256d exactSitePositionsAX = _mm256_add_pd(rotatedSitePositionsAX, xposA);
        const __m256d exactSitePositionsAY = _mm256_add_pd(rotatedSitePositionsAY, yposA);
        const __m256d exactSitePositionsAZ = _mm256_add_pd(rotatedSitePositionsAZ, zposA);

        __m256d sigmaSquared = _mm256_set1_pd(const_sigmaSquared);
        __m256d epsilon24 = _mm256_set1_pd(const_epsilon24);
        __m256d shift6 = applyShift ? _mm256_set1_pd(const_shift6) : _zero;

        for (size_t siteVectorIndex = 0; siteVectorIndex < siteVector.size(); siteVectorIndex += vecLength) {
          const size_t remainder = siteVector.size() - siteVectorIndex;
          const bool remainderCase = remainder < vecLength;
          const __m256i remainderMask = remainderCase ? _masks[remainder - 1] : _mm256_set1_epi64x(-1);

          const __m256i siteIndices =
              load_epi64(remainderCase, &siteVector[siteVectorIndex], remainderMask);

          if constexpr (useMixing) {
            __m256i siteTypeValues = gather_epi64(remainderCase, siteTypesB.data(), siteIndices, remainderMask);
            siteTypeValues =
                _mm256_mul_epu32(siteTypeValues, _mm256_set1_epi64x(3));  // This could maybe be problematic!

            epsilon24 = _mm256_i64gather_pd(mixingPtr, siteTypeValues, 8);
            sigmaSquared = _mm256_i64gather_pd(mixingPtr + 1, siteTypeValues, 8);
            if constexpr (applyShift) {
              shift6 = _mm256_i64gather_pd(mixingPtr + 2, siteTypeValues, 8);
            }
          }

          const __m256d exactSitePositionsBX =
              gather_pd(remainderCase, exactSitePositionBx.data(), siteIndices, remainderMask);
          const __m256d exactSitePositionsBY =
              gather_pd(remainderCase, exactSitePositionBy.data(), siteIndices, remainderMask);
          const __m256d exactSitePositionsBZ =
              gather_pd(remainderCase, exactSitePositionBz.data(), siteIndices, remainderMask);

          const __m256d displacementX = _mm256_sub_pd(exactSitePositionsAX, exactSitePositionsBX);
          const __m256d displacementY = _mm256_sub_pd(exactSitePositionsAY, exactSitePositionsBY);
          const __m256d displacementZ = _mm256_sub_pd(exactSitePositionsAZ, exactSitePositionsBZ);

          const __m256d distanceSquaredX = _mm256_mul_pd(displacementX, displacementX);
          const __m256d distanceSquaredY = _mm256_mul_pd(displacementY, displacementY);
          const __m256d distanceSquaredZ = _mm256_mul_pd(displacementZ, displacementZ);

          const __m256d distanceSquared =
              _mm256_add_pd(distanceSquaredX, _mm256_add_pd(distanceSquaredY, distanceSquaredZ));

          const __m256d invDistSquared = _mm256_div_pd(_one, distanceSquared);
          const __m256d lj2 = _mm256_mul_pd(sigmaSquared, invDistSquared);
          const __m256d lj6 = _mm256_mul_pd(_mm256_mul_pd(lj2, lj2), lj2);
          const __m256d lj12 = _mm256_mul_pd(lj6, lj6);
          const __m256d lj12m6 = _mm256_sub_pd(lj12, lj6);
          const __m256d scalar = _mm256_mul_pd(epsilon24, _mm256_mul_pd(_mm256_add_pd(lj12, lj12m6), invDistSquared));
          const __m256d scalarMultiple =
              _mm256_and_pd(_mm256_castsi256_pd(remainderMask), scalar);  // Idk why this is needed

          // calculate forces
          const __m256d forceX = _mm256_mul_pd(scalarMultiple, displacementX);
          const __m256d forceY = _mm256_mul_pd(scalarMultiple, displacementY);
          const __m256d forceZ = _mm256_mul_pd(scalarMultiple, displacementZ);

          const __m256d torqueAX =
              _mm256_fmsub_pd(rotatedSitePositionsAY, forceZ, _mm256_mul_pd(rotatedSitePositionsAZ, forceY));
          const __m256d torqueAY =
              _mm256_fmsub_pd(rotatedSitePositionsAZ, forceX, _mm256_mul_pd(rotatedSitePositionsAX, forceZ));
          const __m256d torqueAZ =
              _mm256_fmsub_pd(rotatedSitePositionsAX, forceY, _mm256_mul_pd(rotatedSitePositionsAY, forceX));

          forceSumX = _mm256_add_pd(forceSumX, forceX);
          forceSumY = _mm256_add_pd(forceSumY, forceY);
          forceSumZ = _mm256_add_pd(forceSumZ, forceZ);

          torqueSumX = _mm256_add_pd(torqueSumX, torqueAX);
          torqueSumY = _mm256_add_pd(torqueSumY, torqueAY);
          torqueSumZ = _mm256_add_pd(torqueSumZ, torqueAZ);

          // N3L ( total molecular forces + torques to be determined later )
          if constexpr (newton3) {
            __m256d forceSumBX = gather_pd(remainderCase, siteForceBx.data(), siteIndices, remainderMask);
            __m256d forceSumBY = gather_pd(remainderCase, siteForceBy.data(), siteIndices, remainderMask);
            __m256d forceSumBZ = gather_pd(remainderCase, siteForceBz.data(), siteIndices, remainderMask);
            forceSumBX = _mm256_sub_pd(forceSumBX, forceX);
            forceSumBY = _mm256_sub_pd(forceSumBY, forceY);
            forceSumBZ = _mm256_sub_pd(forceSumBZ, forceZ);

            scatter_pd(remainderCase, siteForceBx.data(), siteIndices, remainderMask, forceSumBX);
            scatter_pd(remainderCase, siteForceBy.data(), siteIndices, remainderMask, forceSumBY);
            scatter_pd(remainderCase, siteForceBz.data(), siteIndices, remainderMask, forceSumBZ);
          }

          // globals
          if constexpr (calculateGlobals) {
            const __m256d virialX = _mm256_mul_pd(displacementX, forceX);
            const __m256d virialY = _mm256_mul_pd(displacementY, forceY);
            const __m256d virialZ = _mm256_mul_pd(displacementZ, forceZ);

            const __m256d potentialEnergy6 =
                _mm256_fmadd_pd(epsilon24, lj12m6, shift6);  // FMA may not be supported on all CPUs
            const __m256d potentialEnergy6Masked = _mm256_and_pd(remainderMask, potentialEnergy6);

            __m256i ownedStateA4 = _mm256_set1_epi64x(static_cast<int64_t>(ownedStateA));
            __m256d ownedMaskA = _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateA4),
                                               _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
            __m256d energyFactor = _mm256_blendv_pd(_zero, _one, ownedMaskA);
            if constexpr (newton3) {
              const __m256i ownedStateB =
                  gather_epi64(remainderCase, isSiteOwnedBArr.data(), siteIndices, remainderMask);
              __m256d ownedMaskB = _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateB),
                                                 _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
              energyFactor = _mm256_add_pd(energyFactor, _mm256_blendv_pd(_zero, _one, ownedMaskB));
            }
            potentialEnergySum = _mm256_fmadd_pd(energyFactor, potentialEnergy6Masked, potentialEnergySum);
            virialSumX = _mm256_fmadd_pd(energyFactor, virialX, virialSumX);
            virialSumY = _mm256_fmadd_pd(energyFactor, virialY, virialSumY);
            virialSumZ = _mm256_fmadd_pd(energyFactor, virialZ, virialSumZ);
          }
        }
      }
      fxAptr[molA] += horizontalSum(forceSumX);
      fyAptr[molA] += horizontalSum(forceSumY);
      fzAptr[molA] += horizontalSum(forceSumZ);
      txAptr[molA] += horizontalSum(torqueSumX);
      tyAptr[molA] += horizontalSum(torqueSumY);
      tzAptr[molA] += horizontalSum(torqueSumZ);
    }

    // reduce the forces on individual sites in SoA B to total forces & torques on whole molecules
    if constexpr (useMixing) {
      size_t siteIndex = 0;
      for (size_t mol = 0; mol < soaB.getNumberOfParticles(); ++mol) {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0Bptr[mol], q1Bptr[mol], q2Bptr[mol], q3Bptr[mol]}, _PPLibrary->getSitePositions(typeptrB[mol]));
        for (size_t site = 0; site < _PPLibrary->getNumSites(typeptrB[mol]); ++site) {
          fxBptr[mol] += siteForceBx[siteIndex];
          fyBptr[mol] += siteForceBy[siteIndex];
          fzBptr[mol] += siteForceBz[siteIndex];
          txBptr[mol] += rotatedSitePositions[site][1] * siteForceBz[siteIndex] -
                         rotatedSitePositions[site][2] * siteForceBy[siteIndex];
          tyBptr[mol] += rotatedSitePositions[site][2] * siteForceBx[siteIndex] -
                         rotatedSitePositions[site][0] * siteForceBz[siteIndex];
          tzBptr[mol] += rotatedSitePositions[site][0] * siteForceBy[siteIndex] -
                         rotatedSitePositions[site][1] * siteForceBx[siteIndex];
          ++siteIndex;
        }
      }
    } else {
      size_t siteIndex = 0;
      for (size_t mol = 0; mol < soaB.getNumberOfParticles(); ++mol) {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0Bptr[mol], q1Bptr[mol], q2Bptr[mol], q3Bptr[mol]}, const_unrotatedSitePositions);
        for (size_t site = 0; site < const_unrotatedSitePositions.size(); ++site) {
          fxBptr[mol] += siteForceBx[siteIndex];
          fyBptr[mol] += siteForceBy[siteIndex];
          fzBptr[mol] += siteForceBz[siteIndex];
          txBptr[mol] += rotatedSitePositions[site][1] * siteForceBz[siteIndex] -
                         rotatedSitePositions[site][2] * siteForceBy[siteIndex];
          tyBptr[mol] += rotatedSitePositions[site][2] * siteForceBx[siteIndex] -
                         rotatedSitePositions[site][0] * siteForceBz[siteIndex];
          tzBptr[mol] += rotatedSitePositions[site][0] * siteForceBy[siteIndex] -
                         rotatedSitePositions[site][1] * siteForceBx[siteIndex];
          ++siteIndex;
        }
      }
    }
    if constexpr (calculateGlobals) {
      const auto threadNum = autopas_get_thread_num();
      // SoAFunctorPairImpl obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
      // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.
      const auto newton3Factor = newton3 ? .5 : 1.;

      _aosThreadData[threadNum].potentialEnergySum += horizontalSum(potentialEnergySum) * newton3Factor;
      _aosThreadData[threadNum].virialSum[0] += horizontalSum(virialSumX) * newton3Factor;
      _aosThreadData[threadNum].virialSum[1] += horizontalSum(virialSumY) * newton3Factor;
      _aosThreadData[threadNum].virialSum[2] += horizontalSum(virialSumZ) * newton3Factor;
    }
  }

  // TODO Documentation
  template <bool newton3>
  void SoAFunctorVerletImpl(SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
#ifndef __AVX__
#pragma message "SoAFunctorVerlet functor called without AVX support!"
#endif
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    // Skip if primary particle is dummy
    const auto ownedStatePrime = ownedStatePtr[indexFirst];
    if (ownedStatePrime == OwnershipState::dummy) {
      return;
    }

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict q0ptr = soa.template begin<Particle::AttributeNames::quaternion0>();
    const auto *const __restrict q1ptr = soa.template begin<Particle::AttributeNames::quaternion1>();
    const auto *const __restrict q2ptr = soa.template begin<Particle::AttributeNames::quaternion2>();
    const auto *const __restrict q3ptr = soa.template begin<Particle::AttributeNames::quaternion3>();

    SoAFloatPrecision *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    SoAFloatPrecision *const __restrict txptr = soa.template begin<Particle::AttributeNames::torqueX>();
    SoAFloatPrecision *const __restrict typtr = soa.template begin<Particle::AttributeNames::torqueY>();
    SoAFloatPrecision *const __restrict tzptr = soa.template begin<Particle::AttributeNames::torqueZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();

    __m256d virialSumX = _mm256_setzero_pd();
    __m256d virialSumY = _mm256_setzero_pd();
    __m256d virialSumZ = _mm256_setzero_pd();
    __m256d potentialEnergySum = _mm256_setzero_pd();

    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const SoAFloatPrecision cutoffSquared = _cutoffSquaredAoS;
    const auto const_unrotatedSitePositions = _sitePositionsLJ;

    const auto const_sigmaSquared = _sigmaSquaredAoS;
    const auto const_epsilon24 = _epsilon24AoS;
    const auto const_shift6 = _shift6AoS;

    const size_t neighborListSize = neighborList.size();
    const size_t *const __restrict neighborListPtr = neighborList.data();

    // Count sites
    const size_t siteCountMolPrime =
        useMixing ? _PPLibrary->getNumSites(typeptr[indexFirst]) : const_unrotatedSitePositions.size();

    size_t siteCountNeighbors = 0;  // site count of neighbours of primary molecule
    if constexpr (useMixing) {
      for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
        siteCountNeighbors += _PPLibrary->getNumSites(typeptr[neighborList[neighborMol]]);
      }
    } else {
      siteCountNeighbors = const_unrotatedSitePositions.size() * neighborListSize;
    }

    // initialize site-wise arrays for neighbors
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactNeighborSitePositionsX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactNeighborSitePositionsY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactNeighborSitePositionsZ;

    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteTypesNeighbors;
    std::vector<size_t, autopas::AlignedAllocator<size_t>> isNeighborSiteOwnedArr;

    // we require arrays for forces for sites to maintain SIMD in site-site calculations
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceZ;

    // pre-reserve arrays
    exactNeighborSitePositionsX.reserve(siteCountNeighbors);
    exactNeighborSitePositionsY.reserve(siteCountNeighbors);
    exactNeighborSitePositionsZ.reserve(siteCountNeighbors);

    if constexpr (useMixing) {
      siteTypesNeighbors.reserve(siteCountNeighbors);
    }

    siteForceX.reserve(siteCountNeighbors);
    siteForceY.reserve(siteCountNeighbors);
    siteForceZ.reserve(siteCountNeighbors);

    if constexpr (calculateGlobals) {
      isNeighborSiteOwnedArr.reserve(siteCountNeighbors);
    }

    const auto rotatedSitePositionsPrime =
        useMixing ? autopas::utils::quaternion::rotateVectorOfPositions(
                        {q0ptr[indexFirst], q1ptr[indexFirst], q2ptr[indexFirst], q3ptr[indexFirst]},
                        _PPLibrary->getSitePositions(typeptr[indexFirst]))
                  : autopas::utils::quaternion::rotateVectorOfPositions(
                        {q0ptr[indexFirst], q1ptr[indexFirst], q2ptr[indexFirst], q3ptr[indexFirst]},
                        const_unrotatedSitePositions);

    const auto siteTypesPrime = _PPLibrary->getSiteTypes(typeptr[indexFirst]);  // todo make this work for no mixing

    std::vector<std::array<double, 3>> rotatedSitePositions;
    std::fill_n(siteForceX.begin(), siteCountNeighbors, 0.);
    std::fill_n(siteForceY.begin(), siteCountNeighbors, 0.);
    std::fill_n(siteForceZ.begin(), siteCountNeighbors, 0.);
    for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
      const size_t neighborMolIndex = neighborList[neighborMol];
      if constexpr (useMixing) {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[neighborMolIndex], q1ptr[neighborMolIndex], q2ptr[neighborMolIndex], q3ptr[neighborMolIndex]},
            _PPLibrary->getSitePositions(typeptr[neighborMolIndex]));
      } else {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[neighborMolIndex], q1ptr[neighborMolIndex], q2ptr[neighborMolIndex], q3ptr[neighborMolIndex]},
            const_unrotatedSitePositions);
      }

      size_t localSiteCount =
          useMixing ? _PPLibrary->getNumSites(typeptr[neighborMolIndex]) : const_unrotatedSitePositions.size();
      for (size_t site = 0; site < localSiteCount; ++site) {
        exactNeighborSitePositionsX.push_back(rotatedSitePositions[site][0] + xptr[neighborMolIndex]);
        exactNeighborSitePositionsY.push_back(rotatedSitePositions[site][1] + yptr[neighborMolIndex]);
        exactNeighborSitePositionsZ.push_back(rotatedSitePositions[site][2] + zptr[neighborMolIndex]);
        if constexpr (calculateGlobals) {
          isNeighborSiteOwnedArr.push_back(ownedStatePtr[neighborMolIndex] == OwnershipState::owned);
        }
        if constexpr (useMixing) {
          siteTypesNeighbors.push_back(_PPLibrary->getSiteTypes(typeptr[neighborMolIndex])[site]);
        }
      }
    }

    // -- main force calculation --

    // - calculate mol mask -
    std::vector<double, autopas::AlignedAllocator<double>> molMask;
    molMask.reserve(neighborListSize);

    const __m256d xPosFirst = _mm256_broadcast_sd(&xptr[indexFirst]);
    const __m256d yPosFirst = _mm256_broadcast_sd(&yptr[indexFirst]);
    const __m256d zPosFirst = _mm256_broadcast_sd(&zptr[indexFirst]);

    // WIP vectorized version

//    for (size_t neighborMol = 0; neighborMol < neighborListSize; neighborMol += vecLength) {
//      const size_t rest = neighborListSize - neighborMol;
//      const bool remainderCase = rest < vecLength;
//      __m256i remainderMask = remainderCase ? _masks[rest - 1] : _mm256_set1_epi64x(-1);
//      const __m256i neighborMolIndex = load_epi64(remainderCase, &neighborList[neighborMol], remainderMask);
//
//      const __m256d xposB = gather_pd(remainderCase, xptr, neighborMolIndex, remainderMask);
//      const __m256d yposB = gather_pd(remainderCase, yptr, neighborMolIndex, remainderMask);
//      const __m256d zposB = gather_pd(remainderCase, zptr, neighborMolIndex, remainderMask);
//
//      const __m256d displacementCoMX = _mm256_sub_pd(xPosFirst, xposB);
//      const __m256d displacementCoMY = _mm256_sub_pd(yPosFirst, yposB);
//      const __m256d displacementCoMZ = _mm256_sub_pd(zPosFirst, zposB);
//
//      const __m256d distanceSquaredCoMX = _mm256_mul_pd(displacementCoMX, displacementCoMX);
//      const __m256d distanceSquaredCoMY = _mm256_mul_pd(displacementCoMY, displacementCoMY);
//      const __m256d distanceSquaredCoMZ = _mm256_mul_pd(displacementCoMZ, displacementCoMZ);
//
//      const __m256d distanceSquaredCoM =
//          _mm256_add_pd(_mm256_add_pd(distanceSquaredCoMX, distanceSquaredCoMY), distanceSquaredCoMZ);
//
//      const __m256d cutoffMask = _mm256_cmp_pd(distanceSquaredCoM, _cutoffSquared, _CMP_LE_OS);
//      const __m256i ownedStateB =
//          remainderCase ? _mm256_mask_i64gather_epi64(_zero, ownedStatePtr, neighborMolIndex, remainderMask, 8)
//                        : _mm256_i64gather_pd(ownedStatePtr, neighborMolIndex, 8);
//      const __m256d dummyMask =
//          _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateB), _zero, _CMP_NEQ_OS);  // Assuming that dummy = 0
//      const __m256d totalMask = _mm256_and_pd(cutoffMask, dummyMask);
//
//      if (remainderCase) {
//        _mm256_maskstore_pd((&molMask[neighborMol]), remainderMask, totalMask);
//      } else {
//        _mm256_storeu_pd((&molMask[neighborMol]), totalMask);
//      }
//    }
    unsigned long max_long = 0xFFFFFFFFFFFFFFFF;
    double max_double = *(double *)&max_long;
    for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
      const auto neighborMolIndex = neighborList[neighborMol];  // index of neighbor mol in soa

      const auto ownedState = ownedStatePtr[neighborMolIndex];

      const auto displacementCoMX = xptr[indexFirst] - xptr[neighborMolIndex];
      const auto displacementCoMY = yptr[indexFirst] - yptr[neighborMolIndex];
      const auto displacementCoMZ = zptr[indexFirst] - zptr[neighborMolIndex];

      const auto distanceSquaredCoMX = displacementCoMX * displacementCoMX;
      const auto distanceSquaredCoMY = displacementCoMY * displacementCoMY;
      const auto distanceSquaredCoMZ = displacementCoMZ * displacementCoMZ;

      const auto distanceSquaredCoM = distanceSquaredCoMX + distanceSquaredCoMY + distanceSquaredCoMZ;

      // mask molecules beyond cutoff or if molecule is a dummy
      bool condition = distanceSquaredCoM <= cutoffSquared and ownedState != autopas::OwnershipState::dummy;
      molMask[neighborMol] = condition ? max_double : 0.;
    }

    //     generate mask for each site from molecular mask
    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteVector;
    siteVector.reserve(siteCountNeighbors);

    for (size_t neighborMol = 0, siteIndex = 0; neighborMol < neighborListSize; ++neighborMol) {
      bool condition = static_cast<bool>(molMask[neighborMol]);
      const auto neighborMolIndex = neighborList[neighborMol];
      size_t siteCount = _PPLibrary->getNumSites(typeptr[neighborMolIndex]);  // TODO: change if no mixing is used
      if (condition) {
        for (size_t siteB = 0; siteB < siteCount; ++siteB) {
          siteVector.emplace_back(siteIndex++);
        }
      } else {
        siteIndex += siteCount;
      }
    }

    // sums used for prime mol
    __m256d forceSumX = _zero;
    __m256d forceSumY = _zero;
    __m256d forceSumZ = _zero;
    __m256d torqueSumX = _zero;
    __m256d torqueSumY = _zero;
    __m256d torqueSumZ = _zero;

    __m256d sigmaSquared = _mm256_set1_pd(const_sigmaSquared);
    __m256d epsilon24 = _mm256_set1_pd(const_epsilon24);
    __m256d shift6 = applyShift ? _mm256_set1_pd(const_shift6) : _zero;

    // - actual LJ calculation -

    for (size_t primeSite = 0; primeSite < siteCountMolPrime; ++primeSite) {
      const double *mixingPtr = useMixing ? _PPLibrary->getMixingDataPtr(siteTypesPrime[primeSite], 0) : nullptr;
      const __m256d rotatedPrimeSitePositionX = _mm256_broadcast_sd(&rotatedSitePositionsPrime[primeSite][0]);
      const __m256d rotatedPrimeSitePositionY = _mm256_broadcast_sd(&rotatedSitePositionsPrime[primeSite][1]);
      const __m256d rotatedPrimeSitePositionZ = _mm256_broadcast_sd(&rotatedSitePositionsPrime[primeSite][2]);

      const __m256d exactPrimeSitePositionX = _mm256_add_pd(rotatedPrimeSitePositionX, xPosFirst);
      const __m256d exactPrimeSitePositionY = _mm256_add_pd(rotatedPrimeSitePositionY, yPosFirst);
      const __m256d exactPrimeSitePositionZ = _mm256_add_pd(rotatedPrimeSitePositionZ, zPosFirst);

      for (size_t siteVectorIndex = 0; siteVectorIndex < siteVector.size(); siteVectorIndex += vecLength) {
        const size_t remainder = siteVector.size() - siteVectorIndex;
        const bool remainderCase = remainder < vecLength;
        const __m256i remainderMask = remainderCase ? _masks[remainder - 1] : _mm256_set1_epi64x(-1);

        const __m256i siteIndices =
            load_epi64(remainderCase, &siteVector[siteVectorIndex], remainderMask);

        if constexpr (useMixing) {

          __m256i siteTypeValues = gather_epi64(remainderCase, siteTypesNeighbors.data(), siteIndices, remainderMask);
          siteTypeValues = _mm256_mul_epu32(siteTypeValues, _mm256_set1_epi64x(3));  // This could maybe be problematic!

          epsilon24 = _mm256_i64gather_pd(mixingPtr, siteTypeValues, 8);
          sigmaSquared = _mm256_i64gather_pd(mixingPtr + 1, siteTypeValues, 8);
          if constexpr (applyShift) {
            shift6 = _mm256_i64gather_pd(mixingPtr + 2, siteTypeValues, 8);
          }
        }

        const __m256d exactNeighborSitePositionX =
            gather_pd(remainderCase, exactNeighborSitePositionsX.data(), siteIndices, remainderMask);
        const __m256d exactNeighborSitePositionY =
            gather_pd(remainderCase, exactNeighborSitePositionsY.data(), siteIndices, remainderMask);
        const __m256d exactNeighborSitePositionZ =
            gather_pd(remainderCase, exactNeighborSitePositionsZ.data(), siteIndices, remainderMask);

        const __m256d displacementX = _mm256_sub_pd(exactPrimeSitePositionX, exactNeighborSitePositionX);
        const __m256d displacementY = _mm256_sub_pd(exactPrimeSitePositionY, exactNeighborSitePositionY);
        const __m256d displacementZ = _mm256_sub_pd(exactPrimeSitePositionZ, exactNeighborSitePositionZ);

        const __m256d distanceSquaredX = _mm256_mul_pd(displacementX, displacementX);
        const __m256d distanceSquaredY = _mm256_mul_pd(displacementY, displacementY);
        const __m256d distanceSquaredZ = _mm256_mul_pd(displacementZ, displacementZ);

        const __m256d distanceSquared =
            _mm256_add_pd(distanceSquaredX, _mm256_add_pd(distanceSquaredY, distanceSquaredZ));

        const __m256d invDistSquared = _mm256_div_pd(_one, distanceSquared);
        const __m256d lj2 = _mm256_mul_pd(sigmaSquared, invDistSquared);
        const __m256d lj6 = _mm256_mul_pd(_mm256_mul_pd(lj2, lj2), lj2);
        const __m256d lj12 = _mm256_mul_pd(lj6, lj6);
        const __m256d lj12m6 = _mm256_sub_pd(lj12, lj6);
        const __m256d scalar = _mm256_mul_pd(epsilon24, _mm256_mul_pd(_mm256_add_pd(lj12, lj12m6), invDistSquared));
        const __m256d scalarMultiple =
            _mm256_and_pd(_mm256_castsi256_pd(remainderMask), scalar);  // Idk why this is needed


        // calculate forces
        const __m256d forceX = _mm256_mul_pd(scalarMultiple, displacementX);
        const __m256d forceY = _mm256_mul_pd(scalarMultiple, displacementY);
        const __m256d forceZ = _mm256_mul_pd(scalarMultiple, displacementZ);

        const __m256d torqueAX =
            _mm256_fmsub_pd(rotatedPrimeSitePositionY, forceZ, _mm256_mul_pd(rotatedPrimeSitePositionZ, forceY));
        const __m256d torqueAY =
            _mm256_fmsub_pd(rotatedPrimeSitePositionZ, forceX, _mm256_mul_pd(rotatedPrimeSitePositionX, forceZ));
        const __m256d torqueAZ =
            _mm256_fmsub_pd(rotatedPrimeSitePositionX, forceY, _mm256_mul_pd(rotatedPrimeSitePositionY, forceX));

        forceSumX = _mm256_add_pd(forceSumX, forceX);
        forceSumY = _mm256_add_pd(forceSumY, forceY);
        forceSumZ = _mm256_add_pd(forceSumZ, forceZ);

        torqueSumX = _mm256_add_pd(torqueSumX, torqueAX);
        torqueSumY = _mm256_add_pd(torqueSumY, torqueAY);
        torqueSumZ = _mm256_add_pd(torqueSumZ, torqueAZ);

        // N3L
        if constexpr (newton3) {
          __m256d forceSumBX = gather_pd(remainderCase, siteForceX.data(), siteIndices, remainderMask);
          __m256d forceSumBY = gather_pd(remainderCase, siteForceY.data(), siteIndices, remainderMask);
          __m256d forceSumBZ = gather_pd(remainderCase, siteForceZ.data(), siteIndices, remainderMask);
          forceSumBX = _mm256_sub_pd(forceSumBX, forceX);
          forceSumBY = _mm256_sub_pd(forceSumBY, forceY);
          forceSumBZ = _mm256_sub_pd(forceSumBZ, forceZ);

          scatter_pd(remainderCase, siteForceX.data(), siteIndices, remainderMask, forceSumBX);
          scatter_pd(remainderCase, siteForceY.data(), siteIndices, remainderMask, forceSumBY);
          scatter_pd(remainderCase, siteForceZ.data(), siteIndices, remainderMask, forceSumBZ);
        }

        // calculate globals
        if constexpr (calculateGlobals) {
          const __m256d virialX = _mm256_mul_pd(displacementX, forceX);
          const __m256d virialY = _mm256_mul_pd(displacementY, forceY);
          const __m256d virialZ = _mm256_mul_pd(displacementZ, forceZ);

          const __m256d potentialEnergy6 =
              _mm256_fmadd_pd(epsilon24, lj12m6, shift6);  // FMA may not be supported on all CPUs
          const __m256d potentialEnergy6Masked = _mm256_and_pd(remainderMask, potentialEnergy6);

          __m256i ownedStateP4 = _mm256_set1_epi64x(static_cast<int64_t>(ownedStatePrime));
          __m256d ownedMaskPrime = _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateP4),
                                             _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
          __m256d energyFactor = _mm256_blendv_pd(_zero, _one, ownedMaskPrime);
          if constexpr (newton3) {
            const __m256i ownedStateB =
                gather_epi64(remainderCase, isNeighborSiteOwnedArr.data(), siteIndices, remainderMask);
            __m256d ownedMaskB = _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateB),
                                               _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
            energyFactor = _mm256_add_pd(energyFactor, _mm256_blendv_pd(_zero, _one, ownedMaskB));
          }
          potentialEnergySum = _mm256_fmadd_pd(energyFactor, potentialEnergy6Masked, potentialEnergySum);
          virialSumX = _mm256_fmadd_pd(energyFactor, virialX, virialSumX);
          virialSumY = _mm256_fmadd_pd(energyFactor, virialY, virialSumY);
          virialSumZ = _mm256_fmadd_pd(energyFactor, virialZ, virialSumZ);
        }
      }
    }
    // Add forces to prime mol
    fxptr[indexFirst] += horizontalSum(forceSumX);
    fyptr[indexFirst] += horizontalSum(forceSumY);
    fzptr[indexFirst] += horizontalSum(forceSumZ);
    txptr[indexFirst] += horizontalSum(torqueSumX);
    typtr[indexFirst] += horizontalSum(torqueSumY);
    tzptr[indexFirst] += horizontalSum(torqueSumZ);

    // Reduce forces on individual neighbor sites to molecular forces & torques if newton3=true
    if constexpr (newton3) {
      if constexpr (useMixing) {
        size_t siteIndex = 0;
        for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
          const auto neighborMolIndex = neighborList[neighborMol];
          rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
              {q0ptr[neighborMolIndex], q1ptr[neighborMolIndex], q2ptr[neighborMolIndex], q3ptr[neighborMolIndex]},
              _PPLibrary->getSitePositions(typeptr[neighborMolIndex]));
          for (size_t site = 0; site < _PPLibrary->getNumSites(typeptr[neighborMolIndex]); ++site) {
            fxptr[neighborMolIndex] += siteForceX[siteIndex];
            fyptr[neighborMolIndex] += siteForceY[siteIndex];
            fzptr[neighborMolIndex] += siteForceZ[siteIndex];
            txptr[neighborMolIndex] += rotatedSitePositions[site][1] * siteForceZ[siteIndex] -
                                       rotatedSitePositions[site][2] * siteForceY[siteIndex];
            typtr[neighborMolIndex] += rotatedSitePositions[site][2] * siteForceX[siteIndex] -
                                       rotatedSitePositions[site][0] * siteForceZ[siteIndex];
            tzptr[neighborMolIndex] += rotatedSitePositions[site][0] * siteForceY[siteIndex] -
                                       rotatedSitePositions[site][1] * siteForceX[siteIndex];
            ++siteIndex;
          }
        }
      } else {
        size_t siteIndex = 0;
        for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
          const auto neighborMolIndex = neighborList[neighborMol];
          rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
              {q0ptr[neighborMolIndex], q1ptr[neighborMolIndex], q2ptr[neighborMolIndex], q3ptr[neighborMolIndex]},
              const_unrotatedSitePositions);
          for (size_t site = 0; site < const_unrotatedSitePositions.size(); ++site) {
            fxptr[neighborMolIndex] += siteForceX[siteIndex];
            fyptr[neighborMolIndex] += siteForceY[siteIndex];
            fzptr[neighborMolIndex] += siteForceZ[siteIndex];
            txptr[neighborMolIndex] += rotatedSitePositions[site][1] * siteForceZ[siteIndex] -
                                       rotatedSitePositions[site][2] * siteForceY[siteIndex];
            typtr[neighborMolIndex] += rotatedSitePositions[site][2] * siteForceX[siteIndex] -
                                       rotatedSitePositions[site][0] * siteForceZ[siteIndex];
            txptr[neighborMolIndex] += rotatedSitePositions[site][0] * siteForceY[siteIndex] -
                                       rotatedSitePositions[site][1] * siteForceX[siteIndex];
            ++siteIndex;
          }
        }
      }
    }

    if constexpr (calculateGlobals) {
      const auto threadNum = autopas_get_thread_num();
      // SoAFunctorSingle obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
      // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.
      const auto newton3Factor = newton3 ? .5 : 1.;

      _aosThreadData[threadNum].potentialEnergySum += horizontalSum(potentialEnergySum) * newton3Factor;
      _aosThreadData[threadNum].virialSum[0] += horizontalSum(virialSumX) * newton3Factor;
      _aosThreadData[threadNum].virialSum[1] += horizontalSum(virialSumY) * newton3Factor;
      _aosThreadData[threadNum].virialSum[2] += horizontalSum(virialSumZ) * newton3Factor;
    }
  }

 private:
  /**
   * @brief Utility function to calculate the sum of an AVX register horizontally.
   * @param data register of 4 doubles
   * @return sum of its elements
   */
  inline double horizontalSum(const __m256d &data) {
    __m256d sum = _mm256_hadd_pd(data, data);
    __m256d permuted = _mm256_permute4x64_pd(sum, _MM_PERM_ACBD);
    sum = _mm256_hadd_pd(sum, permuted);
    sum = _mm256_permute4x64_pd(sum, _MM_PERM_BBDD);
    return _mm_cvtsd_f64(_mm256_extractf128_pd(sum, 0));
  }

  static inline __m256d load_pd(bool useMask, const double *const __restrict data, const __m256i &mask) {
    if (useMask) {
      return _mm256_maskload_pd(data, mask);
    } else {
      return _mm256_loadu_pd(data);
    }
  };

  static inline __m256d load_epi64(bool useMask, const size_t *const __restrict data, const __m256i &mask) {
    if (useMask) {
      return _mm256_maskload_epi64(reinterpret_cast<const long long int *>(data), mask);
    } else {
      return _mm256_loadu_si256(reinterpret_cast<const __m256i *>(data));
    }
  }

  static inline __m256d gather_pd(bool useMask, const double *const __restrict data, const __m256i &indices,
                                  const __m256i &mask) {
    if (useMask) {
      return _mm256_mask_i64gather_pd(_mm256_setzero_pd(), data, indices, mask, 8);
    } else {
      return _mm256_i64gather_pd(data, indices, 8);
    }
  }

  static inline __m256i gather_epi64(bool useMask, const size_t *const __restrict data, const __m256i &indices,
                                     const __m256i &mask) {
    if (useMask) {
      return _mm256_mask_i64gather_epi64(_mm256_setzero_si256(), data, indices, mask, 8);
    } else {
      return _mm256_i64gather_epi64(data, indices, 8);
    }
  }

  // Workaround for missing _mm256_i64scatter_pd :(
  static inline void scatter_pd(bool useMask, double *const __restrict ptr, const __m256i &indices, const __m256i &mask,
                                __m256d &values) {
    alignas(32) double temp[4];
    alignas(32) size_t idx[4];

    _mm256_store_pd(temp, values);
    _mm256_store_si256((__m256i *)idx, indices);

    for (size_t i = 0; i < 4; ++i) {
      if (mask[i]) {  // <- This does not work for all compilers >:(
        ptr[idx[i]] = temp[i];
      } else
        break;
    }
  }

 public:
  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
   * TODO: This is just a copy of the AoSThreadData from the LJFunctor. This may need refactoring.
   */
  class AoSThreadData {
   public:
    AoSThreadData() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, __remainingTo64{} {}

    void setZero() {
      virialSum = {0., 0., 0.};
      potentialEnergySum = 0;
    }

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

}  // namespace autopas
