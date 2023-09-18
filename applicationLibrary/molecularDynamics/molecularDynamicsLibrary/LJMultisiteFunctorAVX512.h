/**
 * @file LJMultisiteFunctorAVX512.h
 * @date 15/09/2023
 * @author Q. Behrami
 */
#pragma once

#ifndef __AVX512F__
#pragma message "LJMultisiteFunctorAVX.h included, but AVX is not supported by the compiler."
#else
#include "immintrin.h"
#endif

#include <array>

#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/ExceptionHandler.h"
#include "MultisiteMoleculeLJ.h"
#include "ParticlePropertiesLibrary.h"
#include "autopas/utils/AVXUtils.h"
#include "autopas/utils/Quaternion.h"
#include "autopas/utils/WrapOpenMP.h"

namespace mdLib {

/**
 * A functor to handle Lennard-Jones interactions between two Multisite Molecules.
 * This functor utilizes AVX512 instructions to speed up the computation.
 *
 * @tparam Particle The type of particle.
 * @tparam applyShift Flag for the LJ potential to have a truncated shift.
 * @tparam useMixing Flag for if the functor is to be used with multiple particle types. If set to false, _epsilon and
 * _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off, or both. See FunctorN3Nodes for possible
 * values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 * @tparam useMasks if true, use a 0/1 mask (i.e. calculate all forces for all molecules considered then multiply by
 * 0 if molecules are beyond the cutoff and 1 if molecules are within cutoff). If false, use a Gather/Scatter approach
 * i.e. calculate forces only for
 */
template <class Particle, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool relevantForTuning = true, bool useMasks = true, size_t vecLength = 8>
class LJMultisiteFunctorAVX512
    : public autopas::Functor<Particle, LJMultisiteFunctorAVX512<Particle, applyShift, useMixing, useNewton3,
                                                              calculateGlobals, relevantForTuning, useMasks, vecLength>> {
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

  /**
   * @brief Whether to use site masks or site indices
   *
   * @details If set to true, the sitemask contains either 0 or 1 for all relevant sites to indicate if a site meets the
   * cutoff and ownership condition. All sites are load contiguously, their respective force is calculated, and a
   * bitwise AND-operation is applied in the end. This method may include lots of unnecessary force calculations, but
   * most memory operations are contiguous. It may lead to better cache performance due to improved spatial locality.
   *
   * @details If set to false, the sitemask contains the indices of all sites that meet the necessary criteria. Since
   * these indices are not necessarily sequential, the force calculation requires gather intrinsics instead of regular
   * load operations. This method only calculates forces for owned sites within the cutoff radius but uses less
   * efficient memory operations. It may result in poorer cache performance due to potentially scattered memory access
   * patterns.
   */
//  constexpr static bool useRegularSiteMasks = false;

  /**
   * @brief How to evaluate the cutoff condition
   *
   *  @details If set to true, the cutoff condition is evaluated between the centers of masses between the particles
   * If set to false, the cutoff conditions is evaluated between the center of mass of the first particle and the site
   * positions of the second particle
   */
//  constexpr static bool useCTC = true;

#ifdef __AVX512F__
  const __m512d _cutoffSquared{};
  const __m512d _zero{_mm512_set1_pd(0.)};
  const __m512d _one{_mm512_set1_pd(1.)};
  const __m512i _three{_mm512_set1_epi64(3)}; // todo remove this
  /**
   * Masks for the remainder cases
   */
  std::array<__mmask8,8> _remainderMasks{
//      __mmask8(255), // = 11111111
//      __mmask8(254), // = 11111110
//      __mmask8(252), // = 11111100
//      __mmask8(248), // = 11111000
//      __mmask8(240), // = 11110000
//      __mmask8(224), // = 11100000
//      __mmask8(192), // = 11000000
//      __mmask8(128), // = 10000000
      __mmask8(255), // = 11111111
      __mmask8(127), // = 01111111
      __mmask8(63),  // = 00111111
      __mmask8(31),  // = 00011111
      __mmask8(15),  // = 00001111
      __mmask8(7),   // = 00000111
      __mmask8(3),   // = 00000011
      __mmask8(1),   // = 00000001
  };
  const __m256i _ownedStateOwnedMM256i{_mm256_set1_epi64x(static_cast<int64_t>(autopas::OwnershipState::owned))};
#endif

  // Vectors used to speed up the global verlet functor
  std::vector<double, autopas::AlignedAllocator<double>> _exactSitePositionsX;
  std::vector<double, autopas::AlignedAllocator<double>> _exactSitePositionsY;
  std::vector<double, autopas::AlignedAllocator<double>> _exactSitePositionsZ;

  std::vector<size_t, autopas::AlignedAllocator<size_t>> _siteOwnershipVector;
  std::vector<size_t, autopas::AlignedAllocator<size_t>> _siteTypesVector;

  std::vector<size_t, autopas::AlignedAllocator<size_t>> _molToSiteMap;
  std::vector<size_t, autopas::AlignedAllocator<size_t>> _siteToMolMap;

 public:
  /**
   * Deleted default constructor
   */
  LJMultisiteFunctorAVX512() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy unused, only there to make the signature different from the public constructor.
   */
  explicit LJMultisiteFunctorAVX512(double cutoff, void * /*dummy*/)
#ifdef __AVX512F__
      : autopas::Functor<Particle, LJMultisiteFunctorAVX512<Particle, applyShift, useMixing, useNewton3, calculateGlobals,
                                                relevantForTuning, useMasks, vecLength>>(cutoff),
        _cutoffSquared{_mm512_set1_pd(cutoff * cutoff)},
        _cutoffSquaredAoS(cutoff * cutoff),
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false} {
    if (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
  }
#else
      : autopas::Functor<Particle, LJMultisiteFunctorAVX512<Particle, applyShift, useMixing, useNewton3, calculateGlobals,
                                                relevantForTuning, useMasks, vecLength>>(cutoff) {
    autopas::utils::ExceptionHandler::exception("AutoPas was compiled without AVX support!");
  }
#endif

 public:
  /**
   * Constructor for Functor with particle mixing disabled. setParticleProperties() must be called.
   * @note Only to be used with mixing == false
   * @param cutoff
   */
  explicit LJMultisiteFunctorAVX512(double cutoff) : LJMultisiteFunctorAVX512(cutoff, nullptr) {
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
  explicit LJMultisiteFunctorAVX512(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJMultisiteFunctorAVX512(cutoff, nullptr) {
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
        autopas::utils::quaternion::rotateVectorOfPositions(particleA.getQuaternion(), unrotatedSitePositionsA);
    const auto rotatedSitePositionsB =
        autopas::utils::quaternion::rotateVectorOfPositions(particleB.getQuaternion(), unrotatedSitePositionsB);

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
          const auto virial = newton3 ? autopas::utils::ArrayMath::mulScalar(autopas::utils::ArrayMath::mul(displacement, force), 0.5)
                                      : autopas::utils::ArrayMath::mul(displacement, force);

          const auto threadNum = autopas::autopas_get_thread_num();

          if (particleA.isOwned()) {
            _aosThreadData[threadNum].potentialEnergySum += potentialEnergy6;
            _aosThreadData[threadNum].virialSum = autopas::utils::ArrayMath::add(_aosThreadData[threadNum].virialSum, virial);
          }
          if (newton3 and particleB.isOwned()) {
            _aosThreadData[threadNum].potentialEnergySum += potentialEnergy6;
            _aosThreadData[threadNum].virialSum = autopas::utils::ArrayMath::add(_aosThreadData[threadNum].virialSum, virial);
          }
        }
      }
    }
  }

  /**
   * @copydoc Functor::SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3)
   * @note This functor ignores the newton3 value, as we do not expect any benefit from disabling newton3.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    if (newton3) {
      SoAFunctorSingleImpl<true>(soa);
    } else {
      SoAFunctorSingleImpl<false>(soa);
    }
  }

  /**
   * @copydoc Functor::SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3)
   */
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2, const bool newton3) final {
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
    // Do nothing Todo
  }

 private:
  /**
   * SoAFunctorSingle Implementation.
   */
  template <bool newton3>
  void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {
#ifndef __AVX512F__
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

    // count number of sites in SoA
    size_t siteCount = 0;
    if constexpr (useMixing) {
      for (size_t mol = 0; mol < soa.getNumberOfParticles(); ++mol) {
        siteCount += _PPLibrary->getNumSites(typeptr[mol]);
      }
    } else {
      siteCount = _sitePositionsLJ.size() * soa.getNumberOfParticles();
    }

    // Accumulators for global values
    double potentialEnergyAccumulator = 0;
    std::array<double, 3> virialAccumulator = {0., 0., 0.};

    // ------------------------------ Setup auxiliary vectors -----------------------------
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedSitePositionX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedSitePositionY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedSitePositionZ;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionZ;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceX(siteCount, 0);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceY(siteCount, 0);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceZ(siteCount, 0);

    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteTypes;
    std::vector<size_t, autopas::AlignedAllocator<size_t>> isSiteOwned;

    // pre-reserve site std::vectors
    rotatedSitePositionX.reserve(siteCount);
    rotatedSitePositionY.reserve(siteCount);
    rotatedSitePositionZ.reserve(siteCount);
    exactSitePositionX.reserve(siteCount);
    exactSitePositionY.reserve(siteCount);
    exactSitePositionZ.reserve(siteCount);

    if constexpr (useMixing) {
      siteTypes.reserve(siteCount);
    }

    isSiteOwned.reserve(siteCount);

    // Fill site-wise std::vectors for SIMD

    for (size_t mol = 0; mol < soa.getNumberOfParticles(); ++mol) {
      // todo fix this
      std::vector<std::array<double, 3>> rotatedSitePositions;
      if constexpr (useMixing) {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, _PPLibrary->getSitePositions(typeptr[mol]));
      } else {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, _sitePositionsLJ);
      }

      size_t localSiteCount = useMixing ? _PPLibrary->getNumSites(typeptr[mol]) : _sitePositionsLJ.size();

      for (size_t site = 0; site < localSiteCount; ++site) {
        rotatedSitePositionX.push_back(rotatedSitePositions[site][0]);
        rotatedSitePositionY.push_back(rotatedSitePositions[site][1]);
        rotatedSitePositionZ.push_back(rotatedSitePositions[site][2]);
        exactSitePositionX.push_back(rotatedSitePositions[site][0] + xptr[mol]);
        exactSitePositionY.push_back(rotatedSitePositions[site][1] + yptr[mol]);
        exactSitePositionZ.push_back(rotatedSitePositions[site][2] + zptr[mol]);
        isSiteOwned.push_back(ownedStatePtr[mol] == autopas::OwnershipState::owned);
        if constexpr (useMixing) {
          siteTypes.push_back(_PPLibrary->getSiteTypes(typeptr[mol])[site]);
        }
      }
    }

    // ------------------------------ Main force calculation loop -----------------------------
    auto force_start = std::chrono::high_resolution_clock::now();
    size_t siteIndexMolA = 0;  // index of first site in molA
    for (size_t molA = 0; molA < soa.getNumberOfParticles(); ++molA) {
      const size_t noSitesInMolA = useMixing ? _PPLibrary->getNumSites(typeptr[molA])
                                             : _sitePositionsLJ.size();  // Number of sites in molecule A

      const auto ownedStateA = ownedStatePtr[molA];
      if (ownedStateA == autopas::OwnershipState::dummy) {
        siteIndexMolA += noSitesInMolA;
        continue;
      }

      const size_t siteIndexMolB = siteIndexMolA + noSitesInMolA;  // index of first site in molB
      const size_t noSitesB = (siteCount - siteIndexMolB);         // Number of sites in molecules that A interacts with

      // ------------------- Build site mask -------------------


      std::array<double, 3> centerOfMass{xptr[molA], yptr[molA], zptr[molA]};

//
//      std::vector<size_t, autopas::AlignedAllocator<size_t>> indices(soa.getNumberOfParticles() - molA - 1);
//      std::iota(indices.begin(), indices.end(), molA + 1);

      const auto siteMask = buildSiteMask(xptr, yptr, zptr, typeptr, ownedStatePtr, centerOfMass, molA+1, soa.getNumberOfParticles(), noSitesB);


      // ------------------- Calculate Forces -------------------

      // These values are not needed for the calculation, because there is no torque calculation for this functor
      std::array<double, 3> torqueAccumulator = {0., 0., 0.};

      for (size_t siteA = siteIndexMolA; siteA < siteIndexMolB; ++siteA) {
        const size_t siteTypeA = siteTypes[siteA];
        const std::array<double, 3> rotatedSitePositionA = {rotatedSitePositionX[siteA], rotatedSitePositionY[siteA],
                                                            rotatedSitePositionZ[siteA]};
        const std::array<double, 3> exactSitePositionA = {exactSitePositionX[siteA], exactSitePositionY[siteA],
                                                          exactSitePositionZ[siteA]};

        std::array<double, 3> forceAccumulator = {0., 0., 0.};

        SoAKernelMask<true>(
            siteMask, siteTypes, exactSitePositionX, exactSitePositionY, exactSitePositionZ, siteForceX, siteForceY,
            siteForceZ, isSiteOwned, ownedStateA, siteTypeA, exactSitePositionA, rotatedSitePositionA, forceAccumulator,
            torqueAccumulator, potentialEnergyAccumulator, virialAccumulator, siteIndexMolB);

        // Add forces + torques to molA
        fxptr[molA] += forceAccumulator[0];
        fyptr[molA] += forceAccumulator[1];
        fzptr[molA] += forceAccumulator[2];

        txptr[molA] += torqueAccumulator[0];
        typtr[molA] += torqueAccumulator[1];
        tzptr[molA] += torqueAccumulator[2];
      }
      siteIndexMolA += noSitesInMolA;
    }

    // ------------------------------ Reduction -----------------------------
    // i.e. add the "newton3" site forces to add to the molecular force + torque counters

    if constexpr (useMixing) {
              size_t siteIndex = 0;
              for (size_t mol = 0; mol < soa.getNumberOfParticles(); ++mol) {
        for (size_t site = 0; site < _PPLibrary->getNumSites(typeptr[mol]); ++site) {
          fxptr[mol] += siteForceX[siteIndex];
          fyptr[mol] += siteForceY[siteIndex];
          fzptr[mol] += siteForceZ[siteIndex];
          txptr[mol] += rotatedSitePositionY[siteIndex] * siteForceZ[siteIndex] -
                        rotatedSitePositionZ[siteIndex] * siteForceY[siteIndex];
          typtr[mol] += rotatedSitePositionZ[siteIndex] * siteForceX[siteIndex] -
                        rotatedSitePositionX[siteIndex] * siteForceZ[siteIndex];
          tzptr[mol] += rotatedSitePositionX[siteIndex] * siteForceY[siteIndex] -
                        rotatedSitePositionY[siteIndex] * siteForceX[siteIndex];
          ++siteIndex;
        }
              }
    } else {
              size_t siteIndex = 0;
              for (size_t mol = 0; mol < soa.getNumberOfParticles(); mol++) {
        for (size_t site = 0; site < _sitePositionsLJ.size(); ++site) {
          fxptr[mol] += siteForceX[siteIndex];
          fyptr[mol] += siteForceY[siteIndex];
          fzptr[mol] += siteForceZ[siteIndex];
          txptr[mol] += rotatedSitePositionY[siteIndex] * siteForceZ[siteIndex] -
                        rotatedSitePositionZ[siteIndex] * siteForceY[siteIndex];
          typtr[mol] += rotatedSitePositionZ[siteIndex] * siteForceX[siteIndex] -
                        rotatedSitePositionX[siteIndex] * siteForceZ[siteIndex];
          tzptr[mol] += rotatedSitePositionX[siteIndex] * siteForceY[siteIndex] -
                        rotatedSitePositionY[siteIndex] * siteForceX[siteIndex];
          ++siteIndex;
        }
      }

    }

    if constexpr (calculateGlobals) {
      const auto threadNum = autopas::autopas_get_thread_num();
      // SoAFunctorSingle obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
      // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.
      const auto newton3Factor = newton3 ? .5 : 1.;

      _aosThreadData[threadNum].potentialEnergySum += potentialEnergyAccumulator * newton3Factor;
      _aosThreadData[threadNum].virialSum[0] += virialAccumulator[0] * newton3Factor;
      _aosThreadData[threadNum].virialSum[1] += virialAccumulator[1] * newton3Factor;
      _aosThreadData[threadNum].virialSum[2] += virialAccumulator[2] * newton3Factor;
    }
  }

  /**
   * Implementation of SoAFunctorPair(soa1, soa2, newton3)
   */
  template <bool newton3>
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soaA, autopas::SoAView<SoAArraysType> soaB) {
#ifndef __AVX512F__
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

    // count number of sites in SoA
    size_t siteCountB = 0;
    if constexpr (useMixing) {
      for (size_t mol = 0; mol < soaB.getNumberOfParticles(); ++mol) {
        siteCountB += _PPLibrary->getNumSites(typeptrB[mol]);
      }
    } else {
      siteCountB = _sitePositionsLJ.size() * soaB.getNumberOfParticles();
    }

    // Accumulators for global values
    double potentialEnergyAccumulator = 0;
    std::array<double, 3> virialAccumulator = {0., 0., 0.};

    // ------------------------------ Setup auxiliary vectors -----------------------------

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionBx;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionBy;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionBz;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceBx(siteCountB, 0.);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceBy(siteCountB, 0.);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceBz(siteCountB, 0.);

    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteTypesB;
    std::vector<size_t, autopas::AlignedAllocator<size_t>> isSiteOwnedBArr;

    // pre-reserve std::vectors
    exactSitePositionBx.reserve(siteCountB);
    exactSitePositionBy.reserve(siteCountB);
    exactSitePositionBz.reserve(siteCountB);

    if constexpr (useMixing) {
      siteTypesB.reserve(siteCountB);
    }

    isSiteOwnedBArr.reserve(siteCountB);

    // Fill site-wise std::vectors for SIMD
    std::vector<std::array<double, 3>> rotatedSitePositions;
    for (size_t mol = 0; mol < soaB.getNumberOfParticles(); ++mol) {
      if constexpr (useMixing) {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0Bptr[mol], q1Bptr[mol], q2Bptr[mol], q3Bptr[mol]}, _PPLibrary->getSitePositions(typeptrB[mol]));
      } else {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0Bptr[mol], q1Bptr[mol], q2Bptr[mol], q3Bptr[mol]}, _sitePositionsLJ);
      }

      size_t localSiteCount = useMixing ? _PPLibrary->getNumSites(typeptrB[mol]) : _sitePositionsLJ.size();
      for (size_t site = 0; site < localSiteCount; ++site) {
        exactSitePositionBx.push_back(rotatedSitePositions[site][0] + xBptr[mol]);
        exactSitePositionBy.push_back(rotatedSitePositions[site][1] + yBptr[mol]);
        exactSitePositionBz.push_back(rotatedSitePositions[site][2] + zBptr[mol]);

        isSiteOwnedBArr.push_back(ownedStatePtrB[mol] == autopas::OwnershipState::owned);

        if constexpr (useMixing) {
          siteTypesB.push_back(_PPLibrary->getSiteTypes(typeptrB[mol])[site]);
        }
      }
    }

    // ------------------------------ Main force calculation loop -----------------------------

    for (size_t molA = 0; molA < soaA.getNumberOfParticles(); ++molA) {
      const auto ownedStateA = ownedStatePtrA[molA];
      if (ownedStateA == autopas::OwnershipState::dummy) {
        continue;
      }

      const auto noSitesInMolA = useMixing ? _PPLibrary->getNumSites(typeptrA[molA]) : _sitePositionsLJ.size();
      const auto unrotatedSitePositionsA = useMixing ? _PPLibrary->getSitePositions(typeptrA[molA]) : _sitePositionsLJ;

      const auto rotatedSitePositionsA = autopas::utils::quaternion::rotateVectorOfPositions(
          {q0Aptr[molA], q1Aptr[molA], q2Aptr[molA], q3Aptr[molA]}, unrotatedSitePositionsA);

      // ------------------- Build site mask -------------------

      std::array<double, 3> centerOfMass{xAptr[molA], yAptr[molA], zAptr[molA]};

      std::vector<size_t, autopas::AlignedAllocator<size_t>> indices(soaB.getNumberOfParticles());
      std::iota(indices.begin(), indices.end(), 0);

      const auto siteMask =
          buildSiteMask(xBptr, yBptr, zBptr, typeptrB, ownedStatePtrB, centerOfMass, 0, soaB.getNumberOfParticles(), siteCountB);


      // ------------------- Calculate Forces -------------------

      std::array<double, 3> forceAccumulator = {0., 0., 0.};
      std::array<double, 3> torqueAccumulator = {0., 0., 0.};

      for (size_t siteA = 0; siteA < noSitesInMolA; ++siteA) {
        const size_t siteTypeA = _PPLibrary->getSiteTypes(typeptrA[molA])[siteA];
        const std::array<double, 3> rotatedSitePositionA = rotatedSitePositionsA[siteA];
        const std::array<double, 3> exactSitePositionA =
            autopas::utils::ArrayMath::add(rotatedSitePositionA, centerOfMass);

        SoAKernelMask<newton3>(
            siteMask, siteTypesB, exactSitePositionBx, exactSitePositionBy, exactSitePositionBz, siteForceBx,
            siteForceBy, siteForceBz, isSiteOwnedBArr, ownedStateA, siteTypeA, exactSitePositionA, rotatedSitePositionA,
            forceAccumulator, torqueAccumulator, potentialEnergyAccumulator, virialAccumulator);
      }
      // Add forces and torques to mol A
      fxAptr[molA] += forceAccumulator[0];
      fyAptr[molA] += forceAccumulator[1];
      fzAptr[molA] += forceAccumulator[2];
      txAptr[molA] += torqueAccumulator[0];
      tyAptr[molA] += torqueAccumulator[1];
      tzAptr[molA] += torqueAccumulator[2];
    }

    // ------------------------------ Reduction -----------------------------

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
            {q0Bptr[mol], q1Bptr[mol], q2Bptr[mol], q3Bptr[mol]}, _sitePositionsLJ);
        for (size_t site = 0; site < _sitePositionsLJ.size(); ++site) {
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
      const auto threadNum = autopas::autopas_get_thread_num();
      // SoAFunctorPairImpl obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
      // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.
      const auto newton3Factor = newton3 ? .5 : 1.;

      _aosThreadData[threadNum].potentialEnergySum += potentialEnergyAccumulator * newton3Factor;
      _aosThreadData[threadNum].virialSum[0] += virialAccumulator[0] * newton3Factor;
      _aosThreadData[threadNum].virialSum[1] += virialAccumulator[1] * newton3Factor;
      _aosThreadData[threadNum].virialSum[2] += virialAccumulator[2] * newton3Factor;
    }
  }


  /**
   * @brief Force kernel for the SoA version of the LJFunctor using masks.
   * @tparam newton3 Whether Newton3 optimizations is used.
   * @tparam calculateTorque Whether torque is calculated.
   * @param siteMask The mask used to determine which sites are calculated.
   * @param siteTypesB A vector containing the types of the sites of the neighbor particles
   * @param exactSitePositionX The exact x coordinates of the sites of the neighbor particles
   * @param exactSitePositionY The exact y coordinates of the sites of the neighbor particles
   * @param exactSitePositionZ The exact z coordinates of the sites of the neighbor particles
   * @param siteForceX The site forces in x direction of the neighbor particles (used for newton3)
   * @param siteForceY The site forces in y direction of the neighbor particles (used for newton3)
   * @param siteForceZ The site forces in z direction of the neighbor particles (used for newton3)
   * @param siteOwnership The ownership state of the neighbor particles
   * @param ownedStateA The ownership state of the primary particle
   * @param siteTypeA The site type of the primary particle
   * @param exactSitePositionA The exact site positions of the primary particle
   * @param rotatedSitePositionA The rotated site positions of the primary particle
   * @param forceAccumulator The force accumulator for the primary particle
   * @param torqueAccumulator The torque accumulator for the primary particle
   * @param potentialEnergyAccumulator The potential energy accumulator for all particles
   * @param virialAccumulator The virial accumulator for all particles
   * @param offset An optional offset for neighbor site vectors (only used by SoAFunctorSingle)
   */
  template <bool newton3>
  inline void SoAKernelMask(const std::vector<__mmask8, autopas::AlignedAllocator<__mmask8>> &siteMaskVec,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &siteTypesB,
                        const std::vector<double, autopas::AlignedAllocator<double>> &exactSitePositionX,
                        const std::vector<double, autopas::AlignedAllocator<double>> &exactSitePositionY,
                        const std::vector<double, autopas::AlignedAllocator<double>> &exactSitePositionZ,
                        std::vector<double, autopas::AlignedAllocator<double>> &siteForceX,
                        std::vector<double, autopas::AlignedAllocator<double>> &siteForceY,
                        std::vector<double, autopas::AlignedAllocator<double>> &siteForceZ,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &siteOwnership,
                        const autopas::OwnershipState ownedStateA, const size_t siteTypeA,
                        const std::array<double, 3> exactSitePositionA,
                        const std::array<double, 3> rotatedSitePositionA, std::array<double, 3> &forceAccumulator,
                        std::array<double, 3> &torqueAccumulator, double &potentialEnergyAccumulator,
                        std::array<double, 3> &virialAccumulator, size_t offset = 0) {
#ifndef __AVX512F__
#pragma message "LJMultisiteFunctorAVX.h included, but AVX is not supported by the compiler."
#else
    // Declare mixing variables
    __m512d sigmaSquared = _mm512_set1_pd(_sigmaSquaredAoS);
    __m512d epsilon24 = _mm512_set1_pd(_epsilon24AoS);
    __m512d shift6 = applyShift ? _mm512_set1_pd(_shift6AoS) : _zero;

    // Broadcast exact site positions of particle A
    const __m512d exactSitePositionsAX = _mm512_set1_pd(exactSitePositionA[0]);
    const __m512d exactSitePositionsAY = _mm512_set1_pd(exactSitePositionA[1]);
    const __m512d exactSitePositionsAZ = _mm512_set1_pd(exactSitePositionA[2]);

    // Broadcast rotated site positions of particle A
    const __m512d rotatedSitePositionsAX = _mm512_set1_pd(rotatedSitePositionA[0]);
    const __m512d rotatedSitePositionsAY = _mm512_set1_pd(rotatedSitePositionA[1]);
    const __m512d rotatedSitePositionsAZ = _mm512_set1_pd(rotatedSitePositionA[2]);

    // sums used for siteA
    __m512d forceSumX = _zero;
    __m512d forceSumY = _zero;
    __m512d forceSumZ = _zero;

    __m512d torqueSumX = _zero;
    __m512d torqueSumY = _zero;
    __m512d torqueSumZ = _zero;

    // Globals
    __m512d potentialEnergySum = _zero;
    __m512d virialSumX = _zero;
    __m512d virialSumY = _zero;
    __m512d virialSumZ = _zero;

    // This variable either contains the local mask or the indices of the sites, depending on the template parameter
    __m512i sites = _mm512_setzero_si512();

    for (size_t siteVectorIndex = offset; siteVectorIndex < siteTypesB.size(); siteVectorIndex += vecLength) {
      const size_t remainder = siteTypesB.size() - siteVectorIndex;
      const bool remainderCase = remainder < vecLength;
      const __mmask8 remainderMask = _remainderMasks[8-remainder];
      const __mmask8 siteMask = siteMaskVec[(siteVectorIndex-offset) / vecLength];

      // Load the mixing parameters
      if constexpr (useMixing) {
        const double *const __restrict mixingPtr = useMixing ? _PPLibrary->getMixingDataPtr(siteTypeA, 0) : nullptr;


        const __m512i siteTypeIndices = remainderCase ? _mm512_maskz_loadu_epi64(remainderMask, &siteTypesB[siteVectorIndex])
            : _mm512_loadu_epi64(&siteTypesB[siteVectorIndex]);
//        const __m512i siteTypeIndicesOffset = _mm512_mullox_epi64(siteTypeIndices, _mm256_set1_epi64x(3));  // This could maybe be problematic!

        // todo this is really silly - just generate different pointers for different mixing data types
        const __m512i siteTypeIndicesScaled = _mm512_mullox_epi64(siteTypeIndices, _three);


        // I think this is correct but potentially not
        if (remainderCase) {
          epsilon24 = _mm512_mask_i64gather_pd(_zero, remainderMask, siteTypeIndicesScaled, mixingPtr, 8);
          sigmaSquared =
              _mm512_mask_i64gather_pd(_one, remainderMask, siteTypeIndicesScaled, mixingPtr + 1, 8);  // one used as "filler" in remainder case to avoid dividing by zero
          if constexpr (applyShift) {
            shift6 = _mm512_mask_i64gather_pd(_zero, remainderMask, siteTypeIndicesScaled, mixingPtr + 2, 8);
          }
        } else {
          epsilon24 = _mm512_i64gather_pd(siteTypeIndicesScaled, mixingPtr, 8);
          sigmaSquared = _mm512_i64gather_pd(siteTypeIndicesScaled, mixingPtr + 1, 8);
          if constexpr (applyShift) {
            shift6 = _mm512_i64gather_pd(siteTypeIndicesScaled, mixingPtr + 2, 8);
          }
        }
      }

      // Load the exact site positions of particle B // todo is it faster to just use _mm512_maskz_load_pd?
      const __m512d exactSitePositionsBX = remainderCase ? _mm512_maskz_loadu_pd(remainderMask, &exactSitePositionX[siteVectorIndex])
          : _mm512_loadu_pd(&exactSitePositionX[siteVectorIndex]);
      const __m512d exactSitePositionsBY = remainderCase ? _mm512_maskz_loadu_pd(remainderMask, &exactSitePositionY[siteVectorIndex])
                                                         : _mm512_loadu_pd(&exactSitePositionY[siteVectorIndex]);
      const __m512d exactSitePositionsBZ = remainderCase ? _mm512_maskz_loadu_pd(remainderMask, &exactSitePositionZ[siteVectorIndex])
                                                         : _mm512_loadu_pd(&exactSitePositionZ[siteVectorIndex]);

      // Calculate the Lennard-Jones 12-6 potential

      const __m512d displacementX = _mm512_sub_pd(exactSitePositionsAX, exactSitePositionsBX);
      const __m512d displacementY = _mm512_sub_pd(exactSitePositionsAY, exactSitePositionsBY);
      const __m512d displacementZ = _mm512_sub_pd(exactSitePositionsAZ, exactSitePositionsBZ);

      const __m512d distanceSquaredX = _mm512_mul_pd(displacementX, displacementX);
      const __m512d distanceSquaredY = _mm512_mul_pd(displacementY, displacementY);
      const __m512d distanceSquaredZ = _mm512_mul_pd(displacementZ, displacementZ);

      const __m512d distanceSquared =
          _mm512_add_pd(distanceSquaredX, _mm512_add_pd(distanceSquaredY, distanceSquaredZ));

      const __m512d invDistSquared = _mm512_div_pd(_one, distanceSquared);
      const __m512d lj2 = _mm512_mul_pd(sigmaSquared, invDistSquared);
      const __m512d lj6 = _mm512_mul_pd(_mm512_mul_pd(lj2, lj2), lj2);
      const __m512d lj12 = _mm512_mul_pd(lj6, lj6);
      const __m512d lj12m6 = _mm512_sub_pd(lj12, lj6);

      const __m512d scalar = _mm512_mul_pd(epsilon24, _mm512_mul_pd(_mm512_add_pd(lj12, lj12m6), invDistSquared));

      // Determine forces and mask out sites beyond cutoff
      const __m512d forceX = _mm512_maskz_mul_pd((siteMask), scalar, displacementX);
      const __m512d forceY = _mm512_maskz_mul_pd((siteMask), scalar, displacementY);
      const __m512d forceZ = _mm512_maskz_mul_pd((siteMask), scalar, displacementZ);

      if (remainderCase) {
        forceSumX = _mm512_mask_add_pd(forceSumX, remainderMask, forceSumX, forceX);
        forceSumY = _mm512_mask_add_pd(forceSumY, remainderMask, forceSumY, forceY);
        forceSumZ = _mm512_mask_add_pd(forceSumZ, remainderMask, forceSumZ, forceZ);
      } else {
        forceSumX = _mm512_add_pd( forceSumX, forceX);
        forceSumY = _mm512_add_pd( forceSumY, forceY);
        forceSumZ = _mm512_add_pd( forceSumZ, forceZ);
      }

      // calculate torques

        const __m512d torqueAX =
            _mm512_fmsub_pd(rotatedSitePositionsAY, forceZ, _mm512_mul_pd(rotatedSitePositionsAZ, forceY));
        const __m512d torqueAY =
            _mm512_fmsub_pd(rotatedSitePositionsAZ, forceX, _mm512_mul_pd(rotatedSitePositionsAX, forceZ));
        const __m512d torqueAZ =
            _mm512_fmsub_pd(rotatedSitePositionsAX, forceY, _mm512_mul_pd(rotatedSitePositionsAY, forceX));

        if (remainderCase) {
          torqueSumX = _mm512_mask_add_pd(torqueSumX, remainderMask, torqueSumX, torqueAX);
          torqueSumY = _mm512_mask_add_pd(torqueSumY, remainderMask, torqueSumY, torqueAY);
          torqueSumZ = _mm512_mask_add_pd(torqueSumZ, remainderMask, torqueSumZ, torqueAZ);
        } else {
          torqueSumX = _mm512_add_pd(torqueSumX, torqueAX);
          torqueSumY = _mm512_add_pd(torqueSumY, torqueAY);
          torqueSumZ = _mm512_add_pd(torqueSumZ, torqueAZ);
        }


      // Newton 3 optimization
      if constexpr (newton3) {
        if (remainderCase) {
          const __m512d forceSumBX = _mm512_maskz_loadu_pd(remainderMask, &siteForceX[siteVectorIndex]);
          const __m512d newForceSumBX = _mm512_sub_pd(forceSumBX, forceX);
          _mm512_mask_storeu_pd(&siteForceX[siteVectorIndex], remainderMask, newForceSumBX);

          const __m512d forceSumBY = _mm512_maskz_loadu_pd(remainderMask, &siteForceY[siteVectorIndex]);
          const __m512d newForceSumBY = _mm512_sub_pd(forceSumBY, forceY);
          _mm512_mask_storeu_pd(&siteForceY[siteVectorIndex], remainderMask, newForceSumBY);

          const __m512d forceSumBZ = _mm512_maskz_loadu_pd(remainderMask, &siteForceZ[siteVectorIndex]);
          const __m512d newForceSumBZ = _mm512_sub_pd(forceSumBZ, forceZ);
          _mm512_mask_storeu_pd(&siteForceZ[siteVectorIndex], remainderMask, newForceSumBZ);
        } else {
          const __m512d forceSumBX = _mm512_loadu_pd( &siteForceX[siteVectorIndex]);
          const __m512d newForceSumBX = _mm512_sub_pd(forceSumBX, forceX);
          _mm512_storeu_pd(&siteForceX[siteVectorIndex],  newForceSumBX);

          const __m512d forceSumBY = _mm512_loadu_pd( &siteForceY[siteVectorIndex]);
          const __m512d newForceSumBY = _mm512_sub_pd(forceSumBY, forceY);
          _mm512_storeu_pd(&siteForceY[siteVectorIndex], newForceSumBY);

          const __m512d forceSumBZ = _mm512_loadu_pd( &siteForceZ[siteVectorIndex]);
          const __m512d newForceSumBZ = _mm512_sub_pd(forceSumBZ, forceZ);
          _mm512_storeu_pd(&siteForceZ[siteVectorIndex], newForceSumBZ);
        }

      }

      // globals
      if constexpr (calculateGlobals) {
//        const __m256d virialX = _mm256_mul_pd(displacementX, forceX);
//        const __m256d virialY = _mm256_mul_pd(displacementY, forceY);
//        const __m256d virialZ = _mm256_mul_pd(displacementZ, forceZ);
//
//        const __m256d potentialEnergy6 =
//            _mm256_fmadd_pd(epsilon24, lj12m6, shift6);  // FMA may not be supported on all CPUs
//        __m256d potentialEnergy6Masked = masked ? _mm256_and_pd(castedMask, potentialEnergy6) : potentialEnergy6;
//        potentialEnergy6Masked = _mm256_and_pd(_mm256_castsi256_pd(remainderMask), potentialEnergy6Masked);
//
//        __m256i ownedStateA4 = _mm256_set1_epi64x(static_cast<int64_t>(ownedStateA));
//        __m256d ownedMaskA =
//            _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateA4), _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
//        __m256d energyFactor = _mm256_blendv_pd(_zero, _one, ownedMaskA);
//        if constexpr (newton3) {
//          const __m256i ownedStateB =
//              masked ? autopas::utils::avx::load_epi64(remainderCase, &siteOwnership[offset + siteVectorIndex],
//                                                       remainderMask)
//                     : autopas::utils::avx::gather_epi64(remainderCase, &siteOwnership[offset], sites, remainderMask);
//          __m256d ownedMaskB =
//              _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateB), _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
//          energyFactor = _mm256_add_pd(energyFactor, _mm256_blendv_pd(_zero, _one, ownedMaskB));
//        }
//        potentialEnergySum = _mm256_fmadd_pd(energyFactor, potentialEnergy6Masked, potentialEnergySum);
//        virialSumX = _mm256_fmadd_pd(energyFactor, virialX, virialSumX);
//        virialSumY = _mm256_fmadd_pd(energyFactor, virialY, virialSumY);
//        virialSumZ = _mm256_fmadd_pd(energyFactor, virialZ, virialSumZ);
      }
    }

    // Add up the forces, torques and globals

    forceAccumulator[0] += autopas::utils::avx::horizontalSum(forceSumX);
    forceAccumulator[1] += autopas::utils::avx::horizontalSum(forceSumY);
    forceAccumulator[2] += autopas::utils::avx::horizontalSum(forceSumZ);


    torqueAccumulator[0] += autopas::utils::avx::horizontalSum(torqueSumX);
    torqueAccumulator[1] += autopas::utils::avx::horizontalSum(torqueSumY);
    torqueAccumulator[2] += autopas::utils::avx::horizontalSum(torqueSumZ);


    if constexpr (calculateGlobals) {
      potentialEnergyAccumulator += autopas::utils::avx::horizontalSum(potentialEnergySum);
      virialAccumulator[0] += autopas::utils::avx::horizontalSum(virialSumX);
      virialAccumulator[1] += autopas::utils::avx::horizontalSum(virialSumY);
      virialAccumulator[2] += autopas::utils::avx::horizontalSum(virialSumZ);
    }
#endif
  }


  /**
   * build site mask for SoASingle & SoAPair
   */
  inline std::vector<__mmask8, autopas::AlignedAllocator<__mmask8>> buildSiteMask(
      const double *const __restrict xptr, const double *const __restrict yptr, const double *const __restrict zptr,
      const size_t *const __restrict typeptr, const autopas::OwnershipState *const __restrict ownedStatePtr,
      const std::array<double, 3> &centerOfMass, size_t offset, size_t noMolecules, size_t noSites) {
    std::vector<__mmask8, autopas::AlignedAllocator<__mmask8>> siteMask;

    const size_t sizeOfSiteMask = noSites / 8 + 1; // Every 8 molecules fit into 1 __mmask8 (+ one for remainders)
    siteMask.reserve(sizeOfSiteMask);

    size_t siteIndex = 0;

    __mmask8 maskTemp(0);
    size_t uniqueBitFlipInt = 1; // The integer value, which, when converted into binary, and "OR"ed with maskTemp, flips maskTemp[maskCounter] to 1.
    for (size_t molIndex = offset; molIndex < noMolecules; molIndex++) {
      const size_t siteCount = useMixing ? _PPLibrary->getNumSites(typeptr[molIndex]) : _sitePositionsLJ.size();

      const double xPosB = xptr[molIndex];
      const double yPosB = yptr[molIndex];
      const double zPosB = zptr[molIndex];

      // calculate displacement
      const double displacementCoMX = centerOfMass[0] - xPosB;
      const double displacementCoMY = centerOfMass[1] - yPosB;
      const double displacementCoMZ = centerOfMass[2] - zPosB;

      // calculate distance squared
      const double distanceSquaredCoMX = displacementCoMX * displacementCoMX;
      const double distanceSquaredCoMY = displacementCoMY * displacementCoMY;
      const double distanceSquaredCoMZ = displacementCoMZ * displacementCoMZ;

      const double distanceSquaredCoM = distanceSquaredCoMX + distanceSquaredCoMY + distanceSquaredCoMZ;

      const bool cutoffCondition = distanceSquaredCoM <= _cutoffSquaredAoS;
      const bool dummyCondition = ownedStatePtr[molIndex] != autopas::OwnershipState::dummy;
      const bool condition = cutoffCondition and dummyCondition;

      for (size_t site = 0; site < siteCount; site++) {
        if (condition) {
          maskTemp = _mm512_kor(maskTemp, __mmask8(uniqueBitFlipInt));
        }

        if (uniqueBitFlipInt >= 128) {
          siteMask.emplace_back(maskTemp);
          maskTemp = 0;          // reset maskTemp
          uniqueBitFlipInt = 1;  // reset such that first bit will be flipped next
        } else {
          uniqueBitFlipInt *= 2;  // Shift to the next bit to be modified.
        }
      }
    }
    // emplace back final mask
    siteMask.emplace_back(maskTemp);
    return siteMask;
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
   * @TODO This is maybe not accurate for the vectorized functors
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

 public:
  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
   * TODO: This is just a copy of the AoSThreadData from the LJFunctor. This maybe needs adjustments.
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
};  // namespace autopas


}  // namespace autopas