/**
 * @file LJMultisiteFunctorAVX512_GS.h
 * @date 15/09/2023
 * @author S. Newcome
 * @note Largely derived from LJMultisiteFunctorAVX.h by Q. Behrami which no longer exists.
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
#include "autopas/utils/Quaternion.h"
#include "autopas/utils/WrapOpenMP.h"

namespace mdLib {

/**
 * A functor to handle Lennard-Jones interactions between two Multisite Molecules.
 * This functor utilizes AVX512 instructions to speed up the computation.
 *
 * Cutoff condition is that if two molecule's center-of-masses are within the cutoff, the corresponding site-to-site
 * iteractions are computed. This requires non-continuous access to relevant site-parameter vectors, handled using
 * AVX512's Gather and Scatter operations.
 *
 * @warning no global calculation implemented. Reasoning: At the time of creation, the value of this implementation of
 * multisite molecules (i.e. with AutoPas particles = molecules as opposed to particles = sites) is under question. Calculation
 * of globals is not of priority until we can be sure of which implementation is better and why.
 *
 * @tparam Particle The type of particle.
 * @tparam applyShift Flag for the LJ potential to have a truncated shift.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off, or both. See FunctorN3Nodes for possible
 * values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial). No global calculation
 * implemented so throws exception if true.
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, bool applyShift = false, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool relevantForTuning = true>
class LJMultisiteFunctorAVX512_GS
    : public autopas::Functor<Particle, LJMultisiteFunctorAVX512_GS<Particle, applyShift, useNewton3,
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
   * Particle property library.
   */
  ParticlePropertiesLibrary<SoAFloatPrecision, size_t> *_PPLibrary = nullptr;

  /**
   * Sum of potential energy. Only calculated if calculateGlobals is true.
   */
  double _potentialEnergySum{0};

  /**
   * Sum of the virial. Only calculated if calculateGlobals is true.
   */
  std::array<double, 3> _virialSum{0., 0., 0.};

  /**
   * Defines whether or whether not the global values are already processed
   */
  bool _postProcessed;

#ifdef __AVX512F__
  const __m512d _cutoffSquared{};
  const __m512d _zero{_mm512_set1_pd(0.)};
  const __m512i _zeroI{_mm512_set1_epi64(0)};
  const __m512d _one{_mm512_set1_pd(1.)};
  const __m512i _three{_mm512_set1_epi64(3)};

  /**
   * Masks for the remainder cases to avoid loading data beyond the end of the vectors.
   * Masks are generated with decimal numbers, whose binary equivalent consists of 8 0s or 1s, representing which registers
   * are masked when that mask is used.
   */
  std::array<__mmask8,8> _remainderMasks{
      __mmask8(255), // = 11111111
      __mmask8(127), // = 01111111
      __mmask8(63),  // = 00111111
      __mmask8(31),  // = 00011111
      __mmask8(15),  // = 00001111
      __mmask8(7),   // = 00000111
      __mmask8(3),   // = 00000011
      __mmask8(1),   // = 00000001
  };
#endif


 public:
  /**
   * Deleted default constructor
   */
  LJMultisiteFunctorAVX512_GS() = delete;

  /**
   * Constructor for Functor with particle mixing enabled.
   * Calculating global attributes is done with CoM and overall forces applied
   * @param cutoff
   * @param particlePropertiesLibrary Library used to look up the properties of each type of particle e.g. sigma,
   * epsilon, shift, site positions, site types, moment of inertia.
   */
  LJMultisiteFunctorAVX512_GS(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
#ifdef __AVX512F__
      : autopas::Functor<
            Particle, LJMultisiteFunctorAVX512_GS<Particle, applyShift, useNewton3, calculateGlobals, relevantForTuning>>(
            cutoff),
        _cutoffSquared{_mm512_set1_pd(cutoff * cutoff)},
        _cutoffSquaredAoS(cutoff * cutoff),
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false},
        _PPLibrary{&particlePropertiesLibrary} {
    if (calculateGlobals) {
      autopas::utils::ExceptionHandler::exception("LJMultisiteFunctorAVX512_GS constructed with calculateGlobals=true, but global calculation has not been implemented!");
    }
  }
#else
      : autopas::Functor<Particle, LJMultisiteFunctorAVX512_GS<Particle, applyShift, useNewton3, calculateGlobals,
                                                            relevantForTuning>>(cutoff) {
    autopas::utils::ExceptionHandler::exception("AutoPas was compiled without AVX512 support!");
  }
#endif

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
    const size_t numSitesA = _PPLibrary->getNumSites(particleA.getTypeId());
    const size_t numSitesB = _PPLibrary->getNumSites(particleB.getTypeId());

    // get siteIds
    const std::vector<size_t> siteIdsA = _PPLibrary->getSiteTypes(particleA.getTypeId());
    const std::vector<size_t> siteIdsB = _PPLibrary->getSiteTypes(particleB.getTypeId());

    // get unrotated relative site positions
    const std::vector<std::array<double, 3>> unrotatedSitePositionsA = _PPLibrary->getSitePositions(particleA.getTypeId());
    const std::vector<std::array<double, 3>> unrotatedSitePositionsB = _PPLibrary->getSitePositions(particleB.getTypeId());

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

        const auto sigmaSquared = _PPLibrary->getMixingSigmaSquared(siteIdsA[i], siteIdsB[j]);
        const auto epsilon24 = _PPLibrary->getMixing24Epsilon(siteIdsA[i], siteIdsB[j]);
        const auto shift6 = applyShift ? _PPLibrary->getMixingShift6(siteIdsA[i], siteIdsB[j]) : 0;

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

        // todo handling of globals should go here (see LJMultisiteFunctor.h or 7dc729af)

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
    if (newton3) {
      SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
    }
  }

 private:
  /**
   * SoAFunctorSingle Implementation.
   */
  template <bool newton3>
  void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {
#ifndef __AVX512F__
#pragma message "LJMultisiteFunctorAVX512_GS::SoAFunctorSingleImpl called without AVX512 intrinsics!"
#endif
    if (soa.size() == 0) return;

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
    // todo this could be vectorized
    const size_t siteCount = [&]() {
      size_t siteCountTmp{0};
      for (size_t mol = 0; mol < soa.size(); ++mol) {
        siteCountTmp += _PPLibrary->getNumSites(typeptr[mol]);
      }
      return siteCountTmp;
    }();

    // ------------------------------ Setup auxiliary vectors -----------------------------
    // find out if we waste too much allocating and deallocating these vectors
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


    siteTypes.reserve(siteCount);


    isSiteOwned.reserve(siteCount);

    // Fill site-wise std::vectors for SIMD
    std::vector<std::array<double, 3>> rotatedSitePositions;

    for (size_t mol = 0; mol < soa.size(); ++mol) {
      rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
          {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, _PPLibrary->getSitePositions(typeptr[mol]));

      size_t localSiteCount = _PPLibrary->getNumSites(typeptr[mol]);

      for (size_t site = 0; site < localSiteCount; ++site) {
        rotatedSitePositionX.push_back(rotatedSitePositions[site][0]);
        rotatedSitePositionY.push_back(rotatedSitePositions[site][1]);
        rotatedSitePositionZ.push_back(rotatedSitePositions[site][2]);
        exactSitePositionX.push_back(rotatedSitePositions[site][0] + xptr[mol]);
        exactSitePositionY.push_back(rotatedSitePositions[site][1] + yptr[mol]);
        exactSitePositionZ.push_back(rotatedSitePositions[site][2] + zptr[mol]);
        isSiteOwned.push_back(ownedStatePtr[mol] == autopas::OwnershipState::owned);
        siteTypes.push_back(_PPLibrary->getSiteTypes(typeptr[mol])[site]);
      }
    }

    // ------------------------------ Main force calculation loop -----------------------------

    size_t siteIndexMolA = 0;  // index of site vectors that represent the first site in molA
    for (size_t molA = 0; molA < soa.size(); ++molA) {
      const size_t noSitesInMolA = _PPLibrary->getNumSites(typeptr[molA]);  // Number of sites in molecule A

      const auto ownedStateA = ownedStatePtr[molA];
      if (ownedStateA == autopas::OwnershipState::dummy) {
        siteIndexMolA += noSitesInMolA;
        continue;
      }

      const size_t siteIndexMolB = siteIndexMolA + noSitesInMolA;  // index of first site in molB
      const size_t noSitesB = (siteCount - siteIndexMolB);         // Number of sites in molecules that A interacts with

      const std::array<double, 3> centerOfMass{xptr[molA], yptr[molA], zptr[molA]};

      // Build list of indices
      const auto sitePairIndicies = buildSiteInteractionIndices(xptr, yptr, zptr, typeptr, ownedStatePtr, centerOfMass, molA+1, siteIndexMolB, soa.size(),
                                                         noSitesB);
      // Calculate Forces

      for (size_t siteA = siteIndexMolA; siteA < siteIndexMolB; ++siteA) {
        const size_t siteTypeA = siteTypes[siteA];
        const std::array<double, 3> rotatedSitePositionA = {rotatedSitePositionX[siteA], rotatedSitePositionY[siteA],
                                                            rotatedSitePositionZ[siteA]};
        const std::array<double, 3> exactSitePositionA = {exactSitePositionX[siteA], exactSitePositionY[siteA],
                                                          exactSitePositionZ[siteA]};

        std::array<double, 3> forceAccumulator = {0., 0., 0.};
        std::array<double, 3> torqueAccumulator = {0., 0., 0.};

        interactAllWithSiteA<true>(sitePairIndicies, siteTypes, exactSitePositionX, exactSitePositionY,
                                   exactSitePositionZ, siteForceX, siteForceY, siteForceZ, isSiteOwned, ownedStateA,
                                   siteTypeA, exactSitePositionA, rotatedSitePositionA, forceAccumulator,
                                   torqueAccumulator);

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


    size_t siteIndex = 0;
    for (size_t mol = 0; mol < soa.size(); ++mol) {
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


    // todo processing of global accumulators should go here (see LJMultisiteFunctor.h or 7dc729af)

  }

  /**
   * Implementation of SoAFunctorPair(soa1, soa2, newton3)
   */
  template <bool newton3>
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soaA, autopas::SoAView<SoAArraysType> soaB) {
#ifndef __AVX512F__
#pragma message "LJMultisiteFunctorAVX512_GS::SoAFunctorPairImpl called without AVX512 intrinsics!"
#endif
    using namespace autopas::utils::ArrayMath::literals;

    if (soaA.size() == 0 || soaB.size() == 0) return;

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
    // todo this could be vectorized
    const size_t siteCountB = [&]() {
      size_t siteCountTmp{0};
      for (size_t mol = 0; mol < soaB.size(); ++mol) {
        siteCountTmp += _PPLibrary->getNumSites(typeptrB[mol]);
      }
      return siteCountTmp;
    }();


    // ------------------------------ Setup auxiliary vectors -----------------------------
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedSitePositionBx;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedSitePositionBy;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedSitePositionBz;


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

    siteTypesB.reserve(siteCountB);

    isSiteOwnedBArr.reserve(siteCountB);

    // Fill site-wise std::vectors for SIMD

    for (size_t mol = 0; mol < soaB.size(); ++mol) {
      const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
          {q0Bptr[mol], q1Bptr[mol], q2Bptr[mol], q3Bptr[mol]}, _PPLibrary->getSitePositions(typeptrB[mol]));

      size_t localSiteCount = _PPLibrary->getNumSites(typeptrB[mol]);
      for (size_t site = 0; site < localSiteCount; ++site) {
        rotatedSitePositionBx.push_back(rotatedSitePositions[site][0]);
        rotatedSitePositionBy.push_back(rotatedSitePositions[site][1]);
        rotatedSitePositionBz.push_back(rotatedSitePositions[site][2]);

        exactSitePositionBx.push_back(rotatedSitePositions[site][0] + xBptr[mol]);
        exactSitePositionBy.push_back(rotatedSitePositions[site][1] + yBptr[mol]);
        exactSitePositionBz.push_back(rotatedSitePositions[site][2] + zBptr[mol]);

        isSiteOwnedBArr.push_back(ownedStatePtrB[mol] == autopas::OwnershipState::owned);

        siteTypesB.push_back(_PPLibrary->getSiteTypes(typeptrB[mol])[site]);
      }
    }

    // ------------------------------ Main force calculation loop -----------------------------

    for (size_t molA = 0; molA < soaA.size(); ++molA) {
      const auto ownedStateA = ownedStatePtrA[molA];
      if (ownedStateA == autopas::OwnershipState::dummy) {
        continue;
      }

      const auto noSitesInMolA = _PPLibrary->getNumSites(typeptrA[molA]);
      const auto unrotatedSitePositionsA = _PPLibrary->getSitePositions(typeptrA[molA]);

      const std::array<double, 3> centerOfMass{xAptr[molA], yAptr[molA], zAptr[molA]};
      const auto rotatedSitePositionsA = autopas::utils::quaternion::rotateVectorOfPositions(
          {q0Aptr[molA], q1Aptr[molA], q2Aptr[molA], q3Aptr[molA]}, unrotatedSitePositionsA);

      const auto siteTypesA = _PPLibrary->getSiteTypes(typeptrA[molA]);

      const auto sitePairIndicies = buildSiteInteractionIndices(xBptr, yBptr, zBptr, typeptrB, ownedStatePtrB, centerOfMass, 0, 0, soaB.size(),
                                                         siteCountB);

      for (size_t siteA = 0; siteA < noSitesInMolA; ++siteA) {
        const size_t siteTypeA = siteTypesA[siteA];

        const std::array<double, 3> rotatedSitePositionA = rotatedSitePositionsA[siteA];
        const std::array<double, 3> exactSitePositionA = rotatedSitePositionA + centerOfMass;

        std::array<double, 3> forceAccumulator = {0., 0., 0.};
        std::array<double, 3> torqueAccumulator = {0., 0., 0.};

        interactAllWithSiteA<newton3>(sitePairIndicies, siteTypesB, exactSitePositionBx, exactSitePositionBy,
                                      exactSitePositionBz, siteForceBx, siteForceBy, siteForceBz, isSiteOwnedBArr,
                                      ownedStateA, siteTypeA, exactSitePositionA, rotatedSitePositionA,
                                      forceAccumulator, torqueAccumulator);

        // Add forces + torques to molA
        fxAptr[molA] += forceAccumulator[0];
        fyAptr[molA] += forceAccumulator[1];
        fzAptr[molA] += forceAccumulator[2];
        txAptr[molA] += torqueAccumulator[0];
        tyAptr[molA] += torqueAccumulator[1];
        tzAptr[molA] += torqueAccumulator[2];
      }
    }

    // ------------------------------ Reduction -----------------------------

    size_t siteIndex = 0;
    for (size_t mol = 0; mol < soaB.size(); ++mol) {
      for (size_t site = 0; site < _PPLibrary->getNumSites(typeptrB[mol]); ++site) {
        fxBptr[mol] += siteForceBx[siteIndex];
        fyBptr[mol] += siteForceBy[siteIndex];
        fzBptr[mol] += siteForceBz[siteIndex];
        txBptr[mol] += rotatedSitePositionBy[siteIndex] * siteForceBz[siteIndex] -
                       rotatedSitePositionBz[siteIndex] * siteForceBy[siteIndex];
        tyBptr[mol] += rotatedSitePositionBz[siteIndex] * siteForceBx[siteIndex] -
                       rotatedSitePositionBx[siteIndex] * siteForceBz[siteIndex];
        tzBptr[mol] += rotatedSitePositionBx[siteIndex] * siteForceBy[siteIndex] -
                       rotatedSitePositionBy[siteIndex] * siteForceBx[siteIndex];
        ++siteIndex;
      }
    }

    // todo processing of global accumulators should go here (see LJMultisiteFunctor.h or 7dc729af)

  }

  /** @brief Implementation for SoAFunctorVerlet
   */
  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexPrimary,
                                  const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
#ifndef __AVX512F__
#pragma message "SoAFunctorVerlet functor called without AVX512 support!"
#endif
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    // Skip if primary particle is dummy
    const auto ownedStatePrime = ownedStatePtr[indexPrimary];
    if (ownedStatePrime == autopas::OwnershipState::dummy) {
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

    // Get number of molecules in neighbor list
    const size_t molCountNeighbors = neighborList.size();

    // Count sites of neighbors of primary molecule
    const size_t siteCountNeighbors = [&]() {
      size_t siteCountTmp{0};
      for (size_t neighborMol = 0; neighborMol < molCountNeighbors; ++neighborMol) {
        siteCountTmp += _PPLibrary->getNumSites(typeptr[neighborList[neighborMol]]);
      }
      return siteCountTmp;
    }();

    // Count sites of primary molecule
    const size_t siteCountPrimary = _PPLibrary->getNumSites(typeptr[indexPrimary]);

    // Accumulators for global values
    double potentialEnergyAccumulator = 0;
    std::array<double, 3> virialAccumulator = {0., 0., 0.};

    // ------------------------------ Setup auxiliary vectors -----------------------------
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedNeighborSitePositionsX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedNeighborSitePositionsY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedNeighborSitePositionsZ;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactNeighborSitePositionsX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactNeighborSitePositionsY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactNeighborSitePositionsZ;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceX(siteCountNeighbors, 0.);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceY(siteCountNeighbors, 0.);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceZ(siteCountNeighbors, 0.);

    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteTypesNeighbors;
    std::vector<size_t, autopas::AlignedAllocator<size_t>> isNeighborSiteOwnedArr;


    // pre-reserve site std::vectors
    rotatedNeighborSitePositionsX.reserve(siteCountNeighbors);
    rotatedNeighborSitePositionsY.reserve(siteCountNeighbors);
    rotatedNeighborSitePositionsZ.reserve(siteCountNeighbors);
    exactNeighborSitePositionsX.reserve(siteCountNeighbors);
    exactNeighborSitePositionsY.reserve(siteCountNeighbors);
    exactNeighborSitePositionsZ.reserve(siteCountNeighbors);

    siteTypesNeighbors.reserve(siteCountNeighbors);
    isNeighborSiteOwnedArr.reserve(siteCountNeighbors);

    const auto rotatedSitePositionsInPrimaryMol =
        autopas::utils::quaternion::rotateVectorOfPositions(
                        {q0ptr[indexPrimary], q1ptr[indexPrimary], q2ptr[indexPrimary], q3ptr[indexPrimary]},
                        _PPLibrary->getSitePositions(typeptr[indexPrimary]));

    const auto siteTypesInPrimaryMol = _PPLibrary->getSiteTypes(typeptr[indexPrimary]);

    // Fill site-wise std::vectors for SIMD
    for (size_t neighborMol = 0; neighborMol < molCountNeighbors; ++neighborMol) {
      const size_t neighborMolIndex = neighborList[neighborMol];
      const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
          {q0ptr[neighborMolIndex], q1ptr[neighborMolIndex], q2ptr[neighborMolIndex], q3ptr[neighborMolIndex]},
          _PPLibrary->getSitePositions(typeptr[neighborMolIndex]));

      const size_t localSiteCount = _PPLibrary->getNumSites(typeptr[neighborMolIndex]);

      for (size_t site = 0; site < localSiteCount; ++site) {
        rotatedNeighborSitePositionsX.push_back(rotatedSitePositions[site][0]);
        rotatedNeighborSitePositionsY.push_back(rotatedSitePositions[site][1]);
        rotatedNeighborSitePositionsZ.push_back(rotatedSitePositions[site][2]);
        exactNeighborSitePositionsX.push_back(rotatedSitePositions[site][0] + xptr[neighborMolIndex]);
        exactNeighborSitePositionsY.push_back(rotatedSitePositions[site][1] + yptr[neighborMolIndex]);
        exactNeighborSitePositionsZ.push_back(rotatedSitePositions[site][2] + zptr[neighborMolIndex]);
        isNeighborSiteOwnedArr.push_back(ownedStatePtr[neighborMolIndex] == autopas::OwnershipState::owned);
        siteTypesNeighbors.push_back(_PPLibrary->getSiteTypes(typeptr[neighborMolIndex])[site]);
      }
    }

    const std::array<double, 3> centerOfMassPrimary{xptr[indexPrimary], yptr[indexPrimary], zptr[indexPrimary]};


    const auto sitePairIndicies = buildSiteInteractionIndicesVerlet(xptr, yptr, zptr, typeptr, ownedStatePtr, centerOfMassPrimary, neighborList, molCountNeighbors, siteCountNeighbors);


    // ------------------- Calculate Forces -------------------

    std::array<double, 3> forceAccumulator = {0., 0., 0.};
    std::array<double, 3> torqueAccumulator = {0., 0., 0.};

    for (size_t primarySite = 0; primarySite < siteCountPrimary; ++primarySite) {
      const size_t siteTypePrime = siteTypesInPrimaryMol[primarySite];
      const std::array<double, 3> rotatedSitePositionPrime = rotatedSitePositionsInPrimaryMol[primarySite];
      const std::array<double, 3> exactSitePositionPrime =
          autopas::utils::ArrayMath::add(rotatedSitePositionPrime, centerOfMassPrimary);

      interactAllWithSiteA<newton3>(
          sitePairIndicies, siteTypesNeighbors, exactNeighborSitePositionsX, exactNeighborSitePositionsY,
          exactNeighborSitePositionsZ, siteForceX, siteForceY, siteForceZ, isNeighborSiteOwnedArr, ownedStatePrime,
          siteTypePrime, exactSitePositionPrime, rotatedSitePositionPrime, forceAccumulator, torqueAccumulator);
    }
    // Add forces to prime mol
    fxptr[indexPrimary] += forceAccumulator[0];
    fyptr[indexPrimary] += forceAccumulator[1];
    fzptr[indexPrimary] += forceAccumulator[2];
    txptr[indexPrimary] += torqueAccumulator[0];
    typtr[indexPrimary] += torqueAccumulator[1];
    tzptr[indexPrimary] += torqueAccumulator[2];

    // ------------------------------ Reduction -----------------------------
    if constexpr (newton3) {
      size_t siteIndex = 0;
      for (size_t neighborMol = 0; neighborMol < molCountNeighbors; ++neighborMol) {
        const auto neighborMolIndex = neighborList[neighborMol];
        for (size_t site = 0; site < _PPLibrary->getNumSites(typeptr[neighborMolIndex]); ++site) {
          fxptr[neighborMolIndex] += siteForceX[siteIndex];
          fyptr[neighborMolIndex] += siteForceY[siteIndex];
          fzptr[neighborMolIndex] += siteForceZ[siteIndex];
          txptr[neighborMolIndex] += rotatedNeighborSitePositionsY[siteIndex] * siteForceZ[siteIndex] -
                                     rotatedNeighborSitePositionsZ[siteIndex] * siteForceY[siteIndex];
          typtr[neighborMolIndex] += rotatedNeighborSitePositionsZ[siteIndex] * siteForceX[siteIndex] -
                                     rotatedNeighborSitePositionsX[siteIndex] * siteForceZ[siteIndex];
          tzptr[neighborMolIndex] += rotatedNeighborSitePositionsX[siteIndex] * siteForceY[siteIndex] -
                                     rotatedNeighborSitePositionsY[siteIndex] * siteForceX[siteIndex];
          ++siteIndex;
        }
      }
    }

    // todo processing of global accumulators should go here (see LJMultisiteFunctor.h or 7dc729af)

  }

  template <bool newton3>
  inline void interactAllWithSiteA(const std::vector<size_t, autopas::AlignedAllocator<size_t>> &sitePairIndicies,
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
                            std::array<double, 3> &torqueAccumulator, size_t offset = 0) {
#ifndef __AVX512F__
#pragma message "LJMultisiteFunctorAVX51_GS.h included, but AVX512 is not supported by the compiler."
#else
    // Declare mixing variables
    __m512d sigmaSquared = _one;
    __m512d epsilon24 = _zero;
    __m512d shift6 = _zero;

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

    for (size_t siteVectorIndex = 0; siteVectorIndex < sitePairIndicies.size(); siteVectorIndex += 8) {
      const size_t remainder = sitePairIndicies.size() - siteVectorIndex;
      const bool remainderCase = remainder < 8;
      const __mmask8 remainderMask = remainderCase ? _remainderMasks[8-remainder] : __mmask8(255);

      const __m512i sitePairIndicesVec = remainderCase ? _mm512_maskz_loadu_epi64(remainderMask, &sitePairIndicies[siteVectorIndex])
          : _mm512_loadu_epi64(&sitePairIndicies[siteVectorIndex]);


      // gather mixing parameters.
      const double *const __restrict mixingPtr = _PPLibrary->getMixingDataPtr(siteTypeA, 0);

      // Uses base address of mixingPtr (+0/1/2), gathered data is offset by siteTypeIndicesScaled x 64 bits x 8.
      // We need scaled site-types due to the way the mixing data is stored - todo change this
      if (remainderCase) {
        const auto siteTypeIndicies = _mm512_mask_i64gather_epi64(_zeroI, remainderMask, sitePairIndicesVec, &siteTypesB[0], 8);
        const __m512i siteTypeIndicesScaled = _mm512_mullox_epi64(siteTypeIndicies, _three);


        epsilon24 = _mm512_mask_i64gather_pd(_zero, remainderMask, siteTypeIndicesScaled, mixingPtr, 8);
        sigmaSquared =
            _mm512_mask_i64gather_pd(_one, remainderMask, siteTypeIndicesScaled, mixingPtr + 1, 8);  // one used as "filler" in remainder case to avoid dividing by zero
        if constexpr (applyShift) {
          shift6 = _mm512_mask_i64gather_pd(_zero, remainderMask, siteTypeIndicesScaled, mixingPtr + 2, 8);
        }

      } else {
        const auto siteTypeIndicies = _mm512_i64gather_epi64(sitePairIndicesVec, &siteTypesB[0], 8);
        const __m512i siteTypeIndicesScaled = _mm512_mullox_epi64(siteTypeIndicies, _three);

        epsilon24 = _mm512_i64gather_pd(siteTypeIndicesScaled, mixingPtr, 8);
        sigmaSquared = _mm512_i64gather_pd(siteTypeIndicesScaled, mixingPtr + 1, 8);
        if constexpr (applyShift) {
          shift6 = _mm512_i64gather_pd(siteTypeIndicesScaled, mixingPtr + 2, 8);
        }
    }

      // Load the exact site positions of particle B // todo is it faster to just use _mm512_maskz_load_pd?
      const __m512d exactSitePositionsBX = remainderCase ? _mm512_mask_i64gather_pd(_zero, remainderMask, sitePairIndicesVec, &exactSitePositionX[0], 8)
                                                         : _mm512_i64gather_pd(sitePairIndicesVec, &exactSitePositionX[0], 8);
      const __m512d exactSitePositionsBY = remainderCase ? _mm512_mask_i64gather_pd(_zero, remainderMask, sitePairIndicesVec, &exactSitePositionY[0], 8)
                                                         : _mm512_i64gather_pd(sitePairIndicesVec, &exactSitePositionY[0], 8);
      const __m512d exactSitePositionsBZ = remainderCase ? _mm512_mask_i64gather_pd(_zero, remainderMask, sitePairIndicesVec, &exactSitePositionZ[0], 8)
                                                         : _mm512_i64gather_pd(sitePairIndicesVec, &exactSitePositionZ[0], 8);

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

      // Determine forces
      const __m512d forceX = _mm512_mul_pd(scalar, displacementX);
      const __m512d forceY = _mm512_mul_pd(scalar, displacementY);
      const __m512d forceZ = _mm512_mul_pd(scalar, displacementZ);

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
          const auto forceSumBX = _mm512_mask_i64gather_pd(_zero, remainderMask, sitePairIndicesVec, &siteForceX[0], 8);
          const __m512d newForceSumBX = _mm512_sub_pd(forceSumBX, forceX);
          _mm512_mask_i64scatter_pd(&siteForceX[0], remainderMask, sitePairIndicesVec, newForceSumBX, 8);

          const auto forceSumBY = _mm512_mask_i64gather_pd(_zero, remainderMask, sitePairIndicesVec, &siteForceY[0], 8);
          const __m512d newForceSumBY = _mm512_sub_pd(forceSumBY, forceY);
          _mm512_mask_i64scatter_pd(&siteForceY[0], remainderMask, sitePairIndicesVec, newForceSumBY, 8);

          const auto forceSumBZ = _mm512_mask_i64gather_pd(_zero, remainderMask, sitePairIndicesVec, &siteForceZ[0], 8);
          const __m512d newForceSumBZ = _mm512_sub_pd(forceSumBZ, forceZ);
          _mm512_mask_i64scatter_pd(&siteForceZ[0], remainderMask, sitePairIndicesVec, newForceSumBZ, 8);
        } else {
          const auto forceSumBX = _mm512_i64gather_pd(sitePairIndicesVec, &siteForceX[0], 8);
          const __m512d newForceSumBX = _mm512_sub_pd(forceSumBX, forceX);
          _mm512_i64scatter_pd(&siteForceX[0], sitePairIndicesVec, newForceSumBX, 8);

          const auto forceSumBY = _mm512_i64gather_pd(sitePairIndicesVec, &siteForceY[0], 8);
          const __m512d newForceSumBY = _mm512_sub_pd(forceSumBY, forceY);
          _mm512_i64scatter_pd(&siteForceY[0], sitePairIndicesVec, newForceSumBY, 8);

          const auto forceSumBZ = _mm512_i64gather_pd(sitePairIndicesVec, &siteForceZ[0], 8);
          const __m512d newForceSumBZ = _mm512_sub_pd(forceSumBZ, forceZ);
          _mm512_i64scatter_pd(&siteForceZ[0], sitePairIndicesVec, newForceSumBZ, 8);
        }

      }

      // todo handling of globals should go here (see LJMultisiteFunctor.h or 7dc729af)

    }

    // Add up the forces, torques and globals

    forceAccumulator[0] += _mm512_reduce_add_pd(forceSumX);
    forceAccumulator[1] += _mm512_reduce_add_pd(forceSumY);
    forceAccumulator[2] += _mm512_reduce_add_pd(forceSumZ);


    torqueAccumulator[0] += _mm512_reduce_add_pd(torqueSumX);
    torqueAccumulator[1] += _mm512_reduce_add_pd(torqueSumY);
    torqueAccumulator[2] += _mm512_reduce_add_pd(torqueSumZ);


    // todo accumulation of globals should go here (see LJMultisiteFunctor.h or 7dc729af)

#endif
  }

  /**
   * Build site interaction indices for SoASingle & SoAPair.
   */
  inline std::vector<size_t, autopas::AlignedAllocator<size_t>> buildSiteInteractionIndices(
      const double *const __restrict xptr, const double *const __restrict yptr, const double *const __restrict zptr,
      const size_t *const __restrict typeptr, const autopas::OwnershipState *const __restrict ownedStatePtr,
      const std::array<double, 3> &centerOfMass, size_t firstMolIndex, size_t firstSiteIndex, size_t noMolecules, size_t noSites) {

    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteIndices;
    siteIndices.reserve(noSites / 4); // This is a complete guess, which is only really valid for Linked Cells

    size_t siteIndex = firstSiteIndex;

    for (size_t molIndex = firstMolIndex; molIndex < noMolecules; molIndex++) {
      const size_t siteCount = _PPLibrary->getNumSites(typeptr[molIndex]);

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
          siteIndices.push_back(siteIndex);
        }
        siteIndex++;
      }
    }
    siteIndices.shrink_to_fit();
    return siteIndices;
  }

  /**
   * Build site interaction indices for SoAVerlet.
   * @return a list of site indices that satisfy the cutoff criterion.
   * @warning This site indices returned match the site positions generated, i.e. are indicies of the neighbor list's sites
   * not the SoA's sites.
   */
  inline std::vector<size_t, autopas::AlignedAllocator<size_t>> buildSiteInteractionIndicesVerlet(
      const double *const __restrict xptr, const double *const __restrict yptr, const double *const __restrict zptr,
      const size_t *const __restrict typeptr, const autopas::OwnershipState *const __restrict ownedStatePtr,
      const std::array<double, 3> &centerOfMass, const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
      const size_t noNeighborMol, const size_t noNeighborSite) {

    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteIndices;
    siteIndices.reserve(noNeighborSite); // This is just a guess

    size_t siteIndex = 0;

    for (size_t molIndex = 0; molIndex < noNeighborMol; molIndex++) {
      const size_t siteCount = _PPLibrary->getNumSites(typeptr[neighborList[molIndex]]);

      const double xPosB = xptr[neighborList[molIndex]];
      const double yPosB = yptr[neighborList[molIndex]];
      const double zPosB = zptr[neighborList[molIndex]];

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
      const bool dummyCondition = ownedStatePtr[neighborList[molIndex]] != autopas::OwnershipState::dummy;
      const bool condition = cutoffCondition and dummyCondition;

      for (size_t site = 0; site < siteCount; site++) {
        if (condition) {
          siteIndices.push_back(siteIndex);
        }
        siteIndex++;
      }
    }
    siteIndices.shrink_to_fit();
    return siteIndices;
  }

 public:
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
  constexpr static bool getMixing() { return true; }

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