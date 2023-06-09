/**
 * @file LJMultisiteFunctorAVX.h
 * @date 03/04/2022
 * @author Q. Behrami
 */
#pragma once

#ifndef __AVX__
#pragma message "LJMultisiteFunctorAVX.h included, but AVX is not supported by the compiler."
#else
#include "immintrin.h"
#endif

#include <array>

#include "../pairwiseFunctors/Functor.h"
#include "../utils/ExceptionHandler.h"
#include "MultisiteMoleculeLJ.h"
#include "ParticlePropertiesLibrary.h"
#include "autopas/utils/AVXUtils.h"
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
class LJMultisiteFunctorAVX
    : public autopas::Functor<Particle, LJMultisiteFunctorAVX<Particle, applyShift, useMixing, useNewton3,
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
  constexpr static size_t vecLength = 4;

  /*
   * Whether to build vectors for exact site position locally in each function call or once globally.
   * @note Currently only implemented for Verlet functor
   */
  constexpr static bool useLocalVectors = true;

  /*
   * Whether to use site masks or site indices
   * @note Currently only implemented for Verlet functor
   */
  constexpr static bool useSiteMasks = true;

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

  // Vectors used to speed up the verlet functor
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
  LJMultisiteFunctorAVX() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy unused, only there to make the signature different from the public constructor.
   */
  explicit LJMultisiteFunctorAVX(double cutoff, void * /*dummy*/)
#ifdef __AVX__
      : Functor<Particle, LJMultisiteFunctorAVX<Particle, applyShift, useMixing, useNewton3, calculateGlobals,
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
      : Functor<Particle, LJMultisiteFunctorAVX<Particle, applyShift, useMixing, useNewton3, calculateGlobals,
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
  explicit LJMultisiteFunctorAVX(double cutoff) : LJMultisiteFunctorAVX(cutoff, nullptr) {
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
  explicit LJMultisiteFunctorAVX(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJMultisiteFunctorAVX(cutoff, nullptr) {
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
    if constexpr (useLocalVectors) {
      if (newton3) {
        SoAFunctorVerletImpl_local<true>(soa, indexFirst, neighborList);
      } else {
        SoAFunctorVerletImpl_local<false>(soa, indexFirst, neighborList);
      }
    } else {
      if (newton3) {
        SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
      } else {
        SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
      }
    }
  }

 private:
  /**
   * SoAFunctorSingle Implementation.
   */
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

    // count number of sites in SoA
    size_t siteCount = 0;
    if constexpr (useMixing) {
      for (size_t mol = 0; mol < soa.getNumberOfParticles(); ++mol) {
        siteCount += _PPLibrary->getNumSites(typeptr[mol]);
      }
    } else {
      siteCount = _sitePositionsLJ.size() * soa.getNumberOfParticles();
    }

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionZ;

    // we require arrays for forces for sites to maintain SIMD in site-site calculations
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceX(siteCount, 0);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceY(siteCount, 0);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceZ(siteCount, 0);

    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteTypes;
    std::vector<size_t, autopas::AlignedAllocator<size_t>> isSiteOwned;

    // pre-reserve site std::vectors
    exactSitePositionX.reserve(siteCount);
    exactSitePositionY.reserve(siteCount);
    exactSitePositionZ.reserve(siteCount);

    if constexpr (useMixing) {
      siteTypes.reserve(siteCount);
    }

    // this is only needed for vectorization when calculating globals
    if constexpr (calculateGlobals) {
      isSiteOwned.reserve(siteCount);
    }

    // Fill site-wise std::vectors for SIMD
    std::vector<std::array<double, 3>> rotatedSitePositions;
    for (size_t mol = 0; mol < soa.getNumberOfParticles(); ++mol) {
      if constexpr (useMixing) {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, _PPLibrary->getSitePositions(typeptr[mol]));
      } else {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, _sitePositionsLJ);
      }
      size_t localSiteCount = useMixing ? _PPLibrary->getNumSites(typeptr[mol]) : _sitePositionsLJ.size();

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

    double potentialEnergyAccumulator = 0;
    std::array<double, 3> virialAccumulator = {0., 0., 0.};

    // main force calculation loop
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

        const __m256d xposB = autopas::utils::avx::load_pd(remainderCase, &xptr[molB], remainderMask);
        const __m256d yposB = autopas::utils::avx::load_pd(remainderCase, &yptr[molB], remainderMask);
        const __m256d zposB = autopas::utils::avx::load_pd(remainderCase, &zptr[molB], remainderMask);

        const __m256d displacementCoMX = _mm256_sub_pd(xposA, xposB);
        const __m256d displacementCoMY = _mm256_sub_pd(yposA, yposB);
        const __m256d displacementCoMZ = _mm256_sub_pd(zposA, zposB);

        const __m256d distanceSquaredCoMX = _mm256_mul_pd(displacementCoMX, displacementCoMX);
        const __m256d distanceSquaredCoMY = _mm256_mul_pd(displacementCoMY, displacementCoMY);
        const __m256d distanceSquaredCoMZ = _mm256_mul_pd(displacementCoMZ, displacementCoMZ);

        const __m256d distanceSquaredCoM =
            _mm256_add_pd(_mm256_add_pd(distanceSquaredCoMX, distanceSquaredCoMY), distanceSquaredCoMZ);

        const __m256d cutoffMask = _mm256_cmp_pd(distanceSquaredCoM, _cutoffSquared, _CMP_LE_OS);
        const __m256i ownedStateB = autopas::utils::avx::load_epi64(
            remainderCase, reinterpret_cast<const size_t *>(&ownedStatePtr[molB]), remainderMask);
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

      // These are technically not needed for this functor
      const std::array<double, 3> rotatedSitePositionA{0.0, 0.0, 0.0};
      std::array<double, 3> torqueAccumulator = {0., 0., 0.};

      // ------- Main force calculation loop -------
      for (size_t siteA = siteIndexMolA; siteA < siteIndexMolB; ++siteA) {
        const size_t siteTypeA = siteTypes[siteA];
        const std::array<double, 3> exactSitePositionA = {exactSitePositionX[siteA], exactSitePositionY[siteA],
                                                          exactSitePositionZ[siteA]};

        std::array<double, 3> forceAccumulator = {0., 0., 0.};

        SoAKernel<true, false, false>(siteVector, siteTypes, exactSitePositionX, exactSitePositionY, exactSitePositionZ,
                                      siteForceX, siteForceY, siteForceZ, isSiteOwned, ownedStateA, siteTypeA,
                                      exactSitePositionA, rotatedSitePositionA, forceAccumulator, torqueAccumulator,
                                      potentialEnergyAccumulator, virialAccumulator, siteIndexMolB);

        // sum forces on single site in mol A
        siteForceX[siteA] += forceAccumulator[0];
        siteForceY[siteA] += forceAccumulator[1];
        siteForceZ[siteA] += forceAccumulator[2];
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
            {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, _sitePositionsLJ);
        for (size_t site = 0; site < _sitePositionsLJ.size(); ++site) {
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

      _aosThreadData[threadNum].potentialEnergySum += potentialEnergyAccumulator * newton3Factor;
      _aosThreadData[threadNum].virialSum[0] += virialAccumulator[0] * newton3Factor;
      _aosThreadData[threadNum].virialSum[1] += virialAccumulator[1] * newton3Factor;
      _aosThreadData[threadNum].virialSum[2] += virialAccumulator[2] * newton3Factor;
    }
  }

  /**
   * Implementation of SoAFunctorPair(soa1, soa2, newton3)
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

    // count number of sites in both SoAs
    size_t siteCountB = 0;
    if constexpr (useMixing) {
      for (size_t mol = 0; mol < soaB.getNumberOfParticles(); ++mol) {
        siteCountB += _PPLibrary->getNumSites(typeptrB[mol]);
      }
    } else {
      siteCountB = _sitePositionsLJ.size() * soaB.getNumberOfParticles();
    }

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionBx;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionBy;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionBz;

    // we require arrays for forces for sites to maintain SIMD in site-site calculations
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

    // this is only needed for vectorization when calculating globals
    if constexpr (calculateGlobals) {
      isSiteOwnedBArr.reserve(siteCountB);
    }

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
        if constexpr (calculateGlobals) {
          isSiteOwnedBArr.push_back(ownedStatePtrB[mol] == OwnershipState::owned);
        }
        if constexpr (useMixing) {
          siteTypesB.push_back(_PPLibrary->getSiteTypes(typeptrB[mol])[site]);
        }
      }
    }

    double potentialEnergyAccumulator = 0;
    std::array<double, 3> virialAccumulator = {0., 0., 0.};

    // main force calculation loop
    for (size_t molA = 0; molA < soaA.getNumberOfParticles(); ++molA) {
      const auto ownedStateA = ownedStatePtrA[molA];
      if (ownedStateA == autopas::OwnershipState::dummy) {
        continue;
      }

      const auto noSitesInMolA = useMixing ? _PPLibrary->getNumSites(typeptrA[molA]) : _sitePositionsLJ.size();
      const auto unrotatedSitePositionsA = useMixing ? _PPLibrary->getSitePositions(typeptrA[molA]) : _sitePositionsLJ;

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

        const __m256d xposB = autopas::utils::avx::load_pd(remainderCase, &xBptr[molB], remainderMask);
        const __m256d yposB = autopas::utils::avx::load_pd(remainderCase, &yBptr[molB], remainderMask);
        const __m256d zposB = autopas::utils::avx::load_pd(remainderCase, &zBptr[molB], remainderMask);

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
        size_t siteCount = useMixing ? _PPLibrary->getNumSites(typeptrB[molB]) : _sitePositionsLJ.size();
        if (condition) {
          for (size_t siteB = 0; siteB < siteCount; ++siteB) {
            siteVector.emplace_back(siteIndex++);
          }
        } else {
          siteIndex += siteCount;
        }
      }

      std::array<double, 3> forceAccumulator = {0., 0., 0.};
      std::array<double, 3> torqueAccumulator = {0., 0., 0.};

      for (size_t siteA = 0; siteA < noSitesInMolA; ++siteA) {
        const size_t siteTypeA = _PPLibrary->getSiteTypes(typeptrA[molA])[siteA];
        const std::array<double, 3> rotatedSitePositionA = rotatedSitePositionsA[siteA];
        const std::array<double, 3> exactSitePositionA =
            autopas::utils::ArrayMath::add(rotatedSitePositionA, {xAptr[molA], yAptr[molA], zAptr[molA]});

        SoAKernel<newton3, true, false>(
            siteVector, siteTypesB, exactSitePositionBx, exactSitePositionBy, exactSitePositionBz, siteForceBx,
            siteForceBy, siteForceBz, isSiteOwnedBArr, ownedStateA, siteTypeA, exactSitePositionA, rotatedSitePositionA,
            forceAccumulator, torqueAccumulator, potentialEnergyAccumulator, virialAccumulator);
      }
      fxAptr[molA] += forceAccumulator[0];
      fyAptr[molA] += forceAccumulator[1];
      fzAptr[molA] += forceAccumulator[2];
      txAptr[molA] += torqueAccumulator[0];
      tyAptr[molA] += torqueAccumulator[1];
      tzAptr[molA] += torqueAccumulator[2];
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
      const auto threadNum = autopas_get_thread_num();
      // SoAFunctorPairImpl obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
      // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.
      const auto newton3Factor = newton3 ? .5 : 1.;

      _aosThreadData[threadNum].potentialEnergySum += potentialEnergyAccumulator * newton3Factor;
      _aosThreadData[threadNum].virialSum[0] += virialAccumulator[0] * newton3Factor;
      _aosThreadData[threadNum].virialSum[1] += virialAccumulator[1] * newton3Factor;
      _aosThreadData[threadNum].virialSum[2] += virialAccumulator[2] * newton3Factor;
    }
  }

  // Implementation for SoAFunctorVerlet
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

#pragma omp single  // Needed if multiple threads are running
    buildSiteVectors(soa);

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

    // -- main force calculation --

    // - calculate mol mask -
    std::vector<double, autopas::AlignedAllocator<double>> molMask;
    molMask.reserve(neighborListSize);

    const __m256d xPosFirst = _mm256_broadcast_sd(&xptr[indexFirst]);
    const __m256d yPosFirst = _mm256_broadcast_sd(&yptr[indexFirst]);
    const __m256d zPosFirst = _mm256_broadcast_sd(&zptr[indexFirst]);

    for (size_t neighborMol = 0; neighborMol < neighborListSize; neighborMol += vecLength) {
      const size_t rest = neighborListSize - neighborMol;
      const bool remainderCase = rest < vecLength;
      __m256i remainderMask = remainderCase ? _masks[rest - 1] : _mm256_set1_epi64x(-1);
      const __m256i neighborMolIndex =
          autopas::utils::avx::load_epi64(remainderCase, &neighborList[neighborMol], remainderMask);

      const __m256d xposB = autopas::utils::avx::gather_pd(remainderCase, xptr, neighborMolIndex, remainderMask);
      const __m256d yposB = autopas::utils::avx::gather_pd(remainderCase, yptr, neighborMolIndex, remainderMask);
      const __m256d zposB = autopas::utils::avx::gather_pd(remainderCase, zptr, neighborMolIndex, remainderMask);

      const __m256d displacementCoMX = _mm256_sub_pd(xPosFirst, xposB);
      const __m256d displacementCoMY = _mm256_sub_pd(yPosFirst, yposB);
      const __m256d displacementCoMZ = _mm256_sub_pd(zPosFirst, zposB);

      const __m256d distanceSquaredCoMX = _mm256_mul_pd(displacementCoMX, displacementCoMX);
      const __m256d distanceSquaredCoMY = _mm256_mul_pd(displacementCoMY, displacementCoMY);
      const __m256d distanceSquaredCoMZ = _mm256_mul_pd(displacementCoMZ, displacementCoMZ);

      const __m256d distanceSquaredCoM =
          _mm256_add_pd(_mm256_add_pd(distanceSquaredCoMX, distanceSquaredCoMY), distanceSquaredCoMZ);

      const __m256d cutoffMask = _mm256_cmp_pd(distanceSquaredCoM, _cutoffSquared, _CMP_LE_OS);
      const __m256i ownedStateB =
          remainderCase ? _mm256_mask_i64gather_epi64(_zero, ownedStatePtr, neighborMolIndex, remainderMask, 8)
                        : _mm256_i64gather_pd(ownedStatePtr, neighborMolIndex, 8);
      const __m256d dummyMask =
          _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateB), _zero, _CMP_NEQ_OS);  // Assuming that dummy = 0
      const __m256d totalMask = _mm256_and_pd(cutoffMask, dummyMask);

      if (remainderCase) {
        _mm256_maskstore_pd((&molMask[neighborMol]), remainderMask, totalMask);
      } else {
        _mm256_storeu_pd((&molMask[neighborMol]), totalMask);
      }
    }

    //     generate mask for each site from molecular mask
    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteVector;
    siteVector.reserve(siteCountNeighbors);

    for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
      bool condition = static_cast<bool>(molMask[neighborMol]);
      const auto neighborMolIndex = neighborList[neighborMol];
      size_t siteCount =
          useMixing ? _PPLibrary->getNumSites(typeptr[neighborMolIndex]) : const_unrotatedSitePositions.size();
      if (condition) {
        for (size_t siteB = _molToSiteMap[neighborMolIndex]; siteB < _molToSiteMap[neighborMolIndex] + siteCount;
             ++siteB) {
          siteVector.emplace_back(siteB);
        }
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
      const double *mixingPtr =
          useMixing ? _PPLibrary->getMixingDataPtr(_siteTypesVector[_molToSiteMap[indexFirst] + primeSite], 0)
                    : nullptr;

      const __m256d exactPrimeSitePositionX =
          _mm256_broadcast_sd(&_exactSitePositionsX[_molToSiteMap[indexFirst]] + primeSite);
      const __m256d exactPrimeSitePositionY =
          _mm256_broadcast_sd(&_exactSitePositionsY[_molToSiteMap[indexFirst]] + primeSite);
      const __m256d exactPrimeSitePositionZ =
          _mm256_broadcast_sd(&_exactSitePositionsZ[_molToSiteMap[indexFirst]] + primeSite);

      const __m256d rotatedPrimeSitePositionX = _mm256_sub_pd(exactPrimeSitePositionX, xPosFirst);
      const __m256d rotatedPrimeSitePositionY = _mm256_sub_pd(exactPrimeSitePositionY, yPosFirst);
      const __m256d rotatedPrimeSitePositionZ = _mm256_sub_pd(exactPrimeSitePositionZ, zPosFirst);

      for (size_t siteVectorIndex = 0; siteVectorIndex < siteVector.size(); siteVectorIndex += vecLength) {
        const size_t remainder = siteVector.size() - siteVectorIndex;
        const bool remainderCase = remainder < vecLength;
        const __m256i remainderMask = remainderCase ? _masks[remainder - 1] : _mm256_set1_epi64x(-1);

        const __m256i siteIndices =
            autopas::utils::avx::load_epi64(remainderCase, &siteVector[siteVectorIndex], remainderMask);

        if constexpr (useMixing) {
          __m256i siteTypeValues =
              autopas::utils::avx::gather_epi64(remainderCase, _siteTypesVector.data(), siteIndices, remainderMask);
          siteTypeValues = _mm256_mul_epu32(siteTypeValues, _mm256_set1_epi64x(3));  // This could maybe be problematic!

          epsilon24 = _mm256_i64gather_pd(mixingPtr, siteTypeValues, 8);
          sigmaSquared = _mm256_i64gather_pd(mixingPtr + 1, siteTypeValues, 8);
          if constexpr (applyShift) {
            shift6 = _mm256_i64gather_pd(mixingPtr + 2, siteTypeValues, 8);
          }
        }

        const __m256d exactNeighborSitePositionX =
            autopas::utils::avx::gather_pd(remainderCase, _exactSitePositionsX.data(), siteIndices, remainderMask);
        const __m256d exactNeighborSitePositionY =
            autopas::utils::avx::gather_pd(remainderCase, _exactSitePositionsY.data(), siteIndices, remainderMask);
        const __m256d exactNeighborSitePositionZ =
            autopas::utils::avx::gather_pd(remainderCase, _exactSitePositionsZ.data(), siteIndices, remainderMask);

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
          // TODO don't know if vectorized version is possible here
          for (size_t i = 0; i < std::min(vecLength, remainder); i++) {
            size_t site = siteIndices[i];
            size_t mol = _siteToMolMap[site];
            double siteForceX = -forceX[i];
            double siteForceY = -forceY[i];
            double siteForceZ = -forceZ[i];
            fxptr[mol] += siteForceX;
            fyptr[mol] += siteForceY;
            fzptr[mol] += siteForceZ;

            const auto rotatedSitePositionX = _exactSitePositionsX[site] - xptr[mol];
            const auto rotatedSitePositionY = _exactSitePositionsY[site] - yptr[mol];
            const auto rotatedSitePositionZ = _exactSitePositionsZ[site] - zptr[mol];

            txptr[mol] += rotatedSitePositionY * siteForceZ - rotatedSitePositionZ * siteForceY;
            typtr[mol] += rotatedSitePositionZ * siteForceX - rotatedSitePositionX * siteForceZ;
            tzptr[mol] += rotatedSitePositionX * siteForceY - rotatedSitePositionY * siteForceX;
          }
        }

        // calculate globals
        if constexpr (calculateGlobals) {
          const __m256d virialX = _mm256_mul_pd(displacementX, forceX);
          const __m256d virialY = _mm256_mul_pd(displacementY, forceY);
          const __m256d virialZ = _mm256_mul_pd(displacementZ, forceZ);

          const __m256d potentialEnergy6 =
              _mm256_fmadd_pd(epsilon24, lj12m6, shift6);  // FMA may not be supported on all CPUs
          const __m256d potentialEnergy6Masked = _mm256_and_pd(_mm256_castsi256_pd(remainderMask), potentialEnergy6);

          __m256i ownedStateP4 = _mm256_set1_epi64x(static_cast<int64_t>(ownedStatePrime));
          __m256d ownedMaskPrime =
              _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateP4), _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
          __m256d energyFactor = _mm256_blendv_pd(_zero, _one, ownedMaskPrime);
          if constexpr (newton3) {
            const __m256i ownedStateB = autopas::utils::avx::gather_epi64(remainderCase, _siteOwnershipVector.data(),
                                                                          siteIndices, remainderMask);
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
    fxptr[indexFirst] += autopas::utils::avx::horizontalSum(forceSumX);
    fyptr[indexFirst] += autopas::utils::avx::horizontalSum(forceSumY);
    fzptr[indexFirst] += autopas::utils::avx::horizontalSum(forceSumZ);
    txptr[indexFirst] += autopas::utils::avx::horizontalSum(torqueSumX);
    typtr[indexFirst] += autopas::utils::avx::horizontalSum(torqueSumY);
    tzptr[indexFirst] += autopas::utils::avx::horizontalSum(torqueSumZ);

    if constexpr (calculateGlobals) {
      const auto threadNum = autopas_get_thread_num();
      // SoAFunctorSingle obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
      // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.
      const auto newton3Factor = newton3 ? .5 : 1.;

      _aosThreadData[threadNum].potentialEnergySum +=
          autopas::utils::avx::horizontalSum(potentialEnergySum) * newton3Factor;
      _aosThreadData[threadNum].virialSum[0] += autopas::utils::avx::horizontalSum(virialSumX) * newton3Factor;
      _aosThreadData[threadNum].virialSum[1] += autopas::utils::avx::horizontalSum(virialSumY) * newton3Factor;
      _aosThreadData[threadNum].virialSum[2] += autopas::utils::avx::horizontalSum(virialSumZ) * newton3Factor;
    }
  }

  // Implementation for SoAFunctorVerlet, rebuilds exact site position vectors for each functor call
  template <bool newton3>
  void SoAFunctorVerletImpl_local(SoAView<SoAArraysType> soa, const size_t indexFirst,
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

    const size_t neighborListSize = neighborList.size();
    const size_t *const __restrict neighborListPtr = neighborList.data();

    size_t siteCountNeighbors = 0;  // site count of neighbours of primary molecule
    if constexpr (useMixing) {
      for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
        siteCountNeighbors += _PPLibrary->getNumSites(typeptr[neighborList[neighborMol]]);
      }
    } else {
      siteCountNeighbors = _sitePositionsLJ.size() * neighborListSize;
    }

    // Count sites
    const size_t siteCountMolPrime = useMixing ? _PPLibrary->getNumSites(typeptr[indexFirst]) : _sitePositionsLJ.size();

    // initialize site-wise arrays for neighbors
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactNeighborSitePositionsX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactNeighborSitePositionsY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactNeighborSitePositionsZ;

    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteTypesNeighbors;
    std::vector<size_t, autopas::AlignedAllocator<size_t>> isNeighborSiteOwnedArr;

    // we require arrays for forces for sites to maintain SIMD in site-site calculations
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceX(siteCountNeighbors, 0.);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceY(siteCountNeighbors, 0.);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceZ(siteCountNeighbors, 0.);

    // pre-reserve arrays
    exactNeighborSitePositionsX.reserve(siteCountNeighbors);
    exactNeighborSitePositionsY.reserve(siteCountNeighbors);
    exactNeighborSitePositionsZ.reserve(siteCountNeighbors);

    if constexpr (useMixing) {
      siteTypesNeighbors.reserve(siteCountNeighbors);
    }

    if constexpr (calculateGlobals) {
      isNeighborSiteOwnedArr.reserve(siteCountNeighbors);
    }

    const auto rotatedSitePositionsPrime =
        useMixing ? autopas::utils::quaternion::rotateVectorOfPositions(
                        {q0ptr[indexFirst], q1ptr[indexFirst], q2ptr[indexFirst], q3ptr[indexFirst]},
                        _PPLibrary->getSitePositions(typeptr[indexFirst]))
                  : autopas::utils::quaternion::rotateVectorOfPositions(
                        {q0ptr[indexFirst], q1ptr[indexFirst], q2ptr[indexFirst], q3ptr[indexFirst]}, _sitePositionsLJ);

    const auto siteTypesPrime = _PPLibrary->getSiteTypes(typeptr[indexFirst]);  // todo make this work for no mixing

    std::vector<std::array<double, 3>> rotatedSitePositions;
    for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
      const size_t neighborMolIndex = neighborList[neighborMol];
      if constexpr (useMixing) {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[neighborMolIndex], q1ptr[neighborMolIndex], q2ptr[neighborMolIndex], q3ptr[neighborMolIndex]},
            _PPLibrary->getSitePositions(typeptr[neighborMolIndex]));
      } else {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[neighborMolIndex], q1ptr[neighborMolIndex], q2ptr[neighborMolIndex], q3ptr[neighborMolIndex]},
            _sitePositionsLJ);
      }

      size_t localSiteCount = useMixing ? _PPLibrary->getNumSites(typeptr[neighborMolIndex]) : _sitePositionsLJ.size();
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

    double potentialEnergyAccumulator = 0;
    std::array<double, 3> virialAccumulator = {0., 0., 0.};

    // -- main force calculation --

    // - calculate mol mask -
    std::vector<double, autopas::AlignedAllocator<double>> molMask;
    molMask.reserve(neighborListSize);

    const __m256d xPosFirst = _mm256_broadcast_sd(&xptr[indexFirst]);
    const __m256d yPosFirst = _mm256_broadcast_sd(&yptr[indexFirst]);
    const __m256d zPosFirst = _mm256_broadcast_sd(&zptr[indexFirst]);

    // WIP vectorized version

    for (size_t neighborMol = 0; neighborMol < neighborListSize; neighborMol += vecLength) {
      const size_t rest = neighborListSize - neighborMol;
      const bool remainderCase = rest < vecLength;
      __m256i remainderMask = remainderCase ? _masks[rest - 1] : _mm256_set1_epi64x(-1);
      const __m256i neighborMolIndex =
          autopas::utils::avx::load_epi64(remainderCase, &neighborList[neighborMol], remainderMask);

      const __m256d xposB = autopas::utils::avx::gather_pd(remainderCase, xptr, neighborMolIndex, remainderMask);
      const __m256d yposB = autopas::utils::avx::gather_pd(remainderCase, yptr, neighborMolIndex, remainderMask);
      const __m256d zposB = autopas::utils::avx::gather_pd(remainderCase, zptr, neighborMolIndex, remainderMask);

      const __m256d displacementCoMX = _mm256_sub_pd(xPosFirst, xposB);
      const __m256d displacementCoMY = _mm256_sub_pd(yPosFirst, yposB);
      const __m256d displacementCoMZ = _mm256_sub_pd(zPosFirst, zposB);

      const __m256d distanceSquaredCoMX = _mm256_mul_pd(displacementCoMX, displacementCoMX);
      const __m256d distanceSquaredCoMY = _mm256_mul_pd(displacementCoMY, displacementCoMY);
      const __m256d distanceSquaredCoMZ = _mm256_mul_pd(displacementCoMZ, displacementCoMZ);

      const __m256d distanceSquaredCoM =
          _mm256_add_pd(_mm256_add_pd(distanceSquaredCoMX, distanceSquaredCoMY), distanceSquaredCoMZ);

      const __m256d cutoffMask = _mm256_cmp_pd(distanceSquaredCoM, _cutoffSquared, _CMP_LE_OS);
      const __m256i ownedStateB =
          remainderCase
              ? _mm256_mask_i64gather_epi64(_mm256_set1_epi64x(0),
                                            reinterpret_cast<const long long int *>(ownedStatePtr), neighborMolIndex,
                                            remainderMask, 8)
              : _mm256_i64gather_epi64(reinterpret_cast<const long long int *>(ownedStatePtr), neighborMolIndex, 8);
      const __m256d dummyMask =
          _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateB), _zero, _CMP_NEQ_OS);  // Assuming that dummy = 0
      const __m256d totalMask = _mm256_and_pd(cutoffMask, dummyMask);

      if (remainderCase) {
        _mm256_maskstore_pd((&molMask[neighborMol]), remainderMask, totalMask);
      } else {
        _mm256_storeu_pd((&molMask[neighborMol]), totalMask);
      }
    }

    //     generate mask for each site from molecular mask
    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteVector;
    siteVector.reserve(siteCountNeighbors);

    if constexpr (useSiteMasks) {
      for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
        const auto neighborMolIndex = neighborList[neighborMol];  // index of neighbor mol in soa

        const double xPosB = xptr[neighborMolIndex];
        const double yPosB = yptr[neighborMolIndex];
        const double zPosB = zptr[neighborMolIndex];

        // calculate displacement
        const double displacementCoMX = xptr[indexFirst] - xPosB;
        const double displacementCoMY = yptr[indexFirst] - yPosB;
        const double displacementCoMZ = zptr[indexFirst] - zPosB;

        // calculate distance squared
        const double distanceSquaredCoMX = displacementCoMX * displacementCoMX;
        const double distanceSquaredCoMY = displacementCoMY * displacementCoMY;
        const double distanceSquaredCoMZ = displacementCoMZ * displacementCoMZ;

        const double distanceSquaredCoM = distanceSquaredCoMX + distanceSquaredCoMY + distanceSquaredCoMZ;

        bool cutoffCondition = distanceSquaredCoM <= _cutoffSquaredAoS;
        bool dummyCondition = ownedStatePtr[neighborMolIndex] != OwnershipState::dummy;

        size_t mask = cutoffCondition and dummyCondition ? std::numeric_limits<size_t>::max() : 0;

        for (size_t siteB = 0; siteB < _PPLibrary->getNumSites(typeptr[neighborMolIndex]); ++siteB) {
          siteVector.emplace_back(mask);
        }
      }
    } else {
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
    }

    // - actual LJ calculation -

    std::array<double, 3> forceAccumulator = {0., 0., 0.};
    std::array<double, 3> torqueAccumulator = {0., 0., 0.};

    for (size_t primeSite = 0; primeSite < siteCountMolPrime; ++primeSite) {
      const size_t siteTypePrime = siteTypesPrime[primeSite];
      const std::array<double, 3> rotatedSitePositionPrime = rotatedSitePositionsPrime[primeSite];
      const std::array<double, 3> exactSitePositionPrime = autopas::utils::ArrayMath::add(
          rotatedSitePositionPrime, {xptr[indexFirst], yptr[indexFirst], zptr[indexFirst]});

      SoAKernel<newton3, true, useSiteMasks>(
          siteVector, siteTypesNeighbors, exactNeighborSitePositionsX, exactNeighborSitePositionsY,
          exactNeighborSitePositionsZ, siteForceX, siteForceY, siteForceZ, isNeighborSiteOwnedArr, ownedStatePrime,
          siteTypePrime, exactSitePositionPrime, rotatedSitePositionPrime, forceAccumulator, torqueAccumulator,
          potentialEnergyAccumulator, virialAccumulator);
    }
    // Add forces to prime mol
    fxptr[indexFirst] += forceAccumulator[0];
    fyptr[indexFirst] += forceAccumulator[1];
    fzptr[indexFirst] += forceAccumulator[2];
    txptr[indexFirst] += torqueAccumulator[0];
    typtr[indexFirst] += torqueAccumulator[1];
    tzptr[indexFirst] += torqueAccumulator[2];

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
              _sitePositionsLJ);
          for (size_t site = 0; site < _sitePositionsLJ.size(); ++site) {
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
      // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2here.
      const auto newton3Factor = newton3 ? .5 : 1.;

      _aosThreadData[threadNum].potentialEnergySum += potentialEnergyAccumulator * newton3Factor;
      _aosThreadData[threadNum].virialSum[0] += virialAccumulator[0] * newton3Factor;
      _aosThreadData[threadNum].virialSum[1] += virialAccumulator[1] * newton3Factor;
      _aosThreadData[threadNum].virialSum[2] += virialAccumulator[2] * newton3Factor;
    }
  }

  /** Constructs site vectors for the given SoA.
   * @param soa the SoA to construct the site vectors for.
   * @note This function should be called exactly once per iteration
   */
  inline void buildSiteVectors(SoAView<SoAArraysType> soa) {
    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict q0ptr = soa.template begin<Particle::AttributeNames::quaternion0>();
    const auto *const __restrict q1ptr = soa.template begin<Particle::AttributeNames::quaternion1>();
    const auto *const __restrict q2ptr = soa.template begin<Particle::AttributeNames::quaternion2>();
    const auto *const __restrict q3ptr = soa.template begin<Particle::AttributeNames::quaternion3>();

    auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    std::vector<std::array<double, 3>> rotatedSitePositions;
    size_t localSiteCount = 0;
    size_t totalSiteCount = 0;

    for (size_t mol = 0; mol < soa.getNumberOfParticles(); mol++) {
      _molToSiteMap.push_back(totalSiteCount);
      if constexpr (useMixing) {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, _PPLibrary->getSitePositions(typeptr[mol]));
        localSiteCount = _PPLibrary->getNumSites(typeptr[mol]);
      } else {
        rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, _sitePositionsLJ);
        localSiteCount = _sitePositionsLJ.size();
      }
      totalSiteCount += localSiteCount;

      for (size_t site = 0; site < localSiteCount; ++site) {
        _exactSitePositionsX.push_back(rotatedSitePositions[site][0] + xptr[mol]);
        _exactSitePositionsY.push_back(rotatedSitePositions[site][1] + yptr[mol]);
        _exactSitePositionsZ.push_back(rotatedSitePositions[site][2] + zptr[mol]);
        _siteToMolMap.push_back(mol);
        if constexpr (calculateGlobals) {
          _siteOwnershipVector.push_back(ownedStatePtr[mol] == OwnershipState::owned);
        }
        if constexpr (useMixing) {
          _siteTypesVector.push_back(_PPLibrary->getSiteTypes(typeptr[mol])[site]);
        }
      }
    }
  }

  // TODO: DOCUMENTATION
  template <bool newton3, bool calculateTorque, bool masked>
  inline void SoAKernel(const std::vector<size_t, autopas::AlignedAllocator<size_t>> &siteVector,
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
    __m256d sigmaSquared = _mm256_set1_pd(_sigmaSquaredAoS);
    __m256d epsilon24 = _mm256_set1_pd(_epsilon24AoS);
    __m256d shift6 = applyShift ? _mm256_set1_pd(_shift6AoS) : _zero;
    const double *const __restrict mixingPtr = useMixing ? _PPLibrary->getMixingDataPtr(siteTypeA, 0) : nullptr;

    const __m256d exactSitePositionsAX = _mm256_set1_pd(exactSitePositionA[0]);
    const __m256d exactSitePositionsAY = _mm256_set1_pd(exactSitePositionA[1]);
    const __m256d exactSitePositionsAZ = _mm256_set1_pd(exactSitePositionA[2]);

    const __m256d rotatedSitePositionsAX = _mm256_set1_pd(rotatedSitePositionA[0]);
    const __m256d rotatedSitePositionsAY = _mm256_set1_pd(rotatedSitePositionA[1]);
    const __m256d rotatedSitePositionsAZ = _mm256_set1_pd(rotatedSitePositionA[2]);

    // sums used for siteA
    __m256d forceSumX = _zero;
    __m256d forceSumY = _zero;
    __m256d forceSumZ = _zero;

    __m256d torqueSumX = _zero;
    __m256d torqueSumY = _zero;
    __m256d torqueSumZ = _zero;

    // Globals
    __m256d potentialEnergySum = _zero;
    __m256d virialSumX = _zero;
    __m256d virialSumY = _zero;
    __m256d virialSumZ = _zero;

    __m256i sites = _mm256_setzero_si256();

    for (size_t siteVectorIndex = 0; siteVectorIndex < siteVector.size(); siteVectorIndex += vecLength) {
      const size_t remainder = siteVector.size() - siteVectorIndex;
      const bool remainderCase = remainder < vecLength;
      const __m256i remainderMask = remainderCase ? _masks[remainder - 1] : _mm256_set1_epi64x(-1);

      if constexpr (masked) {
        sites = autopas::utils::avx::load_epi64(remainderCase, &siteVector[siteVectorIndex], remainderMask);
        if (_mm256_testz_si256(sites, sites)) {
          continue;
        }
      } else {
        sites = autopas::utils::avx::load_epi64(remainderCase, &siteVector[siteVectorIndex], remainderMask);
      }

      if constexpr (useMixing) {
        __m256i siteTypeIndices =
            masked
                ? autopas::utils::avx::load_epi64(remainderCase, &siteTypesB[siteVectorIndex + offset], remainderMask)
                : autopas::utils::avx::gather_epi64(remainderCase, &siteTypesB[offset], sites, remainderMask);
        siteTypeIndices = _mm256_mul_epu32(siteTypeIndices, _mm256_set1_epi64x(3));  // This could maybe be problematic!

        epsilon24 = _mm256_i64gather_pd(mixingPtr, siteTypeIndices, 8);
        sigmaSquared = _mm256_i64gather_pd(mixingPtr + 1, siteTypeIndices, 8);
        if constexpr (applyShift) {
          shift6 = _mm256_i64gather_pd(mixingPtr + 2, siteTypeIndices, 8);
        }
      }

      __m256d exactSitePositionsBX;
      __m256d exactSitePositionsBY;
      __m256d exactSitePositionsBZ;

      if constexpr (masked) {
        exactSitePositionsBX =
            autopas::utils::avx::load_pd(remainderCase, &exactSitePositionX[siteVectorIndex + offset], remainderMask);
        exactSitePositionsBY =
            autopas::utils::avx::load_pd(remainderCase, &exactSitePositionY[siteVectorIndex + offset], remainderMask);
        exactSitePositionsBZ =
            autopas::utils::avx::load_pd(remainderCase, &exactSitePositionZ[siteVectorIndex + offset], remainderMask);
      } else {
        exactSitePositionsBX =
            autopas::utils::avx::gather_pd(remainderCase, &exactSitePositionX[offset], sites, remainderMask);
        exactSitePositionsBY =
            autopas::utils::avx::gather_pd(remainderCase, &exactSitePositionY[offset], sites, remainderMask);
        exactSitePositionsBZ =
            autopas::utils::avx::gather_pd(remainderCase, &exactSitePositionZ[offset], sites, remainderMask);
      }

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
      const __m256d castedMask = _mm256_castsi256_pd(sites);
      __m256d scalarMultiple = masked ? _mm256_and_pd(castedMask, scalar) : scalar;
      scalarMultiple = _mm256_and_pd(_mm256_castsi256_pd(remainderMask), scalarMultiple);

      // calculate forces
      const __m256d forceX = _mm256_mul_pd(scalarMultiple, displacementX);
      const __m256d forceY = _mm256_mul_pd(scalarMultiple, displacementY);
      const __m256d forceZ = _mm256_mul_pd(scalarMultiple, displacementZ);

      forceSumX = _mm256_add_pd(forceSumX, forceX);
      forceSumY = _mm256_add_pd(forceSumY, forceY);
      forceSumZ = _mm256_add_pd(forceSumZ, forceZ);

      // calculate torques
      if constexpr (calculateTorque) {
        const __m256d torqueAX =
            _mm256_fmsub_pd(rotatedSitePositionsAY, forceZ, _mm256_mul_pd(rotatedSitePositionsAZ, forceY));
        const __m256d torqueAY =
            _mm256_fmsub_pd(rotatedSitePositionsAZ, forceX, _mm256_mul_pd(rotatedSitePositionsAX, forceZ));
        const __m256d torqueAZ =
            _mm256_fmsub_pd(rotatedSitePositionsAX, forceY, _mm256_mul_pd(rotatedSitePositionsAY, forceX));

        torqueSumX = _mm256_add_pd(torqueSumX, torqueAX);
        torqueSumY = _mm256_add_pd(torqueSumY, torqueAY);
        torqueSumZ = _mm256_add_pd(torqueSumZ, torqueAZ);
      }

      // N3L ( total molecular forces + torques to be determined later )
      if constexpr (newton3) {
        if constexpr (masked) {
          __m256d forceSumBX =
              autopas::utils::avx::load_pd(remainderCase, &siteForceX[siteVectorIndex + offset], remainderMask);
          __m256d forceSumBY =
              autopas::utils::avx::load_pd(remainderCase, &siteForceY[siteVectorIndex + offset], remainderMask);
          __m256d forceSumBZ =
              autopas::utils::avx::load_pd(remainderCase, &siteForceZ[siteVectorIndex + offset], remainderMask);

          forceSumBX = _mm256_sub_pd(forceSumBX, forceX);
          forceSumBY = _mm256_sub_pd(forceSumBY, forceY);
          forceSumBZ = _mm256_sub_pd(forceSumBZ, forceZ);

          remainderCase ? _mm256_maskstore_pd(&siteForceX[siteVectorIndex + offset], remainderMask, forceSumBX)
                        : _mm256_storeu_pd(&siteForceX[siteVectorIndex + offset], forceSumBX);
          remainderCase ? _mm256_maskstore_pd(&siteForceY[siteVectorIndex + offset], remainderMask, forceSumBY)
                        : _mm256_storeu_pd(&siteForceY[siteVectorIndex + offset], forceSumBY);
          remainderCase ? _mm256_maskstore_pd(&siteForceZ[siteVectorIndex + offset], remainderMask, forceSumBZ)
                        : _mm256_storeu_pd(&siteForceZ[siteVectorIndex + offset], forceSumBZ);
        } else {
          __m256d forceSumBX = autopas::utils::avx::gather_pd(remainderCase, &siteForceX[offset], sites, remainderMask);
          __m256d forceSumBY = autopas::utils::avx::gather_pd(remainderCase, &siteForceY[offset], sites, remainderMask);
          __m256d forceSumBZ = autopas::utils::avx::gather_pd(remainderCase, &siteForceZ[offset], sites, remainderMask);
          forceSumBX = _mm256_sub_pd(forceSumBX, forceX);
          forceSumBY = _mm256_sub_pd(forceSumBY, forceY);
          forceSumBZ = _mm256_sub_pd(forceSumBZ, forceZ);

          autopas::utils::avx::scatter_pd(remainderCase, &siteForceX[offset], sites, remainderMask, forceSumBX);
          autopas::utils::avx::scatter_pd(remainderCase, &siteForceY[offset], sites, remainderMask, forceSumBY);
          autopas::utils::avx::scatter_pd(remainderCase, &siteForceZ[offset], sites, remainderMask, forceSumBZ);
        }
      }

      // globals
      if constexpr (calculateGlobals) {
        const __m256d virialX = _mm256_mul_pd(displacementX, forceX);
        const __m256d virialY = _mm256_mul_pd(displacementY, forceY);
        const __m256d virialZ = _mm256_mul_pd(displacementZ, forceZ);

        const __m256d potentialEnergy6 =
            _mm256_fmadd_pd(epsilon24, lj12m6, shift6);  // FMA may not be supported on all CPUs
        __m256d potentialEnergy6Masked = _mm256_and_pd(castedMask, potentialEnergy6);
        potentialEnergy6Masked = _mm256_and_pd(_mm256_castsi256_pd(remainderMask), potentialEnergy6Masked);

        __m256i ownedStateA4 = _mm256_set1_epi64x(static_cast<int64_t>(ownedStateA));
        __m256d ownedMaskA =
            _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateA4), _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
        __m256d energyFactor = _mm256_blendv_pd(_zero, _one, ownedMaskA);
        if constexpr (newton3) {
          const __m256i ownedStateB =
              masked ? autopas::utils::avx::load_epi64(remainderCase, &siteOwnership[offset + siteVectorIndex],
                                                       remainderMask)
                     : autopas::utils::avx::gather_epi64(remainderCase, &siteOwnership[offset], sites, remainderMask);
          __m256d ownedMaskB =
              _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateB), _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
          energyFactor = _mm256_add_pd(energyFactor, _mm256_blendv_pd(_zero, _one, ownedMaskB));
        }
        potentialEnergySum = _mm256_fmadd_pd(energyFactor, potentialEnergy6Masked, potentialEnergySum);
        virialSumX = _mm256_fmadd_pd(energyFactor, virialX, virialSumX);
        virialSumY = _mm256_fmadd_pd(energyFactor, virialY, virialSumY);
        virialSumZ = _mm256_fmadd_pd(energyFactor, virialZ, virialSumZ);
      }
    }

    forceAccumulator[0] += utils::avx::horizontalSum(forceSumX);
    forceAccumulator[1] += utils::avx::horizontalSum(forceSumY);
    forceAccumulator[2] += utils::avx::horizontalSum(forceSumZ);

    if constexpr (calculateTorque) {
      torqueAccumulator[0] += utils::avx::horizontalSum(torqueSumX);
      torqueAccumulator[1] += utils::avx::horizontalSum(torqueSumY);
      torqueAccumulator[2] += utils::avx::horizontalSum(torqueSumZ);
    }

    if constexpr (calculateGlobals) {
      potentialEnergyAccumulator += utils::avx::horizontalSum(potentialEnergySum);
      virialAccumulator[0] += utils::avx::horizontalSum(virialSumX);
      virialAccumulator[1] += utils::avx::horizontalSum(virialSumY);
      virialAccumulator[2] += utils::avx::horizontalSum(virialSumZ);
    }
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

 public:
  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
   * TODO: This is just a copy of the AoSThreadData from the LJFunctor. This should probably be refactored.
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

/**
 * Prevent usage of non-multisite molecule
 * @tparam applyShift
 * @tparam useMixing
 * @tparam useNewton3
 * @tparam calculateGlobals
 * @tparam relevantForTuning
 */
template <bool applyShift, bool useMixing, autopas::FunctorN3Modes useNewton3, bool calculateGlobals,
          bool relevantForTuning>
class LJMultisiteFunctorAVX<autopas::MoleculeLJ, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>
    : public autopas::Functor<autopas::MoleculeLJ,
                              LJMultisiteFunctorAVX<autopas::MoleculeLJ, applyShift, useMixing, useNewton3,
                                                    calculateGlobals, relevantForTuning>> {
 public:
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename autopas::MoleculeLJ::SoAArraysType;

  using SoAFloatPrecision = typename autopas::MoleculeLJ::ParticleSoAFloatPrecision;

  /**
   * Delete Default constructor
   */
  LJMultisiteFunctorAVX() = delete;

 private:
  /**
   * Internal (actually used) constructor
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit LJMultisiteFunctorAVX(SoAFloatPrecision cutoff, void * /*dummy*/)
      : autopas::Functor<autopas::MoleculeLJ, LJMultisiteFunctorAVX<autopas::MoleculeLJ, applyShift, useMixing,
                                                                    useNewton3, calculateGlobals, relevantForTuning>>(
            cutoff) {
    autopas::utils::ExceptionHandler::exception(
        "LJMultisiteFunctor can not be used with MoleculeLJ. Use a MultisiteMoleculeLJ instead.");
  }

 public:
  /**
   * Constructor for Functor with particle mixing disabled. setParticleProperties() must be called.
   * @param cutoff
   */
  explicit LJMultisiteFunctorAVX(double cutoff) : LJMultisiteFunctorAVX(cutoff, nullptr) {
    autopas::utils::ExceptionHandler::exception(
        "LJMultisiteFunctor can not be used with MoleculeLJ. Use a MultisiteMoleculeLJ instead.");
  }

  /**
   * Constructor for Functor with particle mixing enabled.
   * @param cutoff
   * @param particlePropertiesLibrary Library used to look up the properties of each type of particle e.g. sigma,
   * epsilon, shift.
   */
  explicit LJMultisiteFunctorAVX(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJMultisiteFunctorAVX(cutoff, nullptr) {
    autopas::utils::ExceptionHandler::exception(
        "LJMultisiteFunctor can not be used with MoleculeLJ. Use a MultisiteMoleculeLJ instead.");
  }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  constexpr static bool getMixing() { return useMixing; }

  double getVirial() { return 0; }

  /**
   * Functor for AoS. Simply loops over the sites of two particles/molecules to calculate force.
   * @param particleA Particle i
   * @param particleB Particle j
   * @param newton3 Flag for if newton3 is used.
   */
  void AoSFunctor(autopas::MoleculeLJ &particleA, autopas::MoleculeLJ &particleB, bool newton3) final {
    autopas::utils::ExceptionHandler::exception(
        "LJMultisiteFunctorAVX can not be used with MoleculeLJ. Use a MultisiteMoleculeLJ instead.");
  }

  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    autopas::utils::ExceptionHandler::exception(
        "LJMultisiteFunctorAVX can not be used with MoleculeLJ. Use a MultisiteMoleculeLJ instead.");
  }

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {
    autopas::utils::ExceptionHandler::exception(
        "LJMultisiteFunctorAVX can not be used with MoleculeLJ. Use a MultisiteMoleculeLJ instead.");
  }

  void setParticleProperties(SoAFloatPrecision epsilon24, SoAFloatPrecision sigmaSquared) {
    autopas::utils::ExceptionHandler::exception(
        "LJMultisiteFunctorAVX can not be used with MoleculeLJ. Use a MultisiteMoleculeLJ instead.");
  }

  static unsigned long getNumFlopsPerKernelCall(bool newton3, size_t numA, size_t numB) { return 0ul; }

  void initTraversal() final {
    autopas::utils::ExceptionHandler::exception(
        "LJMultisiteFunctorAVX can not be used with MoleculeLJ. Use a MultisiteMoleculeLJ instead.");
  }

  void endTraversal(bool newton3) final {
    autopas::utils::ExceptionHandler::exception(
        "LJMultisiteFunctorAVX can not be used with MoleculeLJ. Use a MultisiteMoleculeLJ instead.");
  }
};
}  // namespace autopas