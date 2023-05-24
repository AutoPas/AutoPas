/**
 * @file LJMultisiteFunctor.h
 * @date 21/02/2022
 * @author S. Newcome
 */

#pragma once

#include "MoleculeLJ.h"
#include "MultisiteMoleculeLJ.h"
#include "ParticlePropertiesLibrary.h"
#include "autopas/pairwiseFunctors/Functor.h"
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
 * @tparam useNewton3 Switch for the functor to support newton3 on, off, or both. See FunctorN3Nodes for possible
 * values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool relevantForTuning = true>
class LJMultisiteFunctor
    : public autopas::Functor<Particle, LJMultisiteFunctor<Particle, applyShift, useMixing, useNewton3,
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
  LJMultisiteFunctor() = delete;

 private:
  /**
   * Internal (actually used) constructor
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit LJMultisiteFunctor(SoAFloatPrecision cutoff, void * /*dummy*/)
      : autopas::Functor<Particle, LJMultisiteFunctor<Particle, applyShift, useMixing, useNewton3, calculateGlobals,
                                                      relevantForTuning>>(cutoff),
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
   * Constructor for Functor with particle mixing disabled. setParticleProperties() must be called.
   * @note Only to be used with mixing == false
   * @param cutoff
   */
  explicit LJMultisiteFunctor(double cutoff) : LJMultisiteFunctor(cutoff, nullptr) {
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
  explicit LJMultisiteFunctor(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJMultisiteFunctor(cutoff, nullptr) {
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
            useMixing ? _PPLibrary->getMixingSigmaSquared(siteIdsA[i], siteIdsB[j]) : _sigmaSquared;
        const auto epsilon24 = useMixing ? _PPLibrary->getMixing24Epsilon(siteIdsA[i], siteIdsB[j]) : _epsilon24;
        const auto shift6 =
            applyShift ? (useMixing ? _PPLibrary->getMixingShift6(siteIdsA[i], siteIdsB[j]) : _shift6) : 0;

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
          const auto virial = newton3 ? (displacement * force) * 0.5 : displacement * force;

          const auto threadNum = autopas::autopas_get_thread_num();

          if (particleA.isOwned()) {
            _aosThreadData[threadNum].potentialEnergySum += potentialEnergy6;
            _aosThreadData[threadNum].virialSum = _aosThreadData[threadNum].virialSum + virial;
          }
          if (newton3 and particleB.isOwned()) {
            _aosThreadData[threadNum].potentialEnergySum += potentialEnergy6;
            _aosThreadData[threadNum].virialSum = _aosThreadData[threadNum].virialSum + virial;
          }
        }
      }
    }
  }

  /**
   * @copydoc autopas::Functor::SoAFunctorSingle()
   * This functor will always use a newton3 like traversing of the soa, however, it still needs to know about newton3
   * to use it correctly for the global values.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
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

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> sigmaSquareds;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> epsilon24s;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> shift6s;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionZ;

    // we require arrays for forces for sites to maintain SIMD in site-site calculations
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceZ;

    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteTypes;
    std::vector<char, autopas::AlignedAllocator<char>> isSiteOwned;

    const SoAFloatPrecision const_sigmaSquared = _sigmaSquared;
    const SoAFloatPrecision const_epsilon24 = _epsilon24;
    const SoAFloatPrecision const_shift6 = _shift6;

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
    if constexpr (useMixing) {
      size_t siteIndex = 0;
      for (size_t mol = 0; mol < soa.getNumberOfParticles(); ++mol) {
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, _PPLibrary->getSitePositions(typeptr[mol]));
        const auto siteTypesOfMol = _PPLibrary->getSiteTypes(typeptr[mol]);

        for (size_t site = 0; site < _PPLibrary->getNumSites(typeptr[mol]); ++site) {
          exactSitePositionX[siteIndex] = rotatedSitePositions[site][0] + xptr[mol];
          exactSitePositionY[siteIndex] = rotatedSitePositions[site][1] + yptr[mol];
          exactSitePositionZ[siteIndex] = rotatedSitePositions[site][2] + zptr[mol];
          siteTypes[siteIndex] = siteTypesOfMol[site];
          siteForceX[siteIndex] = 0.;
          siteForceY[siteIndex] = 0.;
          siteForceZ[siteIndex] = 0.;
          if (calculateGlobals) {
            isSiteOwned[siteIndex] = ownedStatePtr[mol] == autopas::OwnershipState::owned;
          }
          ++siteIndex;
        }
      }
    } else {
      size_t siteIndex = 0;
      for (size_t mol = 0; mol < soa.getNumberOfParticles(); mol++) {
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, const_unrotatedSitePositions);
        for (size_t site = 0; site < const_unrotatedSitePositions.size(); ++site) {
          exactSitePositionX[siteIndex] = rotatedSitePositions[site][0] + xptr[mol];
          exactSitePositionY[siteIndex] = rotatedSitePositions[site][1] + yptr[mol];
          exactSitePositionZ[siteIndex] = rotatedSitePositions[site][2] + zptr[mol];
          siteForceX[siteIndex] = 0.;
          siteForceY[siteIndex] = 0.;
          siteForceZ[siteIndex] = 0.;
          if (calculateGlobals) {
            isSiteOwned[siteIndex] = ownedStatePtr[mol] == autopas::OwnershipState::owned;
          }
          ++siteIndex;
        }
      }
    }

    // main force calculation loop
    size_t siteIndexMolA = 0;  // index of first site in molA
    for (size_t molA = 0; molA < soa.getNumberOfParticles(); ++molA) {
      const size_t noSitesInMolA = useMixing ? _PPLibrary->getNumSites(typeptr[molA])
                                             : const_unrotatedSitePositions.size();  // Number of sites in molecule A

      const auto ownedStateA = ownedStatePtr[molA];
      if (ownedStateA == autopas::OwnershipState::dummy) {
        siteIndexMolA += noSitesInMolA;
        continue;
      }

      const size_t siteIndexMolB = siteIndexMolA + noSitesInMolA;                    // index of first site in molB
      const size_t noSitesB = (siteCount - siteIndexMolB);  // Number of sites in molecules that A interacts with

      // create mask over every mol 'above' molA  (char to keep arrays aligned)
      std::vector<char, autopas::AlignedAllocator<char>> molMask;
      molMask.reserve(soa.getNumberOfParticles() - (molA + 1));

#pragma omp simd
      for (size_t molB = molA + 1; molB < soa.getNumberOfParticles(); ++molB) {
        const auto ownedStateB = ownedStatePtr[molB];

        const auto displacementCoMX = xptr[molA] - xptr[molB];
        const auto displacementCoMY = yptr[molA] - yptr[molB];
        const auto displacementCoMZ = zptr[molA] - zptr[molB];

        const auto distanceSquaredCoMX = displacementCoMX * displacementCoMX;
        const auto distanceSquaredCoMY = displacementCoMY * displacementCoMY;
        const auto distanceSquaredCoMZ = displacementCoMZ * displacementCoMZ;

        const auto distanceSquaredCoM = distanceSquaredCoMX + distanceSquaredCoMY + distanceSquaredCoMZ;

        // mask sites of molecules beyond cutoff or if molecule is a dummy
        molMask[molB - (molA + 1)] =
            distanceSquaredCoM <= cutoffSquared and ownedStateB != autopas::OwnershipState::dummy;
      }

      // generate mask for each site in the mols 'above' molA from molecular mask
      std::vector<char, autopas::AlignedAllocator<char>> siteMask;
      siteMask.reserve(noSitesB);

      for (size_t molB = molA + 1; molB < soa.getNumberOfParticles(); ++molB) {
        for (size_t siteB = 0; siteB < _PPLibrary->getNumSites(typeptr[molB]); ++siteB) {
          siteMask.emplace_back(molMask[molB - (molA + 1)]);
        }
      }

      // calculate LJ forces
      for (size_t siteA = siteIndexMolA; siteA < siteIndexMolB; ++siteA) {
        if (useMixing) {
          // preload sigmas, epsilons, and shifts
          sigmaSquareds.reserve(noSitesB);
          epsilon24s.reserve(noSitesB);
          if constexpr (applyShift) {
            shift6s.reserve(noSitesB);
          }

          for (size_t siteB = 0; siteB < siteCount - (siteIndexMolB); ++siteB) {
            const auto mixingData = _PPLibrary->getMixingData(siteTypes[siteA], siteTypes[siteIndexMolB + siteB]);
            sigmaSquareds[siteB] = mixingData.sigmaSquared;
            epsilon24s[siteB] = mixingData.epsilon24;
            if (applyShift) {
              shift6s[siteB] = mixingData.shift6;
            }
          }
        }
        // sums used for siteA
        SoAFloatPrecision forceSumX = 0.;
        SoAFloatPrecision forceSumY = 0.;
        SoAFloatPrecision forceSumZ = 0.;
        SoAFloatPrecision torqueSumX = 0.;
        SoAFloatPrecision torqueSumY = 0.;
        SoAFloatPrecision torqueSumZ = 0.;

#pragma omp simd reduction (+ : forceSumX, forceSumY, forceSumZ, torqueSumX, torqueSumY, torqueSumZ, potentialEnergySum, virialSumX, virialSumY, virialSumZ)
        for (size_t siteB = 0; siteB < noSitesB; ++siteB) {
          const size_t globalSiteBIndex = siteB + siteIndexMolB;

          const SoAFloatPrecision sigmaSquared = useMixing ? sigmaSquareds[siteB] : const_sigmaSquared;
          const SoAFloatPrecision epsilon24 = useMixing ? epsilon24s[siteB] : const_epsilon24;
          const SoAFloatPrecision shift6 = applyShift ? (useMixing ? shift6s[siteB] : const_shift6) : 0;

          const auto isSiteBOwned = !calculateGlobals || isSiteOwned[globalSiteBIndex];

          const auto displacementX = exactSitePositionX[siteA] - exactSitePositionX[globalSiteBIndex];
          const auto displacementY = exactSitePositionY[siteA] - exactSitePositionY[globalSiteBIndex];
          const auto displacementZ = exactSitePositionZ[siteA] - exactSitePositionZ[globalSiteBIndex];

          const auto distanceSquaredX = displacementX * displacementX;
          const auto distanceSquaredY = displacementY * displacementY;
          const auto distanceSquaredZ = displacementZ * displacementZ;

          const auto distanceSquared = distanceSquaredX + distanceSquaredY + distanceSquaredZ;

          const auto invDistSquared = 1. / distanceSquared;
          const auto lj2 = sigmaSquared * invDistSquared;
          const auto lj6 = lj2 * lj2 * lj2;
          const auto lj12 = lj6 * lj6;
          const auto lj12m6 = lj12 - lj6;
          const auto scalarMultiple = siteMask[siteB] ? epsilon24 * (lj12 + lj12m6) * invDistSquared : 0.;

          // calculate forces
          const auto forceX = scalarMultiple * displacementX;
          const auto forceY = scalarMultiple * displacementY;
          const auto forceZ = scalarMultiple * displacementZ;

          forceSumX += forceX;
          forceSumY += forceY;
          forceSumZ += forceZ;

          // newton's third law
          siteForceX[globalSiteBIndex] -= forceX;
          siteForceY[globalSiteBIndex] -= forceY;
          siteForceZ[globalSiteBIndex] -= forceZ;

          if constexpr (calculateGlobals) {
            const auto virialX = displacementX * forceX;
            const auto virialY = displacementY * forceY;
            const auto virialZ = displacementZ * forceZ;
            const auto potentialEnergy6 = siteMask[siteB] ? (epsilon24 * lj12m6 + shift6) : 0.;

            // Add to the potential energy sum for each particle which is owned.
            // This results in obtaining 12 * the potential energy for the SoA.
            const auto ownershipMask =
                (ownedStateA == autopas::OwnershipState::owned ? 1. : 0.) + (isSiteBOwned ? 1. : 0.);
            potentialEnergySum += potentialEnergy6 * ownershipMask;
            virialSumX += virialX * ownershipMask;
            virialSumY += virialY * ownershipMask;
            virialSumZ += virialZ * ownershipMask;
          }
        }
        // sum forces on single site in mol A
        siteForceX[siteA] += forceSumX;
        siteForceY[siteA] += forceSumY;
        siteForceZ[siteA] += forceSumZ;
      }
      siteIndexMolA += noSitesInMolA;
    }

    // reduce the forces on individual sites to forces & torques on whole molecules.
    if constexpr (useMixing) {
      size_t siteIndex = 0;
      for (size_t mol = 0; mol < soa.getNumberOfParticles(); ++mol) {
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
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
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
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
      const auto threadNum = autopas::autopas_get_thread_num();
      // SoAFunctorSingle obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
      // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.
      const auto newton3Factor = newton3 ? .5 : 1.;

      _aosThreadData[threadNum].potentialEnergySum += potentialEnergySum * newton3Factor;
      _aosThreadData[threadNum].virialSum[0] += virialSumX * newton3Factor;
      _aosThreadData[threadNum].virialSum[1] += virialSumY * newton3Factor;
      _aosThreadData[threadNum].virialSum[2] += virialSumZ * newton3Factor;
    }
  }
  /**
   * @copydoc autopas::Functor::SoAFunctorPair()
   */
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {
    if (newton3) {
      SoAFunctorPairImpl<true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1, soa2);
    }
  }

  // clang-format off
  /**
   * @copydoc autopas::Functor::SoAFunctorVerlet()
   */
  // clang-format on
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
   * Sets the particle properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   * @param epsilon24 epsilon * 24
   * @param sigmaSquared sigma^2
   */
  void setParticleProperties(SoAFloatPrecision epsilon24, SoAFloatPrecision sigmaSquared) {
    _epsilon24 = epsilon24;
    _sigmaSquared = sigmaSquared;
    if (applyShift) {
      _shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(_epsilon24, _sigmaSquared, _cutoffSquared);
    } else {
      _shift6 = 0;
    }
  }

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
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soaA, autopas::SoAView<SoAArraysType> soaB) {
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

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    // local redeclarations to help compilers
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> sigmaSquareds;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> epsilon24s;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> shift6s;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionBx;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionBy;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactSitePositionBz;

    // we require arrays for forces for sites to maintain SIMD in site-site calculations
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceBx;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceBy;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceBz;

    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteTypesB;
    std::vector<char, autopas::AlignedAllocator<char>> isSiteOwnedBArr;

    const SoAFloatPrecision const_sigmaSquared = _sigmaSquared;
    const SoAFloatPrecision const_epsilon24 = _epsilon24;
    const SoAFloatPrecision const_shift6 = _shift6;

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

    if constexpr (useMixing) {
      siteTypesB.reserve(siteCountB);
      sigmaSquareds.reserve(siteCountB);
      epsilon24s.reserve(siteCountB);
      if constexpr (applyShift) {
        shift6s.reserve(siteCountB);
      }
    }

    // Fill site-wise std::vectors for SIMD
    if constexpr (useMixing) {
      size_t siteIndex = 0;
      for (size_t mol = 0; mol < soaB.getNumberOfParticles(); ++mol) {
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0Bptr[mol], q1Bptr[mol], q2Bptr[mol], q3Bptr[mol]}, _PPLibrary->getSitePositions(typeptrB[mol]));
        const auto siteTypesOfMol = _PPLibrary->getSiteTypes(typeptrB[mol]);

        for (size_t site = 0; site < _PPLibrary->getNumSites(typeptrB[mol]); ++site) {
          exactSitePositionBx[siteIndex] = rotatedSitePositions[site][0] + xBptr[mol];
          exactSitePositionBy[siteIndex] = rotatedSitePositions[site][1] + yBptr[mol];
          exactSitePositionBz[siteIndex] = rotatedSitePositions[site][2] + zBptr[mol];
          siteTypesB[siteIndex] = siteTypesOfMol[site];
          siteForceBx[siteIndex] = 0.;
          siteForceBy[siteIndex] = 0.;
          siteForceBz[siteIndex] = 0.;
          if (calculateGlobals) {
            isSiteOwnedBArr[siteIndex] = ownedStatePtrB[mol] == autopas::OwnershipState::owned;
          }
          ++siteIndex;
        }
      }
    } else {
      size_t siteIndex = 0;
      for (size_t mol = 0; mol < soaB.getNumberOfParticles(); mol++) {
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0Bptr[mol], q1Bptr[mol], q2Bptr[mol], q3Bptr[mol]}, const_unrotatedSitePositions);
        for (size_t site = 0; site < const_unrotatedSitePositions.size(); ++site) {
          exactSitePositionBx[siteIndex] = rotatedSitePositions[site][0] + xBptr[mol];
          exactSitePositionBy[siteIndex] = rotatedSitePositions[site][1] + yBptr[mol];
          exactSitePositionBz[siteIndex] = rotatedSitePositions[site][2] + zBptr[mol];
          siteForceBx[siteIndex] = 0.;
          siteForceBy[siteIndex] = 0.;
          siteForceBz[siteIndex] = 0.;
          if (calculateGlobals) {
            isSiteOwnedBArr[siteIndex] = ownedStatePtrB[mol] == autopas::OwnershipState::owned;
          }
          ++siteIndex;
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
      std::vector<char, autopas::AlignedAllocator<char>> molMask;
      molMask.reserve(soaB.getNumberOfParticles());

#pragma omp simd
      for (size_t molB = 0; molB < soaB.getNumberOfParticles(); ++molB) {
        const auto ownedStateB = ownedStatePtrB[molB];

        const auto displacementCoMX = xAptr[molA] - xBptr[molB];
        const auto displacementCoMY = yAptr[molA] - yBptr[molB];
        const auto displacementCoMZ = zAptr[molA] - zBptr[molB];

        const auto distanceSquaredCoMX = displacementCoMX * displacementCoMX;
        const auto distanceSquaredCoMY = displacementCoMY * displacementCoMY;
        const auto distanceSquaredCoMZ = displacementCoMZ * displacementCoMZ;

        const auto distanceSquaredCoM = distanceSquaredCoMX + distanceSquaredCoMY + distanceSquaredCoMZ;

        // mask sites of molecules beyond cutoff or if molecule is a dummy
        molMask[molB] = distanceSquaredCoM <= cutoffSquared and ownedStateB != autopas::OwnershipState::dummy;
      }

      // generate mask for each site in cell B from molecular mask
      std::vector<char, autopas::AlignedAllocator<char>> siteMask;
      siteMask.reserve(siteCountB);

      for (size_t molB = 0; molB < soaB.getNumberOfParticles(); ++molB) {
        for (size_t siteB = 0; siteB < _PPLibrary->getNumSites(typeptrB[molB]); ++siteB) {
          siteMask.emplace_back(molMask[molB]);
        }
      }

      // sums used for molA
      SoAFloatPrecision forceSumX = 0.;
      SoAFloatPrecision forceSumY = 0.;
      SoAFloatPrecision forceSumZ = 0.;
      SoAFloatPrecision torqueSumX = 0.;
      SoAFloatPrecision torqueSumY = 0.;
      SoAFloatPrecision torqueSumZ = 0.;

      for (size_t siteA = 0; siteA < noSitesInMolA; ++siteA) {
        if (useMixing) {
          // preload sigmas, epsilons, and shifts
          for (size_t siteB = 0; siteB < siteCountB; ++siteB) {
            const auto mixingData =
                _PPLibrary->getMixingData(_PPLibrary->getSiteTypes(typeptrA[molA])[siteA], siteTypesB[siteB]);
            sigmaSquareds[siteB] = mixingData.sigmaSquared;
            epsilon24s[siteB] = mixingData.epsilon24;
            if (applyShift) {
              shift6s[siteB] = mixingData.shift6;
            }
          }
        }

        const auto rotatedSitePositionAx = rotatedSitePositionsA[siteA][0];
        const auto rotatedSitePositionAy = rotatedSitePositionsA[siteA][1];
        const auto rotatedSitePositionAz = rotatedSitePositionsA[siteA][2];

        const auto exactSitePositionAx = rotatedSitePositionAx + xAptr[molA];
        const auto exactSitePositionAy = rotatedSitePositionAy + yAptr[molA];
        const auto exactSitePositionAz = rotatedSitePositionAz + zAptr[molA];

#pragma omp simd reduction (+ : forceSumX, forceSumY, forceSumZ, torqueSumX, torqueSumY, torqueSumZ, potentialEnergySum, virialSumX, virialSumY, virialSumZ)
        for (size_t siteB = 0; siteB < siteCountB; ++siteB) {
          const SoAFloatPrecision sigmaSquared = useMixing ? sigmaSquareds[siteB] : const_sigmaSquared;
          const SoAFloatPrecision epsilon24 = useMixing ? epsilon24s[siteB] : const_epsilon24;
          const SoAFloatPrecision shift6 = applyShift ? (useMixing ? shift6s[siteB] : const_shift6) : 0;

          const auto isSiteOwnedB = !calculateGlobals || isSiteOwnedBArr[siteB];

          const auto displacementX = exactSitePositionAx - exactSitePositionBx[siteB];
          const auto displacementY = exactSitePositionAy - exactSitePositionBy[siteB];
          const auto displacementZ = exactSitePositionAz - exactSitePositionBz[siteB];

          const auto distanceSquaredX = displacementX * displacementX;
          const auto distanceSquaredY = displacementY * displacementY;
          const auto distanceSquaredZ = displacementZ * displacementZ;

          const auto distanceSquared = distanceSquaredX + distanceSquaredY + distanceSquaredZ;

          const auto invDistSquared = 1. / distanceSquared;
          const auto lj2 = sigmaSquared * invDistSquared;
          const auto lj6 = lj2 * lj2 * lj2;
          const auto lj12 = lj6 * lj6;
          const auto lj12m6 = lj12 - lj6;
          const auto scalarMultiple = siteMask[siteB] ? epsilon24 * (lj12 + lj12m6) * invDistSquared : 0.;

          // calculate forces
          const auto forceX = scalarMultiple * displacementX;
          const auto forceY = scalarMultiple * displacementY;
          const auto forceZ = scalarMultiple * displacementZ;

          forceSumX += forceX;
          forceSumY += forceY;
          forceSumZ += forceZ;

          torqueSumX += rotatedSitePositionAy * forceZ - rotatedSitePositionAz * forceY;
          torqueSumY += rotatedSitePositionAz * forceX - rotatedSitePositionAx * forceZ;
          torqueSumZ += rotatedSitePositionAx * forceY - rotatedSitePositionAy * forceX;

          // N3L ( total molecular forces + torques to be determined later )
          if constexpr (newton3) {
            siteForceBx[siteB] -= forceX;
            siteForceBy[siteB] -= forceY;
            siteForceBz[siteB] -= forceZ;
          }

          // globals
          if constexpr (calculateGlobals) {
            const auto potentialEnergy6 = siteMask[siteB] ? (epsilon24 * lj12m6 + shift6) : 0.;
            const auto virialX = displacementX * forceX;
            const auto virialY = displacementY * forceY;
            const auto virialZ = displacementZ * forceZ;

            // Add to the potential energy sum for each particle which is owned.
            // This results in obtaining 12 * the potential energy for the SoA.
            const auto ownershipFactor =
                newton3 ? (ownedStateA == autopas::OwnershipState::owned ? 1. : 0.) + (isSiteOwnedB ? 1. : 0.)
                        : (ownedStateA == autopas::OwnershipState::owned ? 1. : 0.);
            potentialEnergySum += potentialEnergy6 * ownershipFactor;
            virialSumX += virialX * ownershipFactor;
            virialSumY += virialY * ownershipFactor;
            virialSumZ += virialZ * ownershipFactor;
          }
        }
      }
      fxAptr[molA] += forceSumX;
      fyAptr[molA] += forceSumY;
      fzAptr[molA] += forceSumZ;
      txAptr[molA] += torqueSumX;
      tyAptr[molA] += torqueSumY;
      tzAptr[molA] += torqueSumZ;
    }

    // reduce the forces on individual sites in SoA B to total forces & torques on whole molecules
    if constexpr (useMixing) {
      size_t siteIndex = 0;
      for (size_t mol = 0; mol < soaB.getNumberOfParticles(); ++mol) {
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
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
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
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
      const auto threadNum = autopas::autopas_get_thread_num();
      // SoAFunctorPairImpl obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
      // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.
      const auto newton3Factor = newton3 ? .5 : 1.;

      _aosThreadData[threadNum].potentialEnergySum += potentialEnergySum * newton3Factor;
      _aosThreadData[threadNum].virialSum[0] += virialSumX * newton3Factor;
      _aosThreadData[threadNum].virialSum[1] += virialSumY * newton3Factor;
      _aosThreadData[threadNum].virialSum[2] += virialSumZ * newton3Factor;
    }
  }

  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexPrime,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    // Skip if primary particle is dummy
    const auto ownedStatePrime = ownedStatePtr[indexPrime];
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

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;
    const auto const_unrotatedSitePositions = _sitePositionsLJ;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> sigmaSquareds;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> epsilon24s;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> shift6s;

    const auto const_sigmaSquared = _sigmaSquared;
    const auto const_epsilon24 = _epsilon24;
    const auto const_shift6 = _shift6;

    const size_t neighborListSize = neighborList.size();
    const size_t *const __restrict neighborListPtr = neighborList.data();

    // Count sites
    const size_t siteCountMolPrime =
        useMixing ? _PPLibrary->getNumSites(typeptr[indexPrime]) : const_unrotatedSitePositions.size();

    size_t siteCountNeighbors = 0;  // site count of neighbours of primary molecule
    if constexpr (useMixing) {
      for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
        siteCountNeighbors += _PPLibrary->getNumSites(typeptr[neighborList[neighborMol]]);
      }
    } else {
      siteCountNeighbors = const_unrotatedSitePositions.size() * neighborListSize;
    }

    // initialize site-wise arrays for neighbors
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactNeighborSitePositionX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactNeighborSitePositionY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactNeighborSitePositionZ;

    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteTypesNeighbors;
    std::vector<char, autopas::AlignedAllocator<char>> isNeighborSiteOwnedArr;

    // we require arrays for forces for sites to maintain SIMD in site-site calculations
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceZ;

    // pre-reserve arrays
    exactNeighborSitePositionX.reserve(siteCountNeighbors);
    exactNeighborSitePositionY.reserve(siteCountNeighbors);
    exactNeighborSitePositionZ.reserve(siteCountNeighbors);

    if constexpr (useMixing) {
      siteTypesNeighbors.reserve(siteCountNeighbors);
    }

    siteForceX.reserve(siteCountNeighbors);
    siteForceY.reserve(siteCountNeighbors);
    siteForceZ.reserve(siteCountNeighbors);

    if constexpr (calculateGlobals) {
      isNeighborSiteOwnedArr.reserve(siteCountNeighbors);
    }

    if constexpr (useMixing) {
      sigmaSquareds.reserve(siteCountNeighbors);
      epsilon24s.reserve(siteCountNeighbors);
      if constexpr (applyShift) {
        shift6s.reserve(siteCountNeighbors);
      }
    }

    const auto rotatedSitePositionsPrime =
        useMixing ? autopas::utils::quaternion::rotateVectorOfPositions(
                        {q0ptr[indexPrime], q1ptr[indexPrime], q2ptr[indexPrime], q3ptr[indexPrime]},
                        _PPLibrary->getSitePositions(typeptr[indexPrime]))
                  : autopas::utils::quaternion::rotateVectorOfPositions(
                        {q0ptr[indexPrime], q1ptr[indexPrime], q2ptr[indexPrime], q3ptr[indexPrime]},
                        const_unrotatedSitePositions);

    const auto siteTypesPrime = _PPLibrary->getSiteTypes(typeptr[indexPrime]);  // this is not used if non-mixing

    // generate site-wise arrays for neighbors of primary mol
    if constexpr (useMixing) {
      size_t siteIndex = 0;
      for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
        const auto neighborMolIndex = neighborList[neighborMol];
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[neighborMolIndex], q1ptr[neighborMolIndex], q2ptr[neighborMolIndex], q3ptr[neighborMolIndex]},
            _PPLibrary->getSitePositions(typeptr[neighborMolIndex]));
        const auto siteTypesOfMol = _PPLibrary->getSiteTypes(typeptr[neighborMolIndex]);

        for (size_t site = 0; site < _PPLibrary->getNumSites(typeptr[neighborMolIndex]); ++site) {
          exactNeighborSitePositionX[siteIndex] = rotatedSitePositions[site][0] + xptr[neighborMolIndex];
          exactNeighborSitePositionY[siteIndex] = rotatedSitePositions[site][1] + yptr[neighborMolIndex];
          exactNeighborSitePositionZ[siteIndex] = rotatedSitePositions[site][2] + zptr[neighborMolIndex];
          siteTypesNeighbors[siteIndex] = siteTypesOfMol[site];
          siteForceX[siteIndex] = 0.;
          siteForceY[siteIndex] = 0.;
          siteForceZ[siteIndex] = 0.;
          if (calculateGlobals) {
            isNeighborSiteOwnedArr[siteIndex] = ownedStatePtr[neighborMolIndex] == autopas::OwnershipState::owned;
          }
          ++siteIndex;
        }
      }
    } else {
      size_t siteIndex = 0;
      for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
        const auto neighborMolIndex = neighborList[neighborMol];
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[neighborMolIndex], q1ptr[neighborMolIndex], q2ptr[neighborMolIndex], q3ptr[neighborMolIndex]},
            const_unrotatedSitePositions);
        for (size_t site = 0; site < const_unrotatedSitePositions.size(); ++site) {
          exactNeighborSitePositionX[siteIndex] = rotatedSitePositions[site][0] + xptr[neighborMolIndex];
          exactNeighborSitePositionY[siteIndex] = rotatedSitePositions[site][1] + yptr[neighborMolIndex];
          exactNeighborSitePositionZ[siteIndex] = rotatedSitePositions[site][2] + zptr[neighborMolIndex];
          siteForceX[siteIndex] = 0.;
          siteForceY[siteIndex] = 0.;
          siteForceZ[siteIndex] = 0.;
          if (calculateGlobals) {
            isNeighborSiteOwnedArr[siteIndex] = ownedStatePtr[neighborMolIndex] == autopas::OwnershipState::owned;
          }
          ++siteIndex;
        }
      }
    }

    // -- main force calculation --

    // - calculate mol mask -
    std::vector<char, autopas::AlignedAllocator<char>> molMask;
    molMask.reserve(neighborListSize);

#pragma omp simd
    for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
      const auto neighborMolIndex = neighborList[neighborMol];  // index of neighbor mol in soa

      const auto ownedState = ownedStatePtr[neighborMolIndex];

      const auto displacementCoMX = xptr[indexPrime] - xptr[neighborMolIndex];
      const auto displacementCoMY = yptr[indexPrime] - yptr[neighborMolIndex];
      const auto displacementCoMZ = zptr[indexPrime] - zptr[neighborMolIndex];

      const auto distanceSquaredCoMX = displacementCoMX * displacementCoMX;
      const auto distanceSquaredCoMY = displacementCoMY * displacementCoMY;
      const auto distanceSquaredCoMZ = displacementCoMZ * displacementCoMZ;

      const auto distanceSquaredCoM = distanceSquaredCoMX + distanceSquaredCoMY + distanceSquaredCoMZ;

      // mask molecules beyond cutoff or if molecule is a dummy
      molMask[neighborMol] = distanceSquaredCoM <= cutoffSquared and ownedState != autopas::OwnershipState::dummy;
    }

    // generate mask for each site from molecular mask
    std::vector<char, autopas::AlignedAllocator<char>> siteMask;
    siteMask.reserve(siteCountNeighbors);

    for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
      const auto neighborMolIndex = neighborList[neighborMol];  // index of neighbor mol in soa
      for (size_t siteB = 0; siteB < _PPLibrary->getNumSites(typeptr[neighborMolIndex]); ++siteB) {
        siteMask.emplace_back(molMask[neighborMol]);
      }
    }

    // sums used for prime mol
    SoAFloatPrecision forceSumX = 0.;
    SoAFloatPrecision forceSumY = 0.;
    SoAFloatPrecision forceSumZ = 0.;
    SoAFloatPrecision torqueSumX = 0.;
    SoAFloatPrecision torqueSumY = 0.;
    SoAFloatPrecision torqueSumZ = 0.;

    // - actual LJ calculation -

    for (size_t primeSite = 0; primeSite < siteCountMolPrime; ++primeSite) {
      const auto rotatedPrimeSitePositionX = rotatedSitePositionsPrime[primeSite][0];
      const auto rotatedPrimeSitePositionY = rotatedSitePositionsPrime[primeSite][1];
      const auto rotatedPrimeSitePositionZ = rotatedSitePositionsPrime[primeSite][2];

      const auto exactPrimeSitePositionX = rotatedPrimeSitePositionX + xptr[indexPrime];
      const auto exactPrimeSitePositionY = rotatedPrimeSitePositionY + yptr[indexPrime];
      const auto exactPrimeSitePositionZ = rotatedPrimeSitePositionZ + zptr[indexPrime];

      // generate parameter data for chosen site
      if constexpr (useMixing) {
        const auto primeSiteType = siteTypesPrime[primeSite];

        size_t siteIndex = 0;
        for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
          const auto neighborMolIndex = neighborList[neighborMol];
          const auto siteTypesOfNeighborMol = _PPLibrary->getSiteTypes(typeptr[neighborMolIndex]);

          for (size_t site = 0; site < _PPLibrary->getNumSites(typeptr[neighborMolIndex]); ++site) {
            const auto mixingData = _PPLibrary->getMixingData(primeSiteType, siteTypesOfNeighborMol[site]);
            sigmaSquareds[siteIndex] = mixingData.sigmaSquared;
            epsilon24s[siteIndex] = mixingData.epsilon24;
            if constexpr (applyShift) {
              shift6s[siteIndex] = mixingData.shift6;
            }
            ++siteIndex;
          }
        }
      }

#pragma omp simd reduction(+ : forceSumX, forceSumY, forceSumZ, torqueSumX, torqueSumY, torqueSumZ)
      for (size_t neighborSite = 0; neighborSite < siteCountNeighbors; ++neighborSite) {
        const SoAFloatPrecision sigmaSquared = useMixing ? sigmaSquareds[neighborSite] : const_sigmaSquared;
        const SoAFloatPrecision epsilon24 = useMixing ? epsilon24s[neighborSite] : const_epsilon24;
        const SoAFloatPrecision shift6 = applyShift ? (useMixing ? shift6s[neighborSite] : const_shift6) : 0;

        const bool isNeighborSiteOwned = !calculateGlobals || isNeighborSiteOwnedArr[neighborSite];

        const auto displacementX = exactPrimeSitePositionX - exactNeighborSitePositionX[neighborSite];
        const auto displacementY = exactPrimeSitePositionY - exactNeighborSitePositionY[neighborSite];
        const auto displacementZ = exactPrimeSitePositionZ - exactNeighborSitePositionZ[neighborSite];

        const auto distanceSquaredX = displacementX * displacementX;
        const auto distanceSquaredY = displacementY * displacementY;
        const auto distanceSquaredZ = displacementZ * displacementZ;

        const auto distanceSquared = distanceSquaredX + distanceSquaredY + distanceSquaredZ;

        const auto invDistSquared = 1. / distanceSquared;
        const auto lj2 = sigmaSquared * invDistSquared;
        const auto lj6 = lj2 * lj2 * lj2;
        const auto lj12 = lj6 * lj6;
        const auto lj12m6 = lj12 - lj6;
        const auto scalarMultiple = siteMask[neighborSite] ? epsilon24 * (lj12 + lj12m6) * invDistSquared : 0.;

        // calculate forces
        const auto forceX = scalarMultiple * displacementX;
        const auto forceY = scalarMultiple * displacementY;
        const auto forceZ = scalarMultiple * displacementZ;

        forceSumX += forceX;
        forceSumY += forceY;
        forceSumZ += forceZ;

        torqueSumX += rotatedPrimeSitePositionY * forceZ - rotatedPrimeSitePositionZ * forceY;
        torqueSumY += rotatedPrimeSitePositionZ * forceX - rotatedPrimeSitePositionX * forceZ;
        torqueSumZ += rotatedPrimeSitePositionX * forceY - rotatedPrimeSitePositionY * forceX;

        // N3L
        if (newton3) {
          siteForceX[neighborSite] -= forceX;
          siteForceY[neighborSite] -= forceY;
          siteForceZ[neighborSite] -= forceZ;
        }

        // calculate globals
        if constexpr (calculateGlobals) {
          const auto potentialEnergy6 = siteMask[neighborSite] ? (epsilon24 * lj12m6 + shift6) : 0.;
          const auto virialX = displacementX * forceX;
          const auto virialY = displacementY * forceY;
          const auto virialZ = displacementZ * forceZ;

          // Add to the potential energy sum for each particle which is owned.
          // This results in obtaining 12 * the potential energy for the SoA.
          const auto ownershipFactor =
              newton3 ? (ownedStatePrime == autopas::OwnershipState::owned ? 1. : 0.) + (isNeighborSiteOwned ? 1. : 0.)
                      : (ownedStatePrime == autopas::OwnershipState::owned ? 1. : 0.);
          potentialEnergySum += potentialEnergy6 * ownershipFactor;
          virialSumX += virialX * ownershipFactor;
          virialSumY += virialY * ownershipFactor;
          virialSumZ += virialZ * ownershipFactor;
        }
      }
    }
    // Add forces to prime mol
    fxptr[indexPrime] += forceSumX;
    fyptr[indexPrime] += forceSumY;
    fzptr[indexPrime] += forceSumZ;
    txptr[indexPrime] += torqueSumX;
    typtr[indexPrime] += torqueSumY;
    tzptr[indexPrime] += torqueSumZ;

    // Reduce forces on individual neighbor sites to molecular forces & torques if newton3=true
    if constexpr (newton3) {
      if constexpr (useMixing) {
        size_t siteIndex = 0;
        for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
          const auto neighborMolIndex = neighborList[neighborMol];
          const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
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
          const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
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
      const auto threadNum = autopas::autopas_get_thread_num();
      // SoAFunctorSingle obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
      // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.
      const auto newton3Factor = newton3 ? .5 : 1.;

      _aosThreadData[threadNum].potentialEnergySum += potentialEnergySum * newton3Factor;
      _aosThreadData[threadNum].virialSum[0] += virialSumX * newton3Factor;
      _aosThreadData[threadNum].virialSum[1] += virialSumY * newton3Factor;
      _aosThreadData[threadNum].virialSum[2] += virialSumZ * newton3Factor;
    }
  }

  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
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
}  // namespace mdLib