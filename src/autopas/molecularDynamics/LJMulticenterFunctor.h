/**
 * @file LJMulticenterFunctor.h
 * @date 21/02/2022
 * @author S. Newcome
*/

#pragma once

#include "MulticenteredMoleculeLJ.h"
#include "autopas/molecularDynamics/MoleculeLJ.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Math.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"
#include "autopas/utils/Quaternion.h"
#include "autopas/utils/AlignedAllocator.h"

namespace autopas {

/**
 * A functor to handle Lennard-Jones interactions between two (potentially multicentered) Molecules.
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
class LJMulticenterFunctor
    : public autopas::Functor<Particle, LJMulticenterFunctor<Particle, applyShift, useMixing, useNewton3,
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
  LJMulticenterFunctor() = delete;

 private:
  /**
   * Internal (actually used) constructor
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit LJMulticenterFunctor(SoAFloatPrecision cutoff, void * /*dummy*/)
      : autopas::Functor<Particle, LJMulticenterFunctor<Particle, applyShift, useMixing, useNewton3, calculateGlobals,
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
  explicit LJMulticenterFunctor(double cutoff) : LJMulticenterFunctor(cutoff, nullptr) {
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
  explicit LJMulticenterFunctor(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJMulticenterFunctor(cutoff, nullptr) {
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
   * Functor for AoS. Simply loops over the sites of two particles/molecules to calculate force.
   * @param particleA Particle i
   * @param particleB Particle j
   * @param newton3 Flag for if newton3 is used.
   */
  void AoSFunctor(Particle &particleA, Particle &particleB, bool newton3) final {
    if (particleA.isDummy() or particleB.isDummy()) {
      return;
    }
    // get site properties and positions

    // the following vectors contain mixed parameters in the order of a,b = 0,0; 0,1; ... ; 0,n; 1,0 ; ... ; n,m
    std::vector<double> sigmaSquareds;
    std::vector<double> epsilon24s;
    std::vector<double> shift6s;

    const size_t numSitesA = useMixing ? _PPLibrary->getNumSites(particleA.getTypeId()) : _sitePositionsLJ.size();
    const size_t numSitesB = useMixing ? _PPLibrary->getNumSites(particleB.getTypeId()) : _sitePositionsLJ.size();

    if constexpr (useMixing) {
      sigmaSquareds.reserve(numSitesA * numSitesB);
      epsilon24s.reserve(numSitesA * numSitesB);
      shift6s.reserve(numSitesA * numSitesB);

      const std::vector<size_t> siteIdsA = _PPLibrary->getSiteTypes(particleA.getTypeId());
      const std::vector<size_t> siteIdsB = _PPLibrary->getSiteTypes(particleB.getTypeId());
      for (int i = 0; i < numSitesA; ++i) {
        for (int j = 0; j < numSitesB; ++j) {
          sigmaSquareds.emplace_back(_PPLibrary->mixingSigmaSquare(siteIdsA[i], siteIdsB[j]));
          epsilon24s.emplace_back(_PPLibrary->mixing24Epsilon(siteIdsA[i], siteIdsB[j]));
          if constexpr (applyShift) {
            shift6s.emplace_back(_PPLibrary->mixingShift6(siteIdsA[i], siteIdsB[j]));
          }
        }
      }
    }

    const std::vector<std::array<double, 3>> unrotatedSitePositionsA =
        useMixing ? _PPLibrary->getSitePositions(particleA.getTypeId()) : _sitePositionsLJ;
    const std::vector<std::array<double, 3>> unrotatedSitePositionsB =
        useMixing ? _PPLibrary->getSitePositions(particleB.getTypeId()) : _sitePositionsLJ;

    double lj12m6Sum = 0;

    const auto displacementCoM = autopas::utils::ArrayMath::sub(particleA.getR(), particleB.getR());
    const auto distanceSquaredCoM = autopas::utils::ArrayMath::dot(displacementCoM, displacementCoM);

    // Don't calculate LJ if particleB outside cutoff of particleA
    if (distanceSquaredCoM > _cutoffSquared) {
      return;
    }

    // calculate relative site positions (rotated correctly)
    const auto rotatedSitePositionsA =
        autopas::utils::quaternion::rotateVectorOfPositions(particleA.getQ(), unrotatedSitePositionsA);
    const auto rotatedSitePositionsB =
        autopas::utils::quaternion::rotateVectorOfPositions(particleB.getQ(), unrotatedSitePositionsB);

    size_t ppl_index = 0;

    for (int m = 0; m < numSitesA; m++) {
      for (int n = 0; n < numSitesB; n++) {
        const auto displacement = autopas::utils::ArrayMath::add(
            autopas::utils::ArrayMath::sub(displacementCoM, rotatedSitePositionsB[n]), rotatedSitePositionsA[m]);
        const auto distanceSquared = autopas::utils::ArrayMath::dot(displacement, displacement);

        const auto sigmaSquared = useMixing ? sigmaSquareds[ppl_index] : _sigmaSquared;
        const auto epsilon24 = useMixing ? epsilon24s[ppl_index] : _epsilon24;
        const auto shift6 = useMixing ? shift6s[ppl_index] : _shift6;

        // Calculate potential between sites and thus force
        // Force = 24 * epsilon * (2*(sigma/distance)^12 - (sigma/distance)^6) * (1/distance)^2 * [x_displacement, y_displacement, z_displacement]
        //         {                         scalarMultiple                                   } * {                     displacement             }
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
        particleA.addTorque(autopas::utils::ArrayMath::cross(rotatedSitePositionsA[m], force));
        if (newton3) {
          particleB.subTorque(autopas::utils::ArrayMath::cross(rotatedSitePositionsB[n], force));
        }

        if (calculateGlobals) {
          // in newton3-case, division by 2 is handled here; in non-newton3-case, division is handled in post-processing
          const auto potentialEnergy = newton3 ? 0.5 * (epsilon24 * lj12m6 + shift6) : (epsilon24 * lj12m6 + shift6);
          const auto virial = newton3 ? utils::ArrayMath::mulScalar(utils::ArrayMath::mul(displacement, force),0.5) : utils::ArrayMath::mul(displacement, force);

          const auto threadNum = autopas_get_thread_num();

          if (particleA.isOwned()) {
            _aosThreadData[threadNum].potentialEnergySum += potentialEnergy;
            _aosThreadData[threadNum].virialSum = utils::ArrayMath::add(_aosThreadData[threadNum].virialSum,virial);
          }
          if (newton3 and particleB.isOwned()) {
            _aosThreadData[threadNum].potentialEnergySum += potentialEnergy;
            _aosThreadData[threadNum].virialSum = utils::ArrayMath::add(_aosThreadData[threadNum].virialSum,virial);
          }

        }

        ++ppl_index;
      }
    }

    // calculate globals
    if (calculateGlobals) {
      // todo sort this out - needs to work for mixing + factor differences for multi-centre case
    }
  }

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

    const SoAFloatPrecision const_sigmaSquared = _sigmaSquared;
    const SoAFloatPrecision const_epsilon24 = _epsilon24;
    const SoAFloatPrecision const_shift6 = _shift6;

    const auto const_unrotatedSitePositions = _sitePositionsLJ;

    size_t siteCount = 0;
    if constexpr (useMixing) {
      for (size_t mol = 0; mol < soa.getNumberOfParticles(); ++mol) {
        siteCount += _PPLibrary->getNumSites(typeptr[mol]);
      }
    } else {
      siteCount = const_unrotatedSitePositions.size() * soa.getNumberOfParticles();
    }

    exactSitePositionX.reserve(siteCount);
    exactSitePositionY.reserve(siteCount);
    exactSitePositionZ.reserve(siteCount);

    if constexpr (useMixing) {
      siteTypes.reserve(siteCount);
    }

    siteForceX.reserve((siteCount));
    siteForceY.reserve((siteCount));
    siteForceZ.reserve((siteCount));

    // Generate site-wise arrays for SIMD
    // todo can we simd anything here? / can we do something 'nicer'?
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
          ++siteIndex;
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
      const size_t siteIndexMolB = siteIndexMolA + noSitesInMolA;
      const size_t noSitesB = (siteCount - siteIndexMolB);  // Number of sites in molecules that A interacts with

      // create mask over every mol 'above' molA  (char to keep arrays aligned)
      std::vector<char, autopas::AlignedAllocator<char>> molMask;
      molMask.reserve(soa.getNumberOfParticles() - (molA+1));

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
        molMask[molB - (molA + 1)] = distanceSquaredCoM <= cutoffSquared and ownedStateB != autopas::OwnershipState::dummy;
      }

      // generate mask for each site in the mols 'above' molA from molecular mask
      // todo investigate if a SIMD implementation is possible
      std::vector<char, autopas::AlignedAllocator<char>> siteMask;
      siteMask.reserve(noSitesB);

      for (size_t molB = molA + 1; molB < soa.getNumberOfParticles(); ++molB) {
        for (size_t siteB = 0; siteB < _PPLibrary->getNumSites(typeptr[molB]); ++siteB) {
          siteMask.template emplace_back(molMask[molB - (molA + 1)]);
        }
      }

      // calculate LJ forces

      SoAFloatPrecision shift6 = const_shift6;
      SoAFloatPrecision sigmaSquared = const_sigmaSquared;
      SoAFloatPrecision epsilon24 = const_epsilon24;

      for (size_t siteA = siteIndexMolA; siteA < siteIndexMolB; ++siteA) {
        if (useMixing) {
          // preload sigmas, epsilons, and shifts
          sigmaSquareds.reserve(noSitesB);
          epsilon24s.reserve(noSitesB);
          if constexpr (applyShift) {
            shift6s.reserve(noSitesB);
          }

          for (size_t siteB = 0; siteB < siteCount - (siteIndexMolB); ++siteB) {
            const auto mixingData =
                _PPLibrary->getMixingData(siteTypes[siteA], siteTypes[siteIndexMolB + siteB]);
            sigmaSquareds[siteB] = mixingData.sigmaSquare;
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
          siteForceX[siteB] -= forceX;
          siteForceY[siteB] -= forceY;
          siteForceZ[siteB] -= forceZ;

          if (calculateGlobals) {
            // todo this
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
          txptr[mol] += rotatedSitePositions[site][0] * siteForceY[siteIndex] -
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
          txptr[mol] += rotatedSitePositions[site][0] * siteForceY[siteIndex] -
                        rotatedSitePositions[site][1] * siteForceX[siteIndex];
          ++siteIndex;
        }
      }
    }

    if (calculateGlobals) {
      // todo this
    }
  }
  /**
   * @copydoc Functor::SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3)
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
   * @copydoc Functor::SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst, const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList, bool newton3)
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
   * @copydoc Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 12>{
        Particle::AttributeNames::id,      Particle::AttributeNames::posX,    Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,    Particle::AttributeNames::forceX,  Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ,  Particle::AttributeNames::torqueX, Particle::AttributeNames::torqueY,
        Particle::AttributeNames::torqueZ, Particle::AttributeNames::typeId,  Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 6>{
        Particle::AttributeNames::id,   Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ, Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
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
   * @param numA number of sites in molecule A
   * @param numB number of sites in molecule B
   * @return #FLOPs
   */
  static unsigned long getNumFlopsPerKernelCall(bool newton3, size_t numA, size_t numB) {
    const unsigned long newton3Flops = newton3 ? 0ul : 3ul;
    return numA * numB * (15ul + newton3Flops);
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
   * todo replace this (and in LJFunctor.h) with the more intuitive potential energu
   * @return the potential energy sum
   */
  double getUpot() {
    if (not calculateGlobals) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Trying to get upot even though calculateGlobals is false. If you want this functor to calculate global "
          "values, please specify calculateGlobals to be true.");
    }
    if (not _postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException("Cannot get upot, because endTraversal was not called.");
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

    const SoAFloatPrecision const_sigmaSquared = _sigmaSquared;
    const SoAFloatPrecision const_epsilon24 = _epsilon24;
    const SoAFloatPrecision const_shift6 = _shift6;

    // load unrotated site positions
    const auto const_unrotatedSitePositions = _sitePositionsLJ;

    size_t siteCountA = 0;
    size_t siteCountB = 0;
    if constexpr (useMixing) {
      for (size_t mol = 0; mol < soaA.getNumberOfParticles(); ++mol) {
        siteCountA += _PPLibrary->getNumSites(typeptrA[mol]);
      }
      for (size_t mol = 0; mol < soaB.getNumberOfParticles(); ++mol) {
        siteCountB += _PPLibrary->getNumSites(typeptrB[mol]);
      }
    } else {
      siteCountA = const_unrotatedSitePositions.size() * soaA.getNumberOfParticles();
      siteCountB = const_unrotatedSitePositions.size() * soaB.getNumberOfParticles();
    }

    exactSitePositionBx.reserve(siteCountB);
    exactSitePositionBy.reserve(siteCountB);
    exactSitePositionBz.reserve(siteCountB);

    if constexpr (useMixing) {
      siteTypesB.reserve(siteCountB);
    }

    siteForceBx.reserve(siteCountB);
    siteForceBy.reserve(siteCountB);
    siteForceBz.reserve(siteCountB);

    if constexpr (useMixing) {
      siteTypesB.reserve(siteCountB);
      sigmaSquareds.reserve(siteCountB);
      epsilon24s.reserve(siteCountB);
      if constexpr (applyShift) {
        shift6s.reserve(siteCountB);
      }
    }

    // Generate site-wise arrays for SIMD
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

      const auto noSitesInMolA = useMixing ? _PPLibrary->getNumSites(typeptrA[molA])
                                             : const_unrotatedSitePositions.size();
      const auto unrotatedSitePositionsA = useMixing ? _PPLibrary->getSitePositions(typeptrA[molA]) : const_unrotatedSitePositions;

      const auto rotatedSitePositionsA = autopas::utils::quaternion::rotateVectorOfPositions(
          {q0Aptr[molA], q1Aptr[molA], q2Aptr[molA], q3Aptr[molA]},unrotatedSitePositionsA);

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

      SoAFloatPrecision shift6 = const_shift6;
      SoAFloatPrecision sigmaSquared = const_sigmaSquared;
      SoAFloatPrecision epsilon24 = const_epsilon24;

      // sums used for molA
      SoAFloatPrecision forceSumX = 0.;
      SoAFloatPrecision forceSumY = 0.;
      SoAFloatPrecision forceSumZ = 0.;
      SoAFloatPrecision torqueSumX = 0.;
      SoAFloatPrecision torqueSumY = 0.;
      SoAFloatPrecision torqueSumZ = 0.;

      for (size_t siteA = 0; siteA < noSitesInMolA; ++siteA) {
        // todo compare performance between preloading parameters and just calculating them
        if (useMixing) {
          // preload sigmas, epsilons, and shifts
          for (size_t siteB = 0; siteB < siteCountB; ++siteB) {
            const auto mixingData = _PPLibrary->getMixingData(_PPLibrary->getSiteTypes(typeptrA[molA])[siteA], siteTypesB[siteB]);
            sigmaSquareds[siteB] = mixingData.sigmaSquare;
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
          if (newton3) {
            siteForceBx[siteB] -= forceX;
            siteForceBy[siteB] -= forceY;
            siteForceBz[siteB] -= forceZ;
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

    // reduce the forces on individual sites in cell B to total forces & torques on whole molecules
    if constexpr (useMixing) {
      size_t siteIndex = 0;
      for (size_t mol = 0; mol < soaB.getNumberOfParticles(); ++mol) {
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0Bptr[mol], q1Bptr[mol], q2Bptr[mol], q3Bptr[mol]}, _PPLibrary->getSitePositions(typeptrB[mol]) );
        for (size_t site = 0; site < _PPLibrary->getNumSites(typeptrB[mol]); ++site) {
          fxBptr[mol] += siteForceBx[siteIndex];
          fyBptr[mol] += siteForceBy[siteIndex];
          fzBptr[mol] += siteForceBz[siteIndex];
          txBptr[mol] += rotatedSitePositions[site][1] * siteForceBz[siteIndex] -
                        rotatedSitePositions[site][2] * siteForceBy[siteIndex];
          tyBptr[mol] += rotatedSitePositions[site][2] * siteForceBx[siteIndex] -
                        rotatedSitePositions[site][0] * siteForceBz[siteIndex];
          txBptr[mol] += rotatedSitePositions[site][0] * siteForceBy[siteIndex] -
                        rotatedSitePositions[site][1] * siteForceBx[siteIndex];
          ++siteIndex;
        }
      }
    } else {
      // to do
    }

  }

  template <bool newton3>
  void SoAFunctorVerletImpl(SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    // Skip if initial particle is dummy
    const auto ownedStateInit = ownedStatePtr[indexFirst];
    if (ownedStateInit == OwnershipState::dummy) { return; }

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
    const size_t siteCountMolI = useMixing ? _PPLibrary->getNumSites(typeptr[indexFirst]) : const_unrotatedSitePositions.size();

    size_t siteCountNeighbors = 0; // site count of neighbours of mol I (which we refer to as mols J)
    if constexpr (useMixing) {
      for (size_t molJ = 0; molJ < neighborListSize; ++molJ) {
        siteCountNeighbors += _PPLibrary->getNumSites(typeptr[neighborList[molJ]]);
      }
    } else {
      siteCountNeighbors = const_unrotatedSitePositions.size() * neighborListSize;
    }

    // initalise site-wise arrays for neighbors
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactNeighborSitePositionX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactNeighborSitePositionY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> exactNeighborSitePositionZ;

    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteTypesI;
    std::vector<size_t, autopas::AlignedAllocator<size_t>> siteTypesNeighbors;

    // we require arrays for forces for sites to maintain SIMD in site-site calculations
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceZ;

    // pre-reserve arrays
    exactNeighborSitePositionX.reserve(siteCountNeighbors);
    exactNeighborSitePositionY.reserve(siteCountNeighbors);
    exactNeighborSitePositionZ.reserve(siteCountNeighbors);

    if constexpr (useMixing) {
      siteTypesI.reserve(siteCountMolI);
      siteTypesNeighbors.reserve(siteCountNeighbors);
    }

    siteForceX.reserve(siteCountNeighbors);
    siteForceY.reserve(siteCountNeighbors);
    siteForceZ.reserve(siteCountNeighbors);

    if constexpr (useMixing) {
      sigmaSquareds.reserve(siteCountNeighbors);
      epsilon24s.reserve(siteCountNeighbors);
      if constexpr (applyShift) {
        shift6s.reserve(siteCountNeighbors);
      }
    }

    const auto rotatedSitePositionsI = useMixing ?
                                                 autopas::utils::quaternion::rotateVectorOfPositions(
                                                       {q0ptr[indexFirst], q1ptr[indexFirst], q2ptr[indexFirst], q3ptr[indexFirst]},
                                                     _PPLibrary->getSitePositions(typeptr[indexFirst])) :
                                                 autopas::utils::quaternion::rotateVectorOfPositions(
                                                     {q0ptr[indexFirst], q1ptr[indexFirst], q2ptr[indexFirst], q3ptr[indexFirst]},
                                                     const_unrotatedSitePositions);

    // generate site-wise arrays for neighbours of mol I
    if constexpr (useMixing) {
      size_t siteIndex = 0;
      for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
        const auto neighborMolIndex = neighborList[neighborMol];
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[neighborMolIndex], q1ptr[neighborMolIndex], q2ptr[neighborMolIndex], q3ptr[neighborMolIndex]}, _PPLibrary->getSitePositions(typeptr[neighborMolIndex]));
        const auto siteTypesOfMol = _PPLibrary->getSiteTypes(typeptr[neighborMolIndex]);

        for (size_t site = 0; site < _PPLibrary->getNumSites(typeptr[neighborMolIndex]); ++site) {
          exactNeighborSitePositionX[siteIndex] = rotatedSitePositions[site][0] + xptr[neighborMolIndex];
          exactNeighborSitePositionY[siteIndex] = rotatedSitePositions[site][1] + yptr[neighborMolIndex];
          exactNeighborSitePositionZ[siteIndex] = rotatedSitePositions[site][2] + zptr[neighborMolIndex];
          siteTypesNeighbors[siteIndex] = siteTypesOfMol[site];
          const auto mixingData = _PPLibrary->getMixingData(typeptr[indexFirst],typeptr[neighborMolIndex]);
          sigmaSquareds[siteIndex] = mixingData.sigmaSquare;
          epsilon24s[siteIndex] = mixingData.epsilon24;
          if constexpr (applyShift) {
            shift6s[siteIndex] = mixingData.shift6;
          }
          siteForceX[siteIndex] = 0.;
          siteForceY[siteIndex] = 0.;
          siteForceZ[siteIndex] = 0.;
          ++siteIndex;
        }
      }
    } else {
      size_t siteIndex = 0;
      for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
        const auto neighborMolIndex = neighborList[neighborMol];
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[neighborMolIndex], q1ptr[neighborMolIndex], q2ptr[neighborMolIndex], q3ptr[neighborMolIndex]}, const_unrotatedSitePositions);
        for (size_t site = 0; site < const_unrotatedSitePositions.size(); ++site) {
          exactNeighborSitePositionX[siteIndex] = rotatedSitePositions[site][0] + xptr[neighborMolIndex];
          exactNeighborSitePositionY[siteIndex] = rotatedSitePositions[site][1] + yptr[neighborMolIndex];
          exactNeighborSitePositionZ[siteIndex] = rotatedSitePositions[site][2] + zptr[neighborMolIndex];
          siteForceX[siteIndex] = 0.;
          siteForceY[siteIndex] = 0.;
          siteForceZ[siteIndex] = 0.;
          ++siteIndex;
        }
      }
    }


    // -- main force calculation --

    // - calculate mol mask -
    std::vector<char, autopas::AlignedAllocator<char>> molMask;
    molMask.reserve(neighborListSize);

    // @note I'm not sure this will get any SIMD benefit due to xptr not accessed consecutively
    // todo try this, using a temp array to store CoM positions consecutively
#pragma omp simd
    for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
      const auto neighborMolIndex = neighborList[neighborMol]; // index of neighbor mol in soa

      const auto ownedState = ownedStatePtr[neighborMolIndex];

      const auto displacementCoMX = xptr[indexFirst] - xptr[neighborMolIndex];
      const auto displacementCoMY = yptr[indexFirst] - yptr[neighborMolIndex];
      const auto displacementCoMZ = zptr[indexFirst] - zptr[neighborMolIndex];

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
      const auto neighborMolIndex = neighborList[neighborMol]; // index of neighbor mol in soa
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

    for (size_t primeSite = 0; primeSite < siteCountMolI; ++primeSite) {
      const auto rotatedPrimeSitePositionX = rotatedSitePositionsI[primeSite][0];
      const auto rotatedPrimeSitePositionY = rotatedSitePositionsI[primeSite][1];
      const auto rotatedPrimeSitePositionZ = rotatedSitePositionsI[primeSite][2];

      const auto exactPrimeSitePositionX = rotatedPrimeSitePositionX + xptr[primeSite];
      const auto exactPrimeSitePositionY = rotatedPrimeSitePositionY + yptr[primeSite];
      const auto exactPrimeSitePositionZ = rotatedPrimeSitePositionZ + zptr[primeSite];

#pragma omp simd reduction (+ : forceSumX, forceSumY, forceSumZ, torqueSumX, torqueSumY, torqueSumZ)
      for (size_t neighborSite = 0; neighborSite < siteCountNeighbors; ++neighborSite) {
        const SoAFloatPrecision sigmaSquared = useMixing ? sigmaSquareds[neighborSite] : const_sigmaSquared;
        const SoAFloatPrecision epsilon24 = useMixing ? epsilon24s[neighborSite] : const_epsilon24;
        const SoAFloatPrecision shift6 = applyShift ? (useMixing ? shift6s[neighborSite] : const_shift6) : 0;

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
      }
    }
    // Add forces to prime mol
    fxptr[indexFirst] += forceSumX;
    fyptr[indexFirst] += forceSumY;
    fzptr[indexFirst] += forceSumZ;
    txptr[indexFirst] += torqueSumX;
    typtr[indexFirst] += torqueSumY;
    tzptr[indexFirst] += torqueSumZ;

    // Reduce forces on individual neighbor sites to molecular forces & torques
    if constexpr (useMixing) {
      size_t siteIndex = 0;
      for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
        const auto neighborMolIndex = neighborList[neighborMol];
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[neighborMolIndex], q1ptr[neighborMolIndex], q2ptr[neighborMolIndex], q3ptr[neighborMolIndex]},
            _PPLibrary->getSitePositions(typeptr[neighborMolIndex]) );
        for (size_t site = 0; site < _PPLibrary->getNumSites(typeptr[neighborMolIndex]); ++site) {
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
    } else {
      // todo
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

template <bool applyShift, bool useMixing, autopas::FunctorN3Modes useNewton3, bool calculateGlobals,
          bool relevantForTuning>
class LJMulticenterFunctor<autopas::MoleculeLJ, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>
    : public autopas::Functor<autopas::MoleculeLJ,
                              LJMulticenterFunctor<autopas::MoleculeLJ, applyShift, useMixing, useNewton3,
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
  LJMulticenterFunctor() = delete;

 private:
  /**
   * Internal (actually used) constructor
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit LJMulticenterFunctor(SoAFloatPrecision cutoff, void * /*dummy*/)
      : autopas::Functor<autopas::MoleculeLJ, LJMulticenterFunctor<autopas::MoleculeLJ, applyShift, useMixing,
                                                                   useNewton3, calculateGlobals, relevantForTuning>>(
            cutoff) {
    autopas::utils::ExceptionHandler::exception(
        "LJMulticenterFunctor can not be used with MoleculeLJ. Use a MulticenteredMoleculeLJ instead.");
  }

 public:
  /**
   * Constructor for Functor with particle mixing disabled. setParticleProperties() must be called.
   * @param cutoff
   */
  explicit LJMulticenterFunctor(double cutoff) : LJMulticenterFunctor(cutoff, nullptr) {
    autopas::utils::ExceptionHandler::exception(
        "LJMulticenterFunctor can not be used with MoleculeLJ. Use a MulticenteredMoleculeLJ instead.");
  }

  /**
   * Constructor for Functor with particle mixing enabled.
   * @param cutoff
   * @param particlePropertiesLibrary Library used to look up the properties of each type of particle e.g. sigma,
   * epsilon, shift.
   */
  explicit LJMulticenterFunctor(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJMulticenterFunctor(cutoff, nullptr) {
    autopas::utils::ExceptionHandler::exception(
        "LJMulticenterFunctor can not be used with MoleculeLJ. Use a MulticenteredMoleculeLJ instead.");
  }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  constexpr static bool getMixing() { return useMixing; }

  double getPotentialEnergy() { return 0; }

  double getVirial() { return 0; }

  /**
   * Functor for AoS. Simply loops over the sites of two particles/molecules to calculate force.
   * @param particleA Particle i
   * @param particleB Particle j
   * @param newton3 Flag for if newton3 is used.
   */
  void AoSFunctor(autopas::MoleculeLJ &particleA, autopas::MoleculeLJ &particleB, bool newton3) final {
    autopas::utils::ExceptionHandler::exception(
        "LJMulticenterFunctor can not be used with MoleculeLJ. Use a MulticenteredMoleculeLJ instead.");
  }

  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    autopas::utils::ExceptionHandler::exception(
        "LJMulticenterFunctor can not be used with MoleculeLJ. Use a MulticenteredMoleculeLJ instead.");
  }

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {
    autopas::utils::ExceptionHandler::exception(
        "LJMulticenterFunctor can not be used with MoleculeLJ. Use a MulticenteredMoleculeLJ instead.");
  }

  void setParticleProperties(SoAFloatPrecision epsilon24, SoAFloatPrecision sigmaSquared) {
    autopas::utils::ExceptionHandler::exception(
        "LJMulticenterFunctor can not be used with MoleculeLJ. Use a MulticenteredMoleculeLJ instead.");
  }

  static unsigned long getNumFlopsPerKernelCall(bool newton3, size_t numA, size_t numB) { return 0ul; }

  void initTraversal() final {
    autopas::utils::ExceptionHandler::exception(
        "LJMulticenterFunctor can not be used with MoleculeLJ. Use a MulticenteredMoleculeLJ instead.");
  }

  void endTraversal(bool newton3) final {
    autopas::utils::ExceptionHandler::exception(
        "LJMulticenterFunctor can not be used with MoleculeLJ. Use a MulticenteredMoleculeLJ instead.");
  }
};
}  // namespace autopas