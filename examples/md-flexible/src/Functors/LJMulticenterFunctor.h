/**
 * @file LJMulticenterFunctor.h
 * @date 21/02/2022
 * @author S. Newcome
*/

#pragma once

#include "src/Particles/MulticenteredMoleculeLJ.h"
#include "autopas/molecularDynamics/MoleculeLJ.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"
#include "autopas/utils/Quaternion.h"
#include "autopas/utils/AlignedAllocator.h"

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
template <class Particle, bool applyShift = false, bool useMixing = false, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGlobals = false, bool relevantForTuning = true>
class LJMulticenterFunctor
    : public autopas::Functor<Particle,
    LJMulticenterFunctor<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>> {
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
  const std::vector<std::array<double,3>> _sitePositionsLJ{};

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
  std::array<double,3> _virialSum;

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
    : autopas::Functor <Particle, LJMulticenterFunctor<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(
              cutoff
              ),
          _cutoffSquared{cutoff*cutoff},
          _potentialEnergySum{0.},
          _virialSum{0.,0.,0.},
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

  bool isRelevantForTuning() final {return relevantForTuning;}

  bool allowsNewton3() final { return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both; }

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
    // get site properties
    auto sigmaSquared = _sigmaSquared;
    auto epsilon24 = _epsilon24; // todo potentially rename
    auto shift6 = _shift6;

    // get relative site positions
    auto unrotatedSitePositionsA = _sitePositionsLJ;
    auto unrotatedSitePositionsB = _sitePositionsLJ;

    if constexpr (useMixing) {
      // todo add this functionality
      autopas::utils::ExceptionHandler::exception("LJMulticenterFunctor: Mixing with multicentered molecules not yet implemented");
    }

    double lj12m6Sum = 0;

    const auto displacementCoM = autopas::utils::ArrayMath::sub(particleA.getR(), particleB.getR());
    const auto distanceSquaredCoM = autopas::utils::ArrayMath::dot(displacementCoM,displacementCoM);

    // Don't calculate LJ if particleB outside cutoff of particleA
    if (distanceSquaredCoM > _cutoffSquared) { return; }

    // calculate relative site positions (rotated correctly)
    const auto rotatedSitePositionsA = autopas::utils::quaternion::rotateVectorOfPositions(particleA.getQ(), unrotatedSitePositionsA);
    const auto rotatedSitePositionsB = autopas::utils::quaternion::rotateVectorOfPositions(particleB.getQ(), unrotatedSitePositionsB);

    std::array<double,3> forceTotal{0.,0.,0.};
    std::array<double,3> torqueTotal{0.};

    // calculate number of sites for each molecule
    const size_t numSitesA = rotatedSitePositionsA.size();
    const size_t numSitesB = rotatedSitePositionsB.size();

    for (int m=0; m<numSitesA; m++) {
      for (int n=0; n<numSitesB; n++) {
        // todo check this
        const auto displacement = autopas::utils::ArrayMath::sub(
            autopas::utils::ArrayMath::add(displacementCoM,rotatedSitePositionsB[n]), rotatedSitePositionsA[m]);
        const auto distanceSquared = autopas::utils::ArrayMath::dot(displacement,displacement);

        // Calculate potential between sites and thus force
        // Force = 24 * epsilon * (2*(sigma/distance)^12 - (sigma/distance)^6) * (1/distance)^2 * [x_displacement, y_displacement, z_displacement]
        //         {                         scalarMultiple                                   } * {                displacement                  }
        const auto invDistSquared = 1. / distanceSquared;
        const auto lj2 = sigmaSquared * invDistSquared;
        const auto lj6 = lj2 * lj2 * lj2;
        const auto lj12 = lj6 * lj6;
        const auto lj12m6 =  lj12 - lj6; // = LJ potential / (4x epsilon)
        const auto scalarMultiple = epsilon24 * (lj12 + lj12m6) * invDistSquared;
        const auto force = autopas::utils::ArrayMath::mulScalar(displacement,scalarMultiple);

        // Add force on site to net force
        particleA.addF(force);
        if (newton3) {particleB.subF(force);}

        // Add torque applied by force
        particleA.addTorque(autopas::utils::ArrayMath::cross(rotatedSitePositionsA[m],force));
        if (newton3) {particleB.subTorque(autopas::utils::ArrayMath::cross(rotatedSitePositionsB[n],force));}

        if (calculateGlobals) {
          forceTotal = autopas::utils::ArrayMath::add(forceTotal,force);
          lj12m6Sum += lj12m6;
        }
      }
    }

    // calculate globals
    if (calculateGlobals) {
      const double newton3Multiplier = newton3 ? 0.5 : 1.0;

      const auto virial = autopas::utils::ArrayMath::mulScalar(autopas::utils::ArrayMath::mul(displacementCoM,forceTotal),newton3Multiplier);
      const auto potentialEnergy = (epsilon24 * lj12m6Sum + shift6 * numSitesA * numSitesB) * newton3Multiplier;

      const int threadNum = autopas::autopas_get_thread_num();

      if (particleA.isOwned()) {
        _aosThreadData[threadNum].potentialEnergySum += potentialEnergy;
        _aosThreadData[threadNum].virialSum = autopas::utils::ArrayMath::add(_aosThreadData[threadNum].virialSum, virial);
      }
      if (newton3 and particleB.isOwned()) {
        _aosThreadData[threadNum].potentialEnergySum += potentialEnergy;
        _aosThreadData[threadNum].virialSum = autopas::utils::ArrayMath::add(_aosThreadData[threadNum].virialSum, virial);
      }
    }
  }

  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    if (soa.getNumParticles() == 0) return;

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

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> sigmaSquares;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> epsilon24s;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> shift6s;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedSitePositionX;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedSitePositionY;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedSitePositionZ;

    const SoAFloatPrecision const_sigmaSquared = _sigmaSquared;
    const SoAFloatPrecision const_epsilon24 = _epsilon24;
    const SoAFloatPrecision const_shift6 = _shift6;

    if constexpr (useMixing) {
      // todo add this functionality
      autopas::utils::ExceptionHandler::exception("LJMulticenterFunctor: Mixing with multicentered molecules not yet implemented");
    }

    auto unrotatedSitePositions = _sitePositionsLJ;

    // calculate rotated site positions
    if constexpr (useMixing) {
      // todo add this
    } else {
      rotatedSitePositionX.resize(unrotatedSitePositions.size() * soa.getNumParticles());
      rotatedSitePositionY.resize(unrotatedSitePositions.size() * soa.getNumParticles());
      rotatedSitePositionZ.resize(unrotatedSitePositions.size() * soa.getNumParticles());


      // ideally, if we had an efficient way of extending molecular arrays to full site arrays, we would do:
      // SIMD for all mol
      //    calculate rotational matrix
      // rotMatPerSite = extend_mol_to_sites(rotMatPerMol)
      // SIMD for all sites
      //    calculate rotated relative site positions
      //
      // alas, we are limited by the technology of the time
      // todo implement this
      // we can use this to implement the vectorised CoM-CoM algorithm as described in Nikola's dissertation
      for (size_t mol = 0; mol < soa.getNumParticles(); mol++) {
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, unrotatedSitePositions);
        for (size_t site = 0; site < unrotatedSitePositions.size(); ++site) {
          rotatedSitePositionX[mol * unrotatedSitePositions.size() + site] = rotatedSitePositions[site][0];
          rotatedSitePositionY[mol * unrotatedSitePositions.size() + site] = rotatedSitePositions[site][1];
          rotatedSitePositionZ[mol * unrotatedSitePositions.size() + site] = rotatedSitePositions[site][2];
        }
      }
    }

    // main force calculation loop
    for (size_t molA = 0; molA < soa.getNumParticles(); ++molA) {
      const auto ownedStateA = ownedStatePtr[molA];
      if (ownedStateA == autopas::OwnershipState::dummy) {
        continue;
      }

      size_t noSitesInMolA = 1; // Number of sites in molecule A
      size_t noSitesB = 0; // Number of sites in molecules that A interacts with
      size_t siteIndexA = 0; // Index at which the number of sites in A starts
      if constexpr (useMixing) {

      } else {
        noSitesInMolA = unrotatedSitePositions.size();
        noSitesB = (soa.getNumParticles()-molA)*unrotatedSitePositions.size();
        siteIndexA = molA*unrotatedSitePositions.size();
      }

      // sums used for molA
      SoAFloatPrecision forceSumX = 0.;
      SoAFloatPrecision forceSumY = 0.;
      SoAFloatPrecision forceSumZ = 0.;
      SoAFloatPrecision torqueSumX = 0.;
      SoAFloatPrecision torqueSumY = 0.;
      SoAFloatPrecision torqueSumZ = 0.;

      // preload sigmas, epsilons, and shifts
      if constexpr (useMixing) {
        // todo add this functionality
        autopas::utils::ExceptionHandler::exception("LJMulticenterFunctor: Mixing with multicentered molecules not yet implemented");
      }

      // create mask over every site 'above' molA
      std::vector<bool, autopas::AlignedAllocator<bool>> mask;
      if constexpr (useMixing) {

      } else {
        mask.resize(noSitesB);
        size_t noSitesMolB = unrotatedSitePositions.size(); // uniform #sites

#pragma omp simd
        // mask sites of molecules with CoM outside cutoff
        for (size_t siteB = 0; siteB < noSitesB; ++siteB) {


          size_t molB = siteB / noSitesMolB + molA;

          const auto ownedStateB = ownedStatePtr[molB];

          const auto displacementCoMX = xptr[molA] - xptr[molB];
          const auto displacementCoMY = yptr[molA] - yptr[molB];
          const auto displacementCoMZ = zptr[molA] - zptr[molB];

          const auto distanceSquaredCoMX = displacementCoMX * displacementCoMX;
          const auto distanceSquaredCoMY = displacementCoMY * displacementCoMY;
          const auto distanceSquaredCoMZ = displacementCoMZ * displacementCoMZ;

          const auto distanceSquaredCoM = distanceSquaredCoMX + distanceSquaredCoMY + distanceSquaredCoMZ;

          // mask sites of molecules beyond cutoff or if molecule is a dummy
          mask[siteB] = distanceSquaredCoM <= cutoffSquared and ownedStateB != autopas::OwnershipState::dummy;
        }
      }

      // calculate LJ forces
      SoAFloatPrecision shift6 = _shift6;
      SoAFloatPrecision sigmaSquared = _sigmaSquared;
      SoAFloatPrecision epsilon24 = _epsilon24;
      // todo implement this for mixing
      for (int siteA = 0; siteA < noSitesInMolA; ++siteA) {
#pragma omp simd reduction (+ : forceSumX, forceSumY, forceSumZ, torqueSumX, torqueSumY, torqueSumZ, potentialEnergySum, virialSumX, virialSumY, virialSumZ)
        for (int siteB = 0; siteB < noSitesB; ++siteB) {
          size_t molB = siteB / unrotatedSitePositions.size() + molA;
          size_t siteIndexB = 0;
          if constexpr (useMixing) {

          } else {
            siteIndexB = molB * _sitePositionsLJ.size();
          }

          const auto displacementX = xptr[molA] - xptr[molB] + rotatedSitePositionX[siteIndexA+siteA] - rotatedSitePositionX[siteIndexB+siteB];
          const auto displacementY = yptr[molA] - yptr[molB] + rotatedSitePositionY[siteIndexA+siteA] - rotatedSitePositionY[siteIndexB+siteB];
          const auto displacementZ = zptr[molA] - zptr[molB] + rotatedSitePositionZ[siteIndexA+siteA] - rotatedSitePositionZ[siteIndexB+siteB];

          const auto distanceSquaredX = displacementX * displacementX;
          const auto distanceSquaredY = displacementY * displacementY;
          const auto distanceSquaredZ = displacementZ * displacementZ;

          const auto distanceSquared = distanceSquaredX + distanceSquaredY + distanceSquaredZ;

          const auto invDistSquared = 1. / distanceSquared;
          const auto lj2 = const_sigmaSquared * invDistSquared;
          const auto lj6 = lj2 * lj2 * lj2;
          const auto lj12 = lj6 * lj6;
          const auto lj12m6 = lj12 - lj6;
          const auto scalarMultiple = mask[siteB] ? const_epsilon24 * (lj12 + lj12m6) * invDistSquared : 0.;

          // calculate forces
          const auto forceX = scalarMultiple * displacementX;
          const auto forceY = scalarMultiple * displacementY;
          const auto forceZ = scalarMultiple * displacementZ;

          forceSumX += forceX;
          forceSumY += forceY;
          forceSumZ += forceZ;

          // newton's third law
          fxptr[molB] -= forceX;
          fyptr[molB] -= forceY;
          fzptr[molB] -= forceZ;

          // calculate torques on molA
          const auto torqueAx = rotatedSitePositionY[siteIndexA+siteA] * forceZ - rotatedSitePositionZ[siteIndexA+siteA] * forceY;
          const auto torqueAy = rotatedSitePositionZ[siteIndexA+siteA] * forceX - rotatedSitePositionX[siteIndexA+siteA] * forceZ;
          const auto torqueAz = rotatedSitePositionX[siteIndexA+siteA] * forceY - rotatedSitePositionY[siteIndexA+siteA] * forceX;

          torqueSumX += torqueAx;
          torqueSumY += torqueAy;
          torqueSumZ += torqueAz;

          // N3L / calculate torques on molB
          //const auto torqueBx = rotatedSitePositionZ[siteIndexB+siteB] * forceY - rotatedSitePositionY[siteIndexB+siteB] * forceZ;
          //const auto torqueBy = rotatedSitePositionX[siteIndexB+siteB] * forceZ - rotatedSitePositionZ[siteIndexB+siteB] * forceX;
          //const auto torqueBz = rotatedSitePositionY[siteIndexB+siteB] * forceX - rotatedSitePositionX[siteIndexB+siteB] * forceY;
          const auto torqueBx = rotatedSitePositionY[siteIndexB+siteB] * forceZ - rotatedSitePositionZ[siteIndexB+siteB] * forceY;
          const auto torqueBy = rotatedSitePositionZ[siteIndexB+siteB] * forceX - rotatedSitePositionX[siteIndexB+siteB] * forceZ;
          const auto torqueBz = rotatedSitePositionX[siteIndexB+siteB] * forceY - rotatedSitePositionY[siteIndexB+siteB] * forceX;

          txptr[molB] -= torqueBx;
          typtr[molB] -= torqueBy;
          tzptr[molB] -= torqueAz;

          if (calculateGlobals) {
            // todo this
          }
        }
      }
      // add all forces + torques for molA
      fxptr[molA] += forceSumX;
      fyptr[molA] += forceSumY;
      fzptr[molA] += forceSumZ;

      txptr[molA] += torqueSumX;
      typtr[molA] += torqueSumY;
      tzptr[molA] += torqueSumZ;

    }
    if (calculateGlobals) {
      // todo this
    }
  }
  /**
   * @copydoc Functor::SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3)
   */
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2, const bool newton3) final {
    if (newton3) {
      SoAFunctorPairImpl<true>(soa1,soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1,soa2);
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
      _shift6 = ParticlePropertiesLibrary<double,size_t>::calcShift6(_epsilon24,_sigmaSquared,_cutoffSquared);
    } else {
      _shift6 = 0;
    }
  }

  /**
   * @return useMixing
   */
  constexpr static bool getMixing() { return  useMixing; }

  /**
   * Get the number of flops used per kernel call - i.e. number of flops to calculate kernel *given* the two particles
   * lie within the cutoff (i.e. distance^2 / cutoff has been already been calculated).
   * @param numA number of sites in molecule A
   * @param numB number of sites in molecule B
   * @return #FLOPs
   */
  static unsigned long getNumFlopsPerKernelCall(bool newton3, size_t numA,size_t numB) {
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
    for (size_t i=0; i < _aosThreadData.size(); i++) {
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
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal()."
          );
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
        _virialSum = autopas::utils::ArrayMath::mulScalar(_virialSum,0.5);
      }
      // we have always calculated 6*potential, so we divide by 6 here!
      _potentialEnergySum /= 6.;
      _postProcessed = true;
    }
  }

  /**
   * Get the potential energy.
   * @return the potential energy sum
   */
  double getPotentialEnergy() {
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
    if (soaA.getNumParticles() == 0 || soaB.getNumParticles() == 0) return;

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

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedSitePositionAx;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedSitePositionAy;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedSitePositionAz;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedSitePositionBx;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedSitePositionBy;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> rotatedSitePositionBz;

    const SoAFloatPrecision const_sigmaSquared = _sigmaSquared;
    const SoAFloatPrecision const_epsilon24 = _epsilon24;
    const SoAFloatPrecision const_shift6 = _shift6;

    if constexpr (useMixing) {
      // todo implement this
      autopas::utils::ExceptionHandler::exception("LJMulticenterFunctor: Mixing with multicentered molecules not yet implemented");
    }

    // load unrotated site positions
    auto unrotatedSitePositionsA = _sitePositionsLJ;
    auto unrotatedSitePositionsB = _sitePositionsLJ;

    // calculate rotated site positions
    if constexpr (useMixing) {
      // todo this
    } else {
      rotatedSitePositionAx.resize(unrotatedSitePositionsA.size() * soaA.getNumParticles());
      rotatedSitePositionAy.resize(unrotatedSitePositionsA.size() * soaA.getNumParticles());
      rotatedSitePositionAz.resize(unrotatedSitePositionsA.size() * soaA.getNumParticles());

      for (size_t mol = 0; mol < soaA.getNumParticles(); mol++) {
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0Aptr[mol], q1Aptr[mol],q2Aptr[mol],q3Aptr[mol]}, unrotatedSitePositionsA);
        for (size_t site = 0; site < unrotatedSitePositionsA.size(); ++site) {
          rotatedSitePositionAx[mol * unrotatedSitePositionsA.size() + site] = rotatedSitePositions[site][0];
          rotatedSitePositionAy[mol * unrotatedSitePositionsA.size() + site] = rotatedSitePositions[site][1];
          rotatedSitePositionAz[mol * unrotatedSitePositionsA.size() + site] = rotatedSitePositions[site][2];
        }
      }

      rotatedSitePositionBx.resize(unrotatedSitePositionsB.size() * soaB.getNumParticles());
      rotatedSitePositionBy.resize(unrotatedSitePositionsB.size() * soaB.getNumParticles());
      rotatedSitePositionBz.resize(unrotatedSitePositionsB.size() * soaB.getNumParticles());

      for (size_t mol = 0; mol < soaB.getNumParticles(); mol++) {
        const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
            {q0Bptr[mol], q1Bptr[mol],q2Bptr[mol],q3Bptr[mol]}, unrotatedSitePositionsB);
        for (size_t site = 0; site < unrotatedSitePositionsB.size(); ++site) {
          rotatedSitePositionBx[mol * unrotatedSitePositionsB.size() + site] = rotatedSitePositions[site][0];
          rotatedSitePositionBy[mol * unrotatedSitePositionsB.size() + site] = rotatedSitePositions[site][1];
          rotatedSitePositionBz[mol * unrotatedSitePositionsB.size() + site] = rotatedSitePositions[site][2];
        }
      }

    }

    // main force calculation loop
    for (size_t molA = 0; molA < soaA.getNumParticles(); ++molA) {
      const auto ownedStateA = ownedStatePtrA[molA];
      if (ownedStateA == autopas::OwnershipState::dummy) {
        continue;
      }

      size_t noSitesMolA = 1; // Number of sites in molecule A
      size_t noSitesB = 0; // Number of sites in cell B
      size_t siteIndexA = 0; // Index at which the number of sites in A starts
      if constexpr (useMixing) {
        // todo this
      } else {
        noSitesMolA = unrotatedSitePositionsA.size();
        noSitesB = unrotatedSitePositionsB.size();
        siteIndexA = molA*noSitesMolA;
      }

      // sums used for molA
      SoAFloatPrecision forceSumX = 0.;
      SoAFloatPrecision forceSumY = 0.;
      SoAFloatPrecision forceSumZ = 0.;
      SoAFloatPrecision torqueSumX = 0.;
      SoAFloatPrecision torqueSumY = 0.;
      SoAFloatPrecision torqueSumZ = 0.;

      // preload sigmas, epsilons, and shifts
      if constexpr (useMixing) {
        // todo this
      }

      // create mask over every site in SoA B
      std::vector<bool, autopas::AlignedAllocator<bool>> mask;

      if constexpr (useMixing) {

      } else {
        mask.resize(noSitesB);

#pragma omp simd
        // mask sites of molecules with CoM outside cutoff
        for (size_t siteB = 0; siteB < noSitesB; ++siteB) {
          size_t molB = siteB / unrotatedSitePositionsB.size();

          const auto ownedStateB = ownedStatePtrB[molB];

          const auto displacementCoMX = xAptr[molA] - xBptr[molB];
          const auto displacementCoMY = yAptr[molA] - yBptr[molB];
          const auto displacementCoMZ = zAptr[molA] - zBptr[molB];

          const auto distanceSquaredCoMX = displacementCoMX * displacementCoMX;
          const auto distanceSquaredCoMY = displacementCoMY * displacementCoMY;
          const auto distanceSquaredCoMZ = displacementCoMZ * displacementCoMZ;

          const auto distanceSquaredCoM = distanceSquaredCoMX + distanceSquaredCoMY + distanceSquaredCoMZ;

          // mask sites of molecules beyond cutoff or if molecule is a dummy
          mask[siteB] = distanceSquaredCoM <= cutoffSquared and ownedStateB != autopas::OwnershipState::dummy;
        }
      }

      // calculate LJ potentials
      for (int siteA = 0; siteA < noSitesMolA; ++siteA) {
#pragma omp simd reduction (+ : forceSumX, forceSumY, forceSumZ, torqueSumX, torqueSumY, torqueSumZ, potentialEnergySum, virialSumX, virialSumY, virialSumZ)
        for (int siteB = 0; siteB < noSitesB; ++siteB) {
          size_t molB = siteB / unrotatedSitePositionsB.size();

          const auto displacementX = xAptr[molA] - xBptr[molB] + rotatedSitePositionAx[siteIndexA+siteA] - rotatedSitePositionBx[siteB];
          const auto displacementY = yAptr[molA] - yBptr[molB] + rotatedSitePositionAy[siteIndexA+siteA] - rotatedSitePositionBy[siteB];
          const auto displacementZ = zAptr[molA] - zBptr[molB] + rotatedSitePositionAz[siteIndexA+siteA] - rotatedSitePositionBz[siteB];

          const auto distanceSquaredX = displacementX * displacementX;
          const auto distanceSquaredY = displacementY * displacementY;
          const auto distanceSquaredZ = displacementZ * displacementZ;

          const auto distanceSquared = distanceSquaredX + distanceSquaredY + distanceSquaredZ;

          const auto invDistSquared = 1. / distanceSquared;
          const auto lj2 = const_sigmaSquared * invDistSquared;
          const auto lj6 = lj2 * lj2 * lj2;
          const auto lj12 = lj6 * lj6;
          const auto lj12m6 = lj12 - lj6;
          const auto scalarMultiple = mask[siteB] ? const_epsilon24 * (lj12 + lj12m6) * invDistSquared : 0.;

          // calculate forces
          const auto forceX = scalarMultiple * displacementX;
          const auto forceY = scalarMultiple * displacementY;
          const auto forceZ = scalarMultiple * displacementZ;

          forceSumX += forceX;
          forceSumY += forceY;
          forceSumZ += forceZ;

          if (newton3) {
            fxBptr[molB] -= forceX;
            fyBptr[molB] -= forceY;
            fzBptr[molB] -= forceZ;
          }

          // calculate torques on molA
          const auto torqueAx = rotatedSitePositionAy[siteIndexA+siteA] * forceZ - rotatedSitePositionAz[siteIndexA+siteA] * forceY;
          const auto torqueAy = rotatedSitePositionAz[siteIndexA+siteA] * forceX - rotatedSitePositionAx[siteIndexA+siteA] * forceZ;
          const auto torqueAz = rotatedSitePositionAx[siteIndexA+siteA] * forceY - rotatedSitePositionAy[siteIndexA+siteA] * forceX;

          torqueSumX += torqueAx;
          torqueSumY += torqueAy;
          torqueSumZ += torqueAz;

          if (newton3) {
            const auto torqueBx = rotatedSitePositionBy[siteB] * forceZ - rotatedSitePositionBz[siteB] * forceY;
            const auto torqueBy = rotatedSitePositionBz[siteB] * forceX - rotatedSitePositionBx[siteB] * forceZ;
            const auto torqueBz = rotatedSitePositionBx[siteB] * forceY - rotatedSitePositionBy[siteB] * forceX;

            txBptr[molB] -= torqueBx;
            tyBptr[molB] -= torqueBy;
            tzBptr[molB] -= torqueAz;
          }

          if (calculateGlobals) {
            // todo this
          }

        }
      }
      // add all forces + torques for molA
      fxAptr[molA] += forceSumX;
      fyAptr[molA] += forceSumY;
      fzAptr[molA] += forceSumZ;

      txAptr[molA] += torqueSumX;
      tyAptr[molA] += torqueSumY;
      tzAptr[molA] += torqueSumZ;

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