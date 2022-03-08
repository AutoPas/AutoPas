/**
 * @file LJMulticenterFunctor.h
 * @date 21/02/2022
 * @author S. Newcome
*/

#pragma once

#include "src/Particles/MulticenteredMoleculeLJ.h"
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
template <class MulticenteredMolecule, bool applyShift = false, bool useMixing = false, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGlobals = false, bool relevantForTuning = true>
class LJMulticenterFunctor
    : public autopas::Functor<MulticenteredMolecule,
    LJMulticenterFunctor<MulticenteredMolecule, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>> {
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename MulticenteredMolecule::SoAArraysType;

  /**
   * Precision of SoA entries
   */
  using SoAFloatPrecision = typename MulticenteredMolecule::ParticleSoAFloatPrecision;

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
    : autopas::Functor <MulticenteredMolecule, LJMulticenterFunctor<MulticenteredMolecule, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(
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
  void AoSFunctor(MulticenteredMolecule &particleA, MulticenteredMolecule &particleB, bool newton3) final {
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

    const auto displacementCoM = autopas::utils::ArrayMath::sub(particleA.get_R(), particleB.get_R());
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
        particleA.addT(autopas::utils::ArrayMath::cross(rotatedSitePositionsA[m],force));
        if (newton3) {particleB.subT(autopas::utils::ArrayMath::cross(rotatedSitePositionsB[n],force));}

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
        _aosThreadData[threadNum].potentialSum += potentialEnergy;
        _aosThreadData[threadNum].virialSum = autopas::utils::ArrayMath::add(_aosThreadData[threadNum].virialSum, virial);
      }
      if (newton3 and particleB.isOwned()) {
        _aosThreadData[threadNum].potentialSum += potentialEnergy;
        _aosThreadData[threadNum].virialSum = autopas::utils::ArrayMath::add(_aosThreadData[threadNum].virialSum, virial);
      }
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