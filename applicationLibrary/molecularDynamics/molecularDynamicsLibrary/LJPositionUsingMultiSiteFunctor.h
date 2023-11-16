/**
* @file PositionUsingMultiSiteFunctor.h
* @date 10/10/2023
* @author Johannes Riemenschneider
 */

#pragma once

#include "MoleculeLJ.h"
#include "PositionStoringMultiSiteMolecule.h"
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
* @warning: Whilst this class allows for mixing to be disabled, this feature is not of much value in real applications
* and as such is experimental only and untested!
* @tparam useNewton3 Switch for the functor to support newton3 on, off, or both. See FunctorN3Nodes for possible
* values.
* @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
* @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
*/
template <class Particle, bool applyShift = false, bool useMixing = false,
         autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
         bool relevantForTuning = true>
class LJPositionUsingMultiSiteFunctor
   : public autopas::Functor<Particle, LJPositionUsingMultiSiteFunctor<Particle, applyShift, useMixing, useNewton3,
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
 LJPositionUsingMultiSiteFunctor() = delete;

private:
 /**
  * Internal (actually used) constructor
  * @param cutoff
  * @note param dummy is unused, only there to make the signature different from the public constructor.
  */
 explicit LJPositionUsingMultiSiteFunctor(SoAFloatPrecision cutoff, void * /*dummy*/)
     : autopas::Functor<Particle, LJPositionUsingMultiSiteFunctor<Particle, applyShift, useMixing, useNewton3, calculateGlobals,
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

 std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> flatten2DimVector(const std::vector<SoAFloatPrecision> *const __restrict unflattened_vector, size_t size_unflattened_vector){
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> returnvalue{};
   for(size_t i{0}; i<size_unflattened_vector; i++){
     for(size_t j{0}; j< unflattened_vector[i].size(); j++){
       returnvalue.push_back(unflattened_vector[i][j]);
     }
   }
   return std::move(returnvalue);
 }

public:
 /**
  * Constructor for Functor with particle mixing disabled. setParticleProperties() must be called.
  * @note Only to be used with mixing == false
  * @param cutoff
  */
 explicit LJPositionUsingMultiSiteFunctor(double cutoff) : LJPositionUsingMultiSiteFunctor(cutoff, nullptr) {
   static_assert(not useMixing,
                 "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                 "mixing to false.");
   AutoPasLog(WARN, "Using LJPositionUsingMultiSiteFunctor with mixing disabled is untested!");
 }

 /**
  * Constructor for Functor with particle mixing enabled.
  * Calculating global attributes is done with the center of mass and overall forces applied
  * @param cutoff
  * @param particlePropertiesLibrary Library used to look up the properties of each type of particle e.g. sigma,
  * epsilon, shift.
  */
 explicit LJPositionUsingMultiSiteFunctor(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
     : LJPositionUsingMultiSiteFunctor(cutoff, nullptr) {
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
   const auto paR = particleA.getR();
   const auto pbR = particleB.getR();
   const auto displacementCoM = autopas::utils::ArrayMath::sub(particleA.getR(), particleB.getR());
   const auto distanceSquaredCoM = autopas::utils::ArrayMath::dot(displacementCoM, displacementCoM);

   if (distanceSquaredCoM > _cutoffSquared) {
     return;
   }

   // get number of sites
   //const size_t numSitesA = useMixing ? _PPLibrary->getNumSites(particleA.getTypeId()) : _sitePositionsLJ.size();
   //const size_t numSitesB = useMixing ? _PPLibrary->getNumSites(particleB.getTypeId()) : _sitePositionsLJ.size();
   const size_t numSitesA = particleA.getNumberOfSites();
   const size_t numSitesB = particleB.getNumberOfSites();

   // get siteIds
   const std::vector<size_t> siteIdsA =
       useMixing ? _PPLibrary->getSiteTypes(particleA.getTypeId()) : std::vector<unsigned long>();
   const std::vector<size_t> siteIdsB =
       useMixing ? _PPLibrary->getSiteTypes(particleB.getTypeId()) : std::vector<unsigned long>();

   // // get unrotated relative site positions
   // const std::vector<std::array<double, 3>> unrotatedSitePositionsA =
   //     useMixing ? _PPLibrary->getSitePositions(particleA.getTypeId()) : _sitePositionsLJ;
   // const std::vector<std::array<double, 3>> unrotatedSitePositionsB =
   //     useMixing ? _PPLibrary->getSitePositions(particleB.getTypeId()) : _sitePositionsLJ;

   // // calculate correctly rotated relative site positions
   // const auto rotatedSitePositionsA =
   //     autopas::utils::quaternion::rotateVectorOfPositions(particleA.getQuaternion(), unrotatedSitePositionsA);
   // const auto rotatedSitePositionsB =
   //     autopas::utils::quaternion::rotateVectorOfPositions(particleB.getQuaternion(), unrotatedSitePositionsB);

   for (int i = 0; i < numSitesA; i++) {
     for (int j = 0; j < numSitesB; j++) {
       const auto relSitePosA = particleA.getRelativeSitePosition(i);
       const auto relSitePosB = particleB.getRelativeSitePosition(j);
       const auto absSitePosA = autopas::utils::ArrayMath::add(particleA.getR(), relSitePosA);
       const auto absSitePosB = autopas::utils::ArrayMath::add(particleB.getR(), relSitePosB);
       const auto displacement = autopas::utils::ArrayMath::sub(absSitePosA, absSitePosB);
       //const auto displacement = autopas::utils::ArrayMath::add(
       //    autopas::utils::ArrayMath::sub(displacementCoM, rotatedSitePositionsB[j]), rotatedSitePositionsA[i]);
       const auto distanceSquared = autopas::utils::ArrayMath::dot(displacement, displacement);

       //if(!std::equal(displacement_copy.begin(), displacement_copy.end(), displacement.begin(), [&](auto lhs, auto rhs){
       //      return (((lhs-rhs) < 0.001) and ((-(lhs-rhs))<0.001));
       //    })){
       //  std::cout << "absSitePositions didn't work as they should" << std::endl;
       //  exit(0);
       //}

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
       particleA.addTorque(autopas::utils::ArrayMath::cross(relSitePosA, force));
       if (newton3) {
         particleB.subTorque(autopas::utils::ArrayMath::cross(relSitePosB, force));
       }

       if (calculateGlobals) {
         // Here we calculate the potential energy * 6.
         // For newton3, this potential energy contribution is distributed evenly to the two molecules.
         // For non-newton3, the full potential energy is added to the one molecule.
         // The division by 6 is handled in endTraversal, as well as the division by two needed if newton3 is not used.
         // There is a similar handling of the virial, but without the mutliplication/division by 6.
         const auto potentialEnergy6 = epsilon24 * lj12m6 + shift6;
         const auto virial = displacement * force;

         const auto threadNum = autopas::autopas_get_thread_num();

         if (particleA.isOwned()) {
           if (newton3) {
             _aosThreadData[threadNum].potentialEnergySumN3 += potentialEnergy6 * 0.5;
             _aosThreadData[threadNum].virialSumN3 += virial * 0.5;
           } else {
             // for non-newton3 the division is in the post-processing step.
             _aosThreadData[threadNum].potentialEnergySumNoN3 += potentialEnergy6;
             _aosThreadData[threadNum].virialSumNoN3 += virial;
           }
         }
         // for non-newton3 the second particle will be considered in a separate calculation
         if (newton3 and particleB.isOwned()) {
           _aosThreadData[threadNum].potentialEnergySumN3 += potentialEnergy6 * 0.5;
           _aosThreadData[threadNum].virialSumN3 += virial * 0.5;
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
   //std::cout<< "SoaFunctor single called" << std::endl;
   if (soa.size() == 0) return;

   const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
   const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
   const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

   const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

   //const auto *const __restrict q0ptr = soa.template begin<Particle::AttributeNames::quaternion0>();
   //const auto *const __restrict q1ptr = soa.template begin<Particle::AttributeNames::quaternion1>();
   //const auto *const __restrict q2ptr = soa.template begin<Particle::AttributeNames::quaternion2>();
   //const auto *const __restrict q3ptr = soa.template begin<Particle::AttributeNames::quaternion3>();

#if not defined(MD_FLEXIBLE_ABSOLUTE_POS_USE_FLATTENING)
   const auto *const __restrict relSitePosXptr = soa.template begin<Particle::AttributeNames::relativeSitePositionsX>();
   const auto *const __restrict relSitePosYptr = soa.template begin<Particle::AttributeNames::relativeSitePositionsY>();
   const auto *const __restrict relSitePosZptr = soa.template begin<Particle::AttributeNames::relativeSitePositionsZ>();
#else
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> relSitePosXptrFlattened = flatten2DimVector(soa.template begin<Particle::AttributeNames::relativeSitePositionsX>(), soa.size());
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> relSitePosYptrFlattened = flatten2DimVector(soa.template begin<Particle::AttributeNames::relativeSitePositionsY>(), soa.size());
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> relSitePosZptrFlattened = flatten2DimVector(soa.template begin<Particle::AttributeNames::relativeSitePositionsZ>(), soa.size());
#endif

   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> absoluteSitePositionsX;
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> absoluteSitePositionsY;
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> absoluteSitePositionsZ;

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
     for (size_t mol = 0; mol < soa.size(); ++mol) {
       siteCount += _PPLibrary->getNumSites(typeptr[mol]);
     }
   } else {
     //we would have returned earlier if this array was empty
     std::cout << "Johnny, you still have to rewrite this case!" << std::endl; //@TODO
     exit(0);
     //siteCount = relSitePosXptr[0].size() * soa.size();
   }

   // pre-reserve site std::vectors
   absoluteSitePositionsX.reserve(siteCount);
   absoluteSitePositionsY.reserve(siteCount);
   absoluteSitePositionsZ.reserve(siteCount);

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
#if defined(MD_FLEXIBLE_ABSOLUTE_POS_USE_FLATTENING)
   {
     size_t siteIndex = 0;
     for (size_t mol = 0; mol < soa.size(); ++mol) {
       // const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
       //     {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, _PPLibrary->getSitePositions(typeptr[mol]));
       const auto siteTypesOfMol = _PPLibrary->getSiteTypes(typeptr[mol]);

       for (size_t site = 0; site < siteTypesOfMol.size(); ++site) {
         absoluteSitePositionsX[siteIndex] = relSitePosXptrFlattened[siteIndex] + xptr[mol];
         absoluteSitePositionsY[siteIndex] = relSitePosXptrFlattened[siteIndex] + yptr[mol];
         absoluteSitePositionsZ[siteIndex] = relSitePosXptrFlattened[siteIndex] + zptr[mol];
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
   }
#else
   {
     size_t siteIndex = 0;
     for (size_t mol = 0; mol < soa.size(); ++mol) {
       const auto relativeSitePositionsX = relSitePosXptr[mol];
       const auto relativeSitePositionsY = relSitePosYptr[mol];
       const auto relativeSitePositionsZ = relSitePosZptr[mol];

       // const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
       //     {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, _PPLibrary->getSitePositions(typeptr[mol]));
       const auto siteTypesOfMol = _PPLibrary->getSiteTypes(typeptr[mol]);

       for (size_t site = 0; site < relativeSitePositionsX.size(); ++site) {
         absoluteSitePositionsX[siteIndex] = relativeSitePositionsX[site] + xptr[mol];
         absoluteSitePositionsY[siteIndex] = relativeSitePositionsY[site] + yptr[mol];
         absoluteSitePositionsZ[siteIndex] = relativeSitePositionsZ[site] + zptr[mol];
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
   }
#endif

   // main force calculation loop
   size_t siteIndexMolA = 0;  // index of first site in molA
   for (size_t molA = 0; molA < soa.size(); ++molA) {
#if defined(MD_FLEXIBLE_ABSOLUTE_POS_USE_FLATTENING)
     const size_t noSitesInMolA = useMixing ? _PPLibrary->getNumSites(typeptr[molA])
                                            : const_unrotatedSitePositions.size();  // Number of sites in molecule A
#else
     const size_t noSitesInMolA = relSitePosXptr[molA].size();
#endif

     const auto ownedStateA = ownedStatePtr[molA];
     if (ownedStateA == autopas::OwnershipState::dummy) {
       siteIndexMolA += noSitesInMolA;
       continue;
     }

     const size_t siteIndexMolB = siteIndexMolA + noSitesInMolA;  // index of first site in molB
     const size_t noSitesB = (siteCount - siteIndexMolB);         // Number of sites in molecules that A interacts with

     // create mask over every mol 'above' molA  (char to keep arrays aligned)
     std::vector<char, autopas::AlignedAllocator<char>> molMask;
     molMask.reserve(soa.size() - (molA + 1));

#pragma omp simd
     for (size_t molB = molA + 1; molB < soa.size(); ++molB) {
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

     for (size_t molB = molA + 1; molB < soa.size(); ++molB) {
#if defined(MD_FLEXIBLE_ABSOLUTE_POS_USE_FLATTENING)
       const size_t noSitesInMolB = useMixing ? _PPLibrary->getNumSites(typeptr[molB])
                                              : const_unrotatedSitePositions.size();  // Number of sites in molecule A
#else
       const size_t noSitesInMolB = relSitePosXptr[molB].size();
#endif
       for (size_t siteB = 0; siteB < noSitesInMolB; ++siteB) {
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

         const auto displacementX = absoluteSitePositionsX[siteA] - absoluteSitePositionsX[globalSiteBIndex];
         const auto displacementY = absoluteSitePositionsY[siteA] - absoluteSitePositionsY[globalSiteBIndex];
         const auto displacementZ = absoluteSitePositionsZ[siteA] - absoluteSitePositionsZ[globalSiteBIndex];

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
#if defined(MD_FLEXIBLE_ABSOLUTE_POS_USE_FLATTENING)
   {
     size_t siteIndex = 0;
     for (size_t mol = 0; mol < soa.size(); mol++) {
       //const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
       //    {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, const_unrotatedSitePositions);
       const size_t noSitesInMol = useMixing ? _PPLibrary->getNumSites(typeptr[mol])
                                              : const_unrotatedSitePositions.size();  // Number of sites in molecule A
       for (size_t site = 0; site < noSitesInMol; ++site) {
         fxptr[mol] += siteForceX[siteIndex];
         fyptr[mol] += siteForceY[siteIndex];
         fzptr[mol] += siteForceZ[siteIndex];
         txptr[mol] += relSitePosYptrFlattened[siteIndex] * siteForceZ[siteIndex] -
                       relSitePosZptrFlattened[siteIndex] * siteForceY[siteIndex];
         typtr[mol] += relSitePosZptrFlattened[siteIndex] * siteForceX[siteIndex] -
                       relSitePosXptrFlattened[siteIndex] * siteForceZ[siteIndex];
         tzptr[mol] += relSitePosXptrFlattened[siteIndex] * siteForceY[siteIndex] -
                       relSitePosYptrFlattened[siteIndex] * siteForceX[siteIndex];
         ++siteIndex;
       }
     }
   }
#else
   {
     size_t siteIndex = 0;
     for (size_t mol = 0; mol < soa.size(); mol++) {
       const auto relativeSitePositionsX = relSitePosXptr[mol];
       const auto relativeSitePositionsY = relSitePosYptr[mol];
       const auto relativeSitePositionsZ = relSitePosZptr[mol];
       //const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
       //    {q0ptr[mol], q1ptr[mol], q2ptr[mol], q3ptr[mol]}, const_unrotatedSitePositions);
       for (size_t site = 0; site < relSitePosXptr[mol].size(); ++site) {
         fxptr[mol] += siteForceX[siteIndex];
         fyptr[mol] += siteForceY[siteIndex];
         fzptr[mol] += siteForceZ[siteIndex];
         txptr[mol] += relativeSitePositionsY[site] * siteForceZ[siteIndex] -
                       relativeSitePositionsZ[site] * siteForceY[siteIndex];
         typtr[mol] += relativeSitePositionsZ[site] * siteForceX[siteIndex] -
                       relativeSitePositionsX[site] * siteForceZ[siteIndex];
         tzptr[mol] += relativeSitePositionsX[site] * siteForceY[siteIndex] -
                       relativeSitePositionsY[site] * siteForceX[siteIndex];
         ++siteIndex;
       }
     }
   }
#endif

   if constexpr (calculateGlobals) {
     const auto threadNum = autopas::autopas_get_thread_num();
     // SoAFunctorSingle obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
     // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.
     if (newton3) {
       _aosThreadData[threadNum].potentialEnergySumN3 += potentialEnergySum * 0.5;
       _aosThreadData[threadNum].virialSumN3[0] += virialSumX * 0.5;
       _aosThreadData[threadNum].virialSumN3[1] += virialSumY * 0.5;
       _aosThreadData[threadNum].virialSumN3[2] += virialSumZ * 0.5;
     } else {
       _aosThreadData[threadNum].potentialEnergySumNoN3 += potentialEnergySum;
       _aosThreadData[threadNum].virialSumNoN3[0] += virialSumX;
       _aosThreadData[threadNum].virialSumNoN3[1] += virialSumY;
       _aosThreadData[threadNum].virialSumNoN3[2] += virialSumZ;
     }
   }
 }

 /**
  * @copydoc autopas::Functor::SoAFunctorPair()
  */
 void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                     const bool newton3) final {
   //std::cout<< "SoaFunctorPair called" << std::endl;

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
   //std::cout << "SoaFunctorVerlet called "<< std::endl;
   if (soa.size() == 0 or neighborList.empty()) return;
   if (newton3) {
     SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
   } else {
     SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
   }
 }

 /**
  * Sets the molecule properties constants for this functor.
  *
  * This is only necessary if no particlePropertiesLibrary is used.
  * @param epsilon24 epsilon * 24
  * @param sigmaSquared sigma^2
  * @param sitePositionsLJ vector of 3D relative unrotated untranslated site positions
  */
 void setParticleProperties(SoAFloatPrecision epsilon24, SoAFloatPrecision sigmaSquared,
                            std::vector<std::array<SoAFloatPrecision, 3>> sitePositionsLJ) {
   _epsilon24 = epsilon24;
   _sigmaSquared = sigmaSquared;
   if (applyShift) {
     _shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(_epsilon24, _sigmaSquared, _cutoffSquared);
   } else {
     _shift6 = 0;
   }
   _sitePositionsLJ = sitePositionsLJ;
 }

 /**
  * @copydoc autopas::Functor::getNeededAttr()
  */
 constexpr static auto getNeededAttr() {
   return std::array<typename Particle::AttributeNames, 19>{
       Particle::AttributeNames::id,          Particle::AttributeNames::posX,
       Particle::AttributeNames::posY,        Particle::AttributeNames::posZ,
       Particle::AttributeNames::forceX,      Particle::AttributeNames::forceY,
       Particle::AttributeNames::forceZ,      Particle::AttributeNames::quaternion0,
       Particle::AttributeNames::quaternion1, Particle::AttributeNames::quaternion2,
       Particle::AttributeNames::quaternion3, Particle::AttributeNames::relativeSitePositionsX,
       Particle::AttributeNames::relativeSitePositionsY, Particle::AttributeNames::relativeSitePositionsZ,
       Particle::AttributeNames::torqueX,
       Particle::AttributeNames::torqueY,     Particle::AttributeNames::torqueZ,
       Particle::AttributeNames::typeId,      Particle::AttributeNames::ownershipState};
 }

 /**
  * @copydoc autopas::Functor::getNeededAttr(std::false_type)
  */
 constexpr static auto getNeededAttr(std::false_type) {
   return std::array<typename Particle::AttributeNames, 19>{
       Particle::AttributeNames::id,          Particle::AttributeNames::posX,
       Particle::AttributeNames::posY,        Particle::AttributeNames::posZ,
       Particle::AttributeNames::forceX,      Particle::AttributeNames::forceY,
       Particle::AttributeNames::forceZ,      Particle::AttributeNames::quaternion0,
       Particle::AttributeNames::quaternion1, Particle::AttributeNames::quaternion2,
       Particle::AttributeNames::quaternion3, Particle::AttributeNames::relativeSitePositionsX,
       Particle::AttributeNames::relativeSitePositionsY, Particle::AttributeNames::relativeSitePositionsZ,
       Particle::AttributeNames::torqueX,
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
     // We distinguish between non-newton3 and newton3 functor calls. Newton3 calls are accumulated directly.
     // Non-newton3 calls are accumulated temporarily and later divided by 2.
     double potentialEnergySumNoN3Acc = 0;
     std::array<double, 3> virialSumNoN3Acc = {0, 0, 0};
     for (size_t i = 0; i < _aosThreadData.size(); ++i) {
       potentialEnergySumNoN3Acc += _aosThreadData[i].potentialEnergySumNoN3;
       _potentialEnergySum += _aosThreadData[i].potentialEnergySumN3;

       virialSumNoN3Acc += _aosThreadData[i].virialSumNoN3;
       _virialSum += _aosThreadData[i].virialSumN3;
     }
     // if the newton3 optimization is disabled we have added every energy contribution twice, so we divide by 2
     // here.
     potentialEnergySumNoN3Acc *= 0.5;
     virialSumNoN3Acc *= 0.5;

     _potentialEnergySum += potentialEnergySumNoN3Acc;
     _virialSum += virialSumNoN3Acc;

     // we have always calculated 6*potentialEnergy, so we divide by 6 here!
     _potentialEnergySum /= 6.;
     _postProcessed = true;

     AutoPasLog(TRACE, "Final potential energy {}", _potentialEnergySum);
     AutoPasLog(TRACE, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
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
   if (soaA.size() == 0 || soaB.size() == 0) return;

   const auto *const __restrict xAptr = soaA.template begin<Particle::AttributeNames::posX>();
   const auto *const __restrict yAptr = soaA.template begin<Particle::AttributeNames::posY>();
   const auto *const __restrict zAptr = soaA.template begin<Particle::AttributeNames::posZ>();
   const auto *const __restrict xBptr = soaB.template begin<Particle::AttributeNames::posX>();
   const auto *const __restrict yBptr = soaB.template begin<Particle::AttributeNames::posY>();
   const auto *const __restrict zBptr = soaB.template begin<Particle::AttributeNames::posZ>();

   const auto *const __restrict ownedStatePtrA = soaA.template begin<Particle::AttributeNames::ownershipState>();
   const auto *const __restrict ownedStatePtrB = soaB.template begin<Particle::AttributeNames::ownershipState>();

#if not defined(MD_FLEXIBLE_ABSOLUTE_POS_USE_FLATTENING)
   const auto *const __restrict relSitePosXAptr = soaA.template begin<Particle::AttributeNames::relativeSitePositionsX>();
   const auto *const __restrict relSitePosYAptr = soaA.template begin<Particle::AttributeNames::relativeSitePositionsY>();
   const auto *const __restrict relSitePosZAptr = soaA.template begin<Particle::AttributeNames::relativeSitePositionsZ>();
   const auto *const __restrict relSitePosXBptr = soaB.template begin<Particle::AttributeNames::relativeSitePositionsX>();
   const auto *const __restrict relSitePosYBptr = soaB.template begin<Particle::AttributeNames::relativeSitePositionsY>();
   const auto *const __restrict relSitePosZBptr = soaB.template begin<Particle::AttributeNames::relativeSitePositionsZ>();
#else
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> relSitePosXAptr = flatten2DimVector(soaA.template begin<Particle::AttributeNames::relativeSitePositionsX>(), soaA.size());
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> relSitePosYAptr = flatten2DimVector(soaA.template begin<Particle::AttributeNames::relativeSitePositionsY>(), soaA.size());
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> relSitePosZAptr = flatten2DimVector(soaA.template begin<Particle::AttributeNames::relativeSitePositionsZ>(), soaA.size());
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> relSitePosXBptr = flatten2DimVector(soaB.template begin<Particle::AttributeNames::relativeSitePositionsX>(), soaB.size());
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> relSitePosYBptr = flatten2DimVector(soaB.template begin<Particle::AttributeNames::relativeSitePositionsY>(), soaB.size());
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> relSitePosZBptr = flatten2DimVector(soaB.template begin<Particle::AttributeNames::relativeSitePositionsZ>(), soaB.size());
#endif

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

   const auto const_unrotatedSitePositions = _sitePositionsLJ;

   SoAFloatPrecision potentialEnergySum = 0.;
   SoAFloatPrecision virialSumX = 0.;
   SoAFloatPrecision virialSumY = 0.;
   SoAFloatPrecision virialSumZ = 0.;

   // local redeclarations to help compilers
   const SoAFloatPrecision cutoffSquared = _cutoffSquared;

   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> sigmaSquareds;
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> epsilon24s;
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> shift6s;

   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> absoluteSitePositionBx;
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> absoluteSitePositionBy;
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> absoluteSitePositionBz;

   // we require arrays for forces for sites to maintain SIMD in site-site calculations
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceBx;
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceBy;
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceBz;

   std::vector<size_t, autopas::AlignedAllocator<size_t>> siteTypesB;
   std::vector<char, autopas::AlignedAllocator<char>> isSiteOwnedBArr;

   const SoAFloatPrecision const_sigmaSquared = _sigmaSquared;
   const SoAFloatPrecision const_epsilon24 = _epsilon24;
   const SoAFloatPrecision const_shift6 = _shift6;

   // count number of sites in both SoAs
   size_t siteCountB = 0;
   if constexpr (useMixing) {
     //we would have returned sooner if soaB was empty
#if defined(MD_FLEXIBLE_ABSOLUTE_POS_USE_FLATTENING)
     siteCountB = relSitePosXBptr.size();
#else
     for (size_t mol = 0; mol < soaB.size(); ++mol) {
       //siteCountB += _PPLibrary->getNumSites(typeptrB[mol]);
       siteCountB += relSitePosXBptr[mol].size();
     }
#endif

   } else {
     siteCountB = const_unrotatedSitePositions.size() * soaB.size();
   }

   // pre-reserve std::vectors
   absoluteSitePositionBx.reserve(siteCountB);
   absoluteSitePositionBy.reserve(siteCountB);
   absoluteSitePositionBz.reserve(siteCountB);

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
#if not defined(MD_FLEXIBLE_ABSOLUTE_POS_USE_FLATTENING)
   {
     size_t siteIndex = 0;
     for (size_t mol = 0; mol < soaB.size(); ++mol) {
       const auto relativeSitePositionsX = relSitePosXBptr[mol];
       const auto relativeSitePositionsY = relSitePosYBptr[mol];
       const auto relativeSitePositionsZ = relSitePosZBptr[mol];
       //const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
       //    {q0Bptr[mol], q1Bptr[mol], q2Bptr[mol], q3Bptr[mol]}, _PPLibrary->getSitePositions(typeptrB[mol]));
       const auto siteTypesOfMol = _PPLibrary->getSiteTypes(typeptrB[mol]);

       for (size_t site = 0; site < relativeSitePositionsX.size(); ++site) {
         absoluteSitePositionBx[siteIndex] = relativeSitePositionsX[site] + xBptr[mol];
         absoluteSitePositionBy[siteIndex] = relativeSitePositionsY[site] + yBptr[mol];
         absoluteSitePositionBz[siteIndex] = relativeSitePositionsZ[site] + zBptr[mol];
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
   }
#else
   {
   size_t siteIndex = 0;
   for (size_t mol = 0; mol < soaB.size(); ++mol) {
     //const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
     //    {q0Bptr[mol], q1Bptr[mol], q2Bptr[mol], q3Bptr[mol]}, _PPLibrary->getSitePositions(typeptrB[mol]));
     const auto siteTypesOfMol = _PPLibrary->getSiteTypes(typeptrB[mol]);

     for (size_t site = 0; site < siteTypesOfMol.size(); ++site) {
       absoluteSitePositionBx[siteIndex] = relSitePosXBptr[siteIndex] + xBptr[mol];
       absoluteSitePositionBy[siteIndex] = relSitePosYBptr[siteIndex] + yBptr[mol];
       absoluteSitePositionBz[siteIndex] = relSitePosZBptr[siteIndex] + zBptr[mol];
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
 }
#endif

   // main force calculation loop

   {
#if defined MD_FLEXIBLE_ABSOLUTE_POS_USE_FLATTENING
     size_t siteIndexA{0};
#endif
     for (size_t molA = 0; molA < soaA.size(); ++molA) {
       const auto ownedStateA = ownedStatePtrA[molA];
       if (ownedStateA == autopas::OwnershipState::dummy) {
         continue;
       }

#if defined(MD_FLEXIBLE_ABSOLUTE_POS_USE_FLATTENING)
       const auto noSitesInMolA =
           useMixing ? _PPLibrary->getNumSites(typeptrA[molA]) : const_unrotatedSitePositions.size();
#else
     const auto noSitesInMolA = relSitePosXAptr[molA].size();

     const std::vector<SoAFloatPrecision> relativeSitePositionsXA = relSitePosXAptr[molA];
     const std::vector<SoAFloatPrecision> relativeSitePositionsYA = relSitePosYAptr[molA];
     const std::vector<SoAFloatPrecision> relativeSitePositionsZA = relSitePosZAptr[molA];
#endif

       // create mask over every mol in cell B (char to keep arrays aligned)
       std::vector<char, autopas::AlignedAllocator<char>> molMask;
       molMask.reserve(soaB.size());

#pragma omp simd
       for (size_t molB = 0; molB < soaB.size(); ++molB) {
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

       for (size_t molB = 0; molB < soaB.size(); ++molB) {
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

       for (size_t siteA = 0; siteA < noSitesInMolA; ++siteA
  #if defined(MD_FLEXIBLE_ABSOLUTE_POS_USE_FLATTENING)
    , ++siteIndexA
  #endif
            ) {   //take that clang-tidy    In all seriousness this is only temporary for measurements
         if (useMixing) {
           // preload sigmas, epsilons, and shifts
           for (size_t siteIndexB = 0; siteIndexB < siteCountB; ++siteIndexB) {
             const auto mixingData =
                 _PPLibrary->getMixingData(_PPLibrary->getSiteTypes(typeptrA[molA])[siteA], siteTypesB[siteIndexB]);
             sigmaSquareds[siteIndexB] = mixingData.sigmaSquared;
             epsilon24s[siteIndexB] = mixingData.epsilon24;
             if (applyShift) {
               shift6s[siteIndexB] = mixingData.shift6;
             }
           }
         }
#if defined(MD_FLEXIBLE_ABSOLUTE_POS_USE_FLATTENING)
         const SoAFloatPrecision relativeSitePositionAx = relSitePosXAptr[siteIndexA];
         const SoAFloatPrecision relativeSitePositionAy = relSitePosYAptr[siteIndexA];
         const SoAFloatPrecision relativeSitePositionAz = relSitePosZAptr[siteIndexA];
#else
         const SoAFloatPrecision relativeSitePositionAx = relativeSitePositionsXA[siteA];
         const SoAFloatPrecision relativeSitePositionAy = relativeSitePositionsYA[siteA];
         const SoAFloatPrecision relativeSitePositionAz = relativeSitePositionsZA[siteA];
#endif
         const SoAFloatPrecision absoluteSitePositionAx = relativeSitePositionAx + xAptr[molA];
         const SoAFloatPrecision absoluteSitePositionAy = relativeSitePositionAy + yAptr[molA];
         const SoAFloatPrecision absoluteSitePositionAz = relativeSitePositionAz + zAptr[molA];

#pragma omp simd reduction(+ : forceSumX, forceSumY, forceSumZ, torqueSumX, torqueSumY, torqueSumZ, potentialEnergySum, virialSumX, virialSumY, virialSumZ)
         for (size_t siteIndexB = 0; siteIndexB < siteCountB; ++siteIndexB) {
           const SoAFloatPrecision sigmaSquared = useMixing ? sigmaSquareds[siteIndexB] : const_sigmaSquared;
           const SoAFloatPrecision epsilon24 = useMixing ? epsilon24s[siteIndexB] : const_epsilon24;
           const SoAFloatPrecision shift6 = applyShift ? (useMixing ? shift6s[siteIndexB] : const_shift6) : 0;

           const auto isSiteOwnedB = !calculateGlobals || isSiteOwnedBArr[siteIndexB];

           const auto displacementX = absoluteSitePositionAx - absoluteSitePositionBx[siteIndexB];
           const auto displacementY = absoluteSitePositionAy - absoluteSitePositionBy[siteIndexB];
           const auto displacementZ = absoluteSitePositionAz - absoluteSitePositionBz[siteIndexB];

           const auto distanceSquaredX = displacementX * displacementX;
           const auto distanceSquaredY = displacementY * displacementY;
           const auto distanceSquaredZ = displacementZ * displacementZ;

           const auto distanceSquared = distanceSquaredX + distanceSquaredY + distanceSquaredZ;

           const auto invDistSquared = 1. / distanceSquared;
           const auto lj2 = sigmaSquared * invDistSquared;
           const auto lj6 = lj2 * lj2 * lj2;
           const auto lj12 = lj6 * lj6;
           const auto lj12m6 = lj12 - lj6;
           const auto scalarMultiple = siteMask[siteIndexB] ? epsilon24 * (lj12 + lj12m6) * invDistSquared : 0.;

           // calculate forces
           const auto forceX = scalarMultiple * displacementX;
           const auto forceY = scalarMultiple * displacementY;
           const auto forceZ = scalarMultiple * displacementZ;

           forceSumX += forceX;
           forceSumY += forceY;
           forceSumZ += forceZ;

           torqueSumX += relativeSitePositionAy * forceZ - relativeSitePositionAz * forceY;
           torqueSumY += relativeSitePositionAz * forceX - relativeSitePositionAx * forceZ;
           torqueSumZ += relativeSitePositionAx * forceY - relativeSitePositionAy * forceX;

           // N3L ( total molecular forces + torques to be determined later )
           if constexpr (newton3) {
             siteForceBx[siteIndexB] -= forceX;
             siteForceBy[siteIndexB] -= forceY;
             siteForceBz[siteIndexB] -= forceZ;
           }

           // globals
           if constexpr (calculateGlobals) {
             const auto potentialEnergy6 = siteMask[siteIndexB] ? (epsilon24 * lj12m6 + shift6) : 0.;
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
   }

   // reduce the forces on individual sites in SoA B to total forces & torques on whole molecules
#if defined(MD_FLEXIBLE_ABSOLUTE_POS_USE_FLATTENING)
   {
     size_t siteIndex = 0;
     for (size_t mol = 0; mol < soaB.size(); ++mol) {
       //const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
       //    {q0Bptr[mol], q1Bptr[mol], q2Bptr[mol], q3Bptr[mol]}, _PPLibrary->getSitePositions(typeptrB[mol]));
       const auto numSites = useMixing? _PPLibrary->getNumSites(typeptrB[mol]):const_unrotatedSitePositions.size();

       for (size_t site = 0; site < numSites; ++site) {
         fxBptr[mol] += siteForceBx[siteIndex];
         fyBptr[mol] += siteForceBy[siteIndex];
         fzBptr[mol] += siteForceBz[siteIndex];
         txBptr[mol] += relSitePosYBptr[siteIndex] * siteForceBz[siteIndex] -
                        relSitePosZBptr[siteIndex] * siteForceBy[siteIndex];
         tyBptr[mol] += relSitePosZBptr[siteIndex] * siteForceBx[siteIndex] -
                        relSitePosXBptr[siteIndex] * siteForceBz[siteIndex];
         tzBptr[mol] += relSitePosXBptr[siteIndex] * siteForceBy[siteIndex] -
                        relSitePosYBptr[siteIndex] * siteForceBx[siteIndex];
         ++siteIndex;
       }
     }
   }
#else
   {
     size_t siteIndex = 0;
     for (size_t mol = 0; mol < soaB.size(); ++mol) {
       //const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
       //    {q0Bptr[mol], q1Bptr[mol], q2Bptr[mol], q3Bptr[mol]}, _PPLibrary->getSitePositions(typeptrB[mol]));
       const auto relativeSitePositionsBx = relSitePosXBptr[mol];
       const auto relativeSitePositionsBy = relSitePosYBptr[mol];
       const auto relativeSitePositionsBz = relSitePosZBptr[mol];

       for (size_t site = 0; site < relativeSitePositionsBx.size(); ++site) {
         fxBptr[mol] += siteForceBx[siteIndex];
         fyBptr[mol] += siteForceBy[siteIndex];
         fzBptr[mol] += siteForceBz[siteIndex];
         txBptr[mol] += relativeSitePositionsBy[site] * siteForceBz[siteIndex] -
                        relativeSitePositionsBz[site] * siteForceBy[siteIndex];
         tyBptr[mol] += relativeSitePositionsBz[site] * siteForceBx[siteIndex] -
                        relativeSitePositionsBx[site] * siteForceBz[siteIndex];
         tzBptr[mol] += relativeSitePositionsBx[site] * siteForceBy[siteIndex] -
                        relativeSitePositionsBy[site] * siteForceBx[siteIndex];
         ++siteIndex;
       }
     }
   }
#endif
   if constexpr (calculateGlobals) {
     const auto threadNum = autopas::autopas_get_thread_num();
     // SoAFunctorPairImpl obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
     // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.
     if constexpr (newton3) {
       _aosThreadData[threadNum].potentialEnergySumN3 += potentialEnergySum * 0.5;
       _aosThreadData[threadNum].virialSumN3[0] += virialSumX * 0.5;
       _aosThreadData[threadNum].virialSumN3[1] += virialSumY * 0.5;
       _aosThreadData[threadNum].virialSumN3[2] += virialSumZ * 0.5;
     } else {
       _aosThreadData[threadNum].potentialEnergySumNoN3 += potentialEnergySum;
       _aosThreadData[threadNum].virialSumNoN3[0] += virialSumX;
       _aosThreadData[threadNum].virialSumNoN3[1] += virialSumY;
       _aosThreadData[threadNum].virialSumNoN3[2] += virialSumZ;
     }
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

   const size_t neighborListSize = neighborList.size();
   const size_t *const __restrict neighborListPtr = neighborList.data();
   if(neighborListSize == 0){
     return;
   }


   const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
   const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
   const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

   //const auto *const __restrict q0ptr = soa.template begin<Particle::AttributeNames::quaternion0>();
   //const auto *const __restrict q1ptr = soa.template begin<Particle::AttributeNames::quaternion1>();
   //const auto *const __restrict q2ptr = soa.template begin<Particle::AttributeNames::quaternion2>();
   //const auto *const __restrict q3ptr = soa.template begin<Particle::AttributeNames::quaternion3>();
   const auto *const __restrict relSitePosXptr = soa.template begin<Particle::AttributeNames::relativeSitePositionsX>();
   const auto *const __restrict relSitePosYptr = soa.template begin<Particle::AttributeNames::relativeSitePositionsY>();
   const auto *const __restrict relSitePosZptr = soa.template begin<Particle::AttributeNames::relativeSitePositionsZ>();

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
   //const auto const_unrotatedSitePositions = _sitePositionsLJ;

   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> sigmaSquareds;
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> epsilon24s;
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> shift6s;

   const auto const_sigmaSquared = _sigmaSquared;
   const auto const_epsilon24 = _epsilon24;
   const auto const_shift6 = _shift6;

   // Count sites
   //const size_t siteCountMolPrime =
   //    useMixing ? _PPLibrary->getNumSites(typeptr[indexPrime]) : const_unrotatedSitePositions.size();
   const size_t siteCountMolPrime = relSitePosXptr[indexPrime].size();

   size_t siteCountNeighbors = 0;  // site count of neighbours of primary molecule
   if constexpr (useMixing) {
     for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
       //siteCountNeighbors += _PPLibrary->getNumSites(typeptr[neighborList[neighborMol]]);
       siteCountNeighbors += relSitePosXptr[neighborList[neighborMol]].size();
     }
   } else {
     //siteCountNeighbors = const_unrotatedSitePositions.size() * neighborListSize;
     //would have returned sooner if list was empty
     siteCountNeighbors = relSitePosXptr[neighborList[0]].size() * neighborListSize;
   }

   // initialize site-wise arrays for neighbors
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> absoluteNeighborSitePositionX;
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> absoluteNeighborSitePositionY;
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> absoluteNeighborSitePositionZ;

   std::vector<size_t, autopas::AlignedAllocator<size_t>> siteTypesNeighbors;
   std::vector<char, autopas::AlignedAllocator<char>> isNeighborSiteOwnedArr;

   // we require arrays for forces for sites to maintain SIMD in site-site calculations
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceX;
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceY;
   std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> siteForceZ;

   // pre-reserve arrays
   absoluteNeighborSitePositionX.reserve(siteCountNeighbors);
   absoluteNeighborSitePositionY.reserve(siteCountNeighbors);
   absoluteNeighborSitePositionZ.reserve(siteCountNeighbors);

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

   //const auto rotatedSitePositionsPrime =
   //    useMixing ? autopas::utils::quaternion::rotateVectorOfPositions(
   //                    {q0ptr[indexPrime], q1ptr[indexPrime], q2ptr[indexPrime], q3ptr[indexPrime]},
   //                    _PPLibrary->getSitePositions(typeptr[indexPrime]))
   //              : autopas::utils::quaternion::rotateVectorOfPositions(
   //                    {q0ptr[indexPrime], q1ptr[indexPrime], q2ptr[indexPrime], q3ptr[indexPrime]},
   //                    const_unrotatedSitePositions);

   const auto relativeSitePositionsXPrime = relSitePosXptr[indexPrime];
   const auto relativeSitePositionsYPrime = relSitePosYptr[indexPrime];
   const auto relativeSitePositionsZPrime = relSitePosZptr[indexPrime];

   const auto siteTypesPrime = _PPLibrary->getSiteTypes(typeptr[indexPrime]);  // this is not used if non-mixing

   // generate site-wise arrays for neighbors of primary mol
   {
     size_t siteIndex = 0;
     for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
       const auto neighborMolIndex = neighborList[neighborMol];
       //const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
       //    {q0ptr[neighborMolIndex], q1ptr[neighborMolIndex], q2ptr[neighborMolIndex], q3ptr[neighborMolIndex]},
       //    _PPLibrary->getSitePositions(typeptr[neighborMolIndex]));
       const auto relativeSitePositionsX = relSitePosXptr[neighborMolIndex];
       const auto relativeSitePositionsY = relSitePosYptr[neighborMolIndex];
       const auto relativeSitePositionsZ = relSitePosZptr[neighborMolIndex];
       const auto siteTypesOfMol = _PPLibrary->getSiteTypes(typeptr[neighborMolIndex]);

       for (size_t site = 0; site < relativeSitePositionsX.size(); ++site) {
         absoluteNeighborSitePositionX[siteIndex] = relativeSitePositionsX[site] + xptr[neighborMolIndex];
         absoluteNeighborSitePositionY[siteIndex] = relativeSitePositionsY[site] + yptr[neighborMolIndex];
         absoluteNeighborSitePositionZ[siteIndex] = relativeSitePositionsZ[site] + zptr[neighborMolIndex];
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
     //const auto rotatedPrimeSitePositionX = relativeSitePositionsXPrime[primeSite];
     //const auto rotatedPrimeSitePositionY = relativeSitePositionsYPrime[primeSite];
     //const auto rotatedPrimeSitePositionZ = relativeSitePositionsZPrime[primeSite];

     const auto absolutePrimeSitePositionX = relativeSitePositionsXPrime[primeSite] + xptr[indexPrime];
     const auto absolutePrimeSitePositionY = relativeSitePositionsYPrime[primeSite] + yptr[indexPrime];
     const auto absolutePrimeSitePositionZ = relativeSitePositionsZPrime[primeSite] + zptr[indexPrime];

     // generate parameter data for chosen site
     if constexpr (useMixing) {
       const auto primeSiteType = siteTypesPrime[primeSite];

       size_t siteIndex = 0;
       for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
         const auto neighborMolIndex = neighborList[neighborMol];
         const auto siteTypesOfNeighborMol = _PPLibrary->getSiteTypes(typeptr[neighborMolIndex]);

         for (size_t site = 0; site < relSitePosXptr[neighborMolIndex].size(); ++site) {
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

       const auto displacementX = absolutePrimeSitePositionX - absoluteNeighborSitePositionX[neighborSite];
       const auto displacementY = absolutePrimeSitePositionY - absoluteNeighborSitePositionY[neighborSite];
       const auto displacementZ = absolutePrimeSitePositionZ - absoluteNeighborSitePositionZ[neighborSite];

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

       torqueSumX += relativeSitePositionsYPrime[primeSite] * forceZ - relativeSitePositionsZPrime[primeSite] * forceY;
       torqueSumY += relativeSitePositionsZPrime[primeSite] * forceX - relativeSitePositionsXPrime[primeSite] * forceZ;
       torqueSumZ += relativeSitePositionsXPrime[primeSite] * forceY - relativeSitePositionsYPrime[primeSite] * forceX;

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
     {
       size_t siteIndex = 0;
       for (size_t neighborMol = 0; neighborMol < neighborListSize; ++neighborMol) {
         const auto neighborMolIndex = neighborList[neighborMol];
         //const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(
         //    {q0ptr[neighborMolIndex], q1ptr[neighborMolIndex], q2ptr[neighborMolIndex], q3ptr[neighborMolIndex]},
         //    _PPLibrary->getSitePositions(typeptr[neighborMolIndex]));
         const auto relativeSitePositionsX = relSitePosXptr[neighborMolIndex];
         const auto relativeSitePositionsY = relSitePosYptr[neighborMolIndex];
         const auto relativeSitePositionsZ = relSitePosZptr[neighborMolIndex];

         for (size_t site = 0; site < _PPLibrary->getNumSites(typeptr[neighborMolIndex]); ++site) {
           fxptr[neighborMolIndex] += siteForceX[siteIndex];
           fyptr[neighborMolIndex] += siteForceY[siteIndex];
           fzptr[neighborMolIndex] += siteForceZ[siteIndex];
           txptr[neighborMolIndex] += relativeSitePositionsY[site] * siteForceZ[siteIndex] -
                                      relativeSitePositionsZ[site] * siteForceY[siteIndex];
           typtr[neighborMolIndex] += relativeSitePositionsZ[site] * siteForceX[siteIndex] -
                                      relativeSitePositionsX[site] * siteForceZ[siteIndex];
           tzptr[neighborMolIndex] += relativeSitePositionsX[site] * siteForceY[siteIndex] -
                                      relativeSitePositionsY[site] * siteForceX[siteIndex];
           ++siteIndex;
         }
       }
     }
   }

   if constexpr (calculateGlobals) {
     const auto threadNum = autopas::autopas_get_thread_num();
     // SoAFunctorSingle obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
     // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.
     if (newton3) {
       _aosThreadData[threadNum].potentialEnergySumN3 += potentialEnergySum * 0.5;
       _aosThreadData[threadNum].virialSumN3[0] += virialSumX * 0.5;
       _aosThreadData[threadNum].virialSumN3[1] += virialSumY * 0.5;
       _aosThreadData[threadNum].virialSumN3[2] += virialSumZ * 0.5;
     } else {
       _aosThreadData[threadNum].potentialEnergySumNoN3 += potentialEnergySum;
       _aosThreadData[threadNum].virialSumNoN3[0] += virialSumX;
       _aosThreadData[threadNum].virialSumNoN3[1] += virialSumY;
       _aosThreadData[threadNum].virialSumNoN3[2] += virialSumZ;
     }
   }
 }

 /**
  * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
  */
 class AoSThreadData {
  public:
   AoSThreadData()
       : virialSumNoN3{0., 0., 0.},
         virialSumN3{0., 0., 0.},
         potentialEnergySumNoN3{0.},
         potentialEnergySumN3{0.},
         __remainingTo64{} {}
   void setZero() {
     virialSumNoN3 = {0., 0., 0.};
     virialSumN3 = {0., 0., 0.};
     potentialEnergySumNoN3 = 0.;
     potentialEnergySumN3 = 0.;
   }

   // variables
   std::array<double, 3> virialSumNoN3;
   std::array<double, 3> virialSumN3;
   double potentialEnergySumNoN3;
   double potentialEnergySumN3;

  private:
   // dummy parameter to get the right size (64 bytes)
   double __remainingTo64[(64 - 8 * sizeof(double)) / sizeof(double)];
 };

 static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size (should be multiple of 64)");

 /**
  * Thread buffer for AoS
  */
 std::vector<AoSThreadData> _aosThreadData;
};
}  // namespace mdLib