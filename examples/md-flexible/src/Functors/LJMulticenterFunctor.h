/**
 * @file LJMulticenterFunctor.h
 * @date 21/02/2022
 * @author S. Newcome
*/

#pragma once

#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopas/pairwiseFunctors/Functor.h"

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
          bool calculatedGlobals = false, bool relevantForTuning = true>
class LJMulticenterFunctor
    : public autopas::Functor<Particle,
    LJMulticenterFunctor<Particle, applyShift, useMixing, useNewton3, calculatedGlobals, relevantForTuning>> {
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
    * What is this? not constant as may be reset through PPL.
    */
   double _epsilon24;

   /**
    * sigma^2. Not constant as may be reset through PPL.
    */
   double _sigmaSquared;

   /**
    * What is this? Not constant as may be reset through PPL.
    */
   double _shift6 = 0;

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
          : autopas::Functor <Particle, LJMulticenterFunctor<Particle, applyShift, useMixing, useNewton3, calculatedGlobals, relevantForTuning>>(
              cutoff
              ),
          _cutoffSquared{cutoff*cutoff},
          _potentialEnergySum{0.},
          _virialSum{0.,0.,0.},
          _aosThreadData(),
          _postProcessed{false}, {
       if constexpr (calculateGlobals) {
         _aosThreadData.resize(autopas_get_max_threads());
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
        * @param p_i Particle i
        * @param p_j Particle j
        * @param newton3 Flag for if newton3 is used.
        */
        void AoSFunctor(Particle &p_i, Particle &p_j, bool newton3) final {
          if (p_i.isDummy() or p_j.isDumm()) {
            return;
          }
          auto sigmasquare = _sigma
        }


};