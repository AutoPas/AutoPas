/**
 * @file Triwise InterpolantFunctor
 * 
 * @date 27 May 2025
 * @author Luis Gall
 */

#pragma once

#include "autopas/baseFunctors/TriwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace mdLib {

template <class TriwiseKernel, class Particle_T, bool useMixing = false, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
    bool calculateGlobals = false, bool countFLOPs = false>
class TriwiseInterpolantFunctor : public autopas::TriwiseFunctor<Particle_T, TriwiseInterpolantFunctor<TriwiseKernel, Particle_T, useMixing, useNewton3, calculateGlobals, countFLOPs>>
{
    
public:
    TriwiseInterpolantFunctor() = delete;

private:

    explicit TriwiseInterpolantFunctor(double cutoff, void * /*dummy*/) : autopas::TriwiseFunctor<Particle_T, TriwiseInterpolantFunctor<TriwiseKernel, Particle_T, useMixing, useNewton3, calculateGlobals, countFLOPs>>(),
        _cutoffSquared{cutoff * cutoff},
            _potentialEnergySum{0.},
            _virialSum{0., 0., 0.},
            _aosThreadDataGlobals(),
            _postProcessed{false} {
        if constexpr (calculateGlobals) {
        _aosThreadDataGlobals.resize(autopas::autopas_get_max_threads());
        }
        if constexpr (countFLOPs) {
        _aosThreadDataFLOPs.resize(autopas::autopas_get_max_threads());
        }
    }

public:

    explicit TriwiseInterpolantFunctor(double cutoff) : TriwiseInterpolantFunctor(cutoff, nullptr) {

    }

    std::string getName() final { return "TriwiseInterpolantFunctor"; }

    bool isRelevantForTuning() final { return true; }

    bool allowsNewton3() final {
        return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
    }

    bool allowsNonNewton3() final {
        return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
    }

    void AoSFunctor(Particle_T &i, Particle_T &j, Particle_T &k, bool newton3) final {
        using namespace autopas::utils::ArrayMath::literals;

        // TODO: implement
    }

    constexpr static auto getNeededAttr() {
        return std::array<typename Particle_T::AttributeNames, 9>{Particle_T::AttributeNames::id,
                                                                  Particle_T::AttributeNames::posX,
                                                                  Particle_T::AttributeNames::posY,
                                                                  Particle_T::AttributeNames::posZ,
                                                                  Particle_T::AttributeNames::forceX,
                                                                  Particle_T::AttributeNames::forceY,
                                                                  Particle_T::AttributeNames::forceZ,
                                                                  Particle_T::AttributeNames::typeId,
                                                                  Particle_T::AttributeNames::ownershipState};
      }

      constexpr static auto getNeededAttr(std::false_type) {
        return std::array<typename Particle_T::AttributeNames, 6>{
            Particle_T::AttributeNames::id,     Particle_T::AttributeNames::posX,
            Particle_T::AttributeNames::posY,   Particle_T::AttributeNames::posZ,
            Particle_T::AttributeNames::typeId, Particle_T::AttributeNames::ownershipState};
      }

      constexpr static auto getComputedAttr() {
        return std::array<typename Particle_T::AttributeNames, 3>{
            Particle_T::AttributeNames::forceX, Particle_T::AttributeNames::forceY, Particle_T::AttributeNames::forceZ};
      }

      constexpr static bool getMixing() { return useMixing; }

      void initTraversal() final {
        _potentialEnergySum = 0.;
        _virialSum = {0., 0., 0.};
        _postProcessed = false;
        for (size_t i = 0; i < _aosThreadDataGlobals.size(); ++i) {
          _aosThreadDataGlobals[i].setZero();
        }
      }

      void endTraversal(bool newton3) final {
        using namespace autopas::utils::ArrayMath::literals;
    
        if (_postProcessed) {
          throw autopas::utils::ExceptionHandler::AutoPasException(
              "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
        }
        if (calculateGlobals) {
          // Accumulate potential energy and virial values.
          for (const auto &data : _aosThreadDataGlobals) {
            _potentialEnergySum += data.potentialEnergySum;
            _virialSum += data.virialSum;
          }
    
          // For each interaction, we added the full contribution for all three particles. Divide by 3 here, so that each
          // contribution is only counted once per triplet.
          _potentialEnergySum /= 3.;
    
          // Additionally, we have always calculated 3*potentialEnergy, so we divide by 3 again.
          _potentialEnergySum /= 3.;
    
          _postProcessed = true;
    
          // TODO: print interpolation error
          
          AutoPasLog(TRACE, "Final potential energy {}", _potentialEnergySum);
          AutoPasLog(TRACE, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
        }
      }

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

      [[nodiscard]] size_t getNumFLOPs() const override {
        // TODO: implement
        return std::numeric_limits<size_t>::max();
      }

      [[nodiscard]] double getHitRate() const override {
        // TODO: implement
        return std::numeric_limits<double>::quiet_NaN();
      }

      class AoSThreadDataGlobals {
        public:
         AoSThreadDataGlobals() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, __remainingTo64{} {}
         void setZero() {
           virialSum = {0., 0., 0.};
           potentialEnergySum = 0.;
         }
     
         // variables
         std::array<double, 3> virialSum;
         double potentialEnergySum;
     
        private:
         // dummy parameter to get the right size (64 bytes)
         double __remainingTo64[(64 - 4 * sizeof(double)) / sizeof(double)];
       };

       class AoSThreadDataFLOPs {
        public:
         AoSThreadDataFLOPs() : __remainingTo64{} {}
     
         /**
          * Set all counters to zero.
          */
         void setZero() {
           numKernelCallsNoN3 = 0;
           numKernelCallsN3 = 0;
           numDistCalls = 0;
           numGlobalCalcsN3 = 0;
           numGlobalCalcsNoN3 = 0;
         }
     
         /**
          * Number of calls to Lennard-Jones Kernel with newton3 disabled.
          * Used for calculating number of FLOPs and hit rate.
          */
         size_t numKernelCallsNoN3 = 0;
     
         /**
          * Number of calls to Lennard-Jones Kernel with newton3 enabled.
          * Used for calculating number of FLOPs and hit rate.
          */
         size_t numKernelCallsN3 = 0;
     
         /**
          * Number of distance calculations.
          * Used for calculating number of FLOPs and hit rate.
          */
         size_t numDistCalls = 0;
     
         /**
          * Counter for the number of times the globals have been calculated with newton3 enabled.
          */
         size_t numGlobalCalcsN3 = 0;
     
         /**
          * Counter for the number of times the globals have been calculated without newton3 enabled.
          */
         size_t numGlobalCalcsNoN3 = 0;
     
        private:
         /**
          * dummy parameter to get the right size (64 bytes)
          */
         double __remainingTo64[(64 - 5 * sizeof(size_t)) / sizeof(size_t)];
       };
     
       // make sure of the size of AoSThreadDataGlobals
       static_assert(sizeof(AoSThreadDataGlobals) % 64 == 0, "AoSThreadDataGlobals has wrong size");
       static_assert(sizeof(AoSThreadDataFLOPs) % 64 == 0, "AoSThreadDataFLOPs has wrong size");

       const double _cutoffSquared;
       double _potentialEnergySum;
       std::array<double, 3> _virialSum;

       std::vector<AoSThreadDataGlobals> _aosThreadDataGlobals;
       std::vector<AoSThreadDataFLOPs> _aosThreadDataFLOPs{};

       bool _postProcessed;
    };
}