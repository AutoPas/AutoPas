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

template <class TriwiseKernel, class Particle_T, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
    bool calculateGlobals = false, bool countFLOPs = false>
class TriwiseInterpolantFunctor : public autopas::TriwiseFunctor<Particle_T, TriwiseInterpolantFunctor<TriwiseKernel, Particle_T, useNewton3, calculateGlobals, countFLOPs>>
{
    
public:
    TriwiseInterpolantFunctor() = delete;

private:
    explicit TriwiseInterpolantFunctor(TriwiseKernel f, double cutoff, const std::array<double, 3>& a,
      const std::array<size_t, 3>& nNodes, const std::array<std::vector<double>, 3>& intervalSplits, void * /*dummy*/)
            : autopas::TriwiseFunctor<Particle_T, TriwiseInterpolantFunctor<TriwiseKernel, Particle_T, useNewton3, calculateGlobals, countFLOPs>>(cutoff),
        _cutoffSquared{cutoff * cutoff},
        _numNodes{nNodes},
        _intervalSplits{intervalSplits},
        _a{a},
        _b{std::array<double, 3>{cutoff, cutoff, cutoff}},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadDataGlobals(),
        _postProcessed{false} {

        // TODO: sanitize nodes and splits lists length

        if constexpr (calculateGlobals) {
        _aosThreadDataGlobals.resize(autopas::autopas_get_max_threads());
        }
        if constexpr (countFLOPs) {
        _aosThreadDataFLOPs.resize(autopas::autopas_get_max_threads());
        }

        constructPolynomial();
        evaluateInterpolant(1., 0.2, 0.3);
    }

    // TODO: recursive definition, clenshaw algorithm
    double evaluateChebyshev(double input, int j) {
      return std::cos(j * std::acos(input));
    }

    double evaluateInterpolant(double x, double y, double z) {
      double result = 0.;
      for (int i = 0; i < _numNodes.at(0); ++i) {
        double mappedX = mapToCheb(x, _a.at(0), _b.at(0));
        double val_x = evaluateChebyshev(mappedX, i); 

        for (int j = 0; j < _numNodes.at(1); ++j) {
          double mappedY = mapToCheb(y, _a.at(1), _b.at(1));
          double val_y = evaluateChebyshev(mappedY, j);

          for (int k = 0; k < _numNodes.at(2); ++k) {
            double mappedZ = mapToCheb(z, _a.at(2), _b.at(2));
            double val_z = evaluateChebyshev(mappedZ, k);

            result += _coefficients.at(i).at(j).at(k) * val_x * val_y * val_z;
          }
        }
      }
      return result;
    }

    double mapToInterval(double x, double a, double b) {
      double intermediate = x + 1;
      double newX = a + ((b - a) * intermediate) / 2.;
      return newX;
    }

    double mapToCheb2(double x, double a, double b) {
      double intermediate = x - a;
      double newX = 2. * intermediate / (b-a) - 1.;
      return newX;
    }

    double mapToCheb(double x, double a, double b) {
      return (2 * x - (b+a)) / (b-a);
    }

    double obtainCoefficient(std::array<int, 3>& indices) {
      double c = 0.;

      for (int k0 = 0; k0 < _numNodes.at(0); ++k0) {
        double xNode = std::cos((k0 + 0.5) * PI / _numNodes.at(0));
        double mappedX = mapToInterval(xNode, _a.at(0), _b.at(0));

        for (int k1 = 0; k1 < _numNodes.at(1); ++k1) {
          double yNode = std::cos((k1 + 0.5) * PI / _numNodes.at(1));
          double mappedY = mapToInterval(yNode, _a.at(1), _b.at(2));
          
          for (int k2 = 0; k2 < _numNodes.at(2); ++k2) {
            double zNode = std::cos((k2 + 0.5) * PI / _numNodes.at(2));
            double mappedZ = mapToInterval(zNode, _a.at(2), _b.at(2));

            double value = _kernel.calculateTriplet(mappedX, mappedY, mappedZ);
            double cosinePart = evaluateChebyshev(xNode, indices.at(0)) * evaluateChebyshev(yNode, indices.at(1)) * evaluateChebyshev(zNode, indices.at(2));
            c += value * cosinePart;
          }
        }
      }

      double factor = 1.;

      for (int d = 0; d < 3; ++d) {
        double upper = (indices.at(d) > 0 && indices.at(d) < _numNodes.at(d)) ? 2. : 1.;
        factor *= upper / _numNodes.at(d);
      }

      return c * factor;
    }

    void constructPolynomial() {
      // TODO: port to interval splits
      // Obtain all coefficients
      for (int i = 0; i < _numNodes.at(0); ++i) {
        _coefficients.resize(_numNodes.at(0));
        for (int j = 0; j < _numNodes.at(1); ++j) {
          _coefficients.at(i).resize(_numNodes.at(1));
          for (int k = 0; k < _numNodes.at(2); ++k) {
            _coefficients.at(i).at(j).resize(_numNodes.at(2));
            std::array<int, 3> indices {i, j, k};
            _coefficients.at(i).at(j).at(k) = obtainCoefficient(indices);
          }
        }
      }
    }

public:

    explicit TriwiseInterpolantFunctor(TriwiseKernel f, double cutoff, const std::array<double, 3>& a,
                                      const std::array<size_t, 3>& nNodes, const std::array<std::vector<double>, 3>& intervalSplits)
      : TriwiseInterpolantFunctor(f, cutoff, a, nNodes, intervalSplits, nullptr) {
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

        // Filter unnecessary force computations
        if (i.isDummy() or j.isDummy() or k.isDummy()) {
          return;
        }

        const auto threadnum = autopas::autopas_get_thread_num();

        if constexpr (countFLOPs) {
          ++_aosThreadDataFLOPs[threadnum].numDistCalls;
        }

        auto drIJ = i.getR() - j.getR();
        auto drIK = i.getR() - k.getR();
        auto drJK = j.getR() - k.getR();

        double dIJ2 = autopas::utils::ArrayMath::dot(drIJ, drIJ);
        double dIK2 = autopas::utils::ArrayMath::dot(drIK, drIK);
        double dJK2 = autopas::utils::ArrayMath::dot(drJK, drJK);

        if (dIJ2 > _cutoffSquared or dIK2 > _cutoffSquared or dJK2 > _cutoffSquared) {
          return;
        }

        double dIJ = std::sqrt(dIJ2);
        double dIK = std::sqrt(dIK2);
        double dJK = std::sqrt(dJK2);

        double fac = evaluateInterpolant(dIJ, dIK, dJK);

        // TODO: guard with compile flag
        double real_fac = _kernel.calculateTriplet(dIJ, dIK, dJK);

        double error = real_fac - fac;

        // TODO: add force contributions to the right particles
        // TODO: decide whether to use potential or derivative -> see paper from HSU, ...
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

      constexpr static bool getMixing() { return false; }

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

       /* Interpolation Parameters */
       TriwiseKernel _kernel;
       const std::array<size_t, 3> _numNodes {};
       const std::array<std::vector<double>, 3> _intervalSplits {};
       std::vector<std::vector<std::vector<double>>> _coefficients {};

       const double PI = 2 * std::acos(0.);

       std::array<double, 3> _a {};
       std::array<double, 3> _b {};
    };
}