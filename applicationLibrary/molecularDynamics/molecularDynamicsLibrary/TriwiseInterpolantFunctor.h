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

#include <fftw3.h>
#include <unsupported/Eigen/CXX11/Tensor>

namespace mdLib {

template <class TriwiseKernel, class Particle_T, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
    bool calculateGlobals = false, bool countFLOPs = false>
class TriwiseInterpolantFunctor : public autopas::TriwiseFunctor<Particle_T, TriwiseInterpolantFunctor<TriwiseKernel, Particle_T, useNewton3, calculateGlobals, countFLOPs>>
{
    
public:
    TriwiseInterpolantFunctor() = delete;

private:
    explicit TriwiseInterpolantFunctor(TriwiseKernel f, double cutoff, const std::array<double, 3>& a,
      const std::array<std::vector<size_t>, 3>& nNodes, const std::array<std::vector<double>, 3>& intervalSplits, void * /*dummy*/)
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

        for (size_t dim = 0; dim < 3; ++dim) {
          if (_numNodes.at(dim).size() != intervalSplits.at(dim).size() + 1) {
            throw autopas::utils::ExceptionHandler::AutoPasException("Size of node list unequal to size of intervals + 1 in dimension " + dim);
          }
        }

        // initialize vectors
        // _coefficients = Eigen::Tensor<double, 3, Eigen::RowMajor>(_numNodes.at(0), _numNodes.at(1), _numNodes.at(2));
        _coefficients = Eigen::Tensor<std::array<Eigen::Tensor<double, 3, Eigen::RowMajor>, 3>, 3, Eigen::RowMajor>(_numNodes.at(0).size(), _numNodes.at(1).size(), _numNodes.at(2).size());
        constructPolynomialFFTW();
    }

    // TODO: recursive definition, clenshaw algorithm
    double evaluateChebyshev(double input, int j) {
      return std::cos(j * std::acos(input));
    }

    double evaluateChebyshevFast1D(double input, int n, const std::vector<double>& coeffs) {
      // see PairwiseInterpolantFunctor.h for reference
      double b_n = 0.;
      double b_n1 = 0.;
      double b_n2 = 0.;
      for (int k = n - 1; k > 0; --k) {
        b_n = coeffs.at(k) + 2. * input * b_n1 - b_n2;
        b_n2 = b_n1;
        b_n1 = b_n;
      }

      b_n = 2. * coeffs.at(0) + 2. * input * b_n1 - b_n2;
      return (b_n - b_n2) / 2.;
    }

    double evaluateChebyshevFast3D(double x, double y, double z, int intervalX, int intervalY, int intervalZ, int dim) {
      // This code is oriented on NumPy's implementation of chebeval3d

      size_t nX = _numNodes.at(0).at(intervalX);
      size_t nY = _numNodes.at(1).at(intervalY);
      size_t nZ = _numNodes.at(2).at(intervalZ);

      auto coeffs = _coefficients(intervalX, intervalY, intervalZ).at(dim);


      /* Flatten x dimension */
      auto jkCoeffs = Eigen::MatrixXd(nY, nZ);

      for (size_t j = 0; j < nY; ++j) {
        for (size_t k = 0; k < nZ; ++k) {

          std::vector<double> xCoeffs {};
          xCoeffs.reserve(nX);
          /* This loop can potentially be avoided */
          for (size_t i = 0; i < nX; ++i) {
            xCoeffs.push_back(coeffs(i,j,k));
          }
          jkCoeffs(j,k) = evaluateChebyshevFast1D(x, nX, xCoeffs);
        }
      }

      std::vector<double> zCoeffs {};
      zCoeffs.reserve(nZ);
      for (size_t k = 0; k < nZ; ++k) {
        std::vector<double> yCoeffs {};
        yCoeffs.reserve(nY);

        /* This loop can potentially be avoided */
        for (size_t j = 0; j < nY; ++j) {
          yCoeffs.push_back(jkCoeffs(j, k));
        }
        zCoeffs.push_back(evaluateChebyshevFast1D(y, nY, yCoeffs));
      }

      return evaluateChebyshevFast1D(z, nZ, zCoeffs);
    }

    /*
    double evaluateInterpolant(double x, double y, double z) {
      int intervalX = 0;
      int intervalY = 0;
      int intervalZ = 0;

      double result = 0.;
      for (int i = 0; i < _numNodes.at(0).at(intervalX); ++i) {
        double mappedX = mapToCheb(x, _a.at(0), _b.at(0));
        double val_x = evaluateChebyshev(mappedX, i); 

        for (int j = 0; j < _numNodes.at(1).at(intervalY); ++j) {
          double mappedY = mapToCheb(y, _a.at(1), _b.at(1));
          double val_y = evaluateChebyshev(mappedY, j);

          for (int k = 0; k < _numNodes.at(2).at(intervalZ); ++k) {
            double mappedZ = mapToCheb(z, _a.at(2), _b.at(2));
            double val_z = evaluateChebyshev(mappedZ, k);

            result += _coefficients(intervalX, intervalY, intervalZ)(i,j,k) * val_x * val_y * val_z;
          }
        }
      }
      return result;
    }
    */

    double mapToInterval(double x, double a, double b) {
      double intermediate = x + 1;
      double newX = a + ((b - a) * intermediate) / 2.;
      return newX;
    }

    double mapToCheb(double x, double a, double b) {
      return (2. * (x-a)) / (b-a) - 1.;
    }

    void constructPolynomialFFTW() {

      double aX = _a.at(0);
      for (size_t intervalX = 0; intervalX < _numNodes.at(0).size(); ++intervalX) {
        double bX = _intervalSplits.at(0).size() > intervalX ? _intervalSplits.at(0).at(intervalX) : _b.at(0);
        
        double aY = _a.at(1);
        for (size_t intervalY = 0; intervalY < _numNodes.at(1).size(); ++intervalY) {
          double bY = _intervalSplits.at(1).size() > intervalY ? _intervalSplits.at(1).at(intervalY) : _b.at(1);

          double aZ = _a.at(2);
          for (size_t intervalZ = 0; intervalZ < _numNodes.at(2).size(); ++intervalZ) {
            double bZ = _intervalSplits.at(2).size() > intervalZ ? _intervalSplits.at(2).at(intervalZ) : _b.at(2);

            int nX = _numNodes.at(0).at(intervalX);
            int nY = _numNodes.at(1).at(intervalY);
            int nZ = _numNodes.at(2).at(intervalZ);

            std::array<Eigen::Tensor<double, 3, Eigen::RowMajor>, 3> & intervalCoeffs = _coefficients(intervalX, intervalY, intervalZ);
            intervalCoeffs = std::array<Eigen::Tensor<double, 3, Eigen::RowMajor>, 3>{
              Eigen::Tensor<double, 3, Eigen::RowMajor>(nX, nY, nZ),
              Eigen::Tensor<double, 3, Eigen::RowMajor>(nX, nY, nZ),
              Eigen::Tensor<double, 3, Eigen::RowMajor>(nX, nY, nZ)
            };

            /* Prepare Tensors */
            for (int i = 0; i < nX; ++i) {
              double nodeX = std::cos((2*i+1)*PI/(2.*nX));

              for (int j = 0; j < nY; ++j) {

                double nodeY = std::cos((2*j+1)*PI/(2.*nY));
                for (int k = 0; k < nZ; ++k) {  
                  double nodeZ = std::cos((2*k+1)*PI/(2.*nZ));

                  std::array<double, 3> forces = _kernel.calculateTripletDerivative(
                    mapToInterval(nodeX, aX, bX),
                    mapToInterval(nodeY, aY, bY),
                    mapToInterval(nodeZ, aZ, bZ)
                  );

                  intervalCoeffs.at(0)(i,j,k) = forces.at(0);
                  intervalCoeffs.at(1)(i,j,k) = forces.at(1);
                  intervalCoeffs.at(2)(i,j,k) = forces.at(2);
                }
              }
            }

            for (int dim = 0; dim < 3; ++dim) {
              /* Let FFTW compute the (unscaled) coefficients */
              fftw_plan plan = fftw_plan_r2r_3d(nX, nY, nZ, intervalCoeffs.at(dim).data(), intervalCoeffs.at(dim).data(),
                FFTW_REDFT10, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
              fftw_execute(plan);
              fftw_destroy_plan(plan);
            }

            /* Scale the coefficients according to the Chebyshev interpolation */
            for (int i = 0; i < nX; ++i) {
              for (int j = 0; j < nY; ++j) {
                for (int k = 0; k < nZ; ++k) {
                  for (int dim = 0; dim < 3; ++dim) {
                    if (i == 0) {
                      intervalCoeffs.at(dim)(i,j,k) *= 0.5;
                    }
                    if (j == 0) {
                      intervalCoeffs.at(dim)(i,j,k) *= 0.5;
                    }
                    if (k == 0) {
                      intervalCoeffs.at(dim)(i,j,k) *= 0.5;
                    }
                    intervalCoeffs.at(dim)(i,j,k) *= (1./nX * 1./nY * 1./ nZ);
                  }
                }
              }
            }

            aZ = bZ;
          } // intervalZ
          aY = bY;
        } // intervalY
        aX = bX;
      } // intervalX
    }

public:

    explicit TriwiseInterpolantFunctor(TriwiseKernel f, double cutoff, const std::array<double, 3>& a,
                                      const std::array<std::vector<size_t>, 3>& nNodes, const std::array<std::vector<double>, 3>& intervalSplits)
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

        auto drIJ = j.getR() - i.getR();
        auto drJK = k.getR() - j.getR();
        auto drKI = i.getR() - k.getR();

        const double dIJ2 = autopas::utils::ArrayMath::dot(drIJ, drIJ);
        const double dJK2 = autopas::utils::ArrayMath::dot(drJK, drJK);
        const double dKI2 = autopas::utils::ArrayMath::dot(drKI, drKI);

        if (dIJ2 > _cutoffSquared or dKI2 > _cutoffSquared or dJK2 > _cutoffSquared) {
          return;
        }

        const double dIJ = std::sqrt(dIJ2);
        const double dJK = std::sqrt(dJK2);
        const double dKI = std::sqrt(dKI2);

        int intervalX = 0;
        int intervalY = 0;
        int intervalZ = 0;

        double aX = _a.at(0);
        double aY = _a.at(1);
        double aZ = _a.at(2);

        double bX = _b.at(0);
        double bY = _b.at(1);
        double bZ = _b.at(2);

        for (; intervalX < _intervalSplits.at(0).size(); ++intervalX) {
          if (dIJ < _intervalSplits.at(0).at(intervalX)) {
            bX = _intervalSplits.at(0).at(intervalX);
            break;
          }
          else {
            aX = _intervalSplits.at(0).at(intervalX);
          }
        }

        for (; intervalY < _intervalSplits.at(1).size(); ++intervalY) {
          if (dJK < _intervalSplits.at(1).at(intervalY)) {
            bY = _intervalSplits.at(1).at(intervalY);
            break;
          }
          else {
            aY = _intervalSplits.at(1).at(intervalY);
          }
        }

        for (; intervalZ < _intervalSplits.at(2).size(); ++intervalZ) {
          if (dKI < _intervalSplits.at(2).at(intervalZ)) {
            bZ = _intervalSplits.at(2).at(intervalZ);
            break;
          }
          else {
            aZ = _intervalSplits.at(2).at(intervalZ);
          }
        }

        const double x = mapToCheb(dIJ, aX, bX);
        const double y = mapToCheb(dJK, aY, bY);
        const double z = mapToCheb(dKI, aZ, bZ);

        const double fac_ij = evaluateChebyshevFast3D(x, y, z, intervalX, intervalY, intervalZ, 0);
        const double fac_ki = evaluateChebyshevFast3D(x, y, z, intervalX, intervalY, intervalZ, 2);
        
        const auto forceI = ( drIJ * fac_ij + drKI * fac_ki ) * -1.;
        i.addF(forceI);

        double fac_jk = 0.;
        if (newton3) {
          fac_jk = evaluateChebyshevFast3D(x, y, z, intervalX, intervalY, intervalZ, 1);

          const auto forceJ = drIJ * -fac_ij + drJK * fac_jk;
          const auto forceK = (forceI + forceJ) * (-1.);

          j.addF(forceJ);
          k.addF(forceK);
        }

        // TODO: guard with compile flag
        double errorX = 0.;
        double errorY = 0.;
        double errorZ = 0.;
#if defined(MD_FLEXIBLE_BENCHMARK_INTERPOLANT_ACCURACY)
        const std::array<double, 3> real_fac = _kernel.calculateTripletDerivative(dIJ, dJK, dKI);

        errorX = real_fac.at(0) - fac_ij;
        errorY = real_fac.at(1) - fac_jk;
        errorZ = real_fac.at(2) - fac_ki;
#endif

        if constexpr (calculateGlobals) {
          
          _aosThreadDataGlobals[threadnum].absErrorSum += {errorX, errorY, errorZ};

        }

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
        _absErrorSum = {0., 0., 0.};
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
            _absErrorSum += data.absErrorSum;
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
          AutoPasLog(DEBUG, "Final absolute error   {} {} {}", _absErrorSum.at(0), _absErrorSum.at(1), _absErrorSum.at(2));
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
         AoSThreadDataGlobals() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, absErrorSum{0.,0.,0.}, __remainingTo64{} {}
         void setZero() {
           virialSum = {0., 0., 0.};
           absErrorSum = {0.,0.,0.};
           potentialEnergySum = 0.;
         }
     
         // variables
         std::array<double, 3> absErrorSum;
         std::array<double, 3> virialSum;
         double potentialEnergySum;
     
        private:
         // dummy parameter to get the right size (64 bytes)
         double __remainingTo64[(64 - 7 * sizeof(double)) / sizeof(double)];
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
       std::array<double, 3> _absErrorSum;

       std::vector<AoSThreadDataGlobals> _aosThreadDataGlobals;
       std::vector<AoSThreadDataFLOPs> _aosThreadDataFLOPs{};

       bool _postProcessed;

       /* Interpolation Parameters */
       TriwiseKernel _kernel;
       const std::array<std::vector<size_t>, 3> _numNodes {};
       const std::array<std::vector<double>, 3> _intervalSplits {};
       Eigen::Tensor<std::array<Eigen::Tensor<double, 3, Eigen::RowMajor>, 3>, 3, Eigen::RowMajor> _coefficients {};

       const double PI = 2 * std::acos(0.);

       std::array<double, 3> _a {};
       std::array<double, 3> _b {};
    };
}