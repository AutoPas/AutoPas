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

#if defined(MD_FLEXIBLE_INTERPOLANT_VECTORIZATION)
#include "immintrin.h"
#endif

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
        _aosThreadDataFLOPs(),
        _postProcessed{false} {

        if constexpr (calculateGlobals) {
        _aosThreadDataGlobals.resize(autopas::autopas_get_max_threads());
        }
        if constexpr (countFLOPs) {
          _aosThreadDataFLOPs.resize(autopas::autopas_get_max_threads());
          _distanceStatistics.resize(autopas::autopas_get_max_threads());
        }

        for (size_t dim = 0; dim < 3; ++dim) {
          if (_numNodes.at(dim).size() != intervalSplits.at(dim).size() + 1 && _numNodes.at(dim).size() != 0) {
            throw autopas::utils::ExceptionHandler::AutoPasException("Size of node list unequal to size of intervals + 1 in dimension " + dim);
          }
        }

        // initialize vectors
        _coeffs = Eigen::Tensor<std::vector<std::vector<std::vector<double>>>, 3, Eigen::RowMajor>
          (_numNodes.at(0).size(), _numNodes.at(1).size(), _numNodes.at(2).size());

#if defined(MD_FLEXIBLE_INTERPOLANT_VECTORIZATION)
        _coeffsVec = Eigen::Tensor<Eigen::Tensor<double, 3, Eigen::RowMajor>, 3, Eigen::RowMajor>
          (_numNodes.at(0).size(), _numNodes.at(1).size(), _numNodes.at(2).size());
#endif

        constructPolyFFTW();
    }

#if defined(MD_FLEXIBLE_INTERPOLANT_VECTORIZATION)
    double evalChebFast1DVector(double input, int n, Eigen::VectorXd& coeffs) {
      double b_n = 0.;
      double b_n1 = 0.;
      double b_n2 = 0.;

      for (int k = n - 1; k > 0; --k) {

        b_n = coeffs(k) + 2. * input * b_n1 - b_n2;
        b_n2 = b_n1;
        b_n1 = b_n;
      }

      b_n = 2. * coeffs(0) + 2. * input * b_n1 - b_n2;
      return (b_n - b_n2) / 2.;
    }

    __m512d evalChebFast1DVector(double input, int n, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& coeffs, size_t i, size_t size) {
      __m512d inputVec = _mm512_set1_pd(input);
      __m512d c;

      __m512d b_n = _mm512_setzero_pd();
      __m512d b_n1 = _mm512_setzero_pd();
      __m512d b_n2 = _mm512_setzero_pd();

      for (int j = n - 1; j > 0; --j) {

        /* Full Load */
        if (size >= 8) {
          c = _mm512_load_pd(&coeffs(j, i));
        }
        /* Masked Load */
        else {
          __mmask8 mask = (1 << size) - 1;
          c = _mm512_maskz_load_pd(mask, &coeffs(j, i));
        }

        b_n = c + 2. * input * b_n1 - b_n2;
        b_n2 = b_n1;
        b_n1 = b_n;
      }
      if (size >= 8) {
        c = _mm512_load_pd(&coeffs(0, i));
      }
      /* Masked Load */
      else {
        __mmask8 mask = (1 << size) - 1;
        c = _mm512_maskz_load_pd(mask, &coeffs(0, i));
      }

      b_n = 2. * c + 2. * input * b_n1 - b_n2;
      return (b_n - b_n2) / 2.;
    }

    __m512d evalChebFast1DVector(double input, int n, Eigen::Tensor<double, 3, Eigen::RowMajor>& coeffs, int i, int j, int size) {

      __m512d inputV =  _mm512_set1_pd(input);
      __m512d c;

      __m512d b_n = _mm512_setzero_pd();
      __m512d b_n1 = _mm512_setzero_pd();
      __m512d b_n2 = _mm512_setzero_pd();

      for (int k = n - 1; k > 0; --k) {

        /* Full Load */
        if (size >= 8) {
          c = _mm512_load_pd(&coeffs(i, k, j));
        }
        /* Masked Load */
        else {
          __mmask8 mask = (1 << size) - 1;
          c = _mm512_maskz_load_pd(mask, &coeffs(i, k, j));
        }

        b_n = c + 2. * input * b_n1 - b_n2;
        b_n2 = b_n1;
        b_n1 = b_n;
      }
      if (size >= 8) {
        c = _mm512_load_pd(&coeffs(i, 0, j));
      }
      /* Masked Load */
      else {
        __mmask8 mask = (1 << size) - 1;
        c = _mm512_maskz_load_pd(mask, &coeffs(i, 0, j));
      }

      b_n = 2. * c + 2. * input * b_n1 - b_n2;
      return (b_n - b_n2) / 2.;
    }

    double evalChebFast3DVector(double x, double y, double z, int intervalX, int intervalY, int intervalZ) {
      size_t nX = _numNodes.at(0).at(intervalX);
      size_t nY = _numNodes.at(1).at(intervalY);
      size_t nZ = _numNodes.at(2).at(intervalZ);

      Eigen::Tensor<double, 3, Eigen::RowMajor> coeffs = _coeffsVec(intervalX, intervalY, intervalZ);

      /* Flatten z dimension */
      size_t newNx = nX + (8 - nX % 8)%8;
      size_t newNy = nY + (8 - nY % 8)%8;
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> xyCoeffs (newNy, newNx);
      for (size_t j = 0; j < nY; j+=8) {
        for (size_t i = 0; i < nX; ++i) {
        
          size_t size = nY - j;
          __m512d values = evalChebFast1DVector(z, nZ, coeffs, i, j, size);

          alignas(64) double helper [8];
          _mm512_store_pd(helper, values);

          for (int jHelp = 0; jHelp < size; ++jHelp) {
            xyCoeffs(jHelp + j, i) = helper[jHelp];
          }
        }
      }

      /* Flatten y dimension */
      Eigen::VectorXd xCoeffs (nX);
      for (size_t i = 0; i < nX; i+=8) {
        int size = nX - i;
        __m512d values = evalChebFast1DVector(y, nY, xyCoeffs, i, size);
        
        alignas(64) double helper [8];
        _mm512_store_pd(helper, values);

        for (int iHelp = 0; iHelp < size; ++iHelp) {
          xCoeffs(iHelp + i) = helper[iHelp];
        }
      }

      /* Flatten x dimension */
      return evalChebFast1DVector(x, nX, xCoeffs);
    }
#endif

    double evalChebFast1D(double input, int n, const std::vector<double>& coeffs) {
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


    double evalChebFast3D(double x, double y, double z, int intervalX, int intervalY, int intervalZ) {
      // This code is oriented on NumPy's implementation of chebeval3d
      size_t nX = _numNodes.at(0).at(intervalX);
      size_t nY = _numNodes.at(1).at(intervalY);
      size_t nZ = _numNodes.at(2).at(intervalZ);

      auto coeffs = _coeffs(intervalX, intervalY, intervalZ);

      /* Flatten z dimension */
      auto xyCoeffs = std::vector<std::vector<double>> {};
      xyCoeffs.reserve(nX);
      for (size_t i = 0; i < nX; ++i) {
        std::vector<double> yCoeffs {};
        yCoeffs.reserve(nY);
        for (size_t j = 0; j < nY; ++j) {
          
          double value = evalChebFast1D(z, nZ, coeffs.at(i).at(j));
          yCoeffs.push_back(value);
        }
        xyCoeffs.push_back(yCoeffs);
      }

      /* Flatten y dimension */
      std::vector<double> xCoeffs {};
      xCoeffs.reserve(nX);
      for (size_t i = 0; i < nX; ++i) {
        xCoeffs.push_back(evalChebFast1D(y, nY, xyCoeffs.at(i)));
      }

      /* Flatten x dimension */
      return evalChebFast1D(x, nX, xCoeffs);
    }

    double mapToInterval(double x, double a, double b) {
      double intermediate = x + 1;
      double newX = a + ((b - a) * intermediate) / 2.;
      return newX;
    }

    double mapToCheb(double x, double a, double b) {
      return (2. * (x-a)) / (b-a) - 1.;
    }

    void constructPolyFFTW() {

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

            double forcesXYZ [nX][nY][nZ];

            /* Prepare Tensors, initialize vectors */
            for (int i = 0; i < nX; ++i) {
              double nodeX = std::cos((2*i+1)*PI/(2.*nX));

              for (int j = 0; j < nY; ++j) {

                double nodeY = std::cos((2*j+1)*PI/(2.*nY));

                for (int k = 0; k < nZ; ++k) {  
                  double nodeZ = std::cos((2*k+1)*PI/(2.*nZ));

                  std::array<double, 3> forces = _kernel.calculateTripletDerivative(
                    mapToInterval(nodeX, aX*aX, bX*bX),
                    mapToInterval(nodeY, aY*aY, bY*bY),
                    mapToInterval(nodeZ, aZ*aZ, bZ*bZ)
                  );

                    forcesXYZ[i][j][k] = forces.at(0);
                }
              }
            }

            /* Let FFTW compute the (unscaled) coefficients */
            fftw_plan plan_test = fftw_plan_r2r_3d(nX, nY, nZ, &forcesXYZ[0][0][0], &forcesXYZ[0][0][0],
                FFTW_REDFT10, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
            fftw_execute(plan_test);
            fftw_destroy_plan(plan_test);
            

            std::vector<std::vector<std::vector<double>>> &coeffs = _coeffs(intervalX, intervalY, intervalZ);
            coeffs = std::vector<std::vector<std::vector<double>>>(nX);

#if defined(MD_FLEXIBLE_INTERPOLANT_VECTORIZATION)
            // Extend Coefficients Tensor such that it is aligned
            size_t newNx = nX + (8 - nX % 8)%8;
            size_t newNy = nY + (8 - nY % 8)%8;
            size_t newNz = nZ + (8 - nZ % 8)%8;

            Eigen::Tensor<double, 3, Eigen::RowMajor> &coeffsVec = _coeffsVec(intervalX, intervalY, intervalZ);
            coeffsVec = Eigen::Tensor<double, 3, Eigen::RowMajor>(newNx, newNz, newNy);
#endif

            /* Scale the coefficients according to the Chebyshev interpolation */
              for (int i = 0; i < nX; ++i) {
                std::vector<std::vector<double>> coeffsYZ {};
                coeffsYZ.reserve(nY);
                for (int j = 0; j < nY; ++j) {
                  std::vector<double> coeffsZ {};
                  coeffsZ.reserve(nZ);
                  for (int k = 0; k < nZ; ++k) {                    
                    if (i == 0) {
                      forcesXYZ[i][j][k] *= 0.5;
                    }
                    if (j == 0) {
                      forcesXYZ[i][j][k] *= 0.5;
                    }
                    if (k == 0) {
                      forcesXYZ[i][j][k] *= 0.5;
                    }
                    forcesXYZ[i][j][k] *= (1./nX * 1./nY * 1./ nZ);
                    coeffsZ.push_back(forcesXYZ[i][j][k]);
#if defined(MD_FLEXIBLE_INTERPOLANT_VECTORIZATION)
                    coeffsVec(i, k, j) = forcesXYZ[i][j][k];
#endif
                  }
                  coeffsYZ.emplace_back(coeffsZ);
                }
                coeffs.at(i) = coeffsYZ;
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

        int intervalX = 0;
        int intervalY = 0;
        int intervalZ = 0;

        double aX = _a.at(0);
        double aY = _a.at(1);
        double aZ = _a.at(2);

        double bX = _b.at(0);
        double bY = _b.at(1);
        double bZ = _b.at(2);

#if defined(MD_FLEXIBLE_BENCHMARK_INTERPOLANT_ACCURACY)
        if constexpr (countFLOPs) {
          _distanceStatistics[threadnum].push_back(dIJ2);
          _distanceStatistics[threadnum].push_back(dKI2);
          _distanceStatistics[threadnum].push_back(dJK2);
        }
#endif

        for (; intervalX < _intervalSplits.at(0).size(); ++intervalX) {
          const double split = _intervalSplits.at(0).at(intervalX);
          if (dIJ2 < split*split) {
            bX = split;
            break;
          }
          else {
            aX = split;
          }
        }

        for (; intervalY < _intervalSplits.at(1).size(); ++intervalY) {
          const double split = _intervalSplits.at(1).at(intervalY);
          if (dJK2 < split*split) {
            bY = split;
            break;
          }
          else {
            aY = split;
          }
        }

        for (; intervalZ < _intervalSplits.at(2).size(); ++intervalZ) {
          const double split = _intervalSplits.at(2).at(intervalZ);
          if (dKI2 < split*split) {
            bZ = split;
            break;
          }
          else {
            aZ = split;
          }
        }

        const double x = mapToCheb(dIJ2, aX*aX, bX*bX);
        const double y = mapToCheb(dJK2, aY*aY, bY*bY);
        const double z = mapToCheb(dKI2, aZ*aZ, bZ*bZ);

#if defined (MD_FLEXIBLE_INTERPOLANT_VECTORIZATION)
        const double fac_ij = evalChebFast3DVector(x, y, z, intervalX, intervalY, intervalZ);
        const double fac_ki = evalChebFast3DVector(z, x, y, intervalZ, intervalX, intervalY);
#else
        const double fac_ij = evalChebFast3D(x, y, z, intervalX, intervalY, intervalZ);
        const double fac_ki = evalChebFast3D(z, x, y, intervalZ, intervalX, intervalY);
#endif
        const auto forceI = ( drIJ * fac_ij + drKI * fac_ki * -1. ) * -1.;
        i.addF(forceI);

        double fac_jk = 0.;
        double fac_jk_vec = 0.;
        if (newton3) {
#if defined (MD_FLEXIBLE_INTERPOLANT_VECTORIZATION)
          fac_jk = evalChebFast3DVector(x, y, z, intervalX, intervalY, intervalZ);
#else
          fac_jk = evalChebFast3D(y, z, x, intervalY, intervalZ, intervalX);
#endif

          const auto forceJ = drIJ * -fac_ij + drJK * fac_jk;
          const auto forceK = (forceI + forceJ) * (-1.);

          j.addF(forceJ);
          k.addF(forceK);
        }

        if constexpr (countFLOPs) {
          if (newton3) {
            ++_aosThreadDataFLOPs[threadnum].numKernelCallsN3;
          }
          else {
            ++_aosThreadDataFLOPs[threadnum].numKernelCallsNoN3;
          }
        }

        double errorX = 0.;
        double errorY = 0.;
        double errorZ = 0.;

        double relErrorX = 0.;
        double relErrorY = 0.;
        double relErrorZ = 0.;
#if defined(MD_FLEXIBLE_BENCHMARK_INTERPOLANT_ACCURACY)
        const std::array<double, 3> real_fac = _kernel.calculateTripletDerivative(dIJ2, dJK2, dKI2);

        errorX = real_fac.at(0) - fac_ij;
        errorY = real_fac.at(1) - fac_jk;
        errorZ = real_fac.at(2) + fac_ki;

        relErrorX = errorX / real_fac.at(0);
        relErrorY = errorY / real_fac.at(1);
        relErrorZ = errorZ / real_fac.at(2);
#endif

        if constexpr (calculateGlobals) {
          
          _aosThreadDataGlobals[threadnum].absErrorSum += {errorX, errorY, errorZ};
          _aosThreadDataGlobals[threadnum].relErrorSum += {relErrorX, relErrorY, relErrorZ};
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
        _relErrorSum = {0., 0., 0.};
        _postProcessed = false;
        if constexpr (calculateGlobals) {
          for (auto &data : _aosThreadDataGlobals) {
            data.setZero();
          }
        }
        if constexpr (countFLOPs) {
          for (auto &data : _aosThreadDataFLOPs) {
            data.setZero();
          }
        }

        for (auto & stat : _distanceStatistics) {
          stat.clear();
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
            _relErrorSum += data.relErrorSum;
          }

          size_t numKernelCalls = 0;
          for (const auto &data : _aosThreadDataFLOPs) {
            numKernelCalls += data.numKernelCallsN3;
            numKernelCalls += data.numKernelCallsNoN3;
          }
    
          // For each interaction, we added the full contribution for all three particles. Divide by 3 here, so that each
          // contribution is only counted once per triplet.
          _potentialEnergySum /= 3.;
    
          // Additionally, we have always calculated 3*potentialEnergy, so we divide by 3 again.
          _potentialEnergySum /= 3.;
    
          _postProcessed = true;

          size_t number = 0;
          std::vector<double> cuts{};
          for (double cut = _a.at(0); cut <= _b.at(0); cut+= 0.25) {
            cuts.push_back(cut);
          }
          double prev_cut = 0.;
          for (auto cut : cuts) {
            number = 0;
            for (const auto &stat : _distanceStatistics) {
              number += std::count_if(stat.begin(), stat.end(),
                                      [cut, prev_cut](double value) { return value <= cut && value > prev_cut; });
            }
            AutoPasLog(DEBUG, "Distances below {} : {}", cut, number);
            prev_cut = cut;
          }
          
          AutoPasLog(TRACE, "Final potential energy {}", _potentialEnergySum);
          AutoPasLog(TRACE, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
          AutoPasLog(DEBUG, "Final absolute error   {} {} {}", _absErrorSum.at(0), _absErrorSum.at(1), _absErrorSum.at(2));
          AutoPasLog(DEBUG, "Final relative error   {} {} {}", _relErrorSum.at(0), _relErrorSum.at(1), _relErrorSum.at(2));
          AutoPasLog(DEBUG, "Number of kernel calls {}", numKernelCalls);
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
         AoSThreadDataGlobals() : virialSum{0.,0.,0.}, potentialEnergySum{0.}, absErrorSum{0.,0.,0.}, relErrorSum{0.,0.,0.}, __remainingTo64{} {}
         void setZero() {
           virialSum = {0., 0., 0.};
           absErrorSum = {0., 0., 0.};
           relErrorSum = {0., 0., 0.};
           potentialEnergySum = 0.;
         }
     
         // variables
         std::array<double, 3> absErrorSum;
         std::array<double, 3> relErrorSum;
         std::array<double, 3> virialSum;
         double potentialEnergySum;
     
        private:
         // dummy parameter to get the right size (64 bytes)
         double __remainingTo64[(64 - 2 * sizeof(double)) / sizeof(double)];
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
       std::array<double, 3> _relErrorSum;

       std::vector<AoSThreadDataGlobals> _aosThreadDataGlobals {};
       std::vector<AoSThreadDataFLOPs> _aosThreadDataFLOPs{};

       std::vector<std::vector<double>> _distanceStatistics{};

       bool _postProcessed;

       /* Interpolation Parameters */
       TriwiseKernel _kernel;
       const std::array<std::vector<size_t>, 3> _numNodes {};
       const std::array<std::vector<double>, 3> _intervalSplits {};

       Eigen::Tensor<
          std::vector<std::vector<std::vector<double>>>,
          3,
          Eigen::RowMajor> _coeffs  {};

       Eigen::Tensor<
          Eigen::Tensor<double, 3, Eigen::RowMajor>,
          3,
          Eigen::RowMajor> _coeffsVec  {};

       const double PI = 2 * std::acos(0.);

       std::array<double, 3> _a {};
       std::array<double, 3> _b {};
    };
}