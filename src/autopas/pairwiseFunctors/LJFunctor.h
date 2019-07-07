/**
 * @file LJFunctor.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"
#if defined(AUTOPAS_CUDA)
#include "autopas/pairwiseFunctors/LJFunctorCuda.cuh"
#include "autopas/utils/CudaStreamHandler.h"
#else
#include "autopas/utils/ExceptionHandler.h"
#endif

namespace autopas {

/**
 * Newton 3 modes for the LJFunctor.
 */
enum class FunctorN3Modes {
  Newton3Only,
  Newton3Off,
  Both,
};

/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * @tparam Particle The type of particle.
 * @tparam ParticleCell The type of particlecell.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, class ParticleCell, FunctorN3Modes useNewton3 = FunctorN3Modes::Both,
          bool calculateGlobals = false, bool relevantForTuning = true>
class LJFunctor
    : public Functor<Particle, ParticleCell, typename Particle::SoAArraysType, LJFunctor<Particle, ParticleCell>> {
  using SoAArraysType = typename Particle::SoAArraysType;
  using floatPrecision = typename Particle::ParticleFloatingPointType;

 public:
  /**
   * Deleted default constructor
   */
  LJFunctor() = delete;

  /**
   * Constructor, which sets the global values, i.e. cutoff, epsilon, sigma and shift.
   * @param cutoff
   * @param epsilon
   * @param sigma
   * @param shift
   * @param duplicatedCalculation Defines whether duplicated calculations are happening across processes / over the
   * simulation boundary. e.g. eightShell: false, fullShell: true.
   */
  explicit LJFunctor(floatPrecision cutoff, floatPrecision epsilon, floatPrecision sigma, floatPrecision shift,
                     bool duplicatedCalculation = true)
      : Functor<Particle, ParticleCell, SoAArraysType, LJFunctor<Particle, ParticleCell>>(cutoff),
        _cutoffsquare{cutoff * cutoff},
        _epsilon24{epsilon * (floatPrecision)24.0},
        _sigmasquare{sigma * sigma},
        _shift6{shift * (floatPrecision)6.0},
        _upotSum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _duplicatedCalculations{duplicatedCalculation},
        _postProcessed {
    false
  }
#if defined(AUTOPAS_CUDA)
  , _streams(32)
#endif
  {
    if (calculateGlobals) {
      _aosThreadData.resize(autopas_get_max_threads());
    }
#if defined(AUTOPAS_CUDA)
    LJFunctorConstants<floatPrecision> constants(_cutoffsquare, _epsilon24, _sigmasquare, _shift6);
    _cudawrapper.loadConstants(&constants);
#endif
  }

  bool isRelevantForTuning() override { return relevantForTuning; }

  bool allowsNewton3() override {
    return useNewton3 == FunctorN3Modes::Newton3Only or useNewton3 == FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() override {
    return useNewton3 == FunctorN3Modes::Newton3Off or useNewton3 == FunctorN3Modes::Both;
  }

  void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
    auto dr = ArrayMath::sub(i.getR(), j.getR());
    floatPrecision dr2 = ArrayMath::dot(dr, dr);

    if (dr2 > _cutoffsquare) return;

    floatPrecision invdr2 = 1. / dr2;
    floatPrecision lj6 = _sigmasquare * invdr2;
    lj6 = lj6 * lj6 * lj6;
    floatPrecision lj12 = lj6 * lj6;
    floatPrecision lj12m6 = lj12 - lj6;
    floatPrecision fac = _epsilon24 * (lj12 + lj12m6) * invdr2;
    auto f = ArrayMath::mulScalar(dr, fac);
    i.addF(f);
    if (newton3) {
      // only if we use newton 3 here, we want to
      j.subF(f);
    }
    if (calculateGlobals) {
      auto virial = ArrayMath::mul(dr, f);
      floatPrecision upot = _epsilon24 * lj12m6 + _shift6;

      const int threadnum = autopas_get_thread_num();
      if (_duplicatedCalculations) {
        // for non-newton3 the division is in the post-processing step.
        if (newton3) {
          upot *= 0.5;
          virial = ArrayMath::mulScalar(virial, (floatPrecision)0.5);
        }
        if (i.isOwned()) {
          _aosThreadData[threadnum].upotSum += upot;
          _aosThreadData[threadnum].virialSum = ArrayMath::add(_aosThreadData[threadnum].virialSum, virial);
        }
        // for non-newton3 the second particle will be considered in a separate calculation
        if (newton3 and j.isOwned()) {
          _aosThreadData[threadnum].upotSum += upot;
          _aosThreadData[threadnum].virialSum = ArrayMath::add(_aosThreadData[threadnum].virialSum, virial);
        }
      } else {
        // for non-newton3 we divide by 2 only in the postprocess step!
        _aosThreadData[threadnum].upotSum += upot;
        _aosThreadData[threadnum].virialSum = ArrayMath::add(_aosThreadData[threadnum].virialSum, virial);
      }
    }
  }

  /**
   * @copydoc Functor::SoAFunctor(SoA<SoAArraysType> &soa, bool newton3)
   * This functor ignores will use a newton3 like traversing of the soa, however, it still needs to know about newton3
   * to use it correctly for the global values.
   */
  void SoAFunctor(SoA<SoAArraysType> &soa, bool newton3) override {
    if (soa.getNumParticles() == 0) return;

    const auto *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict__ ownedPtr = soa.template begin<Particle::AttributeNames::owned>();

    auto *const __restrict__ fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict__ fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict__ fzptr = soa.template begin<Particle::AttributeNames::forceZ>();
    // the local redeclaration of the following values helps the auto-generation of various compilers.
    const floatPrecision cutoffsquare = _cutoffsquare, epsilon24 = _epsilon24, sigmasquare = _sigmasquare,
                         shift6 = _shift6;

    if (calculateGlobals) {
      // Checks if the cell is a halo cell, if it is, we skip it.
      bool isHaloCell = ownedPtr[0] ? false : true;
      if (isHaloCell) {
        return;
      }
    }

    floatPrecision upotSum = 0.;
    floatPrecision virialSumX = 0.;
    floatPrecision virialSumY = 0.;
    floatPrecision virialSumZ = 0.;

    for (unsigned int i = 0; i < soa.getNumParticles(); ++i) {
      floatPrecision fxacc = 0.;
      floatPrecision fyacc = 0.;
      floatPrecision fzacc = 0.;

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, upotSum, virialSumX, virialSumY, virialSumZ)
      for (unsigned int j = i + 1; j < soa.getNumParticles(); ++j) {
        const floatPrecision drx = xptr[i] - xptr[j];
        const floatPrecision dry = yptr[i] - yptr[j];
        const floatPrecision drz = zptr[i] - zptr[j];

        const floatPrecision drx2 = drx * drx;
        const floatPrecision dry2 = dry * dry;
        const floatPrecision drz2 = drz * drz;

        const floatPrecision dr2 = drx2 + dry2 + drz2;

        const floatPrecision mask = (dr2 > cutoffsquare) ? 0. : 1.;

        const floatPrecision invdr2 = 1. / dr2;
        const floatPrecision lj2 = sigmasquare * invdr2;
        const floatPrecision lj6 = lj2 * lj2 * lj2;
        const floatPrecision lj12 = lj6 * lj6;
        const floatPrecision lj12m6 = lj12 - lj6;
        const floatPrecision fac = epsilon24 * (lj12 + lj12m6) * invdr2 * mask;

        const floatPrecision fx = drx * fac;
        const floatPrecision fy = dry * fac;
        const floatPrecision fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;

        // newton 3
        fxptr[j] -= fx;
        fyptr[j] -= fy;
        fzptr[j] -= fz;

        if (calculateGlobals) {
          const floatPrecision virialx = drx * fx;
          const floatPrecision virialy = dry * fy;
          const floatPrecision virialz = drz * fz;
          const floatPrecision upot = (epsilon24 * lj12m6 + shift6) * mask;

          // these calculations assume that this functor is not called for halo cells!
          upotSum += upot;
          virialSumX += virialx;
          virialSumY += virialy;
          virialSumZ += virialz;
        }
      }

      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;
    }
    if (calculateGlobals) {
      const int threadnum = autopas_get_thread_num();
      // we assume newton3 to be enabled in this functor call, thus we multiply by two if the value of newton3 is false,
      // since for newton3 disabled we divide by two later on.
      _aosThreadData[threadnum].upotSum += upotSum * (newton3 ? 1. : 2.);
      _aosThreadData[threadnum].virialSum[0] += virialSumX * (newton3 ? 1. : 2.);
      _aosThreadData[threadnum].virialSum[1] += virialSumY * (newton3 ? 1. : 2.);
      _aosThreadData[threadnum].virialSum[2] += virialSumZ * (newton3 ? 1. : 2.);
    }
  }

  /**
   * @copydoc Functor::SoAFunctor(SoA<SoAArraysType> &soa1, SoA<SoAArraysType> &soa2, bool newton3)
   */
  void SoAFunctor(SoA<SoAArraysType> &soa1, SoA<SoAArraysType> &soa2, const bool newton3) override {
    if (soa1.getNumParticles() == 0 || soa2.getNumParticles() == 0) return;

    const auto *const __restrict__ x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict__ y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict__ z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict__ x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict__ y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict__ z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict__ ownedPtr1 = soa1.template begin<Particle::AttributeNames::owned>();
    const auto *const __restrict__ ownedPtr2 = soa2.template begin<Particle::AttributeNames::owned>();

    auto *const __restrict__ fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict__ fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict__ fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict__ fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict__ fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict__ fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    bool isHaloCell1 = false;
    bool isHaloCell2 = false;
    // Checks whether the cells are halo cells.
    // This check cannot be done if _lowCorner and _highCorner are not set. So we do this only if calculateGlobals is
    // defined. (as of 23.11.2018)
    if (calculateGlobals) {
      isHaloCell1 = ownedPtr1[0] ? false : true;
      isHaloCell2 = ownedPtr2[0] ? false : true;

      // This if is commented out because the AoS vs SoA test would fail otherwise. Even though it is physically
      // correct!
      /*if(_duplicatedCalculations and isHaloCell1 and isHaloCell2){
        return;
      }*/
    }
    floatPrecision upotSum = 0.;
    floatPrecision virialSumX = 0.;
    floatPrecision virialSumY = 0.;
    floatPrecision virialSumZ = 0.;

    const floatPrecision cutoffsquare = _cutoffsquare, epsilon24 = _epsilon24, sigmasquare = _sigmasquare,
                         shift6 = _shift6;
    for (unsigned int i = 0; i < soa1.getNumParticles(); ++i) {
      floatPrecision fxacc = 0;
      floatPrecision fyacc = 0;
      floatPrecision fzacc = 0;

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, upotSum, virialSumX, virialSumY, virialSumZ)
      for (unsigned int j = 0; j < soa2.getNumParticles(); ++j) {
        const floatPrecision drx = x1ptr[i] - x2ptr[j];
        const floatPrecision dry = y1ptr[i] - y2ptr[j];
        const floatPrecision drz = z1ptr[i] - z2ptr[j];

        const floatPrecision drx2 = drx * drx;
        const floatPrecision dry2 = dry * dry;
        const floatPrecision drz2 = drz * drz;

        const floatPrecision dr2 = drx2 + dry2 + drz2;

        const floatPrecision mask = (dr2 > cutoffsquare) ? 0. : 1.;

        const floatPrecision invdr2 = 1. / dr2;
        const floatPrecision lj2 = sigmasquare * invdr2;
        const floatPrecision lj6 = lj2 * lj2 * lj2;
        const floatPrecision lj12 = lj6 * lj6;
        const floatPrecision lj12m6 = lj12 - lj6;
        const floatPrecision fac = epsilon24 * (lj12 + lj12m6) * invdr2 * mask;

        const floatPrecision fx = drx * fac;
        const floatPrecision fy = dry * fac;
        const floatPrecision fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;
        if (newton3) {
          fx2ptr[j] -= fx;
          fy2ptr[j] -= fy;
          fz2ptr[j] -= fz;
        }

        if (calculateGlobals) {
          floatPrecision virialx = drx * fx;
          floatPrecision virialy = dry * fy;
          floatPrecision virialz = drz * fz;
          floatPrecision upot = (epsilon24 * lj12m6 + shift6) * mask;

          upotSum += upot;
          virialSumX += virialx;
          virialSumY += virialy;
          virialSumZ += virialz;
        }
      }
      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }
    if (calculateGlobals) {
      floatPrecision energyfactor = 1.;
      if (_duplicatedCalculations) {
        // if we have duplicated calculations, i.e., we calculate interactions multiple times, we have to take care
        // that we do not add the energy multiple times!
        energyfactor = isHaloCell1 ? 0. : 1.;
        if (newton3) {
          energyfactor += isHaloCell2 ? 0. : 1.;
          energyfactor *= 0.5;  // we count the energies partly to one of the two cells!
        }
      }

      const int threadnum = autopas_get_thread_num();

      _aosThreadData[threadnum].upotSum += upotSum * energyfactor;
      _aosThreadData[threadnum].virialSum[0] += virialSumX * energyfactor;
      _aosThreadData[threadnum].virialSum[1] += virialSumY * energyfactor;
      _aosThreadData[threadnum].virialSum[2] += virialSumZ * energyfactor;
    }
  }

  // clang-format off
  /**
   * @copydoc Functor::SoAFunctor(SoA<SoAArraysType> &soa, const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborList, size_t iFrom, size_t iTo, bool newton3)
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors and specify iFrom and iTo accordingly.
   */
  // clang-format on
  void SoAFunctor(SoA<SoAArraysType> &soa,
                  const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborList, size_t iFrom,
                  size_t iTo, const bool newton3) override {
    if (newton3) {
      if (_duplicatedCalculations) {
        SoAFunctorImpl<true, true>(soa, neighborList, iFrom, iTo);
      } else {
        SoAFunctorImpl<true, false>(soa, neighborList, iFrom, iTo);
      }
    } else {
      if (_duplicatedCalculations) {
        SoAFunctorImpl<false, true>(soa, neighborList, iFrom, iTo);
      } else {
        SoAFunctorImpl<false, false>(soa, neighborList, iFrom, iTo);
      }
    }
  }

  /**
   * @brief Functor using Cuda on SoA in device Memory
   *
   * This Functor calculates the pair-wise interactions between particles in the device_handle on the GPU
   *
   * @param device_handle soa in device memory
   * @param newton3 defines whether or whether not to use newton
   */
  void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle, bool newton3) override {
#if defined(AUTOPAS_CUDA)
    const size_t size = device_handle.template get<Particle::AttributeNames::posX>().size();
    cudaStream_t stream = _streams.getStreamRandom();

    auto cudaSoA = this->createFunctorCudaSoA(device_handle);
    if (newton3) {
      _cudawrapper.SoAFunctorN3Wrapper(cudaSoA.get(), stream);
    } else {
      _cudawrapper.SoAFunctorNoN3Wrapper(cudaSoA.get(), stream);
    }

#else
    utils::ExceptionHandler::exception("AutoPas was compiled without CUDA support!");
#endif
  }

  /**
   * @brief Functor using Cuda on SoAs in device Memory
   *
   * This Functor calculates the pair-wise interactions between particles in the device_handle1 and device_handle2 on
   * the GPU
   *
   * @param device_handle1 first soa in device memory
   * @param device_handle2 second soa in device memory
   * @param newton3 defines whether or whether not to use newton
   */
  void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle1,
                   CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle2, bool newton3) override {
#if defined(AUTOPAS_CUDA)
    const size_t size1 = device_handle1.template get<Particle::AttributeNames::posX>().size();
    const size_t size2 = device_handle2.template get<Particle::AttributeNames::posX>().size();
    cudaStream_t stream = _streams.getStreamRandom();
    auto cudaSoA1 = this->createFunctorCudaSoA(device_handle1);
    auto cudaSoA2 = this->createFunctorCudaSoA(device_handle2);

    if (newton3) {
      if (size1 > size2) {
        _cudawrapper.SoAFunctorN3PairWrapper(cudaSoA1.get(), cudaSoA2.get(), stream);
      } else {
        _cudawrapper.SoAFunctorN3PairWrapper(cudaSoA2.get(), cudaSoA1.get(), stream);
      }
    } else {
      _cudawrapper.SoAFunctorNoN3PairWrapper(cudaSoA1.get(), cudaSoA2.get(), stream);
    }
#else
    utils::ExceptionHandler::exception("AutoPas was compiled without CUDA support!");
#endif
  }
#if defined(AUTOPAS_CUDA)
  CudaWrapperInterface<floatPrecision> *getCudaWrapper() override { return &_cudawrapper; }

  std::unique_ptr<FunctorCudaSoA<floatPrecision>> createFunctorCudaSoA(
      CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) override {
    return std::make_unique<LJFunctorCudaSoA<floatPrecision>>(
        device_handle.template get<Particle::AttributeNames::posX>().size(),
        device_handle.template get<Particle::AttributeNames::posX>().get(),
        device_handle.template get<Particle::AttributeNames::posY>().get(),
        device_handle.template get<Particle::AttributeNames::posZ>().get(),
        device_handle.template get<Particle::AttributeNames::forceX>().get(),
        device_handle.template get<Particle::AttributeNames::forceY>().get(),
        device_handle.template get<Particle::AttributeNames::forceZ>().get());
  }
#endif

  /**
   * @copydoc Functor::deviceSoALoader
   */
  void deviceSoALoader(::autopas::SoA<SoAArraysType> &soa,
                       CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) override {
#if defined(AUTOPAS_CUDA)

    size_t size = soa.getNumParticles();
    if (size == 0) return;

    device_handle.template get<Particle::AttributeNames::posX>().copyHostToDevice(
        size, soa.template begin<Particle::AttributeNames::posX>());
    device_handle.template get<Particle::AttributeNames::posY>().copyHostToDevice(
        size, soa.template begin<Particle::AttributeNames::posY>());
    device_handle.template get<Particle::AttributeNames::posZ>().copyHostToDevice(
        size, soa.template begin<Particle::AttributeNames::posZ>());

    device_handle.template get<Particle::AttributeNames::forceX>().copyHostToDevice(
        size, soa.template begin<Particle::AttributeNames::forceX>());
    device_handle.template get<Particle::AttributeNames::forceY>().copyHostToDevice(
        size, soa.template begin<Particle::AttributeNames::forceY>());
    device_handle.template get<Particle::AttributeNames::forceZ>().copyHostToDevice(
        size, soa.template begin<Particle::AttributeNames::forceZ>());
#else
    utils::ExceptionHandler::exception("AutoPas was compiled without CUDA support!");
#endif
  }

  /**
   * @copydoc Functor::deviceSoAExtractor
   */
  void deviceSoAExtractor(::autopas::SoA<SoAArraysType> &soa,
                          CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) override {
#if defined(AUTOPAS_CUDA)

    size_t size = soa.getNumParticles();
    if (size == 0) return;

    device_handle.template get<Particle::AttributeNames::forceX>().copyDeviceToHost(
        size, soa.template begin<Particle::AttributeNames::forceX>());
    device_handle.template get<Particle::AttributeNames::forceY>().copyDeviceToHost(
        size, soa.template begin<Particle::AttributeNames::forceY>());
    device_handle.template get<Particle::AttributeNames::forceZ>().copyDeviceToHost(
        size, soa.template begin<Particle::AttributeNames::forceZ>());

#else
    utils::ExceptionHandler::exception("AutoPas was compiled without CUDA support!");
#endif
  }

  /**
   * @copydoc Functor::getNeededAttr()
   */
  constexpr static const std::array<typename Particle::AttributeNames, 8> getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 8>{
        Particle::AttributeNames::id,     Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,   Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ, Particle::AttributeNames::owned};
  }

  /**
   * @copydoc Functor::getNeededAttr(std::false_type)
   */
  constexpr static const std::array<typename Particle::AttributeNames, 5> getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 5>{
        Particle::AttributeNames::id, Particle::AttributeNames::posX, Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ, Particle::AttributeNames::owned};
  }

  /**
   * @copydoc Functor::getComputedAttr()
   */
  constexpr static const std::array<typename Particle::AttributeNames, 3> getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 3>{
        Particle::AttributeNames::forceX, Particle::AttributeNames::forceY, Particle::AttributeNames::forceZ};
  }

  /**
   * Get the number of flops used per kernel call. This should count the
   * floating point operations needed for two particles that lie within a cutoff
   * radius.
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall() {
    // Kernel: 12 = 1 (inverse R squared) + 8 (compute scale) + 3 (apply
    // scale) sum Forces: 6 (forces) kernel total = 12 + 6 = 18
    return 18ul;
  }

  /**
   * Reset the global values.
   * Will set the global values to zero to prepare for the next iteration.
   */
  void initTraversal() override {
    _upotSum = 0.;
    _virialSum = {0., 0., 0.};
    _postProcessed = false;
    for (size_t i = 0; i < _aosThreadData.size(); ++i) {
      _aosThreadData[i].setZero();
    }
  }

  /**
   * Postprocesses global values, e.g. upot and virial
   * @param newton3
   */
  void endTraversal(bool newton3) override {
    if (_postProcessed) {
      throw utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
    }
    for (size_t i = 0; i < _aosThreadData.size(); ++i) {
      _upotSum += _aosThreadData[i].upotSum;
      _virialSum = ArrayMath::add(_virialSum, _aosThreadData[i].virialSum);
    }
    if (not newton3) {
      // if the newton3 optimization is disabled we have added every energy contribution twice, so we divide by 2
      // here.
      _upotSum *= 0.5;
      _virialSum = ArrayMath::mulScalar(_virialSum, (floatPrecision)0.5);
    }
    // we have always calculated 6*upot, so we divide by 6 here!
    _upotSum /= 6.;
    _postProcessed = true;
  }

  /**
   * Get the potential Energy
   * @return the potential Energy
   */
  floatPrecision getUpot() {
    if (not calculateGlobals) {
      throw utils::ExceptionHandler::AutoPasException(
          "Trying to get upot even though calculateGlobals is false. If you want this functor to calculate global "
          "values, please specify calculateGlobals to be true.");
    }
    if (not _postProcessed) {
      throw utils::ExceptionHandler::AutoPasException("Cannot get upot, because endTraversal was not called.");
    }
    return _upotSum;
  }

  /**
   * Get the virial
   * @return the virial
   */
  floatPrecision getVirial() {
    if (not calculateGlobals) {
      throw utils::ExceptionHandler::AutoPasException(
          "Trying to get virial even though calculateGlobals is false. If you want this functor to calculate global "
          "values, please specify calculateGlobals to be true.");
    }
    if (not _postProcessed) {
      throw utils::ExceptionHandler::AutoPasException("Cannot get virial, because endTraversal was not called.");
    }
    return _virialSum[0] + _virialSum[1] + _virialSum[2];
  }

 private:
  template <bool newton3, bool duplicatedCalculations>
  void SoAFunctorImpl(SoA<SoAArraysType> &soa,
                      const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborList,
                      size_t iFrom, size_t iTo) {
    if (soa.getNumParticles() == 0) return;

    const auto *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

    auto *const __restrict__ fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict__ fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict__ fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict__ ownedPtr = soa.template begin<Particle::AttributeNames::owned>();

    const floatPrecision cutoffsquare = _cutoffsquare, epsilon24 = _epsilon24, sigmasquare = _sigmasquare,
                         shift6 = _shift6;

    floatPrecision upotSum = 0.;
    floatPrecision virialSumX = 0.;
    floatPrecision virialSumY = 0.;
    floatPrecision virialSumZ = 0.;

    for (size_t i = iFrom; i < iTo; ++i) {
      floatPrecision fxacc = 0;
      floatPrecision fyacc = 0;
      floatPrecision fzacc = 0;
      const size_t listSizeI = neighborList[i].size();
      const size_t *const __restrict__ currentList = neighborList[i].data();

      // checks whether particle 1 is in the domain box, unused if _duplicatedCalculations is false!
      floatPrecision inbox1Mul = 0.;
      if (duplicatedCalculations) {  // only for duplicated calculations we need this value
        inbox1Mul = ownedPtr[i];
        if (newton3) {
          inbox1Mul *= 0.5;
        }
      }

      // this is a magic number, that should correspond to at least
      // vectorization width*N have testet multiple sizes:
      // 4: does not give a speedup, slower than original AoSFunctor
      // 8: small speedup compared to AoS
      // 12: highest speedup compared to Aos
      // 16: smaller speedup
      // in theory this is a variable, we could auto-tune over...
#ifdef __AVX512F__
      // use a multiple of 8 for avx
      constexpr size_t vecsize = 16;
#else
      // for everything else 12 is faster
      constexpr size_t vecsize = 12;
#endif
      size_t joff = 0;

      // if the size of the verlet list is larger than the given size vecsize,
      // we will use a vectorized version.
      if (listSizeI >= vecsize) {
        alignas(64) std::array<floatPrecision, vecsize> xtmp, ytmp, ztmp, xArr, yArr, zArr, fxArr, fyArr, fzArr,
            ownedArr;
        // broadcast of the position of particle i
        for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
          xtmp[tmpj] = xptr[i];
          ytmp[tmpj] = yptr[i];
          ztmp[tmpj] = zptr[i];
        }
        // loop over the verlet list from 0 to x*vecsize
        for (; joff < listSizeI - vecsize + 1; joff += vecsize) {
          // in each iteration we calculate the interactions of particle i with
          // vecsize particles in the neighborlist of particle i starting at
          // particle joff

          // gather position of particle j
#pragma omp simd safelen(vecsize)
          for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
            xArr[tmpj] = xptr[currentList[joff + tmpj]];
            yArr[tmpj] = yptr[currentList[joff + tmpj]];
            zArr[tmpj] = zptr[currentList[joff + tmpj]];
            ownedArr[tmpj] = ownedPtr[currentList[joff + tmpj]];
          }
          // do omp simd with reduction of the interaction
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, upotSum, virialSumX, virialSumY, virialSumZ) safelen(vecsize)
          for (size_t j = 0; j < vecsize; j++) {
            // const size_t j = currentList[jNeighIndex];

            const floatPrecision drx = xtmp[j] - xArr[j];
            const floatPrecision dry = ytmp[j] - yArr[j];
            const floatPrecision drz = ztmp[j] - zArr[j];

            const floatPrecision drx2 = drx * drx;
            const floatPrecision dry2 = dry * dry;
            const floatPrecision drz2 = drz * drz;

            const floatPrecision dr2 = drx2 + dry2 + drz2;

            const floatPrecision mask = (dr2 <= cutoffsquare) ? 1. : 0.;

            const floatPrecision invdr2 = 1. / dr2 * mask;
            const floatPrecision lj2 = sigmasquare * invdr2;
            const floatPrecision lj6 = lj2 * lj2 * lj2;
            const floatPrecision lj12 = lj6 * lj6;
            const floatPrecision lj12m6 = lj12 - lj6;
            const floatPrecision fac = epsilon24 * (lj12 + lj12m6) * invdr2;

            const floatPrecision fx = drx * fac;
            const floatPrecision fy = dry * fac;
            const floatPrecision fz = drz * fac;

            fxacc += fx;
            fyacc += fy;
            fzacc += fz;
            if (newton3) {
              fxArr[j] = fx;
              fyArr[j] = fy;
              fzArr[j] = fz;
            }
            if (calculateGlobals) {
              floatPrecision virialx = drx * fx;
              floatPrecision virialy = dry * fy;
              floatPrecision virialz = drz * fz;
              floatPrecision upot = (epsilon24 * lj12m6 + shift6) * mask;

              if (duplicatedCalculations) {
                // for non-newton3 the division is in the post-processing step.
                floatPrecision inboxMul = inbox1Mul + (newton3 ? ownedArr[j] * .5 : 0.);
                upotSum += upot * inboxMul;
                virialSumX += virialx * inboxMul;
                virialSumY += virialy * inboxMul;
                virialSumZ += virialz * inboxMul;

              } else {
                // for non-newton3 we divide by 2 only in the postprocess step!
                upotSum += upot;
                virialSumX += virialx;
                virialSumY += virialy;
                virialSumZ += virialz;
              }
            }
          }
          // scatter the forces to where they belong, this is only needed for newton3
          if (newton3) {
#pragma omp simd safelen(vecsize)
            for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
              const size_t j = currentList[joff + tmpj];
              fxptr[j] -= fxArr[tmpj];
              fyptr[j] -= fyArr[tmpj];
              fzptr[j] -= fzArr[tmpj];
            }
          }
        }
      }
      // this loop goes over the remainder and uses no optimizations
      for (size_t jNeighIndex = joff; jNeighIndex < listSizeI; ++jNeighIndex) {
        size_t j = neighborList[i][jNeighIndex];
        if (i == j) continue;

        const floatPrecision drx = xptr[i] - xptr[j];
        const floatPrecision dry = yptr[i] - yptr[j];
        const floatPrecision drz = zptr[i] - zptr[j];

        const floatPrecision drx2 = drx * drx;
        const floatPrecision dry2 = dry * dry;
        const floatPrecision drz2 = drz * drz;

        const floatPrecision dr2 = drx2 + dry2 + drz2;

        if (dr2 > cutoffsquare) continue;

        const floatPrecision invdr2 = 1. / dr2;
        const floatPrecision lj2 = sigmasquare * invdr2;
        const floatPrecision lj6 = lj2 * lj2 * lj2;
        const floatPrecision lj12 = lj6 * lj6;
        const floatPrecision lj12m6 = lj12 - lj6;
        const floatPrecision fac = epsilon24 * (lj12 + lj12m6) * invdr2;

        const floatPrecision fx = drx * fac;
        const floatPrecision fy = dry * fac;
        const floatPrecision fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;
        if (newton3) {
          fxptr[j] -= fx;
          fyptr[j] -= fy;
          fzptr[j] -= fz;
        }
        if (calculateGlobals) {
          floatPrecision virialx = drx * fx;
          floatPrecision virialy = dry * fy;
          floatPrecision virialz = drz * fz;
          floatPrecision upot = (epsilon24 * lj12m6 + shift6);

          if (duplicatedCalculations) {
            // for non-newton3 the division is in the post-processing step.
            if (newton3) {
              upot *= 0.5;
              virialx *= 0.5;
              virialy *= 0.5;
              virialz *= 0.5;
            }
            if (inbox1Mul) {
              upotSum += upot;
              virialSumX += virialx;
              virialSumY += virialy;
              virialSumZ += virialz;
            }
            // for non-newton3 the second particle will be considered in a separate calculation
            if (newton3 and ownedPtr[j]) {
              upotSum += upot;
              virialSumX += virialx;
              virialSumY += virialy;
              virialSumZ += virialz;
            }
          } else {
            // for non-newton3 we divide by 2 only in the postprocess step!
            upotSum += upot;
            virialSumX += virialx;
            virialSumY += virialy;
            virialSumZ += virialz;
          }
        }
      }

      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;
    }

    if (calculateGlobals) {
      const int threadnum = autopas_get_thread_num();

      _aosThreadData[threadnum].upotSum += upotSum;
      _aosThreadData[threadnum].virialSum[0] += virialSumX;
      _aosThreadData[threadnum].virialSum[1] += virialSumY;
      _aosThreadData[threadnum].virialSum[2] += virialSumZ;
    }
  }

  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
   */
  class AoSThreadData {
   public:
    AoSThreadData() : virialSum{0., 0., 0.}, upotSum{0.}, __remainingTo64{} {}
    void setZero() {
      virialSum = {0., 0., 0.};
      upotSum = 0.;
    }

    // variables
    std::array<floatPrecision, 3> virialSum;
    floatPrecision upotSum;

   private:
    // dummy parameter to get the right size (64 bytes)
    floatPrecision __remainingTo64[4];
  };
  // make sure of the size of AoSThreadData
  // static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size");

  const floatPrecision _cutoffsquare, _epsilon24, _sigmasquare, _shift6;

  // sum of the potential energy, only calculated if calculateGlobals is true
  floatPrecision _upotSum;
  // sum of the virial, only calculated if calculateGlobals is true
  std::array<floatPrecision, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadData> _aosThreadData;

  // bool that defines whether duplicate calculations are happening
  bool _duplicatedCalculations;

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

#if defined(AUTOPAS_CUDA)
  // contains wrapper functions for cuda calls
  LJFunctorCudaWrapper<floatPrecision> _cudawrapper;

  // Handles all cuda streams
  utils::CudaStreamHandler _streams;
#endif

};  // class LJFunctor
}  // namespace autopas
