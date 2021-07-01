/**
 * @file LJFunctorISPC.h
 *
 * @date 15 June 2021
 * @author humig
 */

#pragma once

#include <array>

#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"
#if defined(AUTOPAS_CUDA)
#include "LJFunctorCudaConstants.cuh"
#include "autopas/molecularDynamics/LJFunctorCuda.cuh"
#include "autopas/molecularDynamics/LJFunctorCudaGlobals.cuh"
#include "autopas/utils/CudaDeviceVector.h"
#include "autopas/utils/CudaStreamHandler.h"
#else
#include "autopas/utils/ExceptionHandler.h"
#endif

extern "C" {

extern void SoAFunctorSingleISPC_N3(int64_t numParticles, const int64_t* ownedStatePtr, double sigmasquare,
                                    double epsilon24, double cutoffsquare, const double* xptr, const double* yptr,
                                    const double* zptr, double* fxptr, double* fyptr, double* fzptr);


extern void SoAFunctorSingleISPC_NoN3(int64_t numParticles, const int64_t* ownedStatePtr, double sigmasquare,
                                    double epsilon24, double cutoffsquare, const double* xptr, const double* yptr,
                                    const double* zptr, double* fxptr, double* fyptr, double* fzptr);

extern void SoAFunctorPairISPC_N3(int64_t numParticles1, int64_t numParticles2,
                               const int64_t* ownedStatePtr1, const int64_t* ownedStatePtr2,
                               double sigmasquare, double epsilon24, double cutoffsquare,
                               const double* x1ptr, const double* x2ptr,
                               const double* y1ptr, const double* y2ptr,
                               const double* z1ptr, const double* z2ptr,
                               double* fx1ptr, double* fx2ptr,
                               double* fy1ptr, double* fy2ptr,
                               double* fz1ptr, double* fz2ptr);

extern void SoAFunctorPairISPC_NoN3(int64_t numParticles1, int64_t numParticles2,
                               const int64_t* ownedStatePtr1, const int64_t* ownedStatePtr2,
                               double sigmasquare, double epsilon24, double cutoffsquare,
                               const double* x1ptr, const double* x2ptr,
                               const double* y1ptr, const double* y2ptr,
                               const double* z1ptr, const double* z2ptr,
                               double* fx1ptr, double* fx2ptr,
                               double* fy1ptr, double* fy2ptr,
                               double* fz1ptr, double* fz2ptr);

extern void SoAFunctorVerletISPC_N3(size_t indexFirst, const size_t* neighborList,
                                 size_t numNeighbors,
                                 double sigmasquare, double epsilon24, double cutoffsquare,
                                 const double* x1ptr, const double* y1ptr,
                                 const double* z1ptr, const int64_t* ownedStates,
                                 double* fx1ptr, double* fy1ptr,
                                 double* fz1ptr);

extern void SoAFunctorVerletISPC_NoN3(size_t indexFirst, const size_t* neighborList,
                                 size_t numNeighbors,
                                 double sigmasquare, double epsilon24, double cutoffsquare,
                                 const double* x1ptr, const double* y1ptr,
                                 const double* z1ptr, const int64_t* ownedStates,
                                 double* fx1ptr, double* fy1ptr,
                                 double* fz1ptr);
}

namespace autopas {

/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * @tparam Particle The type of particle.
 * @tparam applyShift Switch for the lj potential to be truncated shifted.
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, bool applyShift = false, bool useMixing = false,
    FunctorN3Modes useNewton3 = FunctorN3Modes::Both, bool calculateGlobals = false,
    bool relevantForTuning = true>
class LJFunctorISPC
    : public Functor<Particle,
        LJFunctorISPC<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>> {

  static_assert(not applyShift, "Shift not supported yet");
  static_assert(not useMixing, "Mixing not supported yet");
  static_assert(not calculateGlobals, "Globals not supported yet");

  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle::SoAArraysType;

  /**
   * Precision of SoA entries.
   */
  using SoAFloatPrecision = typename Particle::ParticleSoAFloatPrecision;

 public:
  /**
   * Deleted default constructor
   */
  LJFunctorISPC() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit LJFunctorISPC(double cutoff, void * /*dummy*/)
      : Functor<Particle, LJFunctorISPC<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(
      cutoff),
        _cutoffsquare{cutoff * cutoff},
        _upotSum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false} {
    if constexpr (calculateGlobals) {
      _aosThreadData.resize(autopas_get_max_threads());
    }
  }

 public:
  /**
   * Constructor for Functor with mixing disabled. When using this functor it is necessary to call
   * setParticleProperties() to set internal constants because it does not use a particle properties library.
   *
   * @note Only to be used with mixing == false.
   *
   * @param cutoff
   */
  explicit LJFunctorISPC(double cutoff) : LJFunctorISPC(cutoff, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  /**
   * Constructor for Functor with mixing active. This functor takes a ParticlePropertiesLibrary to look up (mixed)
   * properties like sigma, epsilon and shift.
   * @param cutoff
   * @param particlePropertiesLibrary
   */
  explicit LJFunctorISPC(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJFunctorISPC(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final { return useNewton3 == FunctorN3Modes::Newton3Only or useNewton3 == FunctorN3Modes::Both; }

  bool allowsNonNewton3() final {
    return useNewton3 == FunctorN3Modes::Newton3Off or useNewton3 == FunctorN3Modes::Both;
  }

  bool isAppropriateClusterSize(unsigned int clusterSize, DataLayoutOption::Value dataLayout) const final {
    if (dataLayout == DataLayoutOption::cuda) {
#if defined(AUTOPAS_CUDA)
      return _cudawrapper.isAppropriateClusterSize(clusterSize);
#endif
      return false;
    } else {
      return dataLayout == DataLayoutOption::aos;  // LJFunctor does not yet support soa for clusters.
      // The reason for this is that the owned state is not handled correctly, see #396.
    }
  }

  void AoSFunctor(Particle &i, Particle &j, bool newton3) final {
    if (i.isDummy() or j.isDummy()) {
      return;
    }
    auto sigmasquare = _sigmasquare;
    auto epsilon24 = _epsilon24;
    auto shift6 = _shift6;
    if constexpr (useMixing) {
      sigmasquare = _PPLibrary->mixingSigmaSquare(i.getTypeId(), j.getTypeId());
      epsilon24 = _PPLibrary->mixing24Epsilon(i.getTypeId(), j.getTypeId());
      if constexpr (applyShift) {
        shift6 = _PPLibrary->mixingShift6(i.getTypeId(), j.getTypeId());
      }
    }
    auto dr = utils::ArrayMath::sub(i.getR(), j.getR());
    double dr2 = utils::ArrayMath::dot(dr, dr);

    if (dr2 > _cutoffsquare) {
      return;
    }

    double invdr2 = 1. / dr2;
    double lj6 = sigmasquare * invdr2;
    lj6 = lj6 * lj6 * lj6;
    double lj12 = lj6 * lj6;
    double lj12m6 = lj12 - lj6;
    double fac = epsilon24 * (lj12 + lj12m6) * invdr2;
    auto f = utils::ArrayMath::mulScalar(dr, fac);
    i.addF(f);
    if (newton3) {
      // only if we use newton 3 here, we want to
      j.subF(f);
    }
    if (calculateGlobals) {
      auto virial = utils::ArrayMath::mul(dr, f);
      double upot = epsilon24 * lj12m6 + shift6;

      const int threadnum = autopas_get_thread_num();
      // for non-newton3 the division is in the post-processing step.
      if (newton3) {
        upot *= 0.5;
        virial = utils::ArrayMath::mulScalar(virial, (double)0.5);
      }
      if (i.isOwned()) {
        _aosThreadData[threadnum].upotSum += upot;
        _aosThreadData[threadnum].virialSum = utils::ArrayMath::add(_aosThreadData[threadnum].virialSum, virial);
      }
      // for non-newton3 the second particle will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        _aosThreadData[threadnum].upotSum += upot;
        _aosThreadData[threadnum].virialSum = utils::ArrayMath::add(_aosThreadData[threadnum].virialSum, virial);
      }
    }
  }

  /**
   * @copydoc Functor::SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3)
   * This functor ignores will use a newton3 like traversing of the soa, however, it still needs to know about newton3
   * to use it correctly for the global values.
   */
  void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) final {
    if (soa.getNumParticles() == 0) return;

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    SoAFloatPrecision *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();
    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const SoAFloatPrecision cutoffsquare = _cutoffsquare;

    SoAFloatPrecision upotSum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    std::vector<SoAFloatPrecision, AlignedAllocator<SoAFloatPrecision>> sigmaSquares;
    std::vector<SoAFloatPrecision, AlignedAllocator<SoAFloatPrecision>> epsilon24s;
    std::vector<SoAFloatPrecision, AlignedAllocator<SoAFloatPrecision>> shift6s;
    if constexpr (useMixing) {
      // Preload all sigma and epsilons for next vectorized region.
      // Not preloading and directly using the values, will produce worse results.
      sigmaSquares.resize(soa.getNumParticles());
      epsilon24s.resize(soa.getNumParticles());
      // if no mixing or mixing but no shift shift6 is constant therefore we do not need this vector.
      if constexpr (applyShift) {
        shift6s.resize(soa.getNumParticles());
      }
    }

    const SoAFloatPrecision const_shift6 = _shift6;
    const SoAFloatPrecision const_sigmasquare = _sigmasquare;
    const SoAFloatPrecision const_epsilon24 = _epsilon24;

    using OwnershipStateType = std::underlying_type_t<autopas::OwnershipState>;
    auto func = newton3 ? SoAFunctorSingleISPC_N3 : SoAFunctorSingleISPC_NoN3;
    func(soa.getNumParticles(), reinterpret_cast<const OwnershipStateType*>(ownedStatePtr), _sigmasquare, _epsilon24,
         _cutoffsquare, xptr, yptr, zptr, fxptr, fyptr, fzptr);
  }

  /**
   * @copydoc Functor::SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3)
   */
  void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, const bool newton3) final {
    if (newton3) {
      SoAFunctorPairImpl<true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1, soa2);
    }
  }

 private:
  /**
   * Implementation function of SoAFunctorPair(soa1, soa2, newton3)
   *
   * @tparam newton3
   * @param soa1
   * @param soa2
   */
  template <bool newton3>
  void SoAFunctorPairImpl(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2) {
    if (soa1.getNumParticles() == 0 || soa2.getNumParticles() == 0) return;

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();
    [[maybe_unused]] auto *const __restrict typeptr1 = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr2 = soa2.template begin<Particle::AttributeNames::typeId>();

    // Checks whether the cells are halo cells.
    SoAFloatPrecision upotSum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    const SoAFloatPrecision cutoffsquare = _cutoffsquare;
    SoAFloatPrecision shift6 = _shift6;
    SoAFloatPrecision sigmasquare = _sigmasquare;
    SoAFloatPrecision epsilon24 = _epsilon24;

    using OwnershipStateType = std::underlying_type_t<autopas::OwnershipState>;
    auto func = newton3 ? SoAFunctorPairISPC_N3 : SoAFunctorPairISPC_NoN3;
    func(soa1.getNumParticles(), soa2.getNumParticles(),
                       reinterpret_cast<const OwnershipStateType*>(ownedStatePtr1), reinterpret_cast<const OwnershipStateType*>(ownedStatePtr2),
                       sigmasquare, epsilon24, cutoffsquare, x1ptr, x2ptr, y1ptr, y2ptr, z1ptr, z2ptr,
                       fx1ptr, fx2ptr, fy1ptr, fy2ptr, fz1ptr, fz2ptr);
  }

 public:
  // clang-format off
  /**
   * @copydoc Functor::SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst, const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList, bool newton3)
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors!
   */
  // clang-format on
  void SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {
    if (soa.getNumParticles() == 0 or neighborList.empty()) return;
    if (newton3) {
      SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
    }
  }

  /**
   * Functor using Cuda on SoA in device Memory
   *
   * This Functor calculates the pair-wise interactions between particles in the device_handle on the GPU
   *
   * @param device_handle soa in device memory
   * @param newton3 defines whether or whether not to use newton
   */
  void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle, bool newton3) final {
#if defined(AUTOPAS_CUDA)
    const size_t size = device_handle.template get<Particle::AttributeNames::posX>().size();
    if (size == 0) {
      return;
    }
    auto cudaSoA = this->createFunctorCudaSoA(device_handle);
    if (newton3) {
      _cudawrapper.SoAFunctorN3Wrapper(cudaSoA.get(), 0);
    } else {
      _cudawrapper.SoAFunctorNoN3Wrapper(cudaSoA.get(), 0);
    }

#else
    utils::ExceptionHandler::exception("LJFunctor::CudaFunctor: AutoPas was compiled without CUDA support!");
#endif
  }

  /**
   * Sets the particle properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   * If compiled with CUDA this function also loads the values to the GPU.
   *
   * @param epsilon24
   * @param sigmaSquare
   */
  void setParticleProperties(SoAFloatPrecision epsilon24, SoAFloatPrecision sigmaSquare) {
    _epsilon24 = epsilon24;
    _sigmasquare = sigmaSquare;
    if (applyShift) {
      _shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(_epsilon24, _sigmasquare, _cutoffsquare);
    } else {
      _shift6 = 0.;
    }
#if defined(AUTOPAS_CUDA)
    LJFunctorConstants<SoAFloatPrecision> constants(_cutoffsquare, _epsilon24 /* epsilon24 */,
                                                    _sigmasquare /* sigmasquare */, _shift6);
    _cudawrapper.loadConstants(&constants);
#endif
  }

  /**
   * Functor using Cuda on SoAs in device Memory
   *
   * This Functor calculates the pair-wise interactions between particles in the device_handle1 and device_handle2 on
   * the GPU
   *
   * @param device_handle1 first soa in device memory
   * @param device_handle2 second soa in device memory
   * @param newton3 defines whether or whether not to use newton
   */
  void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle1,
                   CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle2, bool newton3) final {
#if defined(AUTOPAS_CUDA)
    const size_t size1 = device_handle1.template get<Particle::AttributeNames::posX>().size();
    const size_t size2 = device_handle2.template get<Particle::AttributeNames::posX>().size();
    if ((size1 == 0) or (size2 == 0)) {
      return;
    }
    auto cudaSoA1 = this->createFunctorCudaSoA(device_handle1);
    auto cudaSoA2 = this->createFunctorCudaSoA(device_handle2);

    if (newton3) {
      if (size1 > size2) {
        _cudawrapper.SoAFunctorN3PairWrapper(cudaSoA1.get(), cudaSoA2.get(), 0);
      } else {
        _cudawrapper.SoAFunctorN3PairWrapper(cudaSoA2.get(), cudaSoA1.get(), 0);
      }
    } else {
      _cudawrapper.SoAFunctorNoN3PairWrapper(cudaSoA1.get(), cudaSoA2.get(), 0);
    }
#else
    utils::ExceptionHandler::exception("AutoPas was compiled without CUDA support!");
#endif
  }
#if defined(AUTOPAS_CUDA)
  CudaWrapperInterface<SoAFloatPrecision> *getCudaWrapper() final { return &_cudawrapper; }

  std::unique_ptr<FunctorCudaSoA<SoAFloatPrecision>> createFunctorCudaSoA(
      CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) final {
    if (calculateGlobals) {
      return std::make_unique<LJFunctorCudaGlobalsSoA<SoAFloatPrecision>>(
          device_handle.template get<Particle::AttributeNames::posX>().size(),
          device_handle.template get<Particle::AttributeNames::posX>().get(),
          device_handle.template get<Particle::AttributeNames::posY>().get(),
          device_handle.template get<Particle::AttributeNames::posZ>().get(),
          device_handle.template get<Particle::AttributeNames::forceX>().get(),
          device_handle.template get<Particle::AttributeNames::forceY>().get(),
          device_handle.template get<Particle::AttributeNames::forceZ>().get(),
          device_handle.template get<Particle::AttributeNames::ownershipState>().get(), _cudaGlobals.get());
    } else {
      return std::make_unique<LJFunctorCudaSoA<SoAFloatPrecision>>(
          device_handle.template get<Particle::AttributeNames::posX>().size(),
          device_handle.template get<Particle::AttributeNames::posX>().get(),
          device_handle.template get<Particle::AttributeNames::posY>().get(),
          device_handle.template get<Particle::AttributeNames::posZ>().get(),
          device_handle.template get<Particle::AttributeNames::forceX>().get(),
          device_handle.template get<Particle::AttributeNames::forceY>().get(),
          device_handle.template get<Particle::AttributeNames::forceZ>().get(),
          device_handle.template get<Particle::AttributeNames::ownershipState>().get());
    }
  }
#endif

  /**
   * @copydoc Functor::deviceSoALoader
   */
  void deviceSoALoader(::autopas::SoA<SoAArraysType> &soa,
                       CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) final {
#if defined(AUTOPAS_CUDA)

    const size_t size = soa.getNumParticles();
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

    device_handle.template get<Particle::AttributeNames::ownershipState>().copyHostToDevice(
        size, soa.template begin<Particle::AttributeNames::ownershipState>());

#else
    utils::ExceptionHandler::exception("LJFunctor::deviceSoALoader: AutoPas was compiled without CUDA support!");
#endif
  }

  /**
   * @copydoc Functor::deviceSoAExtractor
   */
  void deviceSoAExtractor(::autopas::SoA<SoAArraysType> &soa,
                          CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) final {
#if defined(AUTOPAS_CUDA)

    const size_t size = soa.getNumParticles();
    if (size == 0) return;

    device_handle.template get<Particle::AttributeNames::forceX>().copyDeviceToHost(
        size, soa.template begin<Particle::AttributeNames::forceX>());
    device_handle.template get<Particle::AttributeNames::forceY>().copyDeviceToHost(
        size, soa.template begin<Particle::AttributeNames::forceY>());
    device_handle.template get<Particle::AttributeNames::forceZ>().copyDeviceToHost(
        size, soa.template begin<Particle::AttributeNames::forceZ>());

#else
    utils::ExceptionHandler::exception("LJFunctor::deviceSoAExtractor: AutoPas was compiled without CUDA support!");
#endif
  }

  /**
   * @copydoc Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 9>{
        Particle::AttributeNames::id,     Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,   Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ, Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
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
    return std::array<typename Particle::AttributeNames, 3>{
        Particle::AttributeNames::forceX, Particle::AttributeNames::forceY, Particle::AttributeNames::forceZ};
  }

  /**
   *
   * @return useMixing
   */
  constexpr static bool getMixing() { return useMixing; }

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
  void initTraversal() final {
    _upotSum = 0.;
    _virialSum = {0., 0., 0.};
    _postProcessed = false;
    for (size_t i = 0; i < _aosThreadData.size(); ++i) {
      _aosThreadData[i].setZero();
    }
#if defined(AUTOPAS_CUDA)
    if (calculateGlobals) {
      std::array<SoAFloatPrecision, 4> globals{0, 0, 0, 0};
      _cudaGlobals.copyHostToDevice(4, globals.data());
    }
#endif
  }

  /**
   * Postprocesses global values, e.g. upot and virial
   * @param newton3
   */
  void endTraversal(bool newton3) final {
    if (_postProcessed) {
      throw utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
    }
    if (calculateGlobals) {
#if defined(AUTOPAS_CUDA)
      std::array<SoAFloatPrecision, 4> globals{0, 0, 0, 0};
      _cudaGlobals.copyDeviceToHost(4, globals.data());
      _virialSum[0] += globals[0];
      _virialSum[1] += globals[1];
      _virialSum[2] += globals[2];
      _upotSum += globals[3];
#endif
      for (size_t i = 0; i < _aosThreadData.size(); ++i) {
        _upotSum += _aosThreadData[i].upotSum;
        _virialSum = utils::ArrayMath::add(_virialSum, _aosThreadData[i].virialSum);
      }
      if (not newton3) {
        // if the newton3 optimization is disabled we have added every energy contribution twice, so we divide by 2
        // here.
        _upotSum *= 0.5;
        _virialSum = utils::ArrayMath::mulScalar(_virialSum, 0.5);
      }
      // we have always calculated 6*upot, so we divide by 6 here!
      _upotSum /= 6.;
      _postProcessed = true;
    }
  }

  /**
   * Get the potential Energy.
   * @return the potential Energy
   */
  double getUpot() {
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
   * Get the virial.
   * @return
   */
  double getVirial() {
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

  /**
   * Getter for 24*epsilon.
   * @return 24*epsilon
   */
  SoAFloatPrecision getEpsilon24() const { return _epsilon24; }

  /**
   * Getter for the squared sigma.
   * @return squared sigma.
   */
  SoAFloatPrecision getSigmaSquare() const { return _sigmasquare; }

 private:
  template <bool newton3>
  void SoAFunctorVerletImpl(SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();
    [[maybe_unused]] auto *const __restrict typeptr1 = soa.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr2 = soa.template begin<Particle::AttributeNames::typeId>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    const SoAFloatPrecision cutoffsquare = _cutoffsquare;
    SoAFloatPrecision shift6 = _shift6;
    SoAFloatPrecision sigmasquare = _sigmasquare;
    SoAFloatPrecision epsilon24 = _epsilon24;

    if (ownedStatePtr[indexFirst] == OwnershipState::dummy) {
      return;
    }

    auto func = newton3 ? SoAFunctorVerletISPC_N3 : SoAFunctorVerletISPC_NoN3;
    func(indexFirst, neighborList.data(), neighborList.size(), sigmasquare, epsilon24, cutoffsquare,
                         xptr, yptr, zptr, reinterpret_cast<const int64_t*>(ownedStatePtr), fxptr, fyptr, fzptr);
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
    std::array<double, 3> virialSum;
    double upotSum;

   private:
    // dummy parameter to get the right size (64 bytes)
    double __remainingTo64[(64 - 4 * sizeof(double)) / sizeof(double)];
  };
  // make sure of the size of AoSThreadData
  static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size");

  const double _cutoffsquare;
  // not const because they might be reset through PPL
  double _epsilon24, _sigmasquare, _shift6 = 0;

  ParticlePropertiesLibrary<SoAFloatPrecision, size_t> *_PPLibrary = nullptr;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _upotSum;

  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadData> _aosThreadData;

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

#if defined(AUTOPAS_CUDA)
  using CudaWrapperType = typename std::conditional<calculateGlobals, LJFunctorCudaGlobalsWrapper<SoAFloatPrecision>,
                                                    LJFunctorCudaWrapper<SoAFloatPrecision>>::type;
  // contains wrapper functions for cuda calls
  CudaWrapperType _cudawrapper;

  // contains device globals
  utils::CudaDeviceVector<SoAFloatPrecision> _cudaGlobals;

#endif
};
}  // namespace autopas
