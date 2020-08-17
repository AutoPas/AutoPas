/**
 * @file EmptyFunctor.h
 * @author seckler
 * @date 26.03.20
 */

#pragma once

#include <glob.h>
#if defined(AUTOPAS_CUDA)
#include "cuda_runtime.h"
#endif

#include "autopas/cells/ParticleCell.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/Functor.h"
#if defined(AUTOPAS_CUDA)
#include "autopas/utils/CudaSoA.h"
#endif

#ifdef AUTOPAS_CUDA
template <typename floatingPointType>
class EmptyCudaWrapper : public autopas::CudaWrapperInterface<floatingPointType> {
  using floatType = floatingPointType;

 public:
  void setNumThreads(int num_threads) override {}
  void loadConstants(autopas::FunctorCudaConstants<floatType> *constants) override {}
  void SoAFunctorNoN3Wrapper(autopas::FunctorCudaSoA<floatType> *cell1, cudaStream_t stream) override {}
  void SoAFunctorNoN3PairWrapper(autopas::FunctorCudaSoA<floatType> *cell1, autopas::FunctorCudaSoA<floatType> *cell2,
                                 cudaStream_t stream) override {}
  void SoAFunctorN3Wrapper(autopas::FunctorCudaSoA<floatType> *cell1, cudaStream_t stream) override {}
  void SoAFunctorN3PairWrapper(autopas::FunctorCudaSoA<floatType> *cell1, autopas::FunctorCudaSoA<floatType> *cell2,
                               cudaStream_t stream) override {}
  void LinkedCellsTraversalNoN3Wrapper(autopas::FunctorCudaSoA<floatType> *cell1, unsigned int reqThreads,
                                       unsigned int cids_size, unsigned int *cids, unsigned int cellSizes_size,
                                       size_t *cellSizes, cudaStream_t stream) override {}
  void LinkedCellsTraversalN3Wrapper(autopas::FunctorCudaSoA<floatType> *cell1, unsigned int reqThreads,
                                     unsigned int cids_size, unsigned int *cids, unsigned int cellSizes_size,
                                     size_t *cellSizes, cudaStream_t stream) override {}
  void CellVerletTraversalNoN3Wrapper(autopas::FunctorCudaSoA<floatType> *cell1Base, unsigned int ncells,
                                      unsigned int clusterSize, unsigned int others_size, unsigned int *other_ids,
                                      cudaStream_t stream) override {}
  void CellVerletTraversalN3Wrapper(autopas::FunctorCudaSoA<floatType> *cell1Base, unsigned int ncells,
                                    unsigned int clusterSize, unsigned int others_size, unsigned int *other_ids,
                                    cudaStream_t stream) override {}
  void loadLinkedCellsOffsets(unsigned int offsets_size, int *offsets) override {}
  bool isAppropriateClusterSize(unsigned int clusterSize) const override { return true; }
};
#endif

/**
 * Empty Functor, this functor is empty and can be used for testing purposes.
 * It returns that it is applicable for everything.
 */
template <class Particle, class ParticleCell_t, class SoAArraysType = typename Particle::SoAArraysType>
class EmptyFunctor : public autopas::Functor<Particle, ParticleCell_t> {
 private:
  EmptyCudaWrapper<typename Particle::ParticleSoAFloatPrecision> emptyCudaWrapper;

 public:
  /**
   * Default constructor.
   */
  EmptyFunctor() : autopas::Functor<Particle, ParticleCell_t>(0.){};

  void AoSFunctor(Particle &i, Particle &j, bool newton3) override {}

  void SoAFunctorSingle(autopas::SoAView<typename Particle::SoAArraysType> soa, bool newton3) override {}

  void SoAFunctorPair(autopas::SoAView<typename Particle::SoAArraysType> soa,
                      autopas::SoAView<typename Particle::SoAArraysType> soa2, bool newton3) override {}

#ifdef AUTOPAS_CUDA
  autopas::CudaWrapperInterface<typename Particle::ParticleSoAFloatPrecision> *getCudaWrapper() override {
    return &emptyCudaWrapper;
  }
#endif

  void SoAFunctorVerlet(autopas::SoAView<typename Particle::SoAArraysType> soa, size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) override{};

  bool allowsNewton3() override { return true; }

  bool allowsNonNewton3() override { return true; }

  bool isAppropriateClusterSize(unsigned int clusterSize, autopas::DataLayoutOption::Value dataLayout) const override {
    return true;
  }

  bool isRelevantForTuning() override { return true; }

#if defined(AUTOPAS_CUDA)
  void CudaFunctor(autopas::CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle, bool newton3) override {}

  void CudaFunctor(autopas::CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle1,
                   autopas::CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle2, bool newton3) override {}

  void deviceSoALoader(autopas::SoA<SoAArraysType> &soa,
                       autopas::CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) override {}

  void deviceSoAExtractor(autopas::SoA<SoAArraysType> &soa,
                          autopas::CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) override {
  }

#endif
};
