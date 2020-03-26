/**
 * @file EmptyFunctor.h
 * @author seckler
 * @date 26.03.20
 */

#pragma once

#include "autopas/cells/ParticleCell.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"
#if defined(AUTOPAS_CUDA)
#include "autopas/utils/CudaSoA.h"
#endif

template <class Particle, class ParticleCell_t>
class EmptyFunctor : public autopas::Functor<Particle, ParticleCell_t> {
 public:
  EmptyFunctor() : autopas::Functor<Particle, ParticleCell_t>(0.){};

  void AoSFunctor(Particle &i, Particle &j, bool newton3) override {}

  void SoAFunctor(autopas::SoAView<typename Particle::SoAArraysType> soa, bool newton3) override {}

  void SoAFunctor(autopas::SoAView<typename Particle::SoAArraysType> soa,
                  autopas::SoAView<typename Particle::SoAArraysType> soa2, bool newton3) override {}

  void SoAFunctor(autopas::SoAView<typename Particle::SoAArraysType> soa, size_t indexFirst,
                  const std::vector<size_t, autopas::AlignedAllocator<size_t>> &, bool newton3) override{};

  bool allowsNewton3() override { return true; }

  bool allowsNonNewton3() override { return true; }

  bool isAppropriateClusterSize(unsigned int clusterSize, autopas::DataLayoutOption::Value dataLayout) const override {
    return true;
  }

  bool isRelevantForTuning() override { return true; }

#if defined(AUTOPAS_CUDA)
  void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle, bool newton3) override {}

  void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle1,
                   CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle2, bool newton3) override {}

  void deviceSoALoader(::autopas::SoA<SoAArraysType> &soa,
                       CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) override {}
#endif
};