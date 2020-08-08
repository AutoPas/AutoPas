/**
 * @file EmptyFunctor.h
 * @author seckler
 * @date 26.03.20
 */

#pragma once

#include "autopas/cells/ParticleCell.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/Functor.h"
#if defined(AUTOPAS_CUDA)
#include "autopas/utils/CudaSoA.h"
#endif

/**
 * Empty Functor, this functor is empty and can be used for testing purposes.
 * It returns that it is applicable for everything.
 */
template <class Particle, class SoAArraysType = typename Particle::SoAArraysType>
class EmptyFunctor : public autopas::Functor<Particle> {
 public:
  /**
   * Default constructor.
   */
  EmptyFunctor() : autopas::Functor<Particle>(0.){};

  void AoSFunctor(Particle &i, Particle &j, bool newton3) override {}

  /**
   * SoAFunctor for a single Particle.
   * @param soa An autopas::SoAView for the Functor
   * @param newton3 A boolean to indicate whether to allow newton3
   */
  void SoAFunctorSingle(autopas::SoAView<typename Particle::SoAArraysType> soa, bool newton3) override {}

  /**
   * SoAFunctor for a pair of Particle.
   * @param soa An autopas::SoAView for the Functor
   * @param soa2 A second autopas::SoAView for the Functor
   * @param newton3 A boolean to indicate whether to allow newton3
   */
  void SoAFunctorPair(autopas::SoAView<typename Particle::SoAArraysType> soa,
                      autopas::SoAView<typename Particle::SoAArraysType> soa2, bool newton3) override {}

  /**
   * @copydoc autopas::Functor::SoAFunctorVerlet()
   */
  void SoAFunctorVerlet(autopas::SoAView<typename Particle::SoAArraysType> soa, size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) override{};

  bool allowsNewton3() override { return true; }

  bool allowsNonNewton3() override { return true; }

  /**
   * @copydoc autopas::Functor::isAppropriateClusterSize()
   */
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
#endif
};
