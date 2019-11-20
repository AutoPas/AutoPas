/**
 * @file MockFunctor.h
 * @author seckler
 * @date 18.04.18
 */

#pragma once

#include <gmock/gmock.h>
#include "autopas/autopasIncludes.h"
#include "autopas/cells/ParticleCell.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#if defined(AUTOPAS_CUDA)
#include "autopas/utils/CudaSoA.h"
#endif

// gmock does not write overrides, so we suppress that warning here!
#if __GNUC__ >= 5
// Disable GCC 5's -Wsuggest-override warnings in gtest
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-override"
#endif

template <class Particle, class ParticleCell_t>
class MockFunctor : public autopas::Functor<Particle, ParticleCell_t> {
 public:
  MockFunctor() : autopas::Functor<Particle, ParticleCell_t>(0.){};

  // virtual void AoSFunctor(Particle &i, Particle &j, bool newton3)
  MOCK_METHOD3_T(AoSFunctor, void(Particle &i, Particle &j, bool newton3));

  // virtual void SoAFunctor(SoAView &soa, bool newton3)
  MOCK_METHOD2_T(SoAFunctor, void(autopas::SoAView<typename Particle::SoAArraysType> soa, bool newton3));

  // virtual void SoAFunctor(SoAView &soa1, SoAView &soa2, bool newton3)
  MOCK_METHOD3_T(SoAFunctor, void(autopas::SoAView<typename Particle::SoAArraysType> soa,
                                  autopas::SoAView<typename Particle::SoAArraysType> soa2, bool newton3));

  // virtual void SoAFunctor(SoAView &soa, const std::vector<std::vector<size_t,
  // AlignedAllocator<size_t>>> &neighborList, size_t iFrom, size_t iTo, bool newton3)
  MOCK_METHOD5_T(SoAFunctor, void(autopas::SoAView<typename Particle::SoAArraysType> soa,
                                  const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &, size_t,
                                  size_t, bool newton3));

  // virtual void SoALoader(ParticleCell &cell, autopas::SoA &soa, size_t
  // offset=0) {}
  MOCK_METHOD2_T(SoALoader,
                 void(autopas::ParticleCell<Particle> &cell, autopas::SoA<typename Particle::SoAArraysType> &soa));
  MOCK_METHOD3_T(SoALoader, void(autopas::ParticleCell<Particle> &cell,
                                 autopas::SoA<typename Particle::SoAArraysType> &soa, size_t offset));

  MOCK_METHOD3_T(SoALoaderVerlet, void(typename autopas::VerletListHelpers<Particle>::VerletListParticleCellType &cell,
                                       autopas::SoA<typename Particle::SoAArraysType> &soa, size_t offset));

  template <typename /*dummy*/ = void,
            typename = std::enable_if_t<not std::is_same<
                typename autopas::VerletListHelpers<Particle>::VerletListParticleCellType, ParticleCell_t>::value>>
  void SoALoader(typename autopas::VerletListHelpers<Particle>::VerletListParticleCellType &cell,
                 autopas::SoA<typename Particle::SoAArraysType> &soa, size_t offset = 0) {
    SoALoaderVerlet(cell, soa, offset);
  }

  // virtual void SoAExtractor(ParticleCell &cell, autopas::SoA &soa, size_t
  // offset=0) {}
  MOCK_METHOD2_T(SoAExtractor,
                 void(autopas::ParticleCell<Particle> &cell, autopas::SoA<typename Particle::SoAArraysType> &soa));
  MOCK_METHOD3_T(SoAExtractor, void(autopas::ParticleCell<Particle> &cell,
                                    autopas::SoA<typename Particle::SoAArraysType> &soa, size_t offset));

  MOCK_METHOD3_T(SoAExtractorVerlet,
                 void(typename autopas::VerletListHelpers<Particle>::VerletListParticleCellType &cell,
                      autopas::SoA<typename Particle::SoAArraysType> &soa, size_t offset));

  template <typename = std::enable_if_t<not std::is_same<
                typename autopas::VerletListHelpers<Particle>::VerletListParticleCellType, ParticleCell_t>::value>>
  void SoAExtractor(typename autopas::VerletListHelpers<Particle>::VerletListParticleCellType &cell,
                    autopas::SoA<typename Particle::SoAArraysType> &soa, size_t offset = 0) {
    SoAExtractorVerlet(cell, soa, offset);
  }

  // virtual bool allowsNewton3() { return true; }
  MOCK_METHOD0(allowsNewton3, bool());

  // virtual bool allowsNonNewton3() { return false; }
  MOCK_METHOD0(allowsNonNewton3, bool());

  //  bool isRelevantForTuning() { return true; }
  MOCK_METHOD0(isRelevantForTuning, bool());

#if defined(AUTOPAS_CUDA)
  // virtual void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle, bool newton3)
  MOCK_METHOD2_T(CudaFunctor,
                 void(autopas::CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle, bool newton3));

  /*virtual void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle1,
  CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle2, bool newton3)*/
  MOCK_METHOD3_T(CudaFunctor,
                 void(autopas::CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle1,
                      autopas::CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle2, bool newton3));

  //  void deviceSoALoader(::autopas::SoA<SoAArraysType> &soa,
  //                       CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle)
  MOCK_METHOD2_T(deviceSoALoader, void(autopas::SoA<typename Particle::SoAArraysType> &soa,
                                       autopas::CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle));
  MOCK_METHOD2_T(deviceSoAExtractor, void(autopas::SoA<typename Particle::SoAArraysType> &soa,
                                          autopas::CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle));
#endif
};

#if __GNUC__ >= 5
#pragma GCC diagnostic pop
#endif
