/**
 * @file MockFunctor.h
 * @author seckler
 * @date 18.04.18
 */

#pragma once

#include <gmock/gmock.h>

#include "autopas/cells/ParticleCell.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"
#if defined(AUTOPAS_CUDA)
#include "autopas/utils/CudaSoA.h"
#endif

template <class Particle, class ParticleCell_t>
class MockFunctor : public autopas::Functor<Particle, ParticleCell_t> {
 public:
  MockFunctor() : autopas::Functor<Particle, ParticleCell_t>(0.){};

  // virtual void AoSFunctor(Particle &i, Particle &j, bool newton3)
  MOCK_METHOD(void, AoSFunctor, (Particle & i, Particle &j, bool newton3), (override));

  // virtual void SoAFunctorSingle(SoAView &soa, bool newton3)
  MOCK_METHOD(void, SoAFunctorSingle, (autopas::SoAView<typename Particle::SoAArraysType> soa, bool newton3),
              (override));

  // virtual void SoAFunctorPair(SoAView &soa1, SoAView &soa2, bool newton3)
  MOCK_METHOD(void, SoAFunctorPair,
              (autopas::SoAView<typename Particle::SoAArraysType> soa,
               autopas::SoAView<typename Particle::SoAArraysType> soa2, bool newton3),
              (override));

  // virtual void SoAFunctorVerlet(SoAView &soa, const std::vector, (override)<std::vector<size_t,
  // AlignedAllocator<size_t>>> &neighborList, size_t iFrom, size_t iTo, bool newton3)
  MOCK_METHOD(void, SoAFunctorVerlet,
              (autopas::SoAView<typename Particle::SoAArraysType> soa, size_t indexFirst,
               (const std::vector<size_t, autopas::AlignedAllocator<size_t>> &), bool newton3),
              (override));

  // virtual void SoALoader(ParticleCell &cell, autopas::SoA &soa, size_t
  // offset=0) {}
  // no override for the two-input variant, as it only simulates the one with a default argument!
  MOCK_METHOD(void, SoALoader,
              (autopas::ParticleCell<Particle> & cell, autopas::SoA<typename Particle::SoAArraysType> &soa));
  MOCK_METHOD(void, SoALoader,
              (autopas::ParticleCell<Particle> & cell, autopas::SoA<typename Particle::SoAArraysType> &soa,
               size_t offset),
              (override));

  MOCK_METHOD(void, SoALoaderVerlet,
              (typename autopas::VerletListHelpers<Particle>::VerletListParticleCellType & cell,
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
  // no override for the two-input variant, as it only simulates the one with a default argument!
  MOCK_METHOD(void, SoAExtractor,
              (autopas::ParticleCell<Particle> & cell, autopas::SoA<typename Particle::SoAArraysType> &soa));
  MOCK_METHOD(void, SoAExtractor,
              (autopas::ParticleCell<Particle> & cell, autopas::SoA<typename Particle::SoAArraysType> &soa,
               size_t offset),
              (override));

  MOCK_METHOD(void, SoAExtractorVerlet,
              (typename autopas::VerletListHelpers<Particle>::VerletListParticleCellType & cell,
               autopas::SoA<typename Particle::SoAArraysType> &soa, size_t offset));

  template <typename = std::enable_if_t<not std::is_same<
                typename autopas::VerletListHelpers<Particle>::VerletListParticleCellType, ParticleCell_t>::value>>
  void SoAExtractor(typename autopas::VerletListHelpers<Particle>::VerletListParticleCellType &cell,
                    autopas::SoA<typename Particle::SoAArraysType> &soa, size_t offset = 0) {
    SoAExtractorVerlet(cell, soa, offset);
  }

  // virtual bool allowsNewton3() { return true; }
  MOCK_METHOD(bool, allowsNewton3, (), (override));

  // virtual bool allowsNonNewton3() { return false; }
  MOCK_METHOD(bool, allowsNonNewton3, (), (override));

  MOCK_METHOD(bool, isAppropriateClusterSize, (unsigned int clusterSize, autopas::DataLayoutOption::Value dataLayout),
              (const, override));

  //  bool isRelevantForTuning() { return true; }
  MOCK_METHOD(bool, isRelevantForTuning, (), (override));

#if defined(AUTOPAS_CUDA)
  // virtual void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle, bool newton3)
  MOCK_METHOD(void, CudaFunctor,
              (autopas::CudaSoA<typename Particle::CudaDeviceArraysType> & device_handle, bool newton3), (override));

  /*virtual void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle1,
  CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle2, bool newton3)*/
  MOCK_METHOD(void, CudaFunctor,
              (autopas::CudaSoA<typename Particle::CudaDeviceArraysType> & device_handle1,
               autopas::CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle2, bool newton3),
              (override));

  //  void deviceSoALoader(::autopas::SoA<SoAArraysType> &soa,
  //                       CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle)
  MOCK_METHOD(void, deviceSoALoader,
              (autopas::SoA<typename Particle::SoAArraysType> & soa,
               autopas::CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle),
              (override));
  MOCK_METHOD(void, deviceSoAExtractor,
              (autopas::SoA<typename Particle::SoAArraysType> & soa,
               autopas::CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle),
              (override));
#endif
};