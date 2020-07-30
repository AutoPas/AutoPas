/**
 * @file MockFunctor.h
 * @author seckler
 * @date 18.04.18
 */

#pragma once

#include <gmock/gmock.h>
#include <testingHelpers/NonConstructibleParticle.h>

#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"

#if defined(AUTOPAS_CUDA)
#include "autopas/utils/CudaSoA.h"
#endif

template <class Particle>
class MockFunctor : public autopas::Functor<Particle> {
 public:
  MockFunctor() : autopas::Functor<Particle>(0.){};
  // TODO templatisierte mock methoden?
  // TODO Mock method should take template - geht vlt nit -> copy and paste for every combination of ParticleCell and Particle

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

  // void SoALoader<typename ParticleCell>(ParticleCell &cell, autopas::SoA &soa, size_t offset) {}
//  MOCK_METHOD(void, SoALoader,
//              (ParticleCell & cell, autopas::SoA<typename Particle::SoAArraysType> &soa,
//               size_t offset));

    MOCK_METHOD(void, SoALoader,
            ((autopas::FullParticleCell<autopas::ParticleBase<double, unsigned long>>) & cell, (autopas::SoA<autopas::ParticleBase<double, unsigned long>::SoAArraysType>) &soa,
             size_t offset));

//    MOCK_METHOD(void, SoALoader,
//                ((autopas::ReferenceParticleCell<autopas::ParticleBase<double, unsigned long>>) & cell, autopas::SoA<typename Particle::SoAArraysType> &soa,
//                        size_t offset));

    MOCK_METHOD(void, SoALoader,
                (autopas::FullParticleCell<NonConstructibleParticle> & cell, autopas::SoA<NonConstructibleParticle::SoAArraysType> &soa,
                        size_t offset));

//    MOCK_METHOD(void, SoALoader,
//                (autopas::ReferenceParticleCell<NonConstructibleParticle> & cell, autopas::SoA<typename Particle::SoAArraysType> &soa,
//                        size_t offset));


  MOCK_METHOD(void, SoALoaderVerlet,
              ((typename autopas::VerletListHelpers<Particle>::VerletListParticleCellType) & cell,
               autopas::SoA<typename Particle::SoAArraysType> &soa, size_t offset));

  template <typename /*dummy*/ = void>
  void SoALoader(typename autopas::VerletListHelpers<Particle>::VerletListParticleCellType &cell,
                 autopas::SoA<typename Particle::SoAArraysType> &soa, size_t offset) {
    SoALoaderVerlet(cell, soa, offset);
  }

    // TODO templated function mock gmocks
    // TODO try template
    // TODO google
    // TODO stack overflow
    // TODO copy und paste ugly solution
  MOCK_METHOD(void, SoAExtractor,
              ((autopas::FullParticleCell<autopas::ParticleBase<double, unsigned long>>) & cell, autopas::SoA<typename Particle::SoAArraysType> &soa,
               size_t offset));

//    MOCK_METHOD(void, SoAExtractor,
//                ((autopas::ReferenceParticleCell<autopas::ParticleBase<double, unsigned long>>) & cell, autopas::SoA<typename Particle::SoAArraysType> &soa,
//                        size_t offset));

  MOCK_METHOD(void, SoAExtractor,
              (autopas::FullParticleCell<NonConstructibleParticle> & cell, autopas::SoA<typename Particle::SoAArraysType> &soa,
               size_t offset));

//    MOCK_METHOD(void, SoAExtractor,
//                (autopas::ReferenceParticleCell<NonConstructibleParticle> & cell, autopas::SoA<typename Particle::SoAArraysType> &soa,
//                        size_t offset));

  MOCK_METHOD(void, SoAExtractorVerlet,
              ((typename autopas::VerletListHelpers<Particle>::VerletListParticleCellType) & cell,
               autopas::SoA<typename Particle::SoAArraysType> &soa, size_t offset));

  void SoAExtractor(typename autopas::VerletListHelpers<Particle>::VerletListParticleCellType &cell,
                    autopas::SoA<typename Particle::SoAArraysType> &soa, size_t offset) {
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
