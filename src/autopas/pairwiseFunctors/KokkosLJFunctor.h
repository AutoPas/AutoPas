/**
 * @file KokkosLJFunctor.h
 * @author M. Geitner
 * @date 24.06.2019
 */

#pragma once

#include "LJFunctor.h"

#ifdef AUTOPAS_KOKKOS
#include <Kokkos_Core.hpp>
#include "autopas/utils/KokkosTypes.h"
#include "autopas/utils/KokkosHelper.h"
#endif

namespace autopas {

/**
 * Class resembles LJFunctor.h
 * Class can process two KokkosParticles
 */
template <class Particle, class ParticleCell, FunctorN3Modes useNewton3 = FunctorN3Modes::Both>

class KokkosLJFunctor : public Functor<Particle, ParticleCell, typename Particle::SoAArraysType>{
    using SoAArraysType = typename Particle::SoAArraysType;
 private:
  float _cutoff_square, _epsilon24, _sigma_square;
  bool _newton3;

 public:
  KokkosLJFunctor() {
    _cutoff_square = 1;
    _sigma_square = 1;
    _epsilon24 = 1 * 24.0;
    _newton3 = (useNewton3 == FunctorN3Modes::Both || useNewton3 == FunctorN3Modes::Newton3Only);
  };

  /**
   * Constructor, which sets the global values, cutoff, epsilon, sigma and newton3.
   * @param cutoff
   * @param epsilon
   * @param sigma
   * @param newton3
   */
  KokkosLJFunctor(double cutoff, double epsilon, double sigma, bool newton3) {
    _cutoff_square = cutoff * cutoff;
    _sigma_square = sigma * sigma;
    _epsilon24 = epsilon * 24.0;
    _newton3 = newton3;
  }

    bool isRelevantForTuning() override {
      return false;
    }

    ~KokkosLJFunctor() override = default;

    void initTraversal() override {

    }

    void endTraversal(bool newton3) override {

    }

    void AoSFunctor(Particle &i, Particle &j, bool newton3) override {

    }

    void SoAFunctor(SoA<SoAArraysType> &soa, bool newton3) override {

    }

    void SoAFunctor(SoA<SoAArraysType> &soa,
                    const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborList,
                    size_t iFrom, size_t iTo, bool newton3) override {

    }

    void SoAFunctor(SoA<SoAArraysType> &soa1, SoA<SoAArraysType> &soa2, bool newton3) override {

    }

    void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle, bool newton3) override {

    }

    void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle1,
                     CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle2, bool newton3) override {

    }

    void deviceSoALoader(::autopas::SoA<SoAArraysType> &soa,
                         CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) override {

    }

    void deviceSoAExtractor(::autopas::SoA<SoAArraysType> &soa,
                            CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) override {

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

    AUTOPAS_FUNCTOR_SOALOADER(
            cell, soa, offset,
    )

    AUTOPAS_FUNCTOR_SOAEXTRACTOR(
            cell, soa, offset,
    )

    bool allowsNewton3() override {
      return false;
    }

    bool allowsNonNewton3() override {
      return true;
    }



#ifdef AUTOPAS_KOKKOS
  KOKKOS_INLINE_FUNCTION
  void AoSFunctorInline(const Particle &i, const Particle &j) const;
#endif
};

#ifdef AUTOPAS_KOKKOS
    template<class Particle, class ParticleCell, FunctorN3Modes useNewton3>
    KOKKOS_INLINE_FUNCTION
    void KokkosLJFunctor<Particle, ParticleCell, useNewton3>::AoSFunctorInline(const Particle &i, const Particle &j) const {
      KOKKOS_FLOAT dr2 = KokkosHelper::subDot(i.get_r_inline(), j.get_r_inline());
      KOKKOS_FLOAT invdr2 = 1. / dr2;
      KOKKOS_FLOAT lj6 = _sigma_square * invdr2;
      lj6 = lj6 * lj6 * lj6;
      KOKKOS_FLOAT lj12 = lj6 * lj6;
      KOKKOS_FLOAT lj12m6 = lj12 - lj6;
      KOKKOS_FLOAT fac = _epsilon24 * (lj12 + lj12m6) * invdr2;

      KokkosHelper::subDotMulScalarAddF(i.get_r_inline(), j.get_r_inline(), i.get_f_inline(), fac);  // to parallel_for

    }

#endif
}  // namespace autopas