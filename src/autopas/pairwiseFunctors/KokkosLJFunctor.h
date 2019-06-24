/**
 * @file KokkosLJFunctor.h
 * @author M. Geitner
 * @date 24.06.2019
 */

#pragma once

#ifdef AUTOPAS_KOKKOS
#include <Kokkos_Core.hpp>
#endif

namespace autopas {

/**
 * Class resembles LJFunctor.h
 * Class can process two KokkosParticles
 */
class KokkosLJFunctor {
 private:
  float _cutoff_square, _epsilon24, _sigma_square;
  bool _newton3;

 public:
  KokkosLJFunctor() {
    _cutoff_square = 1;
    _sigma_square = 1;
    _epsilon24 = 1 * 24.0;
    _newton3 = false;
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
#ifdef AUTOPAS_KOKKOS

  KOKKOS_INLINE_FUNCTION
  void AoSFunctorInline(const KokkosParticle &i, const KokkosParticle &j) const {
    // precondition, newton3 off
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
};

}  // namespace autopas