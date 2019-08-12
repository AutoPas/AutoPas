/**
 *@file KokkosHelper.h
 *@author M. Geitner
 *@date 24.06.19
 *
 */

#pragma once

#ifdef AUTOPAS_KOKKOS
#include <Kokkos_Core.hpp>
#include "KokkosTypes.h"
#endif

namespace autopas {
class KokkosHelper {
 public:
#ifdef AUTOPAS_KOKKOS
  static void addDir(FloatVectorType &target, FloatVectorType &source) {
    Kokkos::parallel_for(
        KOKKOS_DIM, KOKKOS_LAMBDA(const int i) { target(i) += source(i); });
  }

  static void subDir(FloatVectorType &target, FloatVectorType &source) {
    Kokkos::parallel_for(
        KOKKOS_DIM, KOKKOS_LAMBDA(const int i) { target(i) -= source(i); });
  }

  static FloatVectorType add(FloatVectorType &target, FloatVectorType &source) {
    FloatVectorType c("c", KOKKOS_DIM);
    Kokkos::parallel_for(
        KOKKOS_DIM, KOKKOS_LAMBDA(const int i) { c(i) = target(i) + source(i); });
    return c;
  }

  static FloatVectorType sub(FloatVectorType const &target, FloatVectorType const &source) {
    // printf("KokkosHelper::sub(): Start\n");
    FloatVectorType c("c", KOKKOS_DIM);
    Kokkos::parallel_for(
        KOKKOS_DIM, KOKKOS_LAMBDA(const int i) { c(i) = target(i) - source(i); });
    // printf("KokkosHelper::sub(): End\n");
    return c;
  }



  static FloatVectorType mul(FloatVectorType const &a, FloatVectorType const &b) {
    FloatVectorType c("c", KOKKOS_DIM);
    Kokkos::parallel_for(
        KOKKOS_DIM, KOKKOS_LAMBDA(const int i) { c(i) = a(i) * b(i); });
    return c;
  }

  static FloatVectorType addScalar(FloatVectorType const &a, KOKKOS_FLOAT const &s) {
    FloatVectorType c("c", KOKKOS_DIM);
    Kokkos::parallel_for(
        KOKKOS_DIM, KOKKOS_LAMBDA(const int i) { c(i) = a(i) + s; });
    return c;
  }

  static FloatVectorType mulScalar(FloatVectorType const &a, KOKKOS_FLOAT const &s) {
    FloatVectorType c("c", KOKKOS_DIM);
    Kokkos::parallel_for(
        KOKKOS_DIM, KOKKOS_LAMBDA(const int i) { c(i) = a(i) * s; });
    //KokkosHelper::print(c);
    return c;
  }

  static KOKKOS_FLOAT dot(FloatVectorType const &a, FloatVectorType const &b) {
    KOKKOS_FLOAT res = 0;
    Kokkos::parallel_reduce(
        KOKKOS_DIM, KOKKOS_LAMBDA(int i, double &update) { update += a(i) * b(i); }, res);
    return res;
  }

  static void print(FloatVectorType const &a) {
    std::ostringstream text;
    FloatVectorType::HostMirror h_a = Kokkos::create_mirror_view(a);
    Kokkos::deep_copy(h_a, a);
    std::cout << "KokkosHelper::print(): " << h_a(0) << " | " << h_a(1) << " | " << h_a(2) << "\n";
  }

  // TODO: build two seperate functions
  KOKKOS_INLINE_FUNCTION
  static KOKKOS_FLOAT subDot(FloatVectorType const &target, FloatVectorType const &source) {
    KOKKOS_FLOAT res = 0.0;
    // Kokkos::parallel_for(KOKKOS_DIM, KOKKOS_LAMBDA(const int i){
    KOKKOS_FLOAT temp = 0.0;
    for (int i = 0; i < KOKKOS_DIM; i++) {
      temp = target(i) - source(i);
      res += (temp * temp);
    }

    //});
    // printf("KokkosHelper::sub(): End\n");
    return res;
  }

  /**
   *
   * @param target
   * @param source
   * @param f       modify force here
   * @param s       constant for calculation
   */
  KOKKOS_INLINE_FUNCTION
  static void subDotMulScalarAddF(FloatVectorType const &target, FloatVectorType const &source,
                                  FloatVectorType const &f, KOKKOS_FLOAT const &s) {
    // Kokkos::parallel_for(KOKKOS_DIM, KOKKOS_LAMBDA(const int i){
    KOKKOS_FLOAT temp;
    for (int i = 0; i < KOKKOS_DIM; i++) {
      temp = (target(i) - source(i)) * s;
      f(i) += temp;
    }
  }

    KOKKOS_INLINE_FUNCTION
    static void subDotMulScalarModifyF(FloatVectorType const &target, FloatVectorType const &source,
                                    FloatVectorType const &f_add, FloatVectorType const &f_sub, KOKKOS_FLOAT const &s) {
        // Kokkos::parallel_for(KOKKOS_DIM, KOKKOS_LAMBDA(const int i){
        KOKKOS_FLOAT temp;
        for (int i = 0; i < KOKKOS_DIM; i++) {
            temp = (target(i) - source(i)) * s;
            f_add(i) += temp;
            f_sub(i) -= temp;
        }
    }
#endif
};

}  // namespace autopas
