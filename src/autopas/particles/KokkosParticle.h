/**
 * @file KokkosParticle.h
 * @author M. Geitner
 * @date 24.06.2019
 */

#pragma once
#include "autopas/particles/Particle.h"

namespace autopas {

// template <typename floatType>
/**
 * particle class desinged for kokkos functionalities
 */
class KokkosParticle : public Particle {
 private:
#ifdef AUTOPAS_KOKKOS
  FloatVectorType _r, _v, _f;
#endif
  unsigned long _id;

 public:
  KokkosParticle() : Particle() {
#ifdef AUTOPAS_KOKKOS
    _r = FloatVectorType("_r", 3);
    _v = FloatVectorType("_v", 3);
    _f = FloatVectorType("_f", 3);
#endif
  }

  /**
   * Constructor that accepts array values and store them in the data structure FloatVectorType ( a view) from kokkos
   * @param r
   * @param v
   * @param id
   */
  KokkosParticle(std::array<double, KOKKOS_DIM> r, std::array<double, KOKKOS_DIM> v, unsigned long id) : Particle() {
#ifdef AUTOPAS_KOKKOS
    // create views
    _r = FloatVectorType("_r", 3);
    _v = FloatVectorType("_v", 3);
    _f = FloatVectorType("_f", 3);

    // printf("def particle 1\n");
    // copy values from constructor
    FloatVectorType::HostMirror h_r = Kokkos::create_mirror_view(_r);
    FloatVectorType::HostMirror h_v = Kokkos::create_mirror_view(_v);
    // LongVectorType::HostMirror h_id = Kokkos::create_mirror_view(_id);
    // printf("def particle 2\n");
    for (int i = 0; i < 3; i++) {
      h_r(i) = r[i];
      h_v(i) = v[i];
    }
    // FP_float fx = h_r(0);

    _id = id;
    // target, source
    Kokkos::deep_copy(_r, h_r);
    Kokkos::deep_copy(_v, h_v);
#endif
  }

#ifdef AUTOPAS_KOKKOS
  FloatVectorType get_r() const { return _r; }

  FloatVectorType get_v() { return _v; }
  FloatVectorType get_f() { return _f; }

  /**
   * Set the force acting on the particle
   * @param f force
   */
  void setF(FloatVectorType &f) { _f = f; }
  /**
   * Add a partial force to the force acting on the particle
   * @param f partial force to be added
   */

  void addF(FloatVectorType &f) {
    KokkosHelper::addDir(_f, f);
    //_f = ArrayMath::add(_f, f);
  }

  /**
   * Substract a partial force from the force acting on the particle
   * @param f partial force to be substracted
   */

  void subF(FloatVectorType &f) { KokkosHelper::subDir(_f, f); }

  /**
   * Set the position acting on the particle
   * @param r force
   */

  void setR(FloatVectorType &r) { _r = r; }

  /**
   * Add a partial force to the force acting on the particle
   * @param f partial force to be added
   */

  void addR(FloatVectorType &r) { KokkosHelper::addDir(_r, r); }

  /**
   * Substract a partial force from the force acting on the particle
   * @param f partial force to be substracted
   */

  void subR(FloatVectorType &r) { KokkosHelper::subDir(_r, r); }

  /**
   * Set the position acting on the particle
   * @param r force
   */

  void setV(FloatVectorType &v) { _v = v; }

  /**
   * Add a partial force to the force acting on the particle
   * @param f partial force to be added
   */

  void addV(FloatVectorType &v) { KokkosHelper::addDir(_v, v); }

  std::string toString() override {
    std::ostringstream text;
    FloatVectorType::HostMirror h_r = Kokkos::create_mirror_view(_r);
    Kokkos::deep_copy(h_r, _r);
    FloatVectorType::HostMirror h_v = Kokkos::create_mirror_view(_v);
    FloatVectorType::HostMirror h_f = Kokkos::create_mirror_view(_f);
    Kokkos::deep_copy(h_v, _v);
    Kokkos::deep_copy(h_f, _f);
    // printf("reach toString\n");
    // FP_float fx = h_r(0);

    text << "Particle"
         << "\nID      : " << _id << "\nPosition: " << h_r(0) << " | " << h_r(1) << " | " << h_r(2)

         << "\nVelocity: " << h_v[0] << " | " << h_v[1] << " | " << h_v[2] << "\nForce   : " << h_f[0] << " | "
         << h_f[1] << " | " << h_f[2];
    // printf("reach EndtoString\n");
    return text.str();
  };

  KOKKOS_INLINE_FUNCTION
  FloatVectorType get_r_inline() const { return _r; }

  KOKKOS_INLINE_FUNCTION
  FloatVectorType get_f_inline() const { return _f; }
#endif
};

}  // namespace autopas
