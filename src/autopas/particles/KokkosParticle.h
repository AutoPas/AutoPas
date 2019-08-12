/**
 * @file KokkosParticle.h
 * @author M. Geitner
 * @date 24.06.2019
 */

#pragma once
#include "autopas/particles/Particle.h"
#ifdef AUTOPAS_KOKKOS
#include "Kokkos_Core.hpp"
#include "autopas/utils/KokkosTypes.h"
#include "autopas/utils/KokkosHelper.h"
#endif

namespace autopas {


/**
 * particle class desinged for kokkos functionalities
 */
//template <typename floatType>
typedef double floatType;

class KokkosParticle : public Particle {
 private:
#ifdef AUTOPAS_KOKKOS
  FloatVectorType _r, _v, _f;

#endif
  //unsigned long _id;
  std::array<double, 3> _r_arr;

 public:
  KokkosParticle() : Particle() ,_r_arr({0., 0., 0.}){
#ifdef AUTOPAS_KOKKOS
    _r = FloatVectorType("_r", 3);
    _v = FloatVectorType("_v", 3);
    _f = FloatVectorType("_f", 3);

#endif
  }

  /**
   * Constructor that accepts array values and store them in the data structure FloatVectorType ( a view) from kokkos
   * @param r   position of the particle
   * @param v   velocity of the particle
   * @param id  id of the particle
   */
  KokkosParticle(std::array<double, KOKKOS_DIM> r, std::array<double, KOKKOS_DIM> v, unsigned long id) : Particle(),_r_arr({0., 0., 0.}) {
#ifdef AUTOPAS_KOKKOS
    // create views
    _r = FloatVectorType("_r", 3);
    _v = FloatVectorType("_v", 3);
    _f = FloatVectorType("_f", 3);
    // copy values from constructor
    //FloatVectorType::HostMirror h_r = Kokkos::create_mirror_view(_r);
    FloatVectorType::HostMirror h_v = Kokkos::create_mirror_view(_v);
    for (int i = 0; i < 3; i++) {
      //h_r(i) = r[i];
      h_v(i) = v[i];
    }
    // target, source
    //Kokkos::deep_copy(_r, h_r);
    Kokkos::deep_copy(_v, h_v);
    setR(r);
#endif
    _id = id;
  }

  ~KokkosParticle(){
     //std::cout << "detroy object\n";
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
    //copy values to host space and generate output
    FloatVectorType::HostMirror h_r = Kokkos::create_mirror_view(_r);
    Kokkos::deep_copy(h_r, _r);
    FloatVectorType::HostMirror h_v = Kokkos::create_mirror_view(_v);
    FloatVectorType::HostMirror h_f = Kokkos::create_mirror_view(_f);
    Kokkos::deep_copy(h_v, _v);
    Kokkos::deep_copy(h_f, _f);
    // printf("reach toString\n");

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

  /**
   *
   * @param r set new position of particle
   */
  void setR(const std::array<floatType, 3> &r) override{
    //std::cout<< r[0] << ", " << r[1] << r[2] << "\n";
    FloatVectorType::HostMirror h_r = Kokkos::create_mirror_view(_r);
      _r_arr = r;//save values in a normal array which is used in getR, purpose: autopas needs access to this data and values from the kokkos_view could not be copied in getR() (const qualifier)
    for(unsigned int i = 0; i < 3; i++){
      h_r(i) = _r_arr[i];
    }

    Kokkos::deep_copy(_r, h_r);
    //std::cout << toString() << "\n";
  }

  /**
   * copy data from kokkos view to an array
   * required for placement to correct cell of the
   * @return position of kokkos particle
   */
    const std::array<double, 3> &getR() const override {
    return _r_arr;
    }


#endif

/**
 * type of SoA storage, used for compatibility reasons
 */
typedef autopas::utils::SoAType<size_t, double, double, double, double, double, double>::Type SoAArraysType;

};

}  // namespace autopas
