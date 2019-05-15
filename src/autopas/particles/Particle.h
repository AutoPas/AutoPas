/**
 * @file Particle.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>
#include <sstream>
#include <tuple>
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/CudaSoAType.h"
#include "autopas/utils/SoAStorage.h"
#include "autopas/utils/SoAType.h"

namespace autopas {

/**
 * basic Particle class
 * This class can be used to build your own Particle class. However, you are
 * free to not use it as well.
 */
template <typename floatType>
class ParticleBase {
 public:
  ParticleBase() : _r({0.0, 0.0, 0.0}), _v({0., 0., 0.}), _f({0.0, 0.0, 0.0}), _id(0), _isOwned{true} {}

  /**
   * Constructor of the Particle class
   * @param r position of the particle
   * @param v velocity of the particle
   * @param id id of the particle
   */
  ParticleBase(std::array<floatType, 3> r, std::array<floatType, 3> v, unsigned long id)
      : _r(r), _v(v), _f({0.0, 0.0, 0.0}), _id(id), _isOwned{true} {}

  /**
   * Destructor of ParticleBase class
   */
  virtual ~ParticleBase() = default;

  /**
   * Equality operator for ParticleBase class.
   * @param rhs
   * @return
   */
  bool operator==(const ParticleBase &rhs) const {
    return std::tie(_r, _v, _f, _id) == std::tie(rhs._r, rhs._v, rhs._f, rhs._id);
  }

  /**
   * Not-Equals operator for ParticleBase class.
   * @param rhs
   * @return
   */
  bool operator!=(const ParticleBase &rhs) const { return not(rhs == *this); }

  /**
   * get the force acting on the particle
   * @return force
   */
  const std::array<floatType, 3> &getF() const { return _f; }

  /**
   * Set the force acting on the particle
   * @param f force
   */
  void setF(const std::array<floatType, 3> &f) { _f = f; }

  /**
   * Add a partial force to the force acting on the particle
   * @param f partial force to be added
   */
  void addF(const std::array<floatType, 3> &f) { _f = ArrayMath::add(_f, f); }

  /**
   * Substract a partial force from the force acting on the particle
   * @param f partial force to be substracted
   */
  void subF(const std::array<floatType, 3> &f) { _f = ArrayMath::sub(_f, f); }

  /**
   * Get the id of the particle
   * @return id
   */
  unsigned long getID() const { return _id; }

  /**
   * Set the id of the particle
   * @param id id
   */
  void setID(unsigned long id) { _id = id; }

  /**
   * Get the position of the particle
   * @return current position
   */
  const std::array<floatType, 3> &getR() const { return _r; }

  /**
   * Set the position of the particle
   * @param r new position
   */
  void setR(const std::array<floatType, 3> &r) { _r = r; }

  /**
   * Add a distance vector to the position of the particle
   * @param r vector to be added
   */
  void addR(const std::array<floatType, 3> &r) { _r = ArrayMath::add(_r, r); }

  /**
   * Get the velocity of the particle
   * @return current velocity
   */
  const std::array<floatType, 3> &getV() const { return _v; }

  /**
   * Set the velocity of the particle
   * @param v new velocity
   */
  void setV(const std::array<floatType, 3> &v) { _v = v; }

  /**
   * Add a vector to the current velocity of the particle
   * @param v vector to be added
   */
  void addV(const std::array<floatType, 3> &v) { _v = ArrayMath::add(_v, v); }

  /**
   * Creates a string containing all data of the particle.
   * @return String representation.
   */
  virtual std::string toString() {
    std::ostringstream text;
    // clang-format off
    text << "Particle"
         << "\nID      : " << _id
         << "\nPosition: "
         << _r[0] << " | " << _r[1] << " | " << _r[2]
         << "\nVelocity: "
         << _v[0] << " | " << _v[1] << " | " << _v[2]
         << "\nForce   : "
         << _f[0] << " | " << _f[1] << " | " << _f[2];
    // clang-format on
    return text.str();
  }

  bool isOwned() { return _isOwned; }

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int { id, posX, posY, posZ, forceX, forceY, forceZ };

  /**
   * Floating Point Type used for this particle
   */
  typedef floatType ParticleFloatingPointType;

  /**
   * the type for the soa storage
   */
  typedef
      typename autopas::utils::SoAType<size_t, floatType, floatType, floatType, floatType, floatType, floatType>::Type
          SoAArraysType;

#if defined(AUTOPAS_CUDA)
  /**
   * The type for storage arrays for Cuda
   */
  typedef typename autopas::utils::CudaSoAType<size_t, floatType, floatType, floatType, floatType, floatType,
                                               floatType>::Type CudaDeviceArraysType;
#else
  /**
   * The type for storage arrays for Cuda
   * empty if compiled without Cuda Support
   */
  typedef typename autopas::utils::CudaSoAType<>::Type CudaDeviceArraysType;
#endif

 protected:
  /**
   * Particle position as 3D coordinates.
   */
  std::array<floatType, 3> _r;
  /**
   * Particle velocity as 3D vector.
   */
  std::array<floatType, 3> _v;
  /**
   * Force the particle experiences as 3D vector.
   */
  std::array<floatType, 3> _f;
  /**
   * Particle id.
   */
  unsigned long _id;

  /**
   * Defines whether the particle is owned by the current AutoPas object (aka (MPI-)process)
   */
  bool _isOwned;
};
/**
 * Particle with all variables in 32 bit precision
 */
typedef ParticleBase<float> ParticleFP32;
/**
 * Particle with all variables in 64 bit precision
 */
typedef ParticleBase<double> ParticleFP64;
/**
 * Alias for Particle with all variables in 64 bit precision
 */
typedef ParticleFP64 Particle;

}  // namespace autopas
