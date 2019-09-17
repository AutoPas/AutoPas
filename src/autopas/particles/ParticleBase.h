/**
 * @file ParticleBase.h
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
 * Minimal definition of a basic particle.
 *
 * If a different Particle class should be used with AutoPas this class must be used as a base to build your own
 * Particle class.
 * @tparam floating point type to be used for the SoAs
 */
template <typename floatType, typename idType>
class ParticleBase {
 public:
  ParticleBase() : _r({0.0, 0.0, 0.0}), _v({0., 0., 0.}), _f({0.0, 0.0, 0.0}), _id(0), _isOwned{true} {}

  /**
   * Constructor of the Particle class.
   * @param r Position of the particle.
   * @param v Velocity of the particle.
   * @param id Id of the particle.
   */
  ParticleBase(std::array<double, 3> r, std::array<double, 3> v, idType id)
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
  const std::array<double, 3> &getF() const { return _f; }

  /**
   * Set the force acting on the particle
   * @param f force
   */
  void setF(const std::array<double, 3> &f) { _f = f; }

  /**
   * Add a partial force to the force acting on the particle
   * @param f partial force to be added
   */
  void addF(const std::array<double, 3> &f) { _f = ArrayMath::add(_f, f); }

  /**
   * Substract a partial force from the force acting on the particle
   * @param f partial force to be substracted
   */
  void subF(const std::array<double, 3> &f) { _f = ArrayMath::sub(_f, f); }

  /**
   * Get the id of the particle
   * @return id
   */
  idType getID() const { return _id; }

  /**
   * Set the id of the particle
   * @param id id
   */
  void setID(idType id) { _id = id; }

  /**
   * Get the position of the particle
   * @return current position
   */
  const std::array<double, 3> &getR() const { return _r; }

  /**
   * Set the position of the particle
   * @param r new position
   */
  void setR(const std::array<double, 3> &r) { _r = r; }

  /**
   * Add a distance vector to the position of the particle
   * @param r vector to be added
   */
  void addR(const std::array<double, 3> &r) { _r = ArrayMath::add(_r, r); }

  /**
   * Get the velocity of the particle
   * @return current velocity
   */
  const std::array<double, 3> &getV() const { return _v; }

  /**
   * Set the velocity of the particle
   * @param v new velocity
   */
  void setV(const std::array<double, 3> &v) { _v = v; }

  /**
   * Add a vector to the current velocity of the particle
   * @param v vector to be added
   */
  void addV(const std::array<double, 3> &v) { _v = ArrayMath::add(_v, v); }

  /**
   * Creates a string containing all data of the particle.
   * @return String representation.
   */
  virtual std::string toString() const {
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

  /**
   * Defines whether the particle is owned by the current AutoPas object (aka (MPI-)process)
   * @return true if the particle is owned by the current AutoPas object, false otherwise
   */
  bool isOwned() const { return _isOwned; }

  /**
   * Set the owned state to the given value
   * @param owned
   */
  void setOwned(bool owned) { _isOwned = owned; }

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int { id, posX, posY, posZ, forceX, forceY, forceZ, owned };

  /**
   * Floating Point Type used for this particle
   */
  typedef floatType ParticleSoAFloatPrecision;

  /**
   * Id Type used for this particle
   */
  typedef idType ParticleIdType;

  /**
   * The type for the soa storage.
   * owned is currently used as a floatType to ease calculations within the functors.
   */
  typedef typename autopas::utils::SoAType<idType, floatType, floatType, floatType, floatType, floatType, floatType,
                                           floatType>::Type SoAArraysType;

  /**
   * Getter, which allows access to an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @return Value of the requested attribute.
   * @note The value of owned is return as floating point number (true = 1.0, false = 0.0).
   */
  template <AttributeNames attribute>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() const {
    switch (attribute) {
      case AttributeNames::id:
        return getID();
      case AttributeNames::posX:
        return getR()[0];
      case AttributeNames::posY:
        return getR()[1];
      case AttributeNames::posZ:
        return getR()[2];
      case AttributeNames::forceX:
        return getF()[0];
      case AttributeNames::forceY:
        return getF()[1];
      case AttributeNames::forceZ:
        return getF()[2];
      case AttributeNames::owned:
        return isOwned() ? 1. : 0.;
      default:
        utils::ExceptionHandler::exception("ParticleBase::get: unknown attribute");
        return 0;
    }
  }

  /**
   * Setter, which allows set an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @param value New value of the requested attribute.
   * @note The value of owned is extracted from a floating point number (true = 1.0, false = 0.0).
   */
  template <AttributeNames attribute>
  constexpr void set(typename std::tuple_element<attribute, SoAArraysType>::type::value_type value) {
    switch (attribute) {
      case AttributeNames::id:
        setID(value);
        break;
      case AttributeNames::posX:
        _r[0] = value;
        break;
      case AttributeNames::posY:
        _r[1] = value;
        break;
      case AttributeNames::posZ:
        _r[2] = value;
        break;
      case AttributeNames::forceX:
        _f[0] = value;
        break;
      case AttributeNames::forceY:
        _f[1] = value;
        break;
      case AttributeNames::forceZ:
        _f[2] = value;
        break;
      case AttributeNames::owned:
        setOwned(value == 1.);
        break;
    }
  }

#if defined(AUTOPAS_CUDA)
  /**
   * The type for storage arrays for Cuda
   */
  typedef typename autopas::utils::CudaSoAType<idType, floatType, floatType, floatType, floatType, floatType, floatType,
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
  std::array<double, 3> _r;
  /**
   * Particle velocity as 3D vector.
   */
  std::array<double, 3> _v;
  /**
   * Force the particle experiences as 3D vector.
   */
  std::array<double, 3> _f;
  /**
   * Particle id.
   */
  idType _id;

  /**
   * Defines whether the particle is owned by the current AutoPas object (aka (MPI-)process)
   */
  bool _isOwned;
};

}  // namespace autopas
