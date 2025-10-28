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

#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoAStorage.h"
#include "autopas/utils/SoAType.h"
#include "autopas/utils/markParticleAsDeleted.h"

namespace autopas {

/**
 * Minimal definition of a basic particle.
 *
 * If a different Particle class should be used with AutoPas this class must be used as a base to build your own
 * Particle class.
 * @tparam floatType Floating point type to be used for the SoAs.
 * @tparam idType Integer type to be used for IDs.
 */
template <typename floatType, typename idType>
class ParticleBase {
 public:
  ParticleBase()
      : _r({0.0, 0.0, 0.0}),
        _v({0., 0., 0.}),
        _f({0.0, 0.0, 0.0}),
        _id(0),
        _ownershipState(OwnershipState::owned)
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
        ,
        _rAtRebuild({0.0, 0.0, 0.0})
#endif
  {
  }

  /**
   * Constructor of the Particle class.
   * @param r Position of the particle.
   * @param v Velocity of the particle.
   * @param id Id of the particle.
   * @param ownershipState OwnershipState of the particle (can be either owned, halo, or dummy)
   */
  ParticleBase(const std::array<double, 3> &r, const std::array<double, 3> &v, idType id,
               OwnershipState ownershipState = OwnershipState::owned)
      : _r(r),
        _v(v),
        _f({0.0, 0.0, 0.0}),
        _id(id),
        _ownershipState(ownershipState)
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
        ,
        _rAtRebuild(r)
#endif
  {
  }

  /**
   * Destructor of ParticleBase class
   */
  virtual ~ParticleBase() = default;

 protected:
  /**
   * Particle position as 3D coordinates.
   */
  std::array<floatType, 3> _r;

#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
  /**
   * Particle position during last rebuild as 3D coordinates.
   */
  std::array<floatType, 3> _rAtRebuild;
#endif

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
  idType _id;

  /**
   * Defines the state of the ownership of the particle.
   */
  OwnershipState _ownershipState;

 public:
  /**
   * Stream operator for instances of ParticleBase class.
   * @return String representation.
   */
  template <typename T, typename P>
  friend std::ostream &operator<<(std::ostream &os, const autopas::ParticleBase<T, P> &D);

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
  [[nodiscard]] const std::array<double, 3> &getF() const { return _f; }

  /**
   * Set the force acting on the particle
   * @param f force
   */
  void setF(const std::array<double, 3> &f) { _f = f; }

  /**
   * Add a partial force to the force acting on the particle
   * @param f partial force to be added
   */
  void addF(const std::array<double, 3> &f) {
    using namespace autopas::utils::ArrayMath::literals;
    _f += f;
  }

  /**
   * Substract a partial force from the force acting on the particle
   * @param f partial force to be substracted
   */
  void subF(const std::array<double, 3> &f) {
    using namespace autopas::utils::ArrayMath::literals;
    _f -= f;
  }

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
  [[nodiscard]] const std::array<double, 3> &getR() const { return _r; }

  /**
   * Set the position of the particle
   * @param r new position
   */
  void setR(const std::array<double, 3> &r) { _r = r; }

#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
  /**
   * Get the last rebuild position of the particle
   * @return current rebuild position
   */
  [[nodiscard]] const std::array<double, 3> &getRAtRebuild() const { return _rAtRebuild; }
  /**
   * Set the rebuild position of the particle
   * @param r rebuild position to be set
   */
  void setRAtRebuild(const std::array<double, 3> &r) { _rAtRebuild = r; }

  /**
   * Update the rebuild position of the particle to current position
   */
  void resetRAtRebuild() { this->setRAtRebuild(_r); }

  /**
   * Calculate the distance since the last rebuild.
   * This is used to check if neighbor lists are still valid inside the logic handler
   * @return displacement vector of particle since rebuild
   */
  const std::array<double, 3> calculateDisplacementSinceRebuild() const {
    return utils::ArrayMath::sub(_rAtRebuild, _r);
  }
#endif

  /**
   * Add a distance vector to the position of the particle and check if the distance between the old and new position
   * is less than a given max distance.
   * This max distance usually should be the skin per timestep divided by two.
   *
   * @param r vector to be added
   * @param maxDistSquared The maximum expected movement distance squared.
   * @return true if dot(r - _r) < skinPerTimestepHalvedSquared
   */
  bool setRDistanceCheck(const std::array<double, 3> &r, double maxDistSquared) {
    using namespace autopas::utils::ArrayMath::literals;
    const auto distanceVec = r - _r;
    const double distanceSquared = utils::ArrayMath::dot(distanceVec, distanceVec);
    setR(r);
    const bool distanceIsFine =
        distanceSquared < maxDistSquared or autopas::utils::Math::isNearAbs(maxDistSquared, 0., 1e-12);
    if (not distanceIsFine) {
      AutoPasLog(WARN, "Particle {}: Distance between old and new position is larger than expected: {} > {}", _id,
                 distanceSquared, maxDistSquared);
    }
    return distanceIsFine;
  }

  /**
   * Add a distance vector to the position of the particle
   * @param r vector to be added
   */
  void addR(const std::array<double, 3> &r) {
    using namespace autopas::utils::ArrayMath::literals;
    _r += r;
  }

  /**
   * Add a distance vector to the position of the particle and check if the distance between the old and new position
   * is less than a given max distance.
   * This max distance usually should be the skin per timestep divided by two.
   *
   * @note uses setRDistanceCheck()
   *
   * @param r vector to be added
   * @param maxDistSquared The maximum expected movement distance squared.
   * @return true if dot(r - _r) < skinPerTimestepHalvedSquared
   */
  bool addRDistanceCheck(const std::array<double, 3> &r, double maxDistSquared) {
    using namespace autopas::utils::ArrayMath::literals;
    const auto newR = _r + r;
    return setRDistanceCheck(newR, maxDistSquared);
  }

  /**
   * Get the velocity of the particle
   * @return current velocity
   */
  [[nodiscard]] const std::array<double, 3> &getV() const { return _v; }

  /**
   * Set the velocity of the particle
   * @param v new velocity
   */
  void setV(const std::array<double, 3> &v) { _v = v; }

  /**
   * Add a vector to the current velocity of the particle
   * @param v vector to be added
   */
  void addV(const std::array<double, 3> &v) {
    using namespace autopas::utils::ArrayMath::literals;
    _v += v;
  }

  /**
   * Creates a string containing all data of the particle.
   * @return String representation.
   */
  [[nodiscard]] virtual std::string toString() const {
    std::ostringstream text;
    // clang-format off
    text << "Particle"
         << "\nID      : " << _id
         << "\nPosition: "
         << autopas::utils::ArrayUtils::to_string(_r)
         << "\nVelocity: "
         << autopas::utils::ArrayUtils::to_string(_v)
         << "\nForce   : "
         << autopas::utils::ArrayUtils::to_string(_f)
         << "\nOwnershipState : "
         << _ownershipState;
    // clang-format on
    return text.str();
  }
  /**
   * Defines whether the particle is owned by the current AutoPas object (aka (MPI-)process)
   * @return true if the particle is owned by the current AutoPas object, false otherwise
   */
  [[nodiscard]] bool isOwned() const { return _ownershipState == OwnershipState::owned; }

  /**
   * Defines whether the particle is a halo particle, i.e., not owned by the current AutoPas object (aka (MPI-)process)
   * @return true if the particle is not owned by the current AutoPas object, false otherwise.
   * @note when a
   */
  [[nodiscard]] bool isHalo() const { return _ownershipState == OwnershipState::halo; }

  /**
   * Returns whether the particle is a dummy particle.
   * @return true if the particle is a dummy.
   */
  [[nodiscard]] bool isDummy() const { return _ownershipState == OwnershipState::dummy; }

  /**
   * Returns the particle's ownership state.
   * @return the current OwnershipState
   */
  [[nodiscard]] OwnershipState getOwnershipState() const { return _ownershipState; }

  /**
   * Set the OwnershipState to the given value
   * @param ownershipState
   */
  void setOwnershipState(OwnershipState ownershipState) { _ownershipState = ownershipState; }

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int { ptr, id, posX, posY, posZ, forceX, forceY, forceZ, ownershipState };

  /**
   * Floating Point Type used for this particle
   */
  using ParticleSoAFloatPrecision = floatType;

  /**
   * Id Type used for this particle
   */
  using ParticleIdType = idType;

  /**
   * The type for the soa storage.
   * owned is currently used as a floatType to ease calculations within the functors.
   */
  using SoAArraysType =
      typename autopas::utils::SoAType<ParticleBase<floatType, idType> *, idType /*id*/, floatType /*x*/,
                                       floatType /*y*/, floatType /*z*/, floatType /*fx*/, floatType /*fy*/,
                                       floatType /*fz*/, OwnershipState /*ownershipState*/>::Type;

  /**
   * Non-const getter for the pointer of this object.
   * @tparam attribute Attribute name.
   * @return this.
   */
  template <AttributeNames attribute, std::enable_if_t<attribute == AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() {
    return this;
  }

  /**
   * Getter, which allows access to an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @return Value of the requested attribute.
   * @note The value of owned is return as floating point number (true = 1.0, false = 0.0).
   */
  template <AttributeNames attribute, std::enable_if_t<attribute != AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() const {
    if constexpr (attribute == AttributeNames::id) {
      return getID();
    } else if constexpr (attribute == AttributeNames::posX) {
      return getR()[0];
    } else if constexpr (attribute == AttributeNames::posY) {
      return getR()[1];
    } else if constexpr (attribute == AttributeNames::posZ) {
      return getR()[2];
    } else if constexpr (attribute == AttributeNames::forceX) {
      return getF()[0];
    } else if constexpr (attribute == AttributeNames::forceY) {
      return getF()[1];
    } else if constexpr (attribute == AttributeNames::forceZ) {
      return getF()[2];
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      utils::ExceptionHandler::exception("ParticleBase::get() unknown attribute {}", attribute);
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
    if constexpr (attribute == AttributeNames::id) {
      setID(value);
    } else if constexpr (attribute == AttributeNames::posX) {
      _r[0] = value;
    } else if constexpr (attribute == AttributeNames::posY) {
      _r[1] = value;
    } else if constexpr (attribute == AttributeNames::posZ) {
      _r[2] = value;
    } else if constexpr (attribute == AttributeNames::forceX) {
      _f[0] = value;
    } else if constexpr (attribute == AttributeNames::forceY) {
      _f[1] = value;
    } else if constexpr (attribute == AttributeNames::forceZ) {
      _f[2] = value;
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else {
      utils::ExceptionHandler::exception("MoleculeLJ::set() unknown attribute {}", attribute);
    }
  }

 private:
  /**
   * Marks a particle as deleted.
   * @note: This function should not be used from outside of AutoPas,
   * because the logic handler keeps track of the number of particles. Instead, use AutoPas::deleteParticle(iterator).
   * @note: From within autopas, you might want to use internal::markParticleAsDeleted(Particle &particle)
   */
  void markAsDeleted() {
    // Set ownership as dummy.
    setOwnershipState(OwnershipState::dummy);
  }

  /**
   * Function to access hidden particle.markAsDeleted() to mark it as internal.
   * @tparam ParticleIterator
   */
  template <class T>
  friend void internal::markParticleAsDeleted(T &);
};

/**
 * Stream operator for instances of ParticleBase class.
 * This function enables passing ParticleBase objects to an ostream via `<<`
 * @tparam Floating point type to be used for the SoAs.
 * @param os
 * @param particle
 * @return String representation.
 */
template <typename floatType, typename idType>
std::ostream &operator<<(std::ostream &os, const ParticleBase<floatType, idType> &particle) {
  using utils::ArrayUtils::operator<<;
  os << "Particle"
     << "\nID      : " << particle._id << "\nPosition: " << particle._r << "\nVelocity: " << particle._v
     << "\nForce   : " << particle._f << "\nOwnershipState : " << particle._ownershipState;
  // clang-format on
  return os;
}

}  // namespace autopas
