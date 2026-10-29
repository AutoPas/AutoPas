/**
 * @file KokkosParticleBase.h
 * @date 28.10.2026
 * @author Luis Gall
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS

#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoAType.h"
#include "autopas/utils/markParticleAsDeleted.h"
#include "autopas/utilsKokkos/KokkosSoA.h"

#include "Kokkos_Core.hpp"

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
class KokkosParticleBase {
 public:
  KokkosParticleBase()
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
  KokkosParticleBase(const std::array<double, 3> &r, const std::array<double, 3> &v, idType id,
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
  virtual ~KokkosParticleBase() = default;

 protected:
  /**
   * Particle position as 3D coordinates.
   */
  Kokkos::Array<floatType, 3> _r;

#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
  /**
   * Particle position during last rebuild as 3D coordinates.
   */
  Kokkos::Array<floatType, 3> _rAtRebuild;
#endif

  /**
   * Particle velocity as 3D vector.
   */
  Kokkos::Array<floatType, 3> _v;

  /**
   * Force the particle experiences as 3D vector.
   */
  Kokkos::Array<floatType, 3> _f;

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
  friend std::ostream &operator<<(std::ostream &os, const KokkosParticleBase<T, P> &D);

  /**
   * get the force acting on the particle
   * @return force
   */
  [[nodiscard]] std::array<floatType, 3> getF() const { return {_f[0], _f[1], _f[2]}; }

  /**
   * Set the force acting on the particle
   * @param f force
   */
  void setF(const std::array<floatType, 3> &f) {
    _f[0] = f[0];
    _f[1] = f[1];
    _f[2] = f[2];
  }

  /**
   * Add a partial force to the force acting on the particle
   * @param f partial force to be added
   */
  void addF(const Kokkos::Array<floatType, 3> &f) {
    _f[0] += f[0];
    _f[1] += f[1];
    _f[2] += f[2];
  }

  /**
   * Substract a partial force from the force acting on the particle
   * @param f partial force to be substracted
   */
  void subF(const Kokkos::Array<floatType, 3> &f) {
    _f[0] -= f[0];
    _f[1] -= f[1];
    _f[2] -= f[2];
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
  [[nodiscard]] std::array<floatType, 3> getR() const { return {_r[0], _r[1], _r[2]}; }

  /**
   * Set the position of the particle
   * @param r new position
   */
  void setR(const std::array<floatType, 3> &r) {
    _r[0] = r[0];
    _r[1] = r[1];
    _r[2] = r[2];
  }

#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
  /**
   * Get the last rebuild position of the particle
   * @return current rebuild position
   */
  [[nodiscard]] std::array<floatType, 3> getRAtRebuild() const {
    return {_rAtRebuild[0], _rAtRebuild[1], _rAtRebuild[2]};
  }
  /**
   * Set the rebuild position of the particle
   * @param r rebuild position to be set
   */
  void setRAtRebuild(const std::array<floatType, 3> &r) {
    _rAtRebuild[0] = r[0];
    _rAtRebuild[1] = r[1];
    _rAtRebuild[2] = r[2];
  }

  /**
   * Update the rebuild position of the particle to current position
   */
  void resetRAtRebuild() { this->setRAtRebuild({_r[0], _r[1], _r[2]}); }

  /**
   * Calculate the distance since the last rebuild.
   * This is used to check if neighbor lists are still valid inside the logic handler
   * @return displacement vector of particle since rebuild
   */
  const std::array<floatType, 3> calculateDisplacementSinceRebuild() const {

    return {_rAtRebuild[0]-_r[0], _rAtRebuild[1]-_r[1], _rAtRebuild[2]-_r[2]};
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
  bool setRDistanceCheck(const std::array<floatType, 3> &r, double maxDistSquared) {
    using namespace autopas::utils::ArrayMath::literals;
    const std::array<floatType, 3> distanceVec = {r[0] - _r[0], r[1] - _r[1], r[2] - _r[2]};
    const floatType distanceSquared = utils::ArrayMath::dot(distanceVec, distanceVec);
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
  void addR(const std::array<floatType, 3> &r) {
    _r[0] += r[0];
    _r[1] += r[1];
    _r[2] += r[2];
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
  bool addRDistanceCheck(const std::array<floatType, 3> &r, double maxDistSquared) {
    const auto newR = {_r[0] + r[0], _r[1] + r[1], _r[2] + r[2]};
    return setRDistanceCheck(newR, maxDistSquared);
  }

  /**
   * Get the velocity of the particle
   * @return current velocity
   */
  [[nodiscard]] std::array<floatType, 3> getV() const { return {_v[0], _v[1], _v[2]}; }

  /**
   * Set the velocity of the particle
   * @param v new velocity
   */
  void setV(const std::array<floatType, 3> &v) {
    _v[0] = v[0];
    _v[1] = v[1];
    _v[2] = v[2];
  }

  /**
   * Add a vector to the current velocity of the particle
   * @param v vector to be added
   */
  void addV(const std::array<floatType, 3> &v) {
    _v[0] += v[0];
    _v[1] += v[1];
    _v[2] += v[2];
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
         // todo
         << "\nVelocity: "
         // todo
         << "\nForce   : "
         // todo
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
  enum AttributeNames : int { ptr, id, posX, posY, posZ,
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    rebuildX, rebuildY, rebuildZ,
#endif
    velX, velY, velZ, forceX, forceY, forceZ, ownershipState };

  /**
   * Floating Point Type used for this particle
   */
  using ParticleSoAFloatPrecision = floatType;

  /**
   * Id Type used for this particle
   */
  using ParticleIdType = idType;

  /**
   * The type for the host-side SoA storage. Its column order matches AttributeNames (including ptr) and it determines
   * the types of get() and set().
   */
  using SoAArraysType = typename autopas::utils::SoAType<
      KokkosParticleBase<floatType, idType> * /*ptr*/, idType /*id*/, floatType /*x*/, floatType /*y*/,
      floatType /*z*/,
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
      floatType /*rebuildX*/, floatType /*rebuildY*/, floatType /*rebuildZ*/,
#endif
      floatType /*vx*/, floatType /*vy*/, floatType /*vz*/, floatType /*fx*/, floatType /*fy*/, floatType /*fz*/,
      OwnershipState /*ownershipState*/>::Type;

#ifdef AUTOPAS_ENABLE_KOKKOS
  /**
   * The type for the Kokkos SoA storage. Same columns as SoAArraysType but without ptr.
   */
  using KokkosSoAArraysType = autopas::utilsKokkos::KokkosSoA<
    idType * /*id*/, floatType * /*x*/, floatType * /*y*/,
    floatType * /*z*/,
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    floatType * /*rebuildX*/, floatType * /*rebuildY*/,
    floatType * /*rebuildZ*/,
#endif
    floatType * /*vx*/, floatType * /*vy*/, floatType * /*vz*/,
    floatType * /*fx*/, floatType * /*fy*/, floatType * /*fz*/,
    autopas::OwnershipState * /*ownershipState*/>;

  static_assert(std::tuple_size_v<SoAArraysType> == std::tuple_size_v<typename KokkosSoAArraysType::Type> + 1,
                "KokkosSoAArraysType must have the same columns as SoAArraysType, except for ptr.");
#endif

  /**
   * Non-const accessor used by KokkosAoS, which allows modifying an attribute in place.
   * @tparam attribute Attribute name.
   * @return this for AttributeNames::ptr, otherwise a reference to the requested attribute.
   */
  template <AttributeNames attribute>
  constexpr decltype(auto) operator()() {
    if constexpr (attribute == AttributeNames::ptr) {
      return this;
    } else {
      return getRef<attribute>();
    }
  }

  /**
   * Const accessor used by KokkosAoS.
   * @tparam attribute Attribute name.
   * @return Value of the requested attribute.
   */
  template <AttributeNames attribute>
  constexpr auto operator()() const {
    return get<attribute>();
  }

  /**
   * Getter, which returns a reference to an attribute using the corresponding attribute name (defined in
   * AttributeNames).
   * @tparam attribute Attribute name.
   * @return Reference to the requested attribute.
   */
  template <AttributeNames attribute, std::enable_if_t<attribute != AttributeNames::ptr, bool> = true>
  constexpr decltype(auto) getRef() {
    if constexpr (attribute == AttributeNames::id) {
      return (_id);
    } else if constexpr (attribute == AttributeNames::posX) {
      return (_r[0]);
    } else if constexpr (attribute == AttributeNames::posY) {
      return (_r[1]);
    } else if constexpr (attribute == AttributeNames::posZ) {
      return (_r[2]);
    }
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    else if constexpr (attribute == AttributeNames::rebuildX) {
      return (_rAtRebuild[0]);
    } else if constexpr (attribute == AttributeNames::rebuildY) {
      return (_rAtRebuild[1]);
    } else if constexpr (attribute == AttributeNames::rebuildZ) {
      return (_rAtRebuild[2]);
    }
#endif
    else if constexpr (attribute == AttributeNames::velX) {
      return (_v[0]);
    } else if constexpr (attribute == AttributeNames::velY) {
      return (_v[1]);
    } else if constexpr (attribute == AttributeNames::velZ) {
      return (_v[2]);
    } else if constexpr (attribute == AttributeNames::forceX) {
      return (_f[0]);
    } else if constexpr (attribute == AttributeNames::forceY) {
      return (_f[1]);
    } else if constexpr (attribute == AttributeNames::forceZ) {
      return (_f[2]);
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return (_ownershipState);
    } else {
      utils::ExceptionHandler::exception("KokkosParticleBase::getRef() unknown attribute {}", attribute);
    }
  }

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
    }
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    else if constexpr (attribute == AttributeNames::rebuildX) {
      return getRAtRebuild()[0];
    } else if constexpr (attribute == AttributeNames::rebuildY) {
      return getRAtRebuild()[1];
    } else if constexpr (attribute == AttributeNames::rebuildZ) {
      return getRAtRebuild()[2];
    }
#endif
    else if constexpr (attribute == AttributeNames::velX) {
      return getV()[0];
    } else if constexpr (attribute == AttributeNames::velY) {
      return getV()[1];
    } else if constexpr (attribute == AttributeNames::velZ) {
      return getV()[2];
    } else if constexpr (attribute == AttributeNames::forceX) {
      return getF()[0];
    } else if constexpr (attribute == AttributeNames::forceY) {
      return getF()[1];
    } else if constexpr (attribute == AttributeNames::forceZ) {
      return getF()[2];
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      utils::ExceptionHandler::exception("KokkosParticleBase::get() unknown attribute {}", attribute);
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
    }
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    else if constexpr (attribute == AttributeNames::rebuildX) {
      _rAtRebuild[0] = value;
    } else if constexpr (attribute == AttributeNames::rebuildY) {
      _rAtRebuild[1] = value;
    } else if constexpr (attribute == AttributeNames::rebuildZ) {
      _rAtRebuild[2] = value;
    }
#endif
    else if constexpr (attribute == AttributeNames::velX) {
      _v[0] = value;
    } else if constexpr (attribute == AttributeNames::velY) {
      _v[1] = value;
    } else if constexpr (attribute == AttributeNames::velZ) {
      _v[2] = value;
    } else if constexpr (attribute == AttributeNames::forceX) {
      _f[0] = value;
    } else if constexpr (attribute == AttributeNames::forceY) {
      _f[1] = value;
    } else if constexpr (attribute == AttributeNames::forceZ) {
      _f[2] = value;
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else {
      utils::ExceptionHandler::exception("KokkosParticleBase::set() unknown attribute {}", attribute);
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
   * @tparam T
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
std::ostream &operator<<(std::ostream &os, const KokkosParticleBase<floatType, idType> &particle) {
  using utils::ArrayUtils::operator<<;
  os << "Particle"
     << "\nID      : " << particle._id << "\nPosition: " << particle._r << "\nVelocity: " << particle._v
     << "\nForce   : " << particle._f << "\nOwnershipState : " << particle._ownershipState;
  // clang-format on
  return os;
}

}  // namespace autopas

#endif