/**
 * @file MultisiteMoleculeLJ.h
 * @date 14/02/2022
 * @author S. Newcome
 */

#pragma once

#include "MoleculeLJ.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/particles/ParticleBase.h"

namespace mdLib {
/**
 * Standard multi-site LJ molecules.
 *
 * The molecule is treated as a single particle for the purposes of cutoffs and containers, with a quaternion for
 * angular direction, a 3D vector-array for angular velocity, and a vectors of site positions relative to the CoM and
 * angular direction.
 *
 */
class MultisiteMoleculeLJ : public mdLib::MoleculeLJ {
  using idType = size_t;

 public:
  MultisiteMoleculeLJ() = default;

  /**
   * Constructor of the MultisiteMoleculeLJ Class
   * @param r Position of the particle.
   * @param v Velocity of the particle.
   * @param q Quaternion defining rotation of particle.
   * @param angularVel Rotational velocity of the particle.
   * @param moleculeId Id of the particle.
   * @param typeId Id of the type of the particle. Used in conjunction with ParticlePropertiesLibrary to access
   * molecular information such as site types and relative site positions.
   */
  MultisiteMoleculeLJ(std::array<double, 3> r, std::array<double, 3> v, std::array<double, 4> q,
                      std::array<double, 3> angularVel, unsigned long moleculeId, unsigned long typeId = 0);

  /**
   * Destructor of the MultisiteMoleculeLJ class.
   */
  ~MultisiteMoleculeLJ() override = default;

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int {
    ptr,
    id,
    posX,
    posY,
    posZ,
    velocityX,
    velocityY,
    velocityZ,
    forceX,
    forceY,
    forceZ,
    oldForceX,
    oldForceY,
    oldForceZ,
    quaternion0,
    quaternion1,
    quaternion2,
    quaternion3,
    angularVelX,
    angularVelY,
    angularVelZ,
    torqueX,
    torqueY,
    torqueZ,
    typeId,
    ownershipState
  };

  /**
   * The type for the SoA storage.
   *
   * @note The attribute owned is of type float but treated as a bool.
   * This means it shall always only take values 0.0 (=false) or 1.0 (=true).
   * The reason for this is the easier use of the value in calculations (See LJFunctor "energyFactor")
   */
  // clang-format off
  using SoAArraysType = typename autopas::utils::SoAType<
      MultisiteMoleculeLJ *,
      size_t, // id
      double, // x
      double, // y
      double, // z
      double, // vx
      double, // vy
      double, // vz
      double, // fx
      double, // fy
      double, // fz
      double, // oldFx
      double, // oldFy
      double, // oldFz
      double, // q0
      double, // q1
      double, // q2
      double, // q3
      double, // angVx
      double, // angVy
      double, // angVz
      double, // tx
      double, // ty
      double, // tz
      size_t, // typeid
      autopas::OwnershipState //ownerState
  >::Type;
  // clang-format on

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
   * @note Moving this function to the .cpp leads to undefined references
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
    } else if constexpr (attribute == AttributeNames::velocityX) {
      return getV()[0];
    } else if constexpr (attribute == AttributeNames::velocityY) {
      return getV()[1];
    } else if constexpr (attribute == AttributeNames::velocityZ) {
      return getV()[2];
    } else if constexpr (attribute == AttributeNames::forceX) {
      return getF()[0];
    } else if constexpr (attribute == AttributeNames::forceY) {
      return getF()[1];
    } else if constexpr (attribute == AttributeNames::forceZ) {
      return getF()[2];
    } else if constexpr (attribute == AttributeNames::oldForceX) {
      return getOldF()[0];
    } else if constexpr (attribute == AttributeNames::oldForceY) {
      return getOldF()[1];
    } else if constexpr (attribute == AttributeNames::oldForceZ) {
      return getOldF()[2];
    } else if constexpr (attribute == AttributeNames::quaternion0) {
      return getQ()[0];
    } else if constexpr (attribute == AttributeNames::quaternion1) {
      return getQ()[1];
    } else if constexpr (attribute == AttributeNames::quaternion2) {
      return getQ()[2];
    } else if constexpr (attribute == AttributeNames::quaternion3) {
      return getQ()[3];
    } else if constexpr (attribute == AttributeNames::angularVelX) {
      return getAngularVel()[0];
    } else if constexpr (attribute == AttributeNames::angularVelY) {
      return getAngularVel()[1];
    } else if constexpr (attribute == AttributeNames::angularVelZ) {
      return getAngularVel()[2];
    } else if constexpr (attribute == AttributeNames::torqueX) {
      return getTorque()[0];
    } else if constexpr (attribute == AttributeNames::torqueY) {
      return getTorque()[1];
    } else if constexpr (attribute == AttributeNames::torqueZ) {
      return getTorque()[2];
    } else if constexpr (attribute == AttributeNames::typeId) {
      return getTypeId();
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      autopas::utils::ExceptionHandler::exception("MultisiteMoleculeLJ::get() unknown attribute {}", attribute);
    }
  }

  /**
   * Setter, which allows set an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @param value New value of the requested attribute.
   * @note The value of owned is extracted from a floating point number (true = 1.0, false = 0.0).
   * @note Moving this function to the .cpp leads to undefined references
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
    } else if constexpr (attribute == AttributeNames::velocityX) {
      _v[0] = value;
    } else if constexpr (attribute == AttributeNames::velocityY) {
      _v[1] = value;
    } else if constexpr (attribute == AttributeNames::velocityZ) {
      _v[2] = value;
    } else if constexpr (attribute == AttributeNames::forceX) {
      _f[0] = value;
    } else if constexpr (attribute == AttributeNames::forceY) {
      _f[1] = value;
    } else if constexpr (attribute == AttributeNames::forceZ) {
      _f[2] = value;
    } else if constexpr (attribute == AttributeNames::oldForceX) {
      _oldF[0] = value;
    } else if constexpr (attribute == AttributeNames::oldForceY) {
      _oldF[1] = value;
    } else if constexpr (attribute == AttributeNames::oldForceZ) {
      _oldF[2] = value;
    } else if constexpr (attribute == AttributeNames::quaternion0) {
      _q[0] = value;
    } else if constexpr (attribute == AttributeNames::quaternion1) {
      _q[1] = value;
    } else if constexpr (attribute == AttributeNames::quaternion2) {
      _q[2] = value;
    } else if constexpr (attribute == AttributeNames::quaternion3) {
      _q[3] = value;
    } else if constexpr (attribute == AttributeNames::angularVelX) {
      _angularVel[0] = value;
    } else if constexpr (attribute == AttributeNames::angularVelY) {
      _angularVel[1] = value;
    } else if constexpr (attribute == AttributeNames::angularVelZ) {
      _angularVel[2] = value;
    } else if constexpr (attribute == AttributeNames::torqueX) {
      _torque[0] = value;
    } else if constexpr (attribute == AttributeNames::torqueY) {
      _torque[1] = value;
    } else if constexpr (attribute == AttributeNames::torqueZ) {
      _torque[2] = value;
    } else if constexpr (attribute == AttributeNames::typeId) {
      _typeId = value;
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else {
      autopas::utils::ExceptionHandler::exception("MultisiteMoleculeLJ::set() unknown attribute {}", attribute);
    }
  }

  /**
   * Get the quaternion defining rotation
   * @return quaternion defining rotation
   */
  [[nodiscard]] const std::array<double, 4> &getQ() const;

  /**
   * Set the quaternion defining rotation
   * @param q quaternion defining rotation
   */
  void setQ(const std::array<double, 4> &q);

  /**
   * Get the angular velocity
   * @return angular velocity
   */
  [[nodiscard]] const std::array<double, 3> &getAngularVel() const;

  /**
   * Set the angular velocity
   * @param angularVel
   */
  void setAngularVel(const std::array<double, 3> &angularVel);

  /**
   * Adds given angular velocity to the particle's angular velocity.
   * @param angularVel angular velocity to be added
   */
  void addAngularVel(const std::array<double, 3> &angularVel);

  /**
   * Get the torque.
   * @return torque
   */
  [[nodiscard]] const std::array<double, 3> &getTorque() const;

  /**
   * Set the torque.
   * @param torque
   */
  void setTorque(const std::array<double, 3> &torque);

  /**
   * Adds given torque to the particle's torque.
   * @param torque torque to be added
   */
  void addTorque(const std::array<double, 3> &torque);

  /**
   * Subracts given torque to the particle's torque.
   * @param torque torque to be subtracted
   */
  void subTorque(const std::array<double, 3> &torque);

  /**
   * Creates a string containing all data of the particle.
   * @return String representation.
   */
  [[nodiscard]] std::string toString() const override;

 protected:
  /**
   * Rotational direction of particle as quaternion.
   */
  std::array<double, 4> _q{};

  /**
   * Angular velocity of the particle
   */
  std::array<double, 3> _angularVel{};

  /**
   * Torque applied to particle.
   */
  std::array<double, 3> _torque{};
};

}  // namespace mdLib