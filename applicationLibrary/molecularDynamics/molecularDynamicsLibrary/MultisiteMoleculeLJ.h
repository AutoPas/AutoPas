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
   */
  template <AttributeNames attribute, std::enable_if_t<attribute != AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() const;

  /**
   * Setter, which allows set an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @param value New value of the requested attribute.
   * @note The value of owned is extracted from a floating point number (true = 1.0, false = 0.0).
   */
  template <AttributeNames attribute>
  constexpr void set(typename std::tuple_element<attribute, SoAArraysType>::type::value_type value);

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

 public:
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
};

}  // namespace mdLib