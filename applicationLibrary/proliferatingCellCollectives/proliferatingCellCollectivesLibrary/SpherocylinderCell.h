/**
 * @file SpherocylinderCell.h
 * @date 11/05/2025
 * @author Manuel Lerchner
 */

#pragma once

#include <array>
#include <cmath>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "autopas/particles/ParticleDefinitions.h"
#include "autopas/utils/ExceptionHandler.h"
#include "molecularDynamicsLibrary/MultisiteMoleculeLJ.h"

namespace pccLib {

/**
 * @class SpherocylinderCell
 * @brief Represents a 3D spherocylinder cell for proliferating cell collectives.
 *
 * This class models a spherocylinder (rod-shaped) cell with position, orientation, velocity, and physical properties.
 * It supports growth, division, movement, and overlap detection for use in cell-based simulations.
 */
class SpherocylinderCell : public mdLib::MultisiteMoleculeLJ {
 public:
  /**
   * @brief Default constructor.
   */
  SpherocylinderCell() = default;

  /**
   * @brief Construct a new SpherocylinderCell with given properties.
   * @param position Initial position (3D vector)
   * @param linearVelocity Initial linear velocity (3D vector)
   * @param angularVelocity Initial angular velocity (3D vector)
   * @param quaternion Initial orientation (quaternion, 4D)
   * @param length0 Initial length
   * @param moleculeId Molecule ID
   */
  SpherocylinderCell(const std::array<double, 3> &position, const std::array<double, 3> &linearVelocity,
                     const std::array<double, 3> &angularVelocity, const std::array<double, 4> &quaternion,
                     double length0, unsigned long moleculeId);
  /**
   * @brief Destructor.
   */
  ~SpherocylinderCell() override = default;

  /**
   * @enum AttributeNames
   * @brief Enum for attribute indexing in SoA and functors.
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
    angularVelX,
    angularVelY,
    angularVelZ,
    quaternion0,
    quaternion1,
    quaternion2,
    quaternion3,
    length,
    typeId,
    ownershipState,
    diameter,
    stress,
    forceX,
    forceY,
    forceZ,
    torqueX,
    torqueY,
    torqueZ
  };

  /**
   * @brief Type for Structure of Arrays (SoA) representation.
   */
  using SoAArraysType = typename autopas::utils::SoAType<
      SpherocylinderCell *, size_t /*id*/, double /*x*/, double /*y*/, double /*z*/, double /*vx*/, double /*vy*/,
      double /*vz*/, double /*angularVelX*/, double /*angularVelY*/, double /*angularVelZ*/, double /*quaternion0*/,
      double /*orientation1*/, double /*quaternion2*/, double /*quaternion3*/, double /*length*/, size_t /*typeId*/,
      autopas::OwnershipState /*ownershipState*/, double /*diameter*/, double /*stress*/, double /*forceX*/,
      double /*forceY*/, double /*forceZ*/, double /*torqueX*/, double /*torqueY*/, double /*torqueZ*/>::Type;

  /**
   * @brief Get attribute value by enum (ptr specialization).
   */
  template <AttributeNames attribute, std::enable_if_t<attribute == AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() {
    return this;
  }

  /**
   * @brief Get attribute value by enum (const specialization).
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
    } else if constexpr (attribute == AttributeNames::angularVelX) {
      return getAngularVel()[0];
    } else if constexpr (attribute == AttributeNames::angularVelY) {
      return getAngularVel()[1];
    } else if constexpr (attribute == AttributeNames::angularVelZ) {
      return getAngularVel()[2];
    } else if constexpr (attribute == AttributeNames::quaternion0) {
      return getQuaternion()[0];
    } else if constexpr (attribute == AttributeNames::quaternion1) {
      return getQuaternion()[1];
    } else if constexpr (attribute == AttributeNames::quaternion2) {
      return getQuaternion()[2];
    } else if constexpr (attribute == AttributeNames::quaternion3) {
      return getQuaternion()[3];
    } else if constexpr (attribute == AttributeNames::length) {
      return getLength();
    } else if constexpr (attribute == AttributeNames::typeId) {
      return getTypeId();
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return getOwnershipState();
    } else if constexpr (attribute == AttributeNames::diameter) {
      return getDiameter();
    } else if constexpr (attribute == AttributeNames::stress) {
      return getStress();
    } else if constexpr (attribute == AttributeNames::forceX) {
      return getF()[0];
    } else if constexpr (attribute == AttributeNames::forceY) {
      return getF()[1];
    } else if constexpr (attribute == AttributeNames::forceZ) {
      return getF()[2];
    } else if constexpr (attribute == AttributeNames::torqueX) {
      return getTorque()[0];
    } else if constexpr (attribute == AttributeNames::torqueY) {
      return _torque[1];
    } else if constexpr (attribute == AttributeNames::torqueZ) {
      return _torque[2];
    } else {
      autopas::utils::ExceptionHandler::exception("SpherocylinderCell::get() unknown attribute {}", attribute);
    }
  }

  /**
   * @brief Set attribute value by enum.
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
    } else if constexpr (attribute == AttributeNames::angularVelX) {
      _angularVel[0] = value;
    } else if constexpr (attribute == AttributeNames::angularVelY) {
      _angularVel[1] = value;
    } else if constexpr (attribute == AttributeNames::angularVelZ) {
      _angularVel[2] = value;
    } else if constexpr (attribute == AttributeNames::quaternion0) {
      _q[0] = value;
    } else if constexpr (attribute == AttributeNames::quaternion1) {
      _q[1] = value;
    } else if constexpr (attribute == AttributeNames::quaternion2) {
      _q[2] = value;
    } else if constexpr (attribute == AttributeNames::quaternion3) {
      _q[3] = value;
    } else if constexpr (attribute == AttributeNames::length) {
      _length = value;
    } else if constexpr (attribute == AttributeNames::typeId) {
      _typeId = value;
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else if constexpr (attribute == AttributeNames::diameter) {
      _diameter = value;
    } else if constexpr (attribute == AttributeNames::stress) {
      _stress = value;
    } else if constexpr (attribute == AttributeNames::forceX) {
      _f[0] = value;
    } else if constexpr (attribute == AttributeNames::forceY) {
      _f[1] = value;
    } else if constexpr (attribute == AttributeNames::forceZ) {
      _f[2] = value;
    } else if constexpr (attribute == AttributeNames::torqueX) {
      _torque[0] = value;
    } else if constexpr (attribute == AttributeNames::torqueY) {
      _torque[1] = value;
    } else if constexpr (attribute == AttributeNames::torqueZ) {
      _torque[2] = value;
    } else {
      autopas::utils::ExceptionHandler::exception("SpherocylinderCell::set() unknown attribute {}", attribute);
    }
  }

  /**
   * @brief Check if the cell has a given attribute (for use in functors, SoA, etc.)
   * @tparam attribute Attribute name
   * @return true if the attribute is present
   */
  template <AttributeNames attribute>
  static constexpr bool has() {
    // All attributes in the enum are present
    return true;
  }

  /**
   * @brief Grow the cell according to stress and time step.
   * @param dt Time step
   * @param tao Growth time constant
   * @param lamb Stress sensitivity
   */
  void grow(double dt, double tao, double lamb);

  /**
   * @brief Divide the cell if it exceeds max length, returning the new cell as std::optional.
   * @return New cell (right) if division occurs, std::nullopt otherwise
   */
  std::optional<SpherocylinderCell> divide();

  /**
   * @brief Get the direction vector of the spherocylinder (from orientation).
   * @return 3D direction vector
   */
  std::array<double, 3> getDirectionVector() const;

  /**
   * @brief Get the endpoints of the spherocylinder (for overlap checks).
   * @return Pair of 3D points (end1, end2)
   */
  std::pair<std::array<double, 3>, std::array<double, 3>> getEndpoints() const;

  /**
   * @brief Calculate collision information between this spherocylinder and another.
   * @param other The other spherocylinder to check collision with
   * @return Optional tuple containing (overlap amount, collision normal, contact point) if collision exists
   */
  std::optional<std::tuple<double, std::array<double, 3>, std::array<double, 3>>> getCollisionInfo(
      const SpherocylinderCell &other) const;

  /** @brief Get the stress. */
  double getStress() const { return _stress; }
  /** @brief Set the stress. */
  void setStress(double s) { _stress = s; }

  /** @brief Get the diameter. */
  double getDiameter() const { return _diameter; }
  /** @brief Set the diameter. */
  void setDiameter(double d) { _diameter = d; }

  /** @brief Get the current length of the cell. */
  double getLength() const { return _length; }
  /** @brief Set the length of the cell. */
  void setLength(double l) { _length = l; }

  /** @brief String representation of the cell. */
  std::string toString() const;

 private:
  /** @brief Initial length of the cell */
  double _length0 = 1.0;

  /** @brief Current length of the cell */
  double _length = 1.0;

  /** @brief Maximum length before division (typically 2 * initial length) */
  double _maxLength = 2.0;

  /** @brief Diameter of the spherical caps at the ends of the cell */
  double _diameter = 0.5;

  /** @brief Compressive stress experienced by the cell, affects growth rate */
  double _stress = 0.0;

  /** @brief Type identifier for the cell */
  size_t _typeId = 0;
};

/**
 * @brief Helper: minimum distance between 3D line segments.
 *
 * Calculates the minimum distance between two line segments in 3D space, with optional clamping
 * at endpoints. This is used for collision detection between spherocylinders.
 * @param a0 First endpoint of segment A
 * @param a1 Second endpoint of segment A
 * @param b0 First endpoint of segment B
 * @param b1 Second endpoint of segment B
 * @param clampAll If true, enables all clamping flags (forces points to stay on segments)
 * @param clampA0 If true, clamps result on segment A to a0
 * @param clampA1 If true, clamps result on segment A to a1
 * @param clampB0 If true, clamps result on segment B to b0
 * @param clampB1 If true, clamps result on segment B to b1
 * @return Tuple containing:
 *         - Closest point on segment A
 *         - Closest point on segment B
 *         - Minimum distance between the segments
 */
static std::tuple<std::array<double, 3>, std::array<double, 3>, double> minimumDistanceBetweenLineSegments(
    const std::array<double, 3> &a0, const std::array<double, 3> &a1, const std::array<double, 3> &b0,
    const std::array<double, 3> &b1, bool clampAll = false, bool clampA0 = false, bool clampA1 = false,
    bool clampB0 = false, bool clampB1 = false);

}  // namespace pccLib
