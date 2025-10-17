/**
 * @file GranularDEMParticle.h
 * @author Joon Kim
 * @date 27/03/2025
 */

#pragma once

#include <vector>

#include "autopas/particles/ParticleDefinitions.h"
#include "autopas/utils/ExceptionHandler.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"

namespace demLib {

/**
 * Granular Particle class for the DEMFunctor.
 */

class GranularDEMParticle : public autopas::ParticleBaseFP64 {
 public:
  GranularDEMParticle() = default;

  /**
   * Constructor of granular particle with initialization of typeID.
   * @param pos Position of the particle.
   * @param v Velocity of the particle.
   * @param angularVel Angular velocity of the particle.
   * @param particleId Unique Id of the particle.
   * @param typeId TypeId of the particle.
   */
  GranularDEMParticle(const std::array<double, 3> &pos, const std::array<double, 3> &v,
                      const std::array<double, 3> angularVel, unsigned long particleId, unsigned long typeId = 0,
                      const double temperature = 0.);

  ~GranularDEMParticle() override = default;

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
    angularVelX,
    angularVelY,
    angularVelZ,
    torqueX,
    torqueY,
    torqueZ,
    temperature,
    heatFlux,
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
      GranularDEMParticle *,
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
      double, // angVx
      double, // angVy
      double, // angVz
      double, // tx
      double, // ty
      double, // tz
      double, // temperature
      double, // heatFlux
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
    } else if constexpr (attribute == AttributeNames::temperature) {
      return getTemperature();
    } else if constexpr (attribute == AttributeNames::heatFlux) {
      return getHeatFlux();
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
    } else if constexpr (attribute == AttributeNames::temperature) {
      _temperature = value;
    } else if constexpr (attribute == AttributeNames::heatFlux) {
      _heatFlux = value;
    } else if constexpr (attribute == AttributeNames::typeId) {
      _typeId = value;
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else {
      autopas::utils::ExceptionHandler::exception("MultisiteMoleculeLJ::set() unknown attribute {}", attribute);
    }
  }

  /**
   * Get the old force.
   * @return
   */
  [[nodiscard]] const std::array<double, 3> &getOldF() const;

  /**
   * Set old force.
   * @param oldForce
   */
  void setOldF(const std::array<double, 3> &oldForce);

  /**
   * Get TypeId.
   * @return
   */
  [[nodiscard]] size_t getTypeId() const;

  /**
   * Set the type id of the Molecule.
   * @param typeId
   */
  void setTypeId(size_t typeId);

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
   * Get the temperature.
   * @return temperature
   */
  [[nodiscard]] double getTemperature() const;

  /**
   * Set the temperature
   * @param temperature
   */
  void setTemperature(double temperature);

  /**
   * Adds given temperature to the particle's temperature.
   * @param temperature
   */
  void addTemperature(double temperature);

  /**
   * Get the heat flux.
   * @return heat flux
   */
  [[nodiscard]] double getHeatFlux() const;

  /**
   * Set the heat flux.
   * @param heatFlux
   */
  void setHeatFlux(double heatFlux);

  /**
   * Adds given heat flux to the particle's heat flux.
   * @param heatFlux
   */
  void addHeatFlux(double heatFlux);

  /**
   * Subtracts given heat flux to the particle's heat flux.
   * @param heatFlux
   */
  void subHeatFlux(double heatFlux);

  [[nodiscard]] double getSize() const override;

  /**
   * Set the particle properties library.
   * @param particlePropertiesLibrary Shared pointer to the particle properties library.
   */
  static void setParticlePropertiesLibrary(
      const std::shared_ptr<ParticlePropertiesLibrary<> > &particlePropertiesLibrary) {
    _particlePropertiesLibrary = particlePropertiesLibrary;
  }

  /**
   * Set the cutoff multiplier.
   * @param cutoffMultiplier The new cutoff multiplier value.
   */
  static void setCutoffMultiplier(double cutoffMultiplier) { _cutoffMultiplier = cutoffMultiplier; }

  /**
   * Creates a string containing all data of the particle.
   * @return String representation.
   */
  [[nodiscard]] std::string toString() const override;

 protected:
  /**
   * Particle type id.
   */
  size_t _typeId = 0;

  /**
   * Old Force of the particle experiences as 3D vector.
   */
  std::array<double, 3> _oldF{};

  /**
   * Angular velocity of the particle as 3D vector.
   */
  std::array<double, 3> _angularVel{};

  /**
   * Torque applied to particle.
   */
  std::array<double, 3> _torque{};

  /**
   * Temperature of the particle.
   */
  double _temperature{};

  /**
   * Heatflux (flow of heat energy) of the particle.
   */
  double _heatFlux{};

  /**
 * A static refernce to particlePropertiesLibrary so that sigma of particles can be used in getSize() function.
 */
  static std::shared_ptr<ParticlePropertiesLibrary<> > _particlePropertiesLibrary;
  /**
   * The multiplier the sigma values will be scaled by.
   */
  static double _cutoffMultiplier;
};
}  // namespace demLib