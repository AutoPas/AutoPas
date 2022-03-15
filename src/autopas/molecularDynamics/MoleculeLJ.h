/**
 * @file MoleculeLJ.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <vector>

#include "autopas/particles/Particle.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * Molecule class for the LJFunctor.
 */
class MoleculeLJ : public autopas::Particle {
 public:
  MoleculeLJ() = default;

  /**
   * Constructor of lennard jones molecule with initialization of typeID.
   * @param pos Position of the molecule.
   * @param v Velocitiy of the molecule.
   * @param moleculeId Id of the molecule.
   * @param typeId TypeId of the molecule.
   */
  explicit MoleculeLJ(std::array<double, 3> pos, std::array<double, 3> v, unsigned long moleculeId,
                      unsigned long typeId = 0)
      : Particle(pos, v, moleculeId), _typeId(typeId) {}

  ~MoleculeLJ() = default;

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
  using SoAArraysType = typename autopas::utils::SoAType<
      MoleculeLJ *, size_t /*id*/, double /*x*/, double /*y*/, double /*z*/, double /*vx*/,
      double /*vy*/, double /*vz*/, double /*fx*/, double /*fy*/, double /*fz*/, double /*oldFx*/,
      double /*oldFy*/, double /*oldFz*/, size_t /*typeid*/, OwnershipState /*ownershipState*/>::Type;

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
    } else if constexpr (attribute == AttributeNames::typeId) {
      return getTypeId();
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      utils::ExceptionHandler::exception("MoleculeLJ::get() unknown attribute {}", attribute);
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
    } else if constexpr (attribute == AttributeNames::typeId) {
      setTypeId(value);
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else {
      utils::ExceptionHandler::exception("MoleculeLJ::set() unknown attribute {}", attribute);
    }
  }

  /**
   * Get the old force.
   * @return
   */
  [[nodiscard]] std::array<double, 3> getOldF() const { return _oldF; }

  /**
   * Set old force.
   * @param oldForce
   */
  void setOldF(const std::array<double, 3> &oldForce) { _oldF = oldForce; }

  /**
   * Get TypeId.
   * @return
   */
  [[nodiscard]] size_t getTypeId() const { return _typeId; }

  /**
   * Set the type id of the Molecule.
   * @param typeId
   */
  void setTypeId(size_t typeId) { _typeId = typeId; }

  /**
   * Get the quaternion defining rotation. (Returns nothing of interest + throws exception)
   * @return quaternion defining rotation
   */
  [[nodiscard]] virtual const std::array<double, 4> &getQ() const {
    autopas::utils::ExceptionHandler::exception("Wrong molecule type! MoleculeLJ does not include quaternion");
    return {0.,0.,0.,0.};
  }

  /**
   * Set the quaternion defining rotation. (throws exception)
   * @param q quaternion defining rotation
   */
  virtual void setQ(const std::array<double, 4> &q) {
    autopas::utils::ExceptionHandler::exception("Wrong molecule type! MoleculeLJ does not include quaternion");
  }

  /**
   * Get the angular velocity. (Returns nothing of interest + throws exception)
   * @return angular velocity
   */
  [[nodiscard]] virtual const std::array<double, 3> &getAngularVel() const {
    autopas::utils::ExceptionHandler::exception("Wrong molecule type! MoleculeLJ does not include angular velocity");
    return {0.,0.,0.};
  }

  /**
   * Set the angular velocity. (throws exception)
   * @param angularVelocity
   */
  virtual void setAngularVel(const std::array<double, 3> &angularVel) {
    autopas::utils::ExceptionHandler::exception("Wrong molecule type! MoleculeLJ does not include angular velocity");
  }

  /**
   * Get the torque. (Returns nothing of interest + throws exception)
   * @return torque
   */
  [[nodiscard]] virtual const std::array<double, 3> &getTorque() const {
    autopas::utils::ExceptionHandler::exception("Wrong molecule type! MoleculeLJ does not include torque");
    return {0.,0.,0.};
  }

  /**
   * Set the torque. (throws exception)
   * @param torque
   */
  virtual void setTorque(const std::array<double, 3> &torque) {
    autopas::utils::ExceptionHandler::exception("Wrong molecule type! MoleculeLJ does not include torque");
  }

  /**
    * Adds given torque to the particle's torque. (throws exception)
    * @param torque torque to be added
   */
  void addTorque(const std::array<double, 3> &torque) {
    autopas::utils::ExceptionHandler::exception("Wrong molecule type! MoleculeLJ does not include torque");
  }

  /**
    * Subracts given torque to the particle's torque. (throws exception)
    * @param torque torque to be subtracted
   */
  void subTorque(const std::array<double, 3> &torque) {
    autopas::utils::ExceptionHandler::exception("Wrong molecule type! MoleculeLJ does not include torque");
  }

  /**
   * Returns molecule of type MoleculeLJ, with the same position, velocity, Id, and type Id as this molecule.
   * Throws exception when called (should be used to convert from molecules with more data members to moleculeLJ).
   * @tparam returnedType type of returned
   * @return
   */
  template <class returnedType>
  returnedType returnSimpleMolecule() {
    utils::ExceptionHandler::exception("Converting from MoleculeLJ to MoleculeLJ. This function should not be called.");
    returnedType simpleMolecule;
    simpleMolecule.setR(this->getR());
    simpleMolecule.setV(this->getV());
    simpleMolecule.setID(this->getID());
    simpleMolecule.setTypeId(this->getTypeId());
    return simpleMolecule;
  }

 private:
  /**
   * Particle type id.
   */
  size_t _typeId = 0;

  /**
   * Old Force of the particle experiences as 3D vector.
   */
  std::array<double, 3> _oldF = {0., 0., 0.};
};

}  // namespace autopas
