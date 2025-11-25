/**
 * @file MoleculeLJBase.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <vector>

#include "autopas/particles/ParticleBase.h"
#include "autopas/utils/ExceptionHandler.h"

namespace mdLib {

/**
 * Molecule class for the LJFunctor.
 */
template <typename CalcType, typename AccuType, typename idType>
class MoleculeLJBase : public autopas::ParticleBase<CalcType, AccuType, idType> {
 public:
  MoleculeLJBase() = default;

  /**
   * Constructor of lennard jones molecule with initialization of typeID.
   * @param pos Position of the molecule.
   * @param v Velocity of the molecule.
   * @param moleculeId Unique Id of the molecule.
   * @param typeId TypeId of the molecule.
   */
  MoleculeLJBase(const std::array<CalcType, 3> &pos, const std::array<CalcType, 3> &v, unsigned long moleculeId,
                 unsigned long typeId = 0)
      : autopas::ParticleBase<CalcType, AccuType, idType>(pos, v, moleculeId), _typeId(typeId){};

  ~MoleculeLJBase() override = default;

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
      MoleculeLJBase *, size_t /*id*/, CalcType /*x*/, CalcType /*y*/, CalcType /*z*/, CalcType /*vx*/, CalcType /*vy*/,
      CalcType /*vz*/, AccuType /*fx*/, AccuType /*fy*/, AccuType /*fz*/, AccuType /*oldFx*/, AccuType /*oldFy*/,
      AccuType /*oldFz*/, size_t /*typeid*/, autopas::OwnershipState /*ownershipState*/>::Type;

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
      return this->getID();
    } else if constexpr (attribute == AttributeNames::posX) {
      return this->getR()[0];
    } else if constexpr (attribute == AttributeNames::posY) {
      return this->getR()[1];
    } else if constexpr (attribute == AttributeNames::posZ) {
      return this->getR()[2];
    } else if constexpr (attribute == AttributeNames::velocityX) {
      return this->getV()[0];
    } else if constexpr (attribute == AttributeNames::velocityY) {
      return this->getV()[1];
    } else if constexpr (attribute == AttributeNames::velocityZ) {
      return this->getV()[2];
    } else if constexpr (attribute == AttributeNames::forceX) {
      return this->getF()[0];
    } else if constexpr (attribute == AttributeNames::forceY) {
      return this->getF()[1];
    } else if constexpr (attribute == AttributeNames::forceZ) {
      return this->getF()[2];
    } else if constexpr (attribute == AttributeNames::oldForceX) {
      return this->getOldF()[0];
    } else if constexpr (attribute == AttributeNames::oldForceY) {
      return this->getOldF()[1];
    } else if constexpr (attribute == AttributeNames::oldForceZ) {
      return this->getOldF()[2];
    } else if constexpr (attribute == AttributeNames::typeId) {
      return this->getTypeId();
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      autopas::utils::ExceptionHandler::exception("MoleculeLJBase::get() unknown attribute {}", attribute);
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
      this->setID(value);
    } else if constexpr (attribute == AttributeNames::posX) {
      this->_r[0] = value;
    } else if constexpr (attribute == AttributeNames::posY) {
      this->_r[1] = value;
    } else if constexpr (attribute == AttributeNames::posZ) {
      this->_r[2] = value;
    } else if constexpr (attribute == AttributeNames::velocityX) {
      this->_v[0] = value;
    } else if constexpr (attribute == AttributeNames::velocityY) {
      this->_v[1] = value;
    } else if constexpr (attribute == AttributeNames::velocityZ) {
      this->_v[2] = value;
    } else if constexpr (attribute == AttributeNames::forceX) {
      this->_f[0] = value;
    } else if constexpr (attribute == AttributeNames::forceY) {
      this->_f[1] = value;
    } else if constexpr (attribute == AttributeNames::forceZ) {
      this->_f[2] = value;
    } else if constexpr (attribute == AttributeNames::oldForceX) {
      this->_oldF[0] = value;
    } else if constexpr (attribute == AttributeNames::oldForceY) {
      this->_oldF[1] = value;
    } else if constexpr (attribute == AttributeNames::oldForceZ) {
      this->_oldF[2] = value;
    } else if constexpr (attribute == AttributeNames::typeId) {
      setTypeId(value);
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else {
      autopas::utils::ExceptionHandler::exception("MoleculeLJBase::set() unknown attribute {}", attribute);
    }
  }

  /**
   * Get the old force.
   * @return
   */
  [[nodiscard]] const std::array<AccuType, 3> &getOldF() const { return _oldF; };

  /**
   * Set old force.
   * @param oldForce
   */
  void setOldF(const std::array<AccuType, 3> &oldForce) { _oldF = oldForce; };

  /**
   * Get TypeId.
   * @return
   */
  [[nodiscard]] size_t getTypeId() const { return _typeId; };

  /**
   * Set the type id of the Molecule.
   * @param typeId
   */
  void setTypeId(size_t typeId) { _typeId = typeId; };

  /**
   * Creates a string containing all data of the particle.
   * @return String representation.
   */
  [[nodiscard]] std::string toString() const override {
    using autopas::utils::ArrayUtils::operator<<;
    std::ostringstream text;
    // clang-format off
  text << "MoleculeLJBase"
     << "\nID                 : " << this->_id
     << "\nPosition           : " << this->_r
     << "\nVelocity           : " << this->_v
     << "\nForce              : " << this->_f
     << "\nOld Force          : " << this->_oldF
     << "\nType ID            : " << this->_typeId
     << "\nOwnershipState     : " << this->_ownershipState;
    // clang-format on
    return text.str();
  };

 protected:
  /**
   * Molecule type id. In single-site simulations, this is used as a siteId to look up site attributes in the particle
   * properties library.
   *
   * In multi-site simulations, where a multi-site molecule class inheriting from this class is used, typeId is used as
   * a molId to look up molecular attributes (including siteIds of the sites).
   */
  size_t _typeId = 0;

  /**
   * Old Force of the particle experiences as 3D vector.
   */
  std::array<AccuType, 3> _oldF = {0., 0., 0.};
};

}  // namespace mdLib
