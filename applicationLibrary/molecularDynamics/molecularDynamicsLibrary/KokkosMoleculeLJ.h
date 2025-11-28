/**
 * @file KokkosMoleculeLJ.h
 *
 * @date 06.11.2025
 * @author Luis Gall
 */

#pragma once

#include "autopas/particles/ParticleDefinitions.h"
#include "autopas/utils/ExceptionHandler.h"

#include <Kokkos_Core.hpp>

namespace mdLib {

/**
 * Molecule class for the LJFunctor.
 */
class KokkosMoleculeLJ : public autopas::ParticleBaseFP64 {
 public:
  KokkosMoleculeLJ() = default;

  /**
   * Constructor of lennard jones molecule with initialization of typeID.
   * @param pos Position of the molecule.
   * @param v Velocity of the molecule.
   * @param moleculeId Unique Id of the molecule.
   * @param typeId TypeId of the molecule.
   */
  KokkosMoleculeLJ(const std::array<double, 3> &pos, const std::array<double, 3> &v, unsigned long moleculeId,
             unsigned long typeId = 0) : autopas::ParticleBaseFP64(pos, v, moleculeId), _typeId(typeId) {}

  ~KokkosMoleculeLJ() override = default;

  // TODO: change mem space later when storing SoA elsewhere
  template <class MemSpace, typename... Types>
  using KokkosSoAArraysTypeHelper = std::tuple<Kokkos::View<Types, MemSpace>...>;

  template <class MemSpace>
  using KokkosSoAArraysType = KokkosSoAArraysTypeHelper<MemSpace, MoleculeLJ *, size_t /*id*/, double /*x*/, double /*y*/, double /*z*/,
                                       double /*vx*/, double /*vy*/, double /*vz*/, double /*fx*/, double /*fy*/,
                                       double /*fz*/, double /*oldFx*/, double /*oldFy*/, double /*oldFz*/,
                                       size_t /*typeid*/, autopas::OwnershipState /*ownershipState*/>;

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
    } else if constexpr (attribute == AttributeNames::typeId) {
      return getTypeId();
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      autopas::utils::ExceptionHandler::exception("KokkosMoleculeLJ::get() unknown attribute {}", attribute);
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
    } else if constexpr (attribute == AttributeNames::typeId) {
      setTypeId(value);
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else {
      autopas::utils::ExceptionHandler::exception("KokkosMoleculeLJ::set() unknown attribute {}", attribute);
    }
  }

  /**
   * Get the old force.
   * @return
   */
  const std::array<double, 3> &getOldF() const {
    return _oldF;
  }

  /**
   * Set old force.
   * @param oldForce
   */
  void setOldF(const std::array<double, 3> &oldForce) {
    _oldF = oldForce;
  }

  /**
   * Get TypeId.
   * @return
   */
  size_t getTypeId() const {
    return _typeId;
  }

  /**
   * Set the type id of the Molecule.
   * @param typeId
   */
  void setTypeId(size_t typeId) {
    _typeId = typeId;
  }

  /**
   * Creates a string containing all data of the particle.
   * @return String representation.
   */
  [[nodiscard]] std::string toString() const override {
    using autopas::utils::ArrayUtils::operator<<;
    std::ostringstream text;
    // clang-format off
    text << "KokkosMoleculeLJ"
       << "\nID                 : " << _id
       << "\nPosition           : " << _r
       << "\nVelocity           : " << _v
       << "\nForce              : " << _f
       << "\nOld Force          : " << _oldF
       << "\nType ID            : " << _typeId
       << "\nOwnershipState     : " << _ownershipState;
    // clang-format on
    return text.str();
  }

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
  std::array<double, 3> _oldF = {0., 0., 0.};
};

}  // namespace mdLib
