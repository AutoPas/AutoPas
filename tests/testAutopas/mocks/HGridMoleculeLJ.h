/**
 * @file HGridMoleculeLJ.h
 * @author atacann
 * @date 07.12.2024
 */

#pragma once

#include <gmock/gmock-function-mocker.h>

#include <random>

#include "molecularDynamicsLibrary/MoleculeLJ.h"

class HGridMoleculeLJ : public mdLib::MoleculeLJ {
 public:
  HGridMoleculeLJ(const std::array<double, 3> &pos, const std::array<double, 3> &v, unsigned long moleculeId,
                  unsigned long typeId = 0)
      : mdLib::MoleculeLJ(pos, v, moleculeId, typeId) {
    generateSize();
  }
  HGridMoleculeLJ() : mdLib::MoleculeLJ() { generateSize(); };
  ~HGridMoleculeLJ() override = default;

  // Explicitly define copy constructor and assignment operator
  HGridMoleculeLJ(const HGridMoleculeLJ &other) : mdLib::MoleculeLJ(other) { _size = other._size; }

  HGridMoleculeLJ &operator=(const HGridMoleculeLJ &other) {
    if (this != &other) {
      mdLib::MoleculeLJ::operator=(other);
      _size = other._size;
    }
    return *this;
  }

  void generateSize() {
    std::mt19937 gen(42);
    const int n = 4;
    std::uniform_int_distribution<int> dist(1, n);
    _size = dist(gen) * 1.0 / n;
  }

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

  double getSize() const override { return _size; }

  /**
   * The type for the SoA storage.
   *
   */
  using SoAArraysType =
      typename autopas::utils::SoAType<HGridMoleculeLJ *, size_t /*id*/, double /*x*/, double /*y*/, double /*z*/,
                                       double /*vx*/, double /*vy*/, double /*vz*/, double /*fx*/, double /*fy*/,
                                       double /*fz*/, double /*oldFx*/, double /*oldFy*/, double /*oldFz*/,
                                       size_t /*typeid*/, autopas::OwnershipState /*ownershipState*/>::Type;

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
      autopas::utils::ExceptionHandler::exception("HGridMoleculeLJ::get() unknown attribute {}", attribute);
    }
  }

  /**
   * Setter, which allows set an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @param value New value of the requested attribute.
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
      autopas::utils::ExceptionHandler::exception("HGridMoleculeLJ::set() unknown attribute {}", attribute);
    }
  }

  /**
   * Creates a string containing all data of the particle.
   * @return String representation.
   */
  std::string toString() const override {
    using autopas::utils::ArrayUtils::operator<<;
    std::ostringstream text;
    // clang-format off
    text << MoleculeLJ::toString()
        << "\nSize               : " << _size;
    // clang-format on
    return text.str();
  }

  double _size;
};
