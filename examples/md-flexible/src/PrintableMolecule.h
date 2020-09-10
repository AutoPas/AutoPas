/**
 * @file PrintableMolecule.h
 * @author F. Gratl
 * @date 11/26/18
 */

#pragma once

#include <array>
#include <iomanip>
#include <iostream>

#include "autopas/molecularDynamics/MoleculeLJ.h"

/**
 * Example for a custom particle type derived from a autopas molecule type.
 */
class PrintableMolecule
    : public autopas::MoleculeLJ<> /*apparently c++17 doesnt need <> but doesnt compile without it*/ {
 public:
  /**
   * Empty Constructor.
   */
  PrintableMolecule() : autopas::MoleculeLJ<>() {}

  /**
   * Constructor.
   * @param pos Position
   * @param v Veloctiy
   * @param moleculeId Molecule ID
   * @param typeId Molecule Type ID
   */
  PrintableMolecule(std::array<double, 3> pos, std::array<double, 3> v, unsigned long moleculeId,
                    unsigned int typeId = 0)
      : autopas::MoleculeLJ<>(pos, v, moleculeId, typeId) {}

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int { ptr, id, posX, posY, posZ, forceX, forceY, forceZ, typeId, ownershipState };

  /**
   * The type for the SoA storage.
   */
  using SoAArraysType =
      typename autopas::utils::SoAType<PrintableMolecule *, size_t /*id*/, double /*x*/, double /*y*/, double /*z*/,
                                       double /*fx*/, double /*fy*/, double /*fz*/, size_t /*typeid*/,
                                       autopas::OwnershipState /*ownershipState*/>::Type;

  /**
   * Getter, which allows access to an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @return Value of the requested attribute.
   * @note The value of owned is return as floating point number (true = 1.0, false = 0.0).
   */
  template <AttributeNames attribute>
  constexpr typename std::tuple_element<static_cast<size_t>(attribute), SoAArraysType>::type::value_type get() {
    if constexpr (attribute == AttributeNames::ptr) {
      return this;
    } else if constexpr (attribute == AttributeNames::id) {
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
      autopas::utils::ExceptionHandler::exception("ParticleBase::get() unknown attribute {}", attribute);
    }
  }

  /**
   * Print molecule properties to std out.
   */
  void print() {
    std::cout << "Molecule with position: ";
    for (auto &r : this->getR()) {
      std::cout << std::setw(10) << r << ", ";
    }
    std::cout << "and force: ";

    for (auto &f : this->getF()) {
      std::cout << std::setw(15) << f << ", ";
    }
    std::cout << "ID: " << std::setw(5) << this->getID();
    std::cout << std::endl;
  }
};
