/**
 * @file MoleculeLJ.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <vector>

#include "autopas/particles/Particle.h"

namespace autopas {

/**
 * Molecule class for the LJFunctor.
 */
template <typename floatType = double>
class MoleculeLJ : public Particle {
 public:
  MoleculeLJ() = default;

  /**
   * Constructor of lennard jones molecule with initialization of typeID.
   * @param pos Position of the molecule.
   * @param v Velocitiy of the molecule.
   * @param moleculeId Id of the molecule.
   * @param typeId TypeId of the molecule.
   */
  explicit MoleculeLJ(std::array<floatType, 3> pos, std::array<floatType, 3> v, unsigned long moleculeId,
                      unsigned long typeId = 0)
      : Particle(pos, v, moleculeId), _typeId(typeId) {}

  ~MoleculeLJ() override = default;

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int { id, posX, posY, posZ, forceX, forceY, forceZ, typeId, ownershipState };

  /**
   * The type for the SoA storage.
   *
   * @note The attribute owned is of type float but treated as a bool.
   * This means it shall always only take values 0.0 (=false) or 1.0 (=true).
   * The reason for this is the easier use of the value in calculations (See LJFunctor "energyFactor")
   */
  using SoAArraysType =
      typename autopas::utils::SoAType<size_t /*id*/, floatType /*x*/, floatType /*y*/, floatType /*z*/,
                                       floatType /*fx*/, floatType /*fy*/, floatType /*fz*/, size_t /*typeid*/,
                                       OwnershipState /*ownedState*/>::Type;

#if defined(AUTOPAS_CUDA)
  /**
   * The type for storage arrays for Cuda.
   */
  using CudaDeviceArraysType = typename autopas::utils::CudaSoAType<size_t, floatType, floatType, floatType, floatType,
                                                                    floatType, floatType, size_t, floatType>::Type;
#endif

  /**
   * Getter, which allows access to an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @return Value of the requested attribute.
   * @note The value of owned is return as floating point number (true = 1.0, false = 0.0).
   */
  template <AttributeNames attribute>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() const {
    if constexpr (attribute == AttributeNames::id) {
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
    } else if constexpr (attribute == AttributeNames::forceX) {
      _f[0] = value;
    } else if constexpr (attribute == AttributeNames::forceY) {
      _f[1] = value;
    } else if constexpr (attribute == AttributeNames::forceZ) {
      _f[2] = value;
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
  [[nodiscard]] std::array<double, 3> getOldf() const { return _oldF; }

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
