/**
 * @file MoleculeLJ.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <vector>

#include "autopas/particles/ParticleDefinitions.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/KokkosSoA.h"

namespace mdLib {

/**
 * Molecule class for the LJFunctor.
 */
class MoleculeLJ : public autopas::ParticleBaseFP64 {
 public:
  MoleculeLJ() = default;

  /**
   * Constructor of lennard jones molecule with initialization of typeID.
   * @param pos Position of the molecule.
   * @param v Velocity of the molecule.
   * @param moleculeId Unique Id of the molecule.
   * @param typeId TypeId of the molecule.
   */
  MoleculeLJ(const std::array<ParticleSoAFloatPrecision, 3> &pos, const std::array<ParticleSoAFloatPrecision, 3> &v, unsigned long moleculeId,
             unsigned long typeId = 0);

  ~MoleculeLJ() override = default;

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : size_t {
    ptr,
    id,
    posX,
    posY,
    posZ,
    rebuildX,
    rebuildY,
    rebuildZ,
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
    mass,
    ownershipState
  };

  /**
   * The type for the SoA storage.
   *
   * @note The attribute owned is of type float but treated as a bool.
   * This means it shall always only take values 0.0 (=false) or 1.0 (=true).
   * The reason for this is the easier use of the value in calculations (See LJFunctor "energyFactor")
   */
  using SoAArraysType =
      autopas::utils::SoAType<MoleculeLJ *, size_t /*id*/, ParticleSoAFloatPrecision /*x*/, ParticleSoAFloatPrecision /*y*/, ParticleSoAFloatPrecision /*z*/,
                                       ParticleSoAFloatPrecision /*rebuildX*/, ParticleSoAFloatPrecision /*rebuildY*/, ParticleSoAFloatPrecision /*rebuildZ*/,
                                       ParticleSoAFloatPrecision /*vx*/, ParticleSoAFloatPrecision /*vy*/, ParticleSoAFloatPrecision /*vz*/, ParticleSoAFloatPrecision /*fx*/, ParticleSoAFloatPrecision /*fy*/,
                                       ParticleSoAFloatPrecision /*fz*/, ParticleSoAFloatPrecision /*oldFx*/, ParticleSoAFloatPrecision /*oldFy*/, ParticleSoAFloatPrecision /*oldFz*/,
                                       ParticleSoAFloatPrecision /*mass*/, size_t /*typeid*/, autopas::OwnershipState /*ownershipState*/>::Type;

  using KokkosSoAArraysType = autopas::utils::KokkosSoA<size_t* /*id*/, ParticleSoAFloatPrecision* /*x*/, ParticleSoAFloatPrecision* /*y*/, ParticleSoAFloatPrecision* /*z*/,
                                       ParticleSoAFloatPrecision* /*rebuildX*/, ParticleSoAFloatPrecision* /*rebuildY*/, ParticleSoAFloatPrecision* /*rebuildZ*/,
                                       ParticleSoAFloatPrecision* /*vx*/, ParticleSoAFloatPrecision* /*vy*/, ParticleSoAFloatPrecision* /*vz*/, ParticleSoAFloatPrecision* /*fx*/, ParticleSoAFloatPrecision* /*fy*/,
                                       ParticleSoAFloatPrecision* /*fz*/, ParticleSoAFloatPrecision* /*oldFx*/, ParticleSoAFloatPrecision* /*oldFy*/, ParticleSoAFloatPrecision* /*oldFz*/,
                                       ParticleSoAFloatPrecision* /*mass*/, size_t* /*typeid*/, autopas::OwnershipState* /*ownershipState*/>;

  /**
   * Non-const getter for the pointer of this object.
   * @tparam attribute Attribute name.
   * @return this.
   */

  template <AttributeNames attribute>
  constexpr auto& operator() () {
    auto value = get<attribute>();
    return value;
  }

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
  template <AttributeNames attribute, std::enable_if_t<attribute != ptr, bool> = true>
      constexpr std::tuple_element<attribute, SoAArraysType>::type::value_type get() {
    if constexpr (attribute == id) {
      return _id;
    } else if constexpr (attribute == posX) {
      return _r[0];
    } else if constexpr (attribute == posY) {
      return _r[1];
    } else if constexpr (attribute == posZ) {
      return _r[2];
    } else if constexpr (attribute == rebuildX) {
      return _rAtRebuild[0];
    } else if constexpr (attribute == rebuildY) {
      return _rAtRebuild[1];
    } else if constexpr (attribute == rebuildZ) {
      return _rAtRebuild[2];
    } else if constexpr (attribute == velocityX) {
      return _v[0];
    } else if constexpr (attribute == velocityY) {
      return _v[1];
    } else if constexpr (attribute == velocityZ) {
      return _v[2];
    } else if constexpr (attribute == forceX) {
      return _f[0];
    } else if constexpr (attribute == forceY) {
      return _f[1];
    } else if constexpr (attribute == forceZ) {
      return _f[2];
    } else if constexpr (attribute == oldForceX) {
      return _oldF[0];
    } else if constexpr (attribute == oldForceY) {
      return _oldF[1];
    } else if constexpr (attribute == oldForceZ) {
      return _oldF[2];
    } else if constexpr (attribute == typeId) {
      return _typeId;
    } else if constexpr (attribute == mass) {
      return _mass;
    } else if constexpr (attribute == ownershipState) {
      return _ownershipState;
    } else {
      //autopas::utils::ExceptionHandler::exception("MoleculeLJ::get() unknown attribute {}", attribute);
    }
  }

  template <AttributeNames attribute, std::enable_if_t<attribute != ptr, bool> = true>
    constexpr std::tuple_element<attribute, SoAArraysType>::type::value_type get() const {
        if constexpr (attribute == id) {
            return _id;
        } else if constexpr (attribute == posX) {
            return _r[0];
        } else if constexpr (attribute == posY) {
            return _r[1];
        } else if constexpr (attribute == posZ) {
            return _r[2];
        } else if constexpr (attribute == rebuildX) {
            return _rAtRebuild[0];
        } else if constexpr (attribute == rebuildY) {
            return _rAtRebuild[1];
        } else if constexpr (attribute == rebuildZ) {
            return _rAtRebuild[2];
        } else if constexpr (attribute == velocityX) {
            return _v[0];
        } else if constexpr (attribute == velocityY) {
            return _v[1];
        } else if constexpr (attribute == velocityZ) {
            return _v[2];
        } else if constexpr (attribute == forceX) {
            return _f[0];
        } else if constexpr (attribute == forceY) {
          return _f[1];
        } else if constexpr (attribute == forceZ) {
            return _f[2];
        } else if constexpr (attribute == oldForceX) {
            return _oldF[0];
        } else if constexpr (attribute == oldForceY) {
            return _oldF[1];
        } else if constexpr (attribute == oldForceZ) {
            return _oldF[2];
        } else if constexpr (attribute == typeId) {
            return _typeId;
        } else if constexpr (attribute == mass) {
            return _mass;
        } else if constexpr (attribute == ownershipState) {
            return _ownershipState;
        } else {
            // autopas::utils::ExceptionHandler::exception("MoleculeLJ::get() unknown attribute {}", attribute);
        }
    }

    template <AttributeNames attribute>
    constexpr void set(std::tuple_element<attribute, SoAArraysType>::type::value_type value) {
        if constexpr (attribute == id) {
            _id = value;
        } else if constexpr (attribute == posX) {
            _r[0] = value;
        } else if constexpr (attribute == posY) {
            _r[1] = value;
        } else if constexpr (attribute == posZ) {
            _r[2] = value;
        } else if constexpr (attribute == rebuildX) {
            _rAtRebuild[0] = value;
        } else if constexpr (attribute == rebuildY) {
            _rAtRebuild[1] = value;
        } else if constexpr (attribute == rebuildZ) {
            _rAtRebuild[2] = value;
        } else if constexpr (attribute == velocityX) {
            _v[0] = value;
        } else if constexpr (attribute == velocityY) {
            _v[1] = value;
        } else if constexpr (attribute == velocityZ) {
            _v[2] = value;
        } else if constexpr (attribute == forceX) {
            _f[0] = value;
        } else if constexpr (attribute == forceY) {
            _f[1] = value;
        } else if constexpr (attribute == forceZ) {
            _f[2] = value;
        } else if constexpr (attribute == oldForceX) {
            _oldF[0] = value;
        } else if constexpr (attribute == oldForceY) {
            _oldF[1] = value;
        } else if constexpr (attribute == oldForceZ) {
            _oldF[2] = value;
        } else if constexpr (attribute == typeId) {
            _typeId = value;
        } else if constexpr (attribute == mass) {
            _mass = value;
        } else if constexpr (attribute == ownershipState) {
           _ownershipState = value;
        } else {
            // autopas::utils::ExceptionHandler::exception("MoleculeLJ::set() unknown attribute {}", attribute);
        }
    }

  /**
   * Get the old force.
   * @return
   */
  [[nodiscard]] const std::array<ParticleSoAFloatPrecision, 3> &getOldF() const;

  /**
   * Set old force.
   * @param oldForce
   */
  void setOldF(const std::array<ParticleSoAFloatPrecision, 3> &oldForce);

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
   * Creates a string containing all data of the particle.
   * @return String representation.
   */
  [[nodiscard]] std::string toString() const override;

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
  std::array<ParticleSoAFloatPrecision, 3> _oldF = {0., 0., 0.};
};

}  // namespace mdLib
