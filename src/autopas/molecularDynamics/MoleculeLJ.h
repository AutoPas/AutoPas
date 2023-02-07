/**
 * @file MoleculeLJ.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <vector>

#include "autopas/particles/Particle.h"
#include <jlcxx/jlcxx.hpp>

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
      MoleculeLJ<floatType> *, size_t /*id*/, floatType /*x*/, floatType /*y*/, floatType /*z*/, floatType /*vx*/,
      floatType /*vy*/, floatType /*vz*/, floatType /*fx*/, floatType /*fy*/, floatType /*fz*/, floatType /*oldFx*/,
      floatType /*oldFy*/, floatType /*oldFz*/, size_t /*typeid*/, OwnershipState /*ownershipState*/>::Type;

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
   * adding construcor, getters, setters and add methods for Julia Usage
   */
   
  /**
   * Constructor of lennard jones molecule from julia
   * @param pos Position of the molecule.
   * @param v Velocitiy of the molecule.
   * @param moleculeId Id of the molecule.
   * @param typeId TypeId of the molecule.
   */
   /*
  MoleculeLJ(jlcxx::ArrayRef<double,1> pos, jlcxx::ArrayRef<double,1> v, int moleculeId,
                      int typeId = 0) {
                        std::array<double,3> pos_;
                        std::array<double,3> v_;

                        for(auto i = 0; i < pos.size(); i++) {
                          pos_.at(i) = pos[i];
                        }

                        for(auto i = 0; i < v.size(); i++) {
                          v_.at(i) = v[i];
                        }

                        Particle(pos_, v_, moleculeId);
                        _typeId = typeId;
                    }

  void setPos(jlcxx::ArrayRef<double,1> pos_) {
    ParticleBase::setR({pos_[0], pos_[1], pos_[2]});
  }

  void setV(jlcxx::ArrayRef<double,1> v_) {
    ParticleBase::setV({v_[0], v_[1], v_[2]});
  }

  void setF(jlcxx::ArrayRef<double,1> f_) {
    ParticleBase::setF({f_[0], f_[1], f_[2]});
  }

  void setOldF(jlcxx::ArrayRef<double,1> oldF_) {
    setOldF({oldF_[0], oldF_[1], oldF_[2]});
  }

  jlcxx::ArrayRef<double,1> getPos() {
    // jlcxx::ArrayRef<double,1> pos_{getR().data(), getR().size()};
    return {_r.data(), _r.size()};
  }

  jlcxx::ArrayRef<double,1> getV() {
    return {_v.data(), _v.size()};
  }

  jlcxx::ArrayRef<double,1> getF() {
    return {_f.data(), _f.size()};
  }

  jlcxx::ArrayRef<double,1> getOldF() {
    return {_oldF.data(), _oldF.size()};
  }

  void addPos(jlcxx::ArrayRef<double,1> pos_) {
    ParticleBase::addR({pos_[0], pos_[1], pos_[2]});
  }

  void addV(jlcxx::ArrayRef<double,1> v_) {
    ParticleBase::addV({v_[0], v_[1], v_[2]});
  }

  void addF(jlcxx::ArrayRef<double,1> f_) {
    ParticleBase::addF({f_[0], f_[1], f_[2]});
  }
  */
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
