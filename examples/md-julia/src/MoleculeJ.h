#pragma once
#include "autopas/particles/Particle.h"
#include <jlcxx/jlcxx.hpp>

/**
 * Molecule class for md simulation in Julia
 */

template<typename floatType>
class MoleculeJ : public autopas::Particle {
    public:

    typedef floatType ft;

    /**
     * default constructor
     */
    MoleculeJ() {}

    /**
     * Constructor 
     * @param position position of the molecule
     * @param velocity velocity of the molecule
     * @param moleculeId Id ot the molecule
     * @param typeId indicating the type of the molecule
     */
    MoleculeJ(jlcxx::ArrayRef<floatType,1> pos, jlcxx::ArrayRef<floatType,1> v,
        unsigned long moleculeId, unsigned long typeId = 0) : _typeId(typeId) {
        // : Particle(pos, v, moleculeId),
        std::array<floatType, 3> pos_{pos[0], pos[1], pos[2]};
        std::array<floatType, 3> v_{v[0], v[1], v[2]};
        setPosition(pos);
        setVelocity(v);
        setID(moleculeId);
        // autopas::Particle(pos_, v_, moleculeId);
    }
    
    MoleculeJ(std::array<floatType, 3> pos, std::array<floatType, 3> v, unsigned long moleculeId,
                      unsigned long typeId = 0)
      : autopas::Particle(pos, v, moleculeId), _typeId(typeId) {}


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
      MoleculeJ<floatType> *, size_t /*id*/, floatType /*x*/, floatType /*y*/, floatType /*z*/, floatType /*vx*/,
      floatType /*vy*/, floatType /*vz*/, floatType /*fx*/, floatType /*fy*/, floatType /*fz*/, floatType /*oldFx*/,
      floatType /*oldFy*/, floatType /*oldFz*/, size_t /*typeid*/, autopas::OwnershipState /*ownershipState*/>::Type;
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
      autopas::utils::ExceptionHandler::exception("MoleculeLJ::get() unknown attribute {}", attribute);
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
      autopas::utils::ExceptionHandler::exception("MoleculeLJ::set() unknown attribute {}", attribute);
    }
  }

    /**
     * set position of the molecule
     */
    void setPosition(jlcxx::ArrayRef<floatType,1> pos) {
        ParticleBase::setR({pos[0], pos[1], pos[2]});
    }

    /**
     * set velocity of the molecule
     */
    void setVelocity(jlcxx::ArrayRef<floatType,1> v) {
        ParticleBase::setV({v[0], v[1], v[2]});
    }

    /**
     * set force of the molecule
     */
    void setForce(jlcxx::ArrayRef<floatType,1> f) {
        ParticleBase::setF({f[0], f[1], f[2]});
    }

    /**
     * set oldForce of the molecule
     */
    void setOldF(jlcxx::ArrayRef<floatType,1> oldF) {
        _oldF = {oldF[0], oldF[1], oldF[2]};
    }

    /**
     * get position of the molecule
     */
    jlcxx::ArrayRef<double,1> getPosition() {
        return {_r.data(), _r.size()};
    }

    /**
     * get velocity of the molecule
     */
    jlcxx::ArrayRef<double,1> getVelocity() {
        return {_v.data(), _v.size()};
    }

    /**
     * get force of the molecule
     */
    jlcxx::ArrayRef<double,1> getForce() {
        return {_f.data(), _f.size()};
    }

    /**
     * get oldForce of the molecule
     */
    jlcxx::ArrayRef<floatType,1> getOldF() {
        return {_oldF.data(), _oldF.size()};
    }

    size_t getTypeId() const { return _typeId;}

    /**
     * add position to current position of the molecule
     */
    void addPosition(jlcxx::ArrayRef<floatType,1> pos) {
        ParticleBase::addR({pos[0], pos[1], pos[2]});
    }

    /**
     * add velocity to current velocity of the molecule
     */
    void addVelocity(jlcxx::ArrayRef<floatType,1> v) {
        ParticleBase::addV({v[0], v[1], v[2]});
    }

    /**
     * add force to current force of the molecule
     */
    void addForce(jlcxx::ArrayRef<floatType,1> f) {
        ParticleBase::addF({f[0], f[1], f[2]});
    }

    /**
     * substract force from current force of the molecule
     */
    void subForce(jlcxx::ArrayRef<floatType,1> f) {
        ParticleBase::subF({f[0], f[1], f[2]});
    }

    private:
    /**
     * Particle type id.
     */
    size_t _typeId = 0;
    
    /**
     * Old Force of the particle experiences as 3D vector.
     */
    std::array<floatType,3> _oldF = {0.0, 0.0, 0.0};
};