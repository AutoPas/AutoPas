/**
 * @file TouchableParticle.h
 * @author seckler
 * @date 13.01.20
 */

#pragma once
#include "autopas/particles/Particle.h"

/**
 * Class that extends the default autopas::Particle with a behavior to check how often the particle was touched.
 * This class is useful for all kinds of ParticleIterator tests.
 */
class TouchableParticle : public autopas::Particle {
 public:
  /**
   * Constructor with position and id.
   * @param pos position
   * @param id id of the particle
   */
  TouchableParticle(std::array<double, 3> pos, unsigned long id)
      : autopas::Particle(pos, {0, 0, 0}, id), _numTouched(0){};

  /**
   * Constructor with position, velocity, and id.
   * @param pos position
   * @param velocity velocity
   * @param id id of the particle
   */
  TouchableParticle(std::array<double, 3> pos, std::array<double, 3> velocity, unsigned long id)
      : autopas::Particle(pos, velocity, id), _numTouched(0){};

  /**
   * Default constructor
   */
  TouchableParticle() : TouchableParticle({0., 0., 0.}, 0ul) {}

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int { ptr, id, posX, posY, posZ, forceX, forceY, forceZ, typeId, ownershipState };

  /**
   * The type for the SoA storage.
   */
  using SoAArraysType =
      typename autopas::utils::SoAType<TouchableParticle *, size_t /*id*/, double /*x*/, double /*y*/, double /*z*/,
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
   * Touch the particle.
   * The number of times a particle was touched is saved.
   */
  void touch() { _numTouched++; }

  /**
   * Get the number that indicates how often the particle was touch.
   * @return returns how often the particle was touched.
   */
  unsigned int getNumTouched() const { return _numTouched; }

 private:
  unsigned int _numTouched;
};