/**
 * @file NonConstructibleParticle.h
 * @author F. Gratl
 * @date 14/07/2020
 */

#pragma once

#include "autopas/particles/Particle.h"

/**
 * A particle class without an actual constructor (only copy, etc.).
 */
class NonConstructibleParticle : public autopas::ParticleFP64 {
 public:
  /**
   * Default constructor is deleted.
   */
  NonConstructibleParticle() = delete;

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int { ptr, id, posX, posY, posZ, forceX, forceY, forceZ, ownershipState };

  /**
   * The type for the SoA storage.
   */
  using SoAArraysType =
      typename autopas::utils::SoAType<NonConstructibleParticle *, size_t /*id*/, double /*x*/, double /*y*/,
                                       double /*z*/, double /*fx*/, double /*fy*/, double /*fz*/,
                                       autopas::OwnershipState /*ownershipState*/>::Type;

  /**
   * Getter, which allows access to an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @return Value of the requested attribute.
   * @note The value of owned is return as floating point number (true = 1.0, false = 0.0).
   */
  template <AttributeNames attribute>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() {
    if constexpr (attribute == AttributeNames::ptr) {
      return this;
    } else if constexpr (attribute == AttributeNames::id) {
      return this->getID();
    } else if constexpr (attribute == AttributeNames::posX) {
      return this->getR()[0];
    } else if constexpr (attribute == AttributeNames::posY) {
      return this->getR()[1];
    } else if constexpr (attribute == AttributeNames::posZ) {
      return this->getR()[2];
    } else if constexpr (attribute == AttributeNames::forceX) {
      return this->getF()[0];
    } else if constexpr (attribute == AttributeNames::forceY) {
      return this->getF()[1];
    } else if constexpr (attribute == AttributeNames::forceZ) {
      return this->getF()[2];
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      autopas::utils::ExceptionHandler::exception("ParticleBase::get() unknown attribute {}", attribute);
    }
  }
};