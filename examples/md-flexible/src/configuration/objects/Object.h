/**
 * @file Object.h
 * @author N. Fottner
 * @date 1/8/19
 */
#pragma once

#include <array>
#include <iomanip>
#include <iosfwd>
#include <vector>

#include "autopas/AutoPasDecl.h"
#include "autopas/molecularDynamics/MultisiteMoleculeLJ.h"
#include "autopas/utils/ArrayUtils.h"
#include "src/TypeDefinitions.h"

/**
 * Base class for describing objects made of particles.
 */
class Object {
 public:
  /**
   * Constructor that should be used by inheriting types.
   * @param velocity
   * @param typeId
   */
  Object(const std::array<double, 3> &velocity, unsigned long typeId)
      : _velocity(velocity), _typeId(typeId){}

  virtual ~Object() = default;

  /**
   * Generate the object in the given AutoPas container.
   * @param particles The container to which the new particles will be appended to.
   */
  virtual void generate(std::vector<ParticleType> &particles) const = 0;

  /**
   * Create a particle that acts as blueprint for all particles to be created for the object.
   * @param particleId: Defines the id of the generated dummy particle.
   * @return a particle initialized with default values.
   */
  [[nodiscard]] ParticleType getDummyParticle(const size_t &particleId) const {
    ParticleType particle{};
    particle.setID(particleId);
    particle.setTypeId(_typeId);
    particle.setOwnershipState(autopas::OwnershipState::owned);
    particle.setV(_velocity);
    particle.setF({0.0, 0.0, 0.0});
    particle.setOldF({0.0, 0.0, 0.0});
//#if defined(MD_FLEXIBLE_USE_MULTI_SITE)
    particle.setQ({1.0, 1.0, 0.0, 0.0});
    particle.setAngularVel({0.0, 0.0, 0.0});
    particle.setTorque({0.0, 0.0, 0.0});
//#endif

    return particle;
  }

  /**
   * Getter for Velocity
   * @return velocity
   */
  [[nodiscard]] const std::array<double, 3> &getVelocity() const { return _velocity; }

  /**
   * Getter for typeId of Particles in Objet
   * @return typeId
   */
  [[nodiscard]] unsigned long getTypeId() const { return _typeId; }

  /**
   * Getter for the smallest x,y,z coordinates for Object
   * @return BoxMin of Cube
   */
  [[nodiscard]] virtual std::array<double, 3> getBoxMin() const = 0;

  /**
   * Getter for the highest x,y,z coordinates for Object
   * @return BoxMax of Cube
   */
  [[nodiscard]] virtual std::array<double, 3> getBoxMax() const = 0;

  /**
   * Returns the total amount of Particles in the Object
   * @return ParticlesTotal
   */
  [[nodiscard]] virtual size_t getParticlesTotal() const {
    throw std::runtime_error("Objects::getParticlesTotal() not implemented.");
  };

  /**
   * Getter for ParticleSpacing.
   * Objects that are not based on a grid return 0 since this is the minimal guaranteed spacing.
   * @return particleSpacing
   */
  [[nodiscard]] virtual double getParticleSpacing() const { return 0; }

  /**
   * String description string of the object.
   * @return multiline std::string
   */
  [[nodiscard]] virtual std::string to_string() const {
    std::ostringstream output;
#if defined(MD_FLEXIBLE_USE_MULTI_SITE)
    const auto typeName = "molecule-type";
#else
    const auto typename = "site-type";
#endif
    output << std::setw(_valueOffset) << std::left << "velocity"
           << ":  " << autopas::utils::ArrayUtils::to_string(_velocity) << std::endl;
    output << std::setw(_valueOffset) << std::left << typeName
           << ":  " << _typeId << std::endl;
    return output.str();
  };

  /**
   * Stream operator
   * @param os
   * @param object
   * @return
   */
  friend std::ostream &operator<<(std::ostream &os, const Object &object) {
    os << object.to_string();
    return os;
  }

 protected:
  /**
   * Velocity of every particle in the object.
   */
  std::array<double, 3> _velocity;
  /**
   * Type of every particle in the object. For single-site simulations, this refers directly to the siteId. For multi-site
   * simulations, this refers to the molId.
   */
  unsigned long _typeId;
  /**
   * valueOffset of MDFlexConfig - expected indent
   */
  static constexpr size_t _valueOffset = 33 - 6;
};
