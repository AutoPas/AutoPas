/**
 * @file Objects.h
 * @author N. Fottner
 * @date 1/8/19
 */
#pragma once

#include <array>
#include <iomanip>
#include <iosfwd>
#include <vector>

#include "autopas/AutoPasDecl.h"
#include "autopas/utils/ArrayUtils.h"
#include "src/ParticleAttributes.h"
#include "src/TypeDefinitions.h"

/**
 * Base class for describing objects made of particles.
 */
class Object {
 public:
  /**
   * Type of all particles generated.
   */
  using ParticleType = ::ParticleType;

  /**
   * Constructor that should be used by inheriting types.
   * @param velocity
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass
   */
  Object(const std::array<double, 3> &velocity, unsigned long typeId, double epsilon, double sigma, double mass)
      : _velocity(velocity), _typeId(typeId), _epsilon(epsilon), _sigma(sigma), _mass(mass) {}

  virtual ~Object() = default;

  /**
   * Generate the object in the given AutoPas container.
   * @param particles The container to which the new particles will be appended to.
   */
  virtual void generate(std::vector<ParticleAttributes> &particles) const = 0;

  /**
   * Create a particle that acts as blueprint for all particles to be created for the object.
   * @param autopas
   * @return
   */
  [[nodiscard]] ParticleAttributes getDummyParticle(const size_t &particleId) const {
    ParticleAttributes particle;
    particle.id = particleId;
    particle.typeId = _typeId;
    particle.ownershipState = autopas::OwnershipState::owned;
    particle.velocity = _velocity;
    particle.force = {0.0, 0.0, 0.0};
    particle.oldForce = {0.0, 0.0, 0.0};

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
    output << std::setw(_valueOffset) << std::left << "velocity"
           << ":  " << autopas::utils::ArrayUtils::to_string(_velocity) << std::endl;
    output << std::setw(_valueOffset) << std::left << "particle-type"
           << ":  " << _typeId << std::endl;
    output << std::setw(_valueOffset) << std::left << "particle-epsilon"
           << ":  " << _epsilon << std::endl;
    output << std::setw(_valueOffset) << std::left << "particle-sigma"
           << ":  " << _sigma << std::endl;
    output << std::setw(_valueOffset) << std::left << "particle-mass"
           << ":  " << _mass << std::endl;
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
   * Type of every particle in the object.
   */
  unsigned long _typeId;
  /**
   * Epsilon of every particle in the object.
   */
  double _epsilon;
  /**
   * Sigma of every particle in the object.
   */
  double _sigma;
  /**
   * Mass of every particle in the object.
   */
  double _mass;
  /**
   * valueOffset of MDFlexConfig - expected indent
   */
  static constexpr size_t _valueOffset = 33 - 6;
};
