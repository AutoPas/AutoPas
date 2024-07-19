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
   * @param typeId If single-site mode, this is used to lookup the mixing parameters from the config object during
   * object generation. If multi-site, this is stored directly (and used to lookup these parameters during the
   * simulation).
   */
  Object(const std::array<double, 3> &velocity, unsigned long typeId) : _velocity(velocity), _typeId(typeId) {}

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
    particle.setOwnershipState(autopas::OwnershipState::owned);
    particle.setV(_velocity);
    particle.setF({0.0, 0.0, 0.0});
    particle.setOldF({0.0, 0.0, 0.0});
    particle.setEpsilon(_epsilon);
    particle.setSigma(_sigma);
    particle.setMass(_mass);
#if MD_FLEXIBLE_MODE == MULTISITE
    particle.setTypeId(_typeId);
    particle.setQuaternion({1.0, 0.0, 0.0, 0.0});  // todo: add option for this to be set randomly
    particle.setAngularVel({0.0, 0.0, 0.0});
    particle.setTorque({0.0, 0.0, 0.0});
#endif

    return particle;
  }

  /**
   * Getter for Velocity
   * @return velocity
   */
  [[nodiscard]] const std::array<double, 3> &getVelocity() const { return _velocity; }

  /**
   * Getter for the typeId of the particles in the object
   * @return typeId
   */
  [[nodiscard]] unsigned long getTypeId() const { return _typeId; }

  /**
   * Getter for the epsilon of the particles in the object.
   * @return epsilon
   */
  [[nodiscard]] double getEpsilon() const { return _epsilon; }

  /**
   * Getter for the sigma of the particles in the object.
   * @return sigma
   */
  [[nodiscard]] double getSigma() const { return _sigma; }

  /**
   * Getter for the mass of the particles in the object.
   * @return mass
   */
  [[nodiscard]] double getMass() const { return _mass; }

  /**
   * Setter for the epsilon of the particles in the object.
   * @param epsilon
   */
  void setEpsilon(const double epsilon) { _epsilon = epsilon; }

  /**
   * Setter for the sigma of the particles in the object.
   * @param sigma
   */
  void setSigma(const double &sigma) { _sigma = sigma; }

  /**
   * Setter for the mass of the particles in the object.
   * @param mass
   */
  void setMass(const double &mass) { _mass = mass; }

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
  [[nodiscard]] virtual size_t getParticlesTotal() const = 0;

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
    output << std::setw(_valueOffset) << std::left << "particle-type-id"
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
   * Type of every particle in the object. This is only used for multi-site simulations and this refers to the molId.
   */
  unsigned long _typeId;
  /**
   * Epsilon of every molecule. Only used for single-site simulations.
   */
  double _epsilon;
  /**
   * Sigma of every molecule. Only used for single-site simulations.
   */
  double _sigma;
  /**
   * Mass of every molecule. Only used for single-site simulations.
   */
  double _mass;
  /**
   * valueOffset of MDFlexConfig - expected indent
   */
  static constexpr size_t _valueOffset = 33 - 6;
};
