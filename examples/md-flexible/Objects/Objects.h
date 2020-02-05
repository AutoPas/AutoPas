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

#include "autopas/utils/ArrayUtils.h"

/**
 * Base class for describing objects made of particles.
 */
class Object {
 public:
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
  virtual const std::array<double, 3> getBoxMin() const = 0;

  /**
   * Getter for the highest x,y,z coordinates for Object
   * @return BoxMax of Cube
   */
  virtual const std::array<double, 3> getBoxMax() const = 0;

  /**
   * Returns the total amount of Particles in the Object
   * @return ParticlesTotal
   */
  [[nodiscard]] virtual size_t getParticlesTotal() const = 0;

  /**
   * String description string of the object.
   * @return multiline std::string
   */
  virtual std::string to_string() const {
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

  friend std::ostream &operator<<(std::ostream &os, const Object &object) {
    os << object.to_string();
    return os;
  }

 protected:
  std::array<double, 3> _velocity;
  unsigned long _typeId;
  double _epsilon;
  double _sigma;
  double _mass;

  // valueOffset of MDFlexConfig - expected indent
  static constexpr size_t _valueOffset = 33 - 6;
};