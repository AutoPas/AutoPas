#pragma once
#include <array>
#include <vector>
#include "Generator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"

/**
 * Contains all Objects with their properties and functionalities for their generation in class Generator.h
 * and information prints in the yamlParser class
 */

class Object {
 public:
  virtual ~Object() = default;

  /**
   * Getter for Velocity
   * @return velocity
   */
  [[nodiscard]] const std::array<double, 3> &getVelocity() const { return velocity; }

      /**
       * Getter for typeId of Particles in Objet
       * @return typeId
       */
      [[nodiscard]] unsigned long getTypeId() const {
    return typeId;
  }

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
   * Prints the configuration of the Object
   */
  virtual void printConfig() = 0;

 protected:
  std::array<double, 3> velocity{};
  unsigned long typeId{};
  double epsilon{};
  double sigma{};
  double mass{};
};