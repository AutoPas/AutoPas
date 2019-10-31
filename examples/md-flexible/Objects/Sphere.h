/**
 * @file Sphere.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include "Objects.h"

class Sphere : public Object {
 public:
  /**
   * Constructor.
   * @param velocity
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass
   * @param center
   * @param radius
   * @param particleSpacing
   */
  Sphere(const std::array<double, 3> &velocity, unsigned long typeId, double epsilon, double sigma, double mass,
         const std::array<double, 3> &center, int radius, double particleSpacing)
      : Object(velocity, typeId, epsilon, sigma, mass),
        center(center),
        radius(radius),
        particleSpacing(particleSpacing){}

                /**
                 * Getter for center of Sphere
                 * @return center
                 */
                [[nodiscard]] const std::array<double, 3> &
            getCenter() const {
    return center;
  }

  /**
   * Getter for radius in number of Particles of Sphere
   * @return radius
   */
  [[nodiscard]] int getRadius() const { return radius; }

      /**
       * Getter for particleSpacing
       * @return particleSpacing
       */
      [[nodiscard]] double getParticleSpacing() const {
    return particleSpacing;
  }

  /**
   * Returns the number of particles that will be generated for this Object
   * @return totalNumberOfParticles
   */
  [[nodiscard]] size_t getParticlesTotal() const override {
    // this should look different if the generator for spheres changes
    int counter = 0;
    for (int z = 0; z <= radius; ++z) {
      for (int y = 0; y <= radius; ++y) {
        for (int x = 0; x <= radius; ++x) {
          std::array<double, 3> posDelta = {(double)x, (double)y, (double)z};
          for (int i = -1; i <= 1; i += 2) {
            for (int k = -1; k <= 1; k += 2) {
              for (int l = -1; l <= 1; l += 2) {
                std::array<double, 3> multipliers = {(double)i, (double)k, (double)l};
                std::array<double, 3> posVector = autopas::ArrayMath::add(
                    center,
                    autopas::ArrayMath::mulScalar(autopas::ArrayMath::mul(posDelta, multipliers), particleSpacing));
                double disCheck = autopas::ArrayMath::dot(autopas::ArrayMath::sub(posVector, center),
                                                          autopas::ArrayMath::sub(posVector, center));
                if (disCheck <= (double)(radius + 1) * particleSpacing) {
                  counter++;
                }
                if (z == 0) break;
              }
              if (y == 0) break;
            }
            if (x == 0) break;
          }
        }
      }
    }
    return counter;
  }

  /**
   * Getter for the smallest x,y,z coordinates for Object
   * @return BoxMin of Cube
   */
  const std::array<double, 3> getBoxMin() const override {
    return {center[0] - ((double)radius) * particleSpacing, center[1] - ((double)radius) * particleSpacing,
            center[2] - ((double)radius) * particleSpacing};
  }

  /**
   * Getter for the highest x,y,z coordinates for Object
   * @return BoxMax of Cube
   */
  const std::array<double, 3> getBoxMax() const override {
    return {center[0] + ((double)radius) * particleSpacing, center[1] + ((double)radius) * particleSpacing,
            center[2] + ((double)radius) * particleSpacing};
  }

  /**
   * Prints the Configuration of the current Object
   */
  std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "Center of Sphere"
           << ":  " << autopas::ArrayUtils::to_string(center) << std::endl;
    output << std::setw(_valueOffset) << std::left << "radius in Particles"
           << ":  " << radius << std::endl;
    output << std::setw(_valueOffset) << std::left << "particleSpacing"
           << ":  " << particleSpacing << std::endl;
    output << std::setw(_valueOffset) << std::left << "NumberOfParticles"
           << ":  " << this->getParticlesTotal() << std::endl;
    output << Object::to_string();
    return output.str();
  }

 private:
  /*
   * coordinates of the sphere's center
   */
  std::array<double, 3> center;
  /**
   * radius of the sphere in number of particles
   */
  int radius;
  double particleSpacing;
};