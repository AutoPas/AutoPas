/**
 * @file Sphere.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include "Objects.h"
#include "autopas/utils/ArrayMath.h"

/**
 * Class describing a regular 3D spherical particle grid object.
 */
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
        particleSpacing(particleSpacing) {}

  /**
   * Getter for center of Sphere
   * @return center
   */
  [[nodiscard]] const std::array<double, 3> &getCenter() const { return center; }

  /**
   * Getter for radius in number of Particles of Sphere
   * @return radius
   */
  [[nodiscard]] int getRadius() const { return radius; }

  [[nodiscard]] double getParticleSpacing() const override { return particleSpacing; }

  /**
   * Call f for every point on the sphere where a particle should be.
   * @param f Function called for every point.
   */
  void iteratePositions(const std::function<void(std::array<double, 3>)> &f) const {
    // generate regular grid for 1/8th of the sphere
    for (int z = 0; z <= radius; ++z) {
      for (int y = 0; y <= radius; ++y) {
        for (int x = 0; x <= radius; ++x) {
          // position relative to the center
          std::array<double, 3> relativePos = {(double)x, (double)y, (double)z};
          // mirror to rest of sphere
          for (int i = -1; i <= 1; i += 2) {
            for (int k = -1; k <= 1; k += 2) {
              for (int l = -1; l <= 1; l += 2) {
                std::array<double, 3> mirrorMultipliers = {(double)i, (double)k, (double)l};
                // position mirrored, scaled and absolute
                std::array<double, 3> posVector = autopas::utils::ArrayMath::add(
                    center, autopas::utils::ArrayMath::mulScalar(
                                autopas::utils::ArrayMath::mul(relativePos, mirrorMultipliers), particleSpacing));

                double distFromCentersSquare =
                    autopas::utils::ArrayMath::dot(autopas::utils::ArrayMath::sub(posVector, center),
                                                   autopas::utils::ArrayMath::sub(posVector, center));
                const auto r = (radius + 1) * particleSpacing;
                const auto rSquare = r * r;
                // since the loops create a cubic grid only apply f for positions inside the sphere
                if (distFromCentersSquare <= rSquare) {
                  f(posVector);
                }
                // avoid duplicates
                if (z == 0) break;
              }
              if (y == 0) break;
            }
            if (x == 0) break;
          }
        }
      }
    }
  }

  [[nodiscard]] size_t getParticlesTotal() const override {
    // this should look different if the generator for spheres changes
    int counter = 0;
    iteratePositions([&](auto pos) { ++counter; });
    return counter;
  }

  [[nodiscard]] std::array<double, 3> getBoxMin() const override {
    return {center[0] - ((double)radius) * particleSpacing, center[1] - ((double)radius) * particleSpacing,
            center[2] - ((double)radius) * particleSpacing};
  }

  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    return {center[0] + ((double)radius) * particleSpacing, center[1] + ((double)radius) * particleSpacing,
            center[2] + ((double)radius) * particleSpacing};
  }

  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "center"
           << ":  " << autopas::utils::ArrayUtils::to_string(center) << std::endl;
    output << std::setw(_valueOffset) << std::left << "radius"
           << ":  " << radius << std::endl;
    output << std::setw(_valueOffset) << std::left << "particle-spacing"
           << ":  " << particleSpacing << std::endl;
    output << Object::to_string();
    return output.str();
  }

  void generate(autopas::AutoPas<ParticleType> &autopas) const override {
    ParticleType dummyParticle = getDummyParticle(autopas);
    iteratePositions([&](auto pos) {
      dummyParticle.setR(pos);
      autopas.addParticle(dummyParticle);
      dummyParticle.setID(dummyParticle.getID() + 1);
    });
  }

 private:
  /**
   * coordinates of the sphere's center
   */
  std::array<double, 3> center;
  /**
   * radius of the sphere in number of particles
   */
  int radius;
  double particleSpacing;
};