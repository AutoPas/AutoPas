/**
 * @file Sphere.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include "Object.h"
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
  Sphere(const std::array<double, 3> &velocity, unsigned long typeId, const std::array<double, 3> &center, int radius,
         double particleSpacing)
      : Object(velocity, typeId),
        _center(center),
        _radius(radius),
        _particleSpacing(particleSpacing) {}

  /**
   * Getter for center of Sphere
   * @return center
   */
  [[nodiscard]] const std::array<double, 3> &getCenter() const { return _center; }

  /**
   * Getter for radius in number of Particles of Sphere
   * @return radius
   */
  [[nodiscard]] int getRadius() const { return _radius; }

  /**
   * Returns the amount of space between each particle.
   * @return spacing of the particles.
   */
  [[nodiscard]] double getParticleSpacing() const override { return _particleSpacing; }

  /**
   * Call f for every point on the sphere where a particle should be.
   * @param f Function called for every point.
   */
  void iteratePositions(const std::function<void(std::array<double, 3>)> &f) const {
    using namespace autopas::utils::ArrayMath::literals;

    // generate regular grid for 1/8th of the sphere
    for (int z = 0; z <= _radius; ++z) {
      for (int y = 0; y <= _radius; ++y) {
        for (int x = 0; x <= _radius; ++x) {
          // position relative to the center
          const std::array<double, 3> relativePos = {(double)x, (double)y, (double)z};
          // mirror to rest of sphere
          for (int i = -1; i <= 1; i += 2) {
            for (int k = -1; k <= 1; k += 2) {
              for (int l = -1; l <= 1; l += 2) {
                const std::array<double, 3> mirrorMultipliers = {(double)i, (double)k, (double)l};
                // position mirrored, scaled and absolute
                const std::array<double, 3> posVector =
                    _center + ((relativePos * mirrorMultipliers) * _particleSpacing);

                double distFromCentersSquare = autopas::utils::ArrayMath::dot(posVector - _center, posVector - _center);
                const auto r = (_radius + 1) * _particleSpacing;
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

  /**
   * Returns the total amount of particles which will be / have been generated.
   * @return the total number of particles initialized for the sphere object.
   */
  [[nodiscard]] size_t getParticlesTotal() const override {
    // this should look different if the generator for spheres changes
    int counter = 0;
    iteratePositions([&](auto pos) { ++counter; });
    return counter;
  }

  /**
   * Returns the coordinates of box's the bottom left front corner.
   * @return the bottom left front corner of the sphere's box domain.
   */
  [[nodiscard]] std::array<double, 3> getBoxMin() const override {
    return {_center[0] - ((double)_radius) * _particleSpacing, _center[1] - ((double)_radius) * _particleSpacing,
            _center[2] - ((double)_radius) * _particleSpacing};
  }

  /**
   * Returns the coordinates of box's the top right back corner.
   * @return the top right back corner of the sphere's box domain.
   */
  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    return {_center[0] + ((double)_radius) * _particleSpacing, _center[1] + ((double)_radius) * _particleSpacing,
            _center[2] + ((double)_radius) * _particleSpacing};
  }

  /**
   * Converts the object to a human readable string.
   * @return a human readable string containing information about the sphere object.
   */
  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "center"
           << ":  " << autopas::utils::ArrayUtils::to_string(_center) << std::endl;
    output << std::setw(_valueOffset) << std::left << "radius"
           << ":  " << _radius << std::endl;
    output << std::setw(_valueOffset) << std::left << "particle-spacing"
           << ":  " << _particleSpacing << std::endl;
    output << Object::to_string();
    return output.str();
  }

  /**
   * Generates the particles based on the configuration of the sphere object provided in the yaml file.
   * @param particles The container where the generated particles will be stored.
   */
  void generate(std::vector<ParticleType> &particles) const override {
    ParticleType particle = getDummyParticle(particles.size());
    iteratePositions([&](const auto &pos) {
      particle.setR(pos);
      particles.push_back(particle);
      particle.setID(particle.getID() + 1);
    });
  }

 private:
  /**
   * coordinates of the sphere's center.
   */
  std::array<double, 3> _center;

  /**
   * radius of the sphere in number of particles.
   */
  int _radius;

  /**
   * The amount of space between each particle.
   */
  double _particleSpacing;
};
