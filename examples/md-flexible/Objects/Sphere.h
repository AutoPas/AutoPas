/**
 * @file Sphere.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

class Sphere : public Object {
 public:
  /**
   * Constructor for Sphere that is created in YamlParser and then included into the Simulation via the Generator class
   * @param center
   * @param radius as number of particles
   * @param particleSpacing
   * @param velocity_arg
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass
   */
  Sphere(const std::array<double, 3> &center, int radius, double particleSpacing,
         const std::array<double, 3> &velocity_arg, const unsigned long &typeId_arg, const double &epsilon_arg,
         const double &sigma_arg, const double &mass_arg)
      : center(center), radius(radius), particleSpacing(particleSpacing) {
    velocity = velocity_arg;
    typeId = typeId_arg;
    epsilon = epsilon_arg;
    sigma = sigma_arg;
    mass = mass_arg;
  }

  /**
   * Getter for center of Sphere
   * @return center
   */
  [[nodiscard]] const std::array<double, 3> &getCenter() const { return center; }

      /**
       * Getter for radius in number of Particles of Sphere
       * @return radius
       */
      [[nodiscard]] int getRadius() const {
    return radius;
  }

  /**
   * Getter for particleSpacing
   * @return particleSpacing
   */
  [[nodiscard]] double getParticleSpacing() const { return particleSpacing; }

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
  void printConfig() override {
    using namespace std;
    cout << std::setw(valueOffset) << left << "Center of Sphere"
         << ":  " << autopas::ArrayUtils::to_string(center) << endl;
    cout << std::setw(valueOffset) << left << "radius in Particles"
         << ":  " << radius << endl;
    cout << std::setw(valueOffset) << left << "particleSpacing"
         << ":  " << particleSpacing << endl;
    cout << std::setw(valueOffset) << left << "NumberOfParticles"
         << ":  " << this->getParticlesTotal() << endl;
    cout << std::setw(valueOffset) << left << "Initial velocities"
         << ":  " << autopas::ArrayUtils::to_string(velocity) << endl;
    cout << std::setw(valueOffset) << left << "Particle Properties in Object:" << endl;
    cout << setw(valueOffset) << left << "Particle TypeId"
         << ":  " << typeId << endl;
    cout << setw(valueOffset) << left << "Particles Epsilon"
         << ":  " << epsilon << endl;
    cout << setw(valueOffset) << left << "Particles Sigma"
         << ":  " << sigma << endl;
    cout << setw(valueOffset) << left << "Particles Mass"
         << ":  " << mass << endl
         << endl;
  }

 private:
  static constexpr size_t valueOffset = 32;
  std::array<double, 3> center;
  // radius of the sphere in number of particles
  int radius;
  double particleSpacing;
};