#pragma once
#include <array>
#include <vector>
#include "Generator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"

/**Containes all Objects with their properties and functionalities for their generation in class Generator.h
 * and information prints in the yamlParser class
 * */

class CubeGrid {
 public:
  /**constructor for CubeGrid that is created in YamlParser and then included into the Simulation via the Generator
   * class
   * @param particlesPerDim
   * @param particleSpacing
   * @param velocity
   * @param center
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass
   * */
  CubeGrid(const std::array<size_t, 3> &particlesPerDim, double particleSpacing, const std::array<double, 3> &velocity,
           const std::array<double, 3> &center, const unsigned long &typeId, const double &epsilon, const double &sigma,
           const double &mass)
      : particlesPerDim(particlesPerDim),
        particleSpacing(particleSpacing),
        velocity(velocity),
        particlesTotal(particlesPerDim[0] * particlesPerDim[1] * particlesPerDim[2]),
        center(center),
        typeId(typeId),
        epsilon(epsilon),
        sigma(sigma),
        mass(mass){}

                /**Getter for ParticlesPerDim
                 * @return particlePerDim
                 * */
                [[nodiscard]] const std::array<size_t, 3> &
            getParticlesPerDim() const {
    return particlesPerDim;
  }
  /**GEtter for ParticleSpacing
   * @return particleSpacing
   * */
  [[nodiscard]] double getParticleSpacing() const { return particleSpacing; }
      /**Getter for Velocity
       * @return velocity
       * */
      [[nodiscard]] const std::array<double, 3> &getVelocity() const {
    return velocity;
  }
  /**Getter for total number of Particles for object
   * @return particlesTotal
   * */
  [[nodiscard]] int getParticlesTotal() const { return particlesTotal; }
  /**Getter for the smallest x,y,z coordinates for Object
   * @return BoxMin of Cube
   * */
  std::array<double, 3> getBoxMin() {
    return {center[0] - 0.5 * particlesPerDim[0] * particleSpacing,
            center[1] - 0.5 * particlesPerDim[1] * particleSpacing,
            center[2] - 0.5 * particlesPerDim[2] * particleSpacing};
  }
  /**Getter for the highest x,y,z coordinates for Object
   * @return BoxMax of Cube
   * */
  std::array<double, 3> getBoxMax() {
    return {center[0] + 0.5 * particlesPerDim[0] * particleSpacing,
            center[1] + 0.5 * particlesPerDim[1] * particleSpacing,
            center[2] + 0.5 * particlesPerDim[2] * particleSpacing};
  }
  /**Prints the Configuration of the current Object
   * */
  void printConfig() {
    using namespace std;

    cout << std::setw(valueOffset) << left << "Particles per dimension"
         << ":  " << autopas::ArrayUtils::to_string(particlesPerDim) << endl;
    cout << std::setw(valueOffset) << left << "Particle spacing"
         << ":  " << particleSpacing << endl;
    cout << std::setw(valueOffset) << left << "Number of Particles"
         << ":  " << (particlesPerDim[0] * particlesPerDim[1] * particlesPerDim[2]) << endl;
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
         << ":  " << mass << endl << endl;
  }

  /**Getter for typeId of Particles in Objet
   * @return typeId
   * */
  unsigned long getTypeId() const { return typeId; }

 private:
  static constexpr size_t valueOffset = 32;
  std::array<size_t, 3> particlesPerDim;
  double particleSpacing;
  std::array<double, 3> velocity;
  int particlesTotal;
  std::array<double, 3> center;
  unsigned long typeId;
  double epsilon;
  double sigma;
  double mass;
};

class CubeGauss {
 public:
  /**constructor for CubeGauss that is created in YamlParser and then included into the Simulation via the Generator
   * class
   * @param numParticles
   * @param boxLength
   * @param distributionMean
   * @param distributionStdDev
   * @param velocity
   * @param center
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass
   * */
  CubeGauss(size_t numParticles, const std::array<double, 3> &boxLength, double distributionMean,
            double distributionStdDev, const std::array<double, 3> &velocity, const std::array<double, 3> &center,
            const unsigned long &typeId, const double &epsilon, const double &sigma, const double &mass)
      : numParticles(numParticles),
        boxLength(boxLength),
        distributionMean(distributionMean),
        distributionStdDev(distributionStdDev),
        velocity(velocity),
        center(center),
        typeId(typeId),
        epsilon(epsilon),
        sigma(sigma),
        mass(mass){} /**Getter total number of Particles of Object
                      * @return numParticles
                      * */
            [[nodiscard]] size_t getParticlesTotal() const {
    return numParticles;
  }
  /**Getter for typeId of Particles in Objet
   * @return typeId
   * */
  unsigned long getTypeId() const { return typeId; }

  /**Getter for distribution mean
   * @return distributionMean
   * */
  double getDistributionMean() const { return distributionMean; }
  /**Getter for distributionStdDev
   * @return distributionStdDev
   * */
  [[nodiscard]] double getDistributionStdDev() const { return distributionStdDev; }
      /**Getter for velocities of the Particles generated
       * @return velocity
       * */
      [[nodiscard]] const std::array<double, 3> &getVelocity() const {
    return velocity;
  }
  /**Getter for the smallest x,y,z coordinates for Object
   * @return BoxMin of Cube
   * */
  std::array<double, 3> getBoxMin() {
    return {center[0] - 0.5 * boxLength[0], center[1] - 0.5 * boxLength[1], center[2] - 0.5 * boxLength[2]};
  }
  /**Getter for the highest x,y,z coordinates for Object
   * @return BoxMax of Cube
   * */
  std::array<double, 3> getBoxMax() {
    return {center[0] + 0.5 * boxLength[0], center[1] + 0.5 * boxLength[1], center[2] + 0.5 * boxLength[2]};
  }

  /**Prints the Configuration of the current Object
   * */
  void printConfig() {
    using namespace std;

    cout << std::setw(valueOffset) << left << "Distribution-Mean"
         << ":  " << distributionMean << endl;
    cout << std::setw(valueOffset) << left << "Distribution-StdDev"
         << ":  " << distributionStdDev << endl;
    cout << std::setw(valueOffset) << left << "NumberOfParticles"
         << ":  " << numParticles << endl;
    cout << std::setw(valueOffset) << left << "BoxLength"
         << ":  " << autopas::ArrayUtils::to_string(boxLength) << endl;
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
  size_t numParticles;
  std::array<double, 3> boxLength;
  double distributionMean;
  double distributionStdDev;
  std::array<double, 3> velocity;
  std::array<double, 3> center;
  unsigned long typeId;
  double epsilon;
  double sigma;
  double mass;
};

class CubeUniform {
 public:
  /**constructor for CubeUniform that is created in YamlParser and then included into the Simulation via the Generator
   * class
   * @param numParticles
   * @param boxLength
   * @param velocity
   * @param center
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass*/
  CubeUniform(size_t numParticles, const std::array<double, 3> &boxLength, const std::array<double, 3> &velocity,
              const std::array<double, 3> &center, const unsigned long &typeId, const double &epsilon,
              const double &sigma, const double &mass)
      : numParticles(numParticles),
        boxLength(boxLength),
        velocity(velocity),
        center(center),
        typeId(typeId),
        epsilon(epsilon),
        sigma(sigma),
        mass(mass){}

            /**Getter for total number of Particles in Object
             * @return numParticles
             * */
            [[nodiscard]] size_t getParticlesTotal() const {
    return numParticles;
  }

  [[nodiscard]] const std::array<double, 3> &getVelocity() const { return velocity; }
  /**Getter for the smallest x,y,z coordinates for Object
   * @return BoxMin of Cube
   * */
  std::array<double, 3> getBoxMin() {
    return {center[0] - 0.5 * boxLength[0], center[1] - 0.5 * boxLength[1], center[2] - 0.5 * boxLength[2]};
  }
  /**Getter for the highest x,y,z coordinates for Object
   * @return BoxMax of Cube
   * */
  std::array<double, 3> getBoxMax() {
    return {center[0] + 0.5 * boxLength[0], center[1] + 0.5 * boxLength[1], center[2] + 0.5 * boxLength[2]};
  }
  /**Getter for typeId of Particles in Objet
   * @return typeId
   * */
  unsigned long getTypeId() const { return typeId; }

  /**Prints the Configuration of the current Object
   * */
  void printConfig() {
    using namespace std;

    cout << std::setw(valueOffset) << left << "Center"
         << ":  " << autopas::ArrayUtils::to_string(center) << endl;
    cout << std::setw(valueOffset) << left << "NumberOfParticles"
         << ":  " << numParticles << endl;
    cout << std::setw(valueOffset) << left << "BoxLength"
         << ":  " << autopas::ArrayUtils::to_string(boxLength) << endl;
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
  size_t numParticles;
  std::array<double, 3> boxLength;
  std::array<double, 3> velocity;
  std::array<double, 3> center;
  unsigned long typeId;
  double epsilon;
  double sigma;
  double mass;
};
class Sphere {
 public:
  /**constructor for Sphere that is created in YamlParser and then included into the Simulation via the Generator class
   * @param center
   * @param radius
   * @param particleSpacing
   * @param velocity
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass*/
  Sphere(const std::array<double, 3> &center, int radius, double particleSpacing, const std::array<double, 3> &velocity,
         const unsigned long &typeId, const double &epsilon, const double &sigma, const double &mass)
      : center(center),
        radius(radius),
        particleSpacing(particleSpacing),
        velocity(velocity),
        typeId(typeId),
        epsilon(epsilon),
        sigma(sigma),
        mass(mass){} /**Getter for center of Sphere
                      * @return center
                      * */
                [[nodiscard]] const std::array<double, 3> &
            getCenter() const {
    return center;
  }
  /**Getter for radius of Sphere
   * @return radius
   * */
  [[nodiscard]] int getRadius() const { return radius; }
      /**Getter for particleSpacing
       * @return particleSpacing
       * */
      [[nodiscard]] double getParticleSpacing() const {
    return particleSpacing;
  }
  /**Getter for typeId of Particles in Objet
   * @return typeId
   * */
  unsigned long getTypeId() const { return typeId; }
  /**Getter for initial velocity of Particles
   * @return velocity
   * */
  [[nodiscard]] const std::array<double, 3> &getVelocity() const { return velocity; }
  /**Returns the number of particles that will be generated for this Object
   * @return totalNumberOfParticles
   * */
  int getParticlesTotal() {
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
                double disCheck = Generator::L2Norm(autopas::ArrayMath::sub(posVector, center));
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
  /**Getter for the smallest x,y,z coordinates for Object
   * @return BoxMin of Cube
   * */
  std::array<double, 3> getBoxMin() {
    return {center[0] - ((double)radius) * particleSpacing, center[1] - ((double)radius) * particleSpacing,
            center[2] - ((double)radius) * particleSpacing};
  }
  /**Getter for the highest x,y,z coordinates for Object
   * @return BoxMax of Cube
   * */
  std::array<double, 3> getBoxMax() {
    return {center[0] + ((double)radius) * particleSpacing, center[1] + ((double)radius) * particleSpacing,
            center[2] + ((double)radius) * particleSpacing};
  }
  /**Prints the Configuration of the current Object
   * */
  void printConfig() {
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
  int radius;
  double particleSpacing;
  std::array<double, 3> velocity;
  unsigned long typeId;
  double epsilon;
  double sigma;
  double mass;
};