#pragma once
#include <array>
#include <vector>
#include "Generator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
class CubeGrid {
 public:
  CubeGrid(const std::array<size_t, 3> &particlesPerDim, double particleSpacing, const std::array<double, 3> &velocity,
           const std::array<double, 3> &center)
      : particlesPerDim(particlesPerDim),
        particleSpacing(particleSpacing),
        velocity(velocity),
        particlesTotal(particlesPerDim[0] * particlesPerDim[1] * particlesPerDim[2]),
        center(center) {}
  /**Getter for ParticlesPerDim
   * @return particlePerDim
   * */
  const std::array<size_t, 3> &getParticlesPerDim() const { return particlesPerDim; }
  /**GEtter for ParticleSpacing
   * @return particleSpacing
   * */
  double getParticleSpacing() const { return particleSpacing; }
    /**Getter for Velocity
     * @return velocity
     * */
  const std::array<double, 3> &getVelocity() const { return velocity; }
    /**Getter for total number of Particles for object
     * @return particlesTotal
     * */
  int getParticlesTotal() const { return particlesTotal; }
    /**Getter for the smallest x,y,z coordinates for Object
     * @return BoxMin of Cube
     * */
  std::array<double, 3> getBoxMin() {
    return {center[0] - 0.5 * particlesPerDim[0] * particleSpacing,
            center[1] - 0.5 * particlesPerDim[1] * particleSpacing,
            center[1] - 0.5 * particlesPerDim[1] * particleSpacing};
  }
    /**Getter for the highest x,y,z coordinates for Object
      * @return BoxMax of Cube
      * */
  std::array<double, 3> getBoxMax() {
    return {center[0] + 0.5 * particlesPerDim[0] * particleSpacing,
            center[1] + 0.5 * particlesPerDim[1] * particleSpacing,
            center[1] + 0.5 * particlesPerDim[1] * particleSpacing};
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
         << ":  " << autopas::ArrayUtils::to_string(velocity) << endl
         << endl;
  }

 private:
  static constexpr size_t valueOffset = 32;
  std::array<size_t, 3> particlesPerDim;
  double particleSpacing;
  std::array<double, 3> velocity;
  int particlesTotal;
  std::array<double, 3> center;
};

class CubeGauss {
 public:
  CubeGauss(size_t numParticles, const std::array<double, 3> &boxLength, double distributionMean,
            double distributionStdDev, const std::array<double, 3> &velocity, const std::array<double, 3> &center)
      : numParticles(numParticles),
        boxLength(boxLength),
        distributionMean(distributionMean),
        distributionStdDev(distributionStdDev),
        velocity(velocity),
        center(center) {}
    /**Getter total number of Particles of Object
     * @return numParticles
     * */
  size_t getNumParticles() const { return numParticles; }
  /**Getter for distribution mean
   * @return distributionMean
   * */
  double getDistributionMean() const { return distributionMean; }
    /**Getter for distributionStdDev
     * @return distributionStdDev
     * */
  double getDistributionStdDev() const { return distributionStdDev; }
   /**Getter for velocities of the Particles generated
    * @return velocity
    * */
  const std::array<double, 3> &getVelocity() const { return velocity; }
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
         << ":  " << autopas::ArrayUtils::to_string(boxLength) << endl
         << endl;
    cout << std::setw(valueOffset) << left << "Initial velocities"
         << ":  " << autopas::ArrayUtils::to_string(velocity) << endl
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
};

class CubeUniform {
 public:
  CubeUniform(size_t numParticles, const std::array<double, 3> &boxLength, const std::array<double, 3> &velocity,
              const std::array<double, 3> &center)
      : numParticles(numParticles), boxLength(boxLength), velocity(velocity), center(center) {}

      /**Getter for total number of Particles in Object
       * @return numParticles
       * */
  size_t getNumParticles() const { return numParticles; }

  const std::array<double, 3> &getVelocity() const { return velocity; }
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

    cout << std::setw(valueOffset) << left << "Center"
         << ":  " << autopas::ArrayUtils::to_string(center) << endl;
    cout << std::setw(valueOffset) << left << "NumberOfParticles"
         << ":  " << numParticles << endl;
    cout << std::setw(valueOffset) << left << "BoxLength"
         << ":  " << autopas::ArrayUtils::to_string(boxLength) << endl
         << endl;
    cout << std::setw(valueOffset) << left << "Initial velocities"
         << ":  " << autopas::ArrayUtils::to_string(velocity) << endl
         << endl;
  }

 private:
  static constexpr size_t valueOffset = 32;
  size_t numParticles;
  std::array<double, 3> boxLength;
  std::array<double, 3> velocity;
  std::array<double, 3> center;
};
class Sphere {
 public:
  Sphere(const std::array<double, 3> &center, int radius, double particleSpacing, unsigned long id,
         const std::array<double, 3> &velocity)
      : center(center), radius(radius), particleSpacing(particleSpacing), id(id), velocity(velocity) {}
    /**Getter for center of Sphere
     * @return center
     * */
  const std::array<double, 3> &getCenter() const { return center; }
    /**Getter for radius of Sphere
     * @return radius
     * */
  int getRadius() const { return radius; }
    /**Getter for particleSpacing
     * @return particleSpacing
     * */
  double getParticleSpacing() const { return particleSpacing; }
    /**Getter for initial Id of SPhere(id of first Particle Generated during generation Phase)
     * @return id
     * */
  unsigned long getId() const { return id; }
    /**Getter for initial velocity of Particles
     * @return velocity
     * */
  const std::array<double, 3> &getVelocity() const { return velocity; }
  /**Returns the number of particles that will be generated for this Object
   * @return totalNumberOfParticles
   * */
  //@todo besser implementieren: (anderen Sphere Generator?)
  int particlesTotal() {
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
    //        cout << setw(valueOffset) << left << "first Particle in Sphere"
    //             << ":  " << id << endl;
    cout << std::setw(valueOffset) << left << "NumberOfParticles"
         << ":  " << this->particlesTotal() << endl;
    cout << std::setw(valueOffset) << left << "Initial velocities"
         << ":  " << autopas::ArrayUtils::to_string(velocity) << endl
         << endl;
  }

 private:
  static constexpr size_t valueOffset = 32;
  std::array<double, 3> center;
  int radius;
  double particleSpacing;
  unsigned long id;
  std::array<double, 3> velocity;
};