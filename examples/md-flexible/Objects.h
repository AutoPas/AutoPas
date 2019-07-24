#pragma once
#include <array>
#include <vector>
#include "Generator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
class CubeGrid {
 public:
  CubeGrid(const std::array<size_t, 3> &particlesPerDim, double particleSpacing, const std::array<double, 3> &velocity,const std::array<double,3> &center)
      : particlesPerDim(particlesPerDim),
        particleSpacing(particleSpacing),
        velocity(velocity),
        particlesTotal(particlesPerDim[0] * particlesPerDim[1] * particlesPerDim[2]),center(center) {}

  const array<size_t, 3> &getParticlesPerDim() const { return particlesPerDim; }

  double getParticleSpacing() const { return particleSpacing; }

  const array<double, 3> &getVelocity() const { return velocity; }

  int getParticlesTotal() const { return particlesTotal; }

  std::array<double,3> getBoxMin(){
    return {center[0]-0.5*particlesPerDim[0]*particleSpacing,center[1]-0.5*particlesPerDim[1]*particleSpacing,center[1]-0.5*particlesPerDim[1]*particleSpacing};
  }
    std::array<double,3> getBoxMax(){
        return {center[0]+0.5*particlesPerDim[0]*particleSpacing,center[1]+0.5*particlesPerDim[1]*particleSpacing,center[1]+0.5*particlesPerDim[1]*particleSpacing};
    }


  void printConfig() {
    cout << setw(valueOffset) << left << "Particles per dimension"
         << ":  " << ArrayUtils::to_string(particlesPerDim) << endl;
    cout << setw(valueOffset) << left << "Particle spacing"
         << ":  " << particleSpacing << endl;
      cout << setw(valueOffset) << left << "Number of Particles"
           << ":  " << (particlesPerDim[0] * particlesPerDim[1] * particlesPerDim[2]) << endl;
    cout << setw(valueOffset) << left << "Initial velocities"
         << ":  " << ArrayUtils::to_string(velocity) << endl;
  }

 private:
  static constexpr size_t valueOffset = 32;
  std::array<size_t, 3> particlesPerDim;
  double particleSpacing;
  std::array<double, 3> velocity;
  int particlesTotal;
    std::array<double,3> center;

};

class CubeGauss {
 public:
  CubeGauss(size_t numParticles,const std::array<double, 3> &boxLength, double distributionMean,
            double distributionStdDev, const std::array<double, 3> &velocity,const std::array<double,3> &center)
      :
        numParticles(numParticles),
        boxLength(boxLength),
        distributionMean(distributionMean),
        distributionStdDev(distributionStdDev),
        velocity(velocity),center(center) {}


  size_t getNumParticles() const { return numParticles; }

  double getDistributionMean() const { return distributionMean; }

  double getDistributionStdDev() const { return distributionStdDev; }

  const array<double, 3> &getVelocity() const { return velocity; }

  std::array<double,3> getBoxMin(){
    return {center[0]-0.5*boxLength[0],center[1]-0.5*boxLength[1],center[2]-0.5*boxLength[2]};
  }
    std::array<double,3> getBoxMax(){
        return {center[0]+0.5*boxLength[0],center[1]+0.5*boxLength[1],center[2]+0.5*boxLength[2]};
    }

  void printConfig() {
    cout << setw(valueOffset) << left << "Distribution-Mean"
         << ":  " << distributionMean << endl;
    cout << setw(valueOffset) << left << "Distribution-StdDev"
         << ":  " << distributionStdDev << endl;
    cout << setw(valueOffset) << left << "NumberOfParticles"
         << ":  " << numParticles << endl;
    cout << setw(valueOffset) << left << "Initial velocities"
         << ":  " << ArrayUtils::to_string(velocity) << endl;
  }

 private:
  static constexpr size_t valueOffset = 32;
  size_t numParticles;
    std::array<double, 3> boxLength;
  double distributionMean;
  double distributionStdDev;
  std::array<double, 3> velocity;
  std::array<double,3> center;
};

class CubeUniform {
 public:
  CubeUniform(size_t numParticles, const std::array<double, 3> &boxLength,const std::array<double, 3> &velocity,const std::array<double,3> &center)
      : numParticles(numParticles),boxLength(boxLength), velocity(velocity),center(center) {}


  size_t getNumParticles() const { return numParticles; }

  const array<double, 3> &getVelocity() const { return velocity; }

    std::array<double,3> getBoxMin(){
        return {center[0]-0.5*boxLength[0],center[1]-0.5*boxLength[1],center[2]-0.5*boxLength[2]};
    }

    std::array<double,3> getBoxMax(){
        return {center[0]+0.5*boxLength[0],center[1]+0.5*boxLength[1],center[2]+0.5*boxLength[2]};
    }

  void printConfig() {
    cout << setw(valueOffset) << left << "Center"
         << ":  " << ArrayUtils::to_string(center) << endl;
    cout << setw(valueOffset) << left << "NumberOfParticles"
         << ":  " << numParticles << endl;
    cout << setw(valueOffset) << left << "Initial velocities"
         << ":  " << ArrayUtils::to_string(velocity) << endl;
  }

 private:
  static constexpr size_t valueOffset = 32;
  size_t numParticles;
    std::array<double, 3> boxLength;
  std::array<double, 3> velocity;
    std::array<double,3> center;

};
class Sphere {
 public:
  Sphere(const std::array<double, 3> &center, int radius, double particleSpacing, unsigned long id,
         const std::array<double, 3> &velocity)
      : center(center), radius(radius), particleSpacing(particleSpacing), id(id), velocity(velocity) {}

  const std::array<double, 3> &getCenter() const { return center; }

  int getRadius() const { return radius; }

  double getParticleSpacing() const { return particleSpacing; }

  unsigned long getId() const { return id; }

  const std::array<double, 3> &getVelocity() const { return velocity; }

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
                std::array<double, 3> posVector = ArrayMath::add(
                    center, ArrayMath::mulScalar(ArrayMath::mul(posDelta, multipliers), particleSpacing));
                double disCheck = Generator::L2Norm(ArrayMath::sub(posVector, center));
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

  std::array<double,3> getBoxMin(){
      return {center[0]-((double)radius)*particleSpacing,center[1]-((double)radius)*particleSpacing,center[2]-((double)radius)*particleSpacing};
  }
    std::array<double,3> getBoxMax(){
        return {center[0]+((double)radius)*particleSpacing,center[1]+((double)radius)*particleSpacing,center[2]+((double)radius)*particleSpacing};
    }

  void printConfig() {
    cout << setw(valueOffset) << left << "Center of Sphere"
         << ":  " << ArrayUtils::to_string(center) << endl;
    cout << setw(valueOffset) << left << "radius in Particles"
         << ":  " << radius << endl;
    cout << setw(valueOffset) << left << "particleSpacing"
         << ":  " << particleSpacing << endl;
    //        cout << setw(valueOffset) << left << "first Particle in Sphere"
    //             << ":  " << id << endl;
    cout << setw(valueOffset) << left << "NumberOfParticles"
         << ":  " << this->particlesTotal() << endl;
    cout << setw(valueOffset) << left << "Initial velocities"
         << ":  " << ArrayUtils::to_string(velocity) << endl;
  }

 private:
  static constexpr size_t valueOffset = 32;
  std::array<double, 3> center;
  int radius;
  double particleSpacing;
  unsigned long id;
  std::array<double, 3> velocity;
};