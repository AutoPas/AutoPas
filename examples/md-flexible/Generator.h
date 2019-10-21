#pragma once

#include <vector>
#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "PrintableMolecule.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayMath.h"
/**Class for contructing a container and generating Objects and Shapes filled with Particles
 * */
class Generator {
 public:
  /**Calculated the L2Norm of an array
   * @param array
   * @return L2Norm(std::array<double,3>)
   * */
  static double L2Norm(std::array<double, 3> array) {
    double square_sum = 0;
    for (auto e : array) {
      square_sum += (e * e);
    }
    return sqrt(square_sum);
  }

  /**Generates a Cube filled with Particles with dimensions: @param particlesPerDim
   * @param autopas
   * @param particlesPerDim
   * @param particleSpacing
   * */
  template <class Particle, class ParticleCell>
  static void CubeGrid(autopas::AutoPas<Particle, ParticleCell> &autopas, size_t typeId, size_t id,
                       const std::array<double, 3> &BoxMin, const std::array<size_t, 3> &particlesPerDim,
                       double particleSpacing, const std::array<double, 3> &velocity);

  /**Fills Autopas Object with Particles with Gauss distribution
   * @param autopas
   * @param boxLength
   * @param numParticles
   * @param distributionMean
   * @param distributionStdDev
   * */
  template <class Particle, class ParticleCell>
  static void CubeGauss(autopas::AutoPas<Particle, ParticleCell> &autopas, size_t typeId, size_t id,
                        const std::array<double, 3> &BoxMin, const std::array<double, 3> &BoxMax, size_t numParticles,
                        double distributionMean, double distributionStdDev, const std::array<double, 3> &velocity);

  /**Fills Autopas Object randomly with Particles
   * @param autopas
   * @param boxLength
   * @param numParticles
   * */
  template <class Particle, class ParticleCell>
  static void CubeRandom(autopas::AutoPas<Particle, ParticleCell> &autopas, size_t typeId, size_t id,
                         const std::array<double, 3> &BoxMin, const std::array<double, 3> &BoxMax, size_t numParticles,
                         const std::array<double, 3> &velocity);

  /**Generates a Sphere with @param radius number of Particles with initial @param velocity
   * @param Autopas
   * @param center
   * @param radius
   * @param velocity
   * @param particleSpacing
   * @param id
   * */
  template <class Particle, class ParticleCell>
  static void Sphere(autopas::AutoPas<Particle, ParticleCell> &autopas, const std::array<double, 3> &center, int radius,
                     double particleSpacing, unsigned long id, unsigned long typeId,
                     const std::array<double, 3> &velocity = {0., 0., 0.});
};

template <class Particle, class ParticleCell>
void Generator::CubeGrid(autopas::AutoPas<Particle, ParticleCell> &autopas, size_t typeId, size_t id,
                         const std::array<double, 3> &BoxMin, const std::array<size_t, 3> &particlesPerDim,
                         double particleSpacing, const std::array<double, 3> &velocity) {
  Particle dummyParticle;
  dummyParticle.setV(velocity);
  GridGenerator::fillWithParticles(autopas, particlesPerDim, typeId, id, dummyParticle,
                                   {particleSpacing, particleSpacing, particleSpacing}, BoxMin);
}

template <class Particle, class ParticleCell>
void Generator::CubeGauss(autopas::AutoPas<Particle, ParticleCell> &autopas, size_t typeId, size_t id,
                          const std::array<double, 3> &BoxMin, const std::array<double, 3> &BoxMax, size_t numParticles,
                          double distributionMean, double distributionStdDev, const std::array<double, 3> &velocity) {
  Particle dummyParticle;

  GaussianGenerator::fillWithParticles(autopas, typeId, id, BoxMin, BoxMax, numParticles, dummyParticle,
                                       distributionMean, distributionStdDev);
}

template <class Particle, class ParticleCell>
void Generator::CubeRandom(autopas::AutoPas<Particle, ParticleCell> &autopas, size_t typeId, size_t id,
                           const std::array<double, 3> &BoxMin, const std::array<double, 3> &BoxMax,
                           size_t numParticles, const std::array<double, 3> &velocity) {
  Particle dummyParticle;
  dummyParticle.setV(velocity);
  RandomGenerator::fillWithParticles(autopas, typeId, id, dummyParticle, BoxMin, BoxMax, numParticles);
}

template <class Particle, class ParticleCell>
void Generator::Sphere(autopas::AutoPas<Particle, ParticleCell> &autopas, const std::array<double, 3> &center,
                       int radius, double particleSpacing, unsigned long id, unsigned long typeId,
                       const std::array<double, 3> &velocity) {
  for (int z = 0; z <= radius; ++z) {      // generate circles along the z-axis; uses symmetry of sphere
    for (int y = 0; y <= radius; ++y) {    // generate lines among the y-axis
      for (int x = 0; x <= radius; ++x) {  // generate particles among the x-axis
        std::array<double, 3> posDelta = {(double)x, (double)y, (double)z};           // offset of center as array
        for (int i = -1; i <= 1; i += 2) {                                            // mirror x-coordinate
          for (int k = -1; k <= 1; k += 2) {                                          // mirror y-coordinate
            for (int l = -1; l <= 1; l += 2) {                                        // mirror z-coordinate
              std::array<double, 3> multipliers = {(double)i, (double)k, (double)l};  // multipliers for mirroring
              std::array<double, 3> posVector = autopas::ArrayMath::add(
                  center, autopas::ArrayMath::mulScalar(autopas::ArrayMath::mul(posDelta, multipliers),
                                                        particleSpacing));  // actual coordinates of new particle
              double disCheck = L2Norm(autopas::ArrayMath::sub(posVector, center));
              if (disCheck <= (double)(radius + 1) * particleSpacing) {
                Particle p(posVector, velocity, id, typeId);
                autopas.addParticle(p);
                id++;
              }
              if (z == 0)  // prevent duplicates
                break;
            }
            if (y == 0)  // prevent duplicates
              break;
          }
          if (x == 0)  // prevent duplicates
            break;
        }
      }
    }
  }
}