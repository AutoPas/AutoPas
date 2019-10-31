#pragma once

#include <vector>
#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "Objects/CubeGauss.h"
#include "Objects/CubeGrid.h"
#include "Objects/CubeUniform.h"
#include "Objects/Sphere.h"
#include "PrintableMolecule.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayMath.h"
/**
 * Class for contructing a container and generating Objects and Shapes filled with Particles
 */
class Generator {
 public:
  /**
   * Generates a Cube filled with Particles with dimensions: @param particlesPerDim
   * @param autopas
   * @param particlesPerDim
   * @param particleSpacing
   */
  template <class Particle, class ParticleCell>
  static void CubeGrid(autopas::AutoPas<Particle, ParticleCell> &autopas, size_t typeId, size_t id,
                       const std::array<double, 3> &boxMin, const std::array<size_t, 3> &particlesPerDim,
                       double particleSpacing, const std::array<double, 3> &velocity);

  /**
   * Fills Autopas Object with Particles with Gauss distribution
   * @param autopas
   * @param boxLength
   * @param numParticles
   * @param distributionMean
   * @param distributionStdDev
   */
  template <class Particle, class ParticleCell>
  static void CubeGauss(autopas::AutoPas<Particle, ParticleCell> &autopas, size_t typeId, size_t id,
                        const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, size_t numParticles,
                        double distributionMean, double distributionStdDev, const std::array<double, 3> &velocity);

  /**
   * Fills Autopas Object randomly with Particles
   * @param autopas
   * @param boxLength
   * @param numParticles
   */
  template <class Particle, class ParticleCell>
  static void CubeRandom(autopas::AutoPas<Particle, ParticleCell> &autopas, size_t typeId, size_t id,
                         const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, size_t numParticles,
                         const std::array<double, 3> &velocity);

  /**
   * Generates a Sphere with @param radius number of Particles with initial @param velocity
   * @param Autopas
   * @param center
   * @param radius
   * @param velocity
   * @param particleSpacing
   * @param id
   */
  template <class Particle, class ParticleCell>
  static void Sphere(autopas::AutoPas<Particle, ParticleCell> &autopas, const std::array<double, 3> &center, int radius,
                     double particleSpacing, unsigned long id, unsigned long typeId,
                     const std::array<double, 3> &velocity = {0., 0., 0.});
};

template <class Particle, class ParticleCell>
void Generator::CubeGrid(autopas::AutoPas<Particle, ParticleCell> &autopas, size_t typeId, size_t id,
                         const std::array<double, 3> &boxMin, const std::array<size_t, 3> &particlesPerDim,
                         double particleSpacing, const std::array<double, 3> &velocity) {
  Particle dummyParticle;
  dummyParticle.setV(velocity);
  dummyParticle.setID(id);
  dummyParticle.setTypeId(typeId);
  GridGenerator::fillWithParticles(autopas, particlesPerDim, dummyParticle,
                                   {particleSpacing, particleSpacing, particleSpacing}, boxMin);
}

template <class Particle, class ParticleCell>
void Generator::CubeGauss(autopas::AutoPas<Particle, ParticleCell> &autopas, size_t typeId, size_t id,
                          const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, size_t numParticles,
                          double distributionMean, double distributionStdDev, const std::array<double, 3> &velocity) {
  Particle dummyParticle;
  dummyParticle.setV(velocity);
  dummyParticle.setID(id);
  dummyParticle.setTypeId(typeId);
  GaussianGenerator::fillWithParticles(autopas, boxMin, boxMax, numParticles, dummyParticle, distributionMean,
                                       distributionStdDev);
}

template <class Particle, class ParticleCell>
void Generator::CubeRandom(autopas::AutoPas<Particle, ParticleCell> &autopas, size_t typeId, size_t id,
                           const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                           size_t numParticles, const std::array<double, 3> &velocity) {
  Particle dummyParticle;
  dummyParticle.setV(velocity);
  dummyParticle.setTypeId(typeId);
  dummyParticle.setID(id);
  RandomGenerator::fillWithParticles(autopas, dummyParticle, boxMin, boxMax, numParticles);
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
              double disCheck = sqrt(autopas::ArrayMath::dot(autopas::ArrayMath::sub(posVector, center),
                                                             autopas::ArrayMath::sub(posVector, center)));
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