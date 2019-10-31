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
  static void cubeGrid(autopas::AutoPas<Particle, ParticleCell> &autopas, const CubeGrid &object);

  /**
   * Fills Autopas Object with Particles with Gauss distribution
   * @param autopas
   * @param boxLength
   * @param numParticles
   * @param distributionMean
   * @param distributionStdDev
   */
  template <class Particle, class ParticleCell>
  static void cubeGauss(autopas::AutoPas<Particle, ParticleCell> &autopas, const CubeGauss &object);

  /**
   * Fills Autopas Object randomly with Particles
   * @param autopas
   * @param boxLength
   * @param numParticles
   */
  template <class Particle, class ParticleCell>
  static void cubeRandom(autopas::AutoPas<Particle, ParticleCell> &autopas, const CubeUniform &object);

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
  static void sphere(autopas::AutoPas<Particle, ParticleCell> &autopas, const Sphere &object);
};

template <class Particle, class ParticleCell>
void Generator::cubeGrid(autopas::AutoPas<Particle, ParticleCell> &autopas, const CubeGrid &object) {
  Particle dummyParticle;
  dummyParticle.setV(object.getVelocity());
  dummyParticle.setID(autopas.getNumberOfParticles());
  dummyParticle.setTypeId(object.getTypeId());
  GridGenerator::fillWithParticles(
      autopas, object.getParticlesPerDim(), dummyParticle,
      {object.getParticleSpacing(), object.getParticleSpacing(), object.getParticleSpacing()}, object.getBoxMin());
}

template <class Particle, class ParticleCell>
void Generator::cubeGauss(autopas::AutoPas<Particle, ParticleCell> &autopas, const CubeGauss &object) {
  Particle dummyParticle;
  dummyParticle.setV(object.getVelocity());
  dummyParticle.setID(autopas.getNumberOfParticles());
  dummyParticle.setTypeId(object.getTypeId());
  GaussianGenerator::fillWithParticles(autopas, object.getBoxMin(), object.getBoxMax(), object.getParticlesTotal(),
                                       dummyParticle, object.getDistributionMean(), object.getDistributionMean());
}

template <class Particle, class ParticleCell>
void Generator::cubeRandom(autopas::AutoPas<Particle, ParticleCell> &autopas, const CubeUniform &object) {
  Particle dummyParticle;
  dummyParticle.setV(object.getVelocity());
  dummyParticle.setTypeId(object.getTypeId());
  dummyParticle.setID(autopas.getNumberOfParticles());
  RandomGenerator::fillWithParticles(autopas, dummyParticle, object.getBoxMin(), object.getBoxMax(),
                                     object.getParticlesTotal());
}

template <class Particle, class ParticleCell>
void Generator::sphere(autopas::AutoPas<Particle, ParticleCell> &autopas, const Sphere &object) {
  Particle dummyParticle({0, 0, 0}, object.getVelocity(), autopas.getNumberOfParticles(), object.getTypeId());
  for (int z = 0; z <= object.getRadius(); ++z) {      // generate circles along the z-axis; uses symmetry of sphere
    for (int y = 0; y <= object.getRadius(); ++y) {    // generate lines among the y-axis
      for (int x = 0; x <= object.getRadius(); ++x) {  // generate particles among the x-axis
        std::array<double, 3> posDelta = {(double)x, (double)y, (double)z};           // offset of center as array
        for (int i = -1; i <= 1; i += 2) {                                            // mirror x-coordinate
          for (int k = -1; k <= 1; k += 2) {                                          // mirror y-coordinate
            for (int l = -1; l <= 1; l += 2) {                                        // mirror z-coordinate
              std::array<double, 3> multipliers = {(double)i, (double)k, (double)l};  // multipliers for mirroring
              std::array<double, 3> posVector = autopas::ArrayMath::add(
                  object.getCenter(),
                  autopas::ArrayMath::mulScalar(autopas::ArrayMath::mul(posDelta, multipliers),
                                                object.getParticleSpacing()));  // actual coordinates of new particle
              double disCheck = sqrt(autopas::ArrayMath::dot(autopas::ArrayMath::sub(posVector, object.getCenter()),
                                                             autopas::ArrayMath::sub(posVector, object.getCenter())));
              if (disCheck <= (double)(object.getRadius() + 1) * object.getParticleSpacing()) {
                dummyParticle.setR(posVector);
                autopas.addParticle(dummyParticle);
                dummyParticle.setID(dummyParticle.getID() + 1);
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