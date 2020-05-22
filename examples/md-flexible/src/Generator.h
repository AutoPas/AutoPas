/**
 * @file Generator.h
 * @author N. Fottner
 * @date 1/8/19
 */

#pragma once

#include <vector>

#include "Objects/CubeGauss.h"
#include "Objects/CubeGrid.h"
#include "Objects/CubeUniform.h"
#include "Objects/Sphere.h"
#include "PrintableMolecule.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayMath.h"
#include "autopasTools/generators/GaussianGenerator.h"
#include "autopasTools/generators/GridGenerator.h"
#include "autopasTools/generators/RandomGenerator.h"
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
  autopasTools::generators::GridGenerator::fillWithParticles(
      autopas, object.getParticlesPerDim(), dummyParticle,
      {object.getParticleSpacing(), object.getParticleSpacing(), object.getParticleSpacing()}, object.getBoxMin());
}

template <class Particle, class ParticleCell>
void Generator::cubeGauss(autopas::AutoPas<Particle, ParticleCell> &autopas, const CubeGauss &object) {
  Particle dummyParticle;
  dummyParticle.setV(object.getVelocity());
  dummyParticle.setID(autopas.getNumberOfParticles());
  dummyParticle.setTypeId(object.getTypeId());
  autopasTools::generators::GaussianGenerator::fillWithParticles(
      autopas, object.getBoxMin(), object.getBoxMax(), object.getParticlesTotal(), dummyParticle,
      object.getDistributionMean(), object.getDistributionMean());
}

template <class Particle, class ParticleCell>
void Generator::cubeRandom(autopas::AutoPas<Particle, ParticleCell> &autopas, const CubeUniform &object) {
  Particle dummyParticle;
  dummyParticle.setV(object.getVelocity());
  dummyParticle.setTypeId(object.getTypeId());
  dummyParticle.setID(autopas.getNumberOfParticles());
  autopasTools::generators::RandomGenerator::fillWithParticles(autopas, dummyParticle, object.getBoxMin(),
                                                               object.getBoxMax(), object.getParticlesTotal());
}

template <class Particle, class ParticleCell>
void Generator::sphere(autopas::AutoPas<Particle, ParticleCell> &autopas, const Sphere &object) {
  Particle dummyParticle({0, 0, 0}, object.getVelocity(), autopas.getNumberOfParticles(), object.getTypeId());
  object.iteratePositions([&](auto pos) {
    dummyParticle.setR(pos);
    autopas.addParticle(dummyParticle);
    dummyParticle.setID(dummyParticle.getID() + 1);
  });
}