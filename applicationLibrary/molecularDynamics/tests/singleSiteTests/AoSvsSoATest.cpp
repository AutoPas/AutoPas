/**
 * @file AoSvsSoATest.cpp
 * @author F.Gratl
 * @date 8.02.18
 */

#include "AoSvsSoATest.h"

#include "molecularDynamicsLibrary/LJFunctor.h"

#define PARTICLES_PER_DIM 16

/**
 * Generates a reproducible set of particles
 * @param particles Vector where particles will be stored.
 */
void AoSvsSoATest::generateParticles(std::vector<Molecule> *particles) {
  particles->resize(PARTICLES_PER_DIM * PARTICLES_PER_DIM);

  for (unsigned int i = 0; i < PARTICLES_PER_DIM; ++i) {
    for (unsigned int j = 0; j < PARTICLES_PER_DIM; ++j) {
      particles->at(i * PARTICLES_PER_DIM + j).setID(i * PARTICLES_PER_DIM + j);
      particles->at(i * PARTICLES_PER_DIM + j).setR({(double)i, (double)j, 0});
      particles->at(i * PARTICLES_PER_DIM + j).setF({0, 0, 0});
      particles->at(i * PARTICLES_PER_DIM + j).setV({0, 0, 0});
    }
  }
}

/**
 * Compares the result of the LJFunctor on a AoS and an SoA which contain the
 * same set of particles.
 */
TEST_F(AoSvsSoATest, testAoSvsSoA) {
  auto particlesAoS = std::vector<Molecule>();
  generateParticles(&particlesAoS);
  auto particlesSoA = particlesAoS;
  mdLib::LJFunctor<Molecule> ljFunctor(PARTICLES_PER_DIM * 10);
  ljFunctor.setParticleProperties(1., 1.);
  // AoS
  std::chrono::high_resolution_clock::time_point start, stop;
  start = std::chrono::high_resolution_clock::now();
  for (unsigned int i = 0; i < PARTICLES_PER_DIM * PARTICLES_PER_DIM; ++i) {
    for (unsigned int j = i + 1; j < PARTICLES_PER_DIM * PARTICLES_PER_DIM; ++j) {
      if (i != j) {
        ljFunctor.AoSFunctor(particlesAoS[i], particlesAoS[j], true);
      }
    }
  }
  stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

  std::cout << "AoS : " << duration << " \u03bcs" << std::endl;

  // SoA
  FMCell cell;
  for (auto &&p : particlesSoA) {
    cell.addParticle(p);
  }

  ljFunctor.SoALoader(cell, cell._particleSoABuffer, 0);
  start = std::chrono::high_resolution_clock::now();
  ljFunctor.SoAFunctorSingle(cell._particleSoABuffer, true);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

  std::cout << "SoA : " << duration << " \u03bcs" << std::endl;

  // copy back to particle array
  particlesSoA.clear();

  ljFunctor.SoAExtractor(cell, cell._particleSoABuffer, 0);

  //  ASSERT_EQ(particlesAoS.size(), particlesSoA.size());
  ASSERT_EQ(particlesAoS.size(), cell.numParticles());

  // compare particle vectors
  for (unsigned int i = 0; i < particlesAoS.size(); ++i) {
    ASSERT_NEAR(particlesAoS[i].getF()[0], cell._particles[i].getF()[0], 1.0e-13);
    ASSERT_NEAR(particlesAoS[i].getF()[1], cell._particles[i].getF()[1], 1.0e-13);
    ASSERT_NEAR(particlesAoS[i].getF()[2], cell._particles[i].getF()[2], 1.0e-13);
  }
}