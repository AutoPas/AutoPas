#include "AoSvsSoATest.h"
#include <chrono>

using namespace autopas;

#define PARTICLES_PER_DIM 16

/**
 * @brief Generates a reproducible set of particles
 * @param particles Vector where particles will be stored.
 */
void AoSvsSoATest::generateParticles(
    std::vector<autopas::Particle> *particles) {
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
  auto particlesAoS = *(new std::vector<autopas::Particle>);
  generateParticles(&particlesAoS);
  auto particlesSoA = particlesAoS;

  LJFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>>
      ljFunctor;
  ljFunctor.setGlobals(PARTICLES_PER_DIM * 10, 1, 1, 0);

  // AoS
  std::chrono::high_resolution_clock::time_point start, stop;
  start = std::chrono::high_resolution_clock::now();
  for (unsigned int i = 0; i < PARTICLES_PER_DIM * PARTICLES_PER_DIM; ++i) {
    for (unsigned int j = 0; j < PARTICLES_PER_DIM * PARTICLES_PER_DIM; ++j) {
      if (i != j) {
        ljFunctor.AoSFunctor(particlesAoS[i], particlesAoS[j]);
      }
    }
  }
  stop = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start)
          .count();

  std::cout << "AoS : " << duration << " \u03bcs" << std::endl;

  // SoA
  autopas::FullParticleCell<autopas::Particle> cell;
  for (auto &&p : particlesSoA) {
    cell.addParticle(p);
  }

  ljFunctor.SoALoader(cell, &cell._particleSoABuffer);
  start = std::chrono::high_resolution_clock::now();
  ljFunctor.SoAFunctor(cell._particleSoABuffer, cell._particleSoABuffer);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start)
                 .count();

  std::cout << "SoA : " << duration << " \u03bcs" << std::endl;

  // copy back to particle array
  particlesSoA.clear();

  ljFunctor.SoAExtractor(&cell, &cell._particleSoABuffer);

  //  ASSERT_EQ(particlesAoS.size(), particlesSoA.size());
  ASSERT_EQ(particlesAoS.size(), cell.numParticles());

  // compare particle vectors
  for (unsigned int i = 0; i < particlesAoS.size(); ++i) {
    ASSERT_NEAR(particlesAoS[i].getF()[0], cell._particles[i].getF()[0],
                1.0e-13);
    ASSERT_NEAR(particlesAoS[i].getF()[1], cell._particles[i].getF()[1],
                1.0e-13);
    ASSERT_NEAR(particlesAoS[i].getF()[2], cell._particles[i].getF()[2],
                1.0e-13);
  }
}