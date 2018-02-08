#include "AoSvsSoATest.h"

using namespace autopas;

#define PARTICLES_PER_DIM 8

/**
 * @brief Generates a reproducable set of particles
 * @param particles Vector where particles will be stored.
 */
void AoSvsSoATest::generateParticles(
    std::vector<autopas::Particle> *particles) {
  particles->resize(PARTICLES_PER_DIM * PARTICLES_PER_DIM);

  for (unsigned int i = 0; i < PARTICLES_PER_DIM; ++i) {
    for (unsigned int j = 0; j < PARTICLES_PER_DIM; ++j) {
      particles->at(i * PARTICLES_PER_DIM + j).setID(i * PARTICLES_PER_DIM + j);
      particles->at(i * PARTICLES_PER_DIM + j).setR({(double) i, (double) j, 0});
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

  LJFunctor<autopas::Particle> ljFunctor;
  ljFunctor.setGlobals(PARTICLES_PER_DIM * 10, 1, 1, 0);

  // AoS
  for (unsigned int i = 0; i < PARTICLES_PER_DIM * PARTICLES_PER_DIM; ++i) {
    for (unsigned int j = 0; j < PARTICLES_PER_DIM * PARTICLES_PER_DIM; ++j) {
      if (i != j) {
        ljFunctor.AoSFunctor(particlesAoS[i], particlesAoS[j]);
      }
    }
  }

  std::cout << std::endl << std::endl;
  // SoA
  SoA soa1;

  ljFunctor.SoALoader(particlesSoA, &soa1);
  ljFunctor.SoAFunctor(soa1, soa1);

  // copy back to particle array
  particlesSoA.clear();

  ljFunctor.SoAExtractor(&particlesSoA, &soa1);

  ASSERT_EQ(particlesAoS.size(), particlesSoA.size());

  // compare particle vectors
  for (unsigned int i = 0; i < particlesAoS.size(); ++i) {
    ASSERT_DOUBLE_EQ(particlesAoS[i].getF()[0], particlesSoA[i].getF()[0]);
    ASSERT_DOUBLE_EQ(particlesAoS[i].getF()[1], particlesSoA[i].getF()[1]);
    ASSERT_DOUBLE_EQ(particlesAoS[i].getF()[2], particlesSoA[i].getF()[2]);
  }
}