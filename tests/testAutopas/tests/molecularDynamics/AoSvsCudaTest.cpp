/**
 * @file AoSvsCudaTest.cpp
 * @author jspahl
 * @date 11.02.19
 */

#if defined(AUTOPAS_CUDA)

#include "AoSvsCudaTest.h"

#include "autopas/molecularDynamics/LJFunctor.h"
#include "testingHelpers/commonTypedefs.h"

using namespace autopas;

constexpr size_t PARTICLES_PER_DIM = 16;

/**
 * @brief Generates a reproducible set of particles
 * @param particles Vector where particles will be stored.
 */
void AoSvsCudaTest::generateParticles(std::vector<Molecule> *particles, bool withDeletions) {
  particles->resize(PARTICLES_PER_DIM * PARTICLES_PER_DIM);

  double deletionPercentage = 30;
  constexpr unsigned seed = 42;
  std::mt19937 generator(seed);
  std::uniform_real_distribution<double> uniform0_100(0.0, 100.0);

  for (unsigned int i = 0; i < PARTICLES_PER_DIM; ++i) {
    for (unsigned int j = 0; j < PARTICLES_PER_DIM; ++j) {
      particles->at(i * PARTICLES_PER_DIM + j).setID(i * PARTICLES_PER_DIM + j);
      particles->at(i * PARTICLES_PER_DIM + j).setR({(double)i, (double)j, 0});
      particles->at(i * PARTICLES_PER_DIM + j).setF({0, 0, 0});
      particles->at(i * PARTICLES_PER_DIM + j).setV({0, 0, 0});
      if (withDeletions and uniform0_100(generator) < deletionPercentage) {
        particles->at(i * PARTICLES_PER_DIM + j).markAsDeleted();
      }
    }
  }
}

/**
 * Compares the result of the LJFunctor on a AoS and Cuda which contain the
 * same set of particles.
 */
TEST_P(AoSvsCudaTest, testAoSvsCuda) {
  auto [withDeletions] = GetParam();

  auto particlesAoS = std::vector<Molecule>();
  generateParticles(&particlesAoS, withDeletions);
  auto particlesSoA = particlesAoS;
  double epsilon = 1.0;
  double sigma = 1.0;
  LJFunctor<Molecule, FMCell> ljFunctor(PARTICLES_PER_DIM * 10);
  ljFunctor.setParticleProperties(epsilon * 24.0, sigma * sigma);

  // AoS
  std::chrono::high_resolution_clock::time_point start, stop;
  start = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < PARTICLES_PER_DIM * PARTICLES_PER_DIM; ++i) {
    for (size_t j = i + 1; j < PARTICLES_PER_DIM * PARTICLES_PER_DIM; ++j) {
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
  ljFunctor.deviceSoALoader(cell._particleSoABuffer, cell._particleSoABufferDevice);

  start = std::chrono::high_resolution_clock::now();
  ljFunctor.CudaFunctor(cell._particleSoABufferDevice, false);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

  std::cout << "SoA : " << duration << " \u03bcs" << std::endl;

  // copy back to particle array
  particlesSoA.clear();

  ljFunctor.deviceSoAExtractor(cell._particleSoABuffer, cell._particleSoABufferDevice);
  ljFunctor.SoAExtractor(cell, cell._particleSoABuffer, 0);

  //  ASSERT_EQ(particlesAoS.size(), particlesSoA.size());
  ASSERT_EQ(particlesAoS.size(), cell.numParticles());

  // compare particle vectors
  for (size_t i = 0; i < particlesAoS.size(); ++i) {
    ASSERT_NEAR(particlesAoS[i].getF()[0], cell._particles[i].getF()[0], 1.0e-13);
    ASSERT_NEAR(particlesAoS[i].getF()[1], cell._particles[i].getF()[1], 1.0e-13);
    ASSERT_NEAR(particlesAoS[i].getF()[2], cell._particles[i].getF()[2], 1.0e-13);
  }
}

static auto toString = [](const auto &info) {
  auto [deleteSomeParticles] = info.param;
  std::stringstream resStream;
  resStream << (deleteSomeParticles ? "withDeletions" : "withoutDeletions");
  std::string res = resStream.str();
  std::replace(res.begin(), res.end(), '-', '_');
  std::replace(res.begin(), res.end(), '.', '_');
  return res;
};

INSTANTIATE_TEST_SUITE_P(Generated, AoSvsCudaTest, ::testing::Combine(::testing::Bool()), toString);

#endif  // AUTOPAS_CUDA
