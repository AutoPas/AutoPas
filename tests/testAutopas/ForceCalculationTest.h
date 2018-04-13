#ifndef AUTOPAS_LJFORCECALCULATIONTEST_H
#define AUTOPAS_LJFORCECALCULATIONTEST_H

#include <AutoPas.h>
#include <gtest/gtest.h>
#include <vector>

class ForceCalculationTest : public ::testing::Test {
 public:
  ForceCalculationTest() = default;

  ~ForceCalculationTest() override = default;

  void fillWithParticles(
      AutoPas<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>>
      &autoPas,
      std::vector<size_t> particlesPerDim, double spacing);
};

#endif  // AUTOPAS_LJFORCECALCULATIONTEST_H
