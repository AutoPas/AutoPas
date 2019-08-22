/**
 * @file PeriodicBoundariesTest.h
 * @author N. Fottner
 * @date 2/8/19
 */

#include "SimulationTest.h"

void SimulationTest::initFillWithParticles(
    autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autopas,
    std::array<size_t, 3> particlesPerDim, double particleSpacing, double cutoff) {
  double minimalBoxLength = cutoff + 0.2 /*default skin value*/;
  std::array<double, 3> boxmax = {std::max(particlesPerDim[0] * particleSpacing, minimalBoxLength),
                                  std::max(particlesPerDim[1] * particleSpacing, minimalBoxLength),
                                  std::max(particlesPerDim[2] * particleSpacing, minimalBoxLength)};
  std::array<double, 3> boxmin = {0., 0., 0.};
  autopas.setBoxMin(boxmin);
  autopas.setBoxMax(boxmax);
  autopas.setCutoff(cutoff);
  autopas.init();
  PrintableMolecule dummy;
  GridGenerator::fillWithParticles(autopas, particlesPerDim, 0, 0, dummy,
                                   {particleSpacing, particleSpacing, particleSpacing}, {0., 0., 0.});
}

void SimulationTest::VisualizeSmallSzenario(std::array<size_t, 3> particlesPerDim, double cutoff,
                                            double particleSpacing, double epsilon, double sigma, double mass,
                                            int iterations, double delta_t, const std::string &filename) {
  // initializes all members for a VisualizeSmallSzenario to test behaviour of Particles interacting with eachother
  auto autopas = autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>();
  initFillWithParticles(autopas, particlesPerDim, particleSpacing, cutoff);
  ParticlePropertiesLibrary<double, size_t> _particlePropertiesLibrary;
  _particlePropertiesLibrary.addType(0, epsilon, sigma, mass);
  auto functor = autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>>(
      cutoff, 0.0, _particlePropertiesLibrary);
  TimeDiscretization<decltype(autopas), decltype(_particlePropertiesLibrary)> _timeDiscretization(
      delta_t, _particlePropertiesLibrary);
  ASSERT_TRUE(autopas.getNumberOfParticles() == particlesPerDim[0] * particlesPerDim[1] * particlesPerDim[2]);
  // starting simulation
  for (int i = 0; i < iterations; i++) {
    writeVTKFile(autopas, i, filename);
    _timeDiscretization.CalculateX(autopas);
    autopas.iteratePairwise(&functor);
    _timeDiscretization.CalculateV(autopas);
  }
}

double SimulationTest::distanceBetween2Points(const std::array<double, 3> &iPos, const std::array<double, 3> &jPos) {
  return std::sqrt(std::pow((iPos[0] - jPos[0]), 2.0) + std::pow((iPos[1] + jPos[1]), 2.0) +
                   std::pow((iPos[2] + jPos[2]), 2.0));
}

void SimulationTest::initWithTwoParticlesWithDistance(
    autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autopas, double distance,
    double cutoff) {
  autopas.setBoxMin({0., 0., 0.});
  autopas.setBoxMax({5., 5., 5.});
  autopas.setCutoff(cutoff);
  autopas.init();
  PrintableMolecule p1({0., 0., 0.}, {0., 0., 0.}, 0);
  PrintableMolecule p2({distance, 0., 0.}, {0., 0., 0.}, 1);
  autopas.addParticle(p1);
  autopas.addParticle(p2);
}

TEST_F(SimulationTest, NoMovementBetweenParticles) {
  // tests Particle movement when particleSpacing = 1.12 * sigma and cutoff = 1.5 > particleSpacing
  auto autopas = autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>();
  double cutoff = 1.5;
  double epsilon = 1.0;
  double sigma = 1.0;
  double mass = 1.0;
  double delta_t = 0.002;
  ParticlePropertiesLibrary<double, size_t> _particlePropertiesLibrary;
  _particlePropertiesLibrary.addType(0, epsilon, sigma, mass);
  auto functor = autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>, true>(
      cutoff, 0.0, _particlePropertiesLibrary);
  TimeDiscretization<decltype(autopas), decltype(_particlePropertiesLibrary)> _timeDiscretization(
      delta_t, _particlePropertiesLibrary);
  // as particleSpacing=1.12 * sigma and cutoff > particleSpacing, testing if particlesPositions stay the same
  initWithTwoParticlesWithDistance(autopas, 1.12, 1.5);
  for (int i = 0; i < 10 /*number of iterations*/; i++) {
    //            std::cout << "ITERATION " << i << std::endl;
    //            for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    //                std::cout << iter->toString() << std::endl;
    //            }
    std::cout << "Iteration: " << i << std::endl;
    auto iter = autopas.begin();
    auto p1Position = iter->getR();
    ++iter;
    auto p2Position = iter->getR();
    ASSERT_GE(distanceBetween2Points(p1Position, p2Position), 1.12);
    _timeDiscretization.CalculateX(autopas);
    autopas.iteratePairwise(&functor);
    _timeDiscretization.CalculateV(autopas);
  }
}

// TEST_F(SimulationTest,smallSzenariosVTKVisualization) {
//    // epsilon=1.0,sigma=1.0,particleSpacing=1.12
//
//    //cutoff=1.
//    VisualizeSmallSzenario({2, 2, 1}, /*particlesPerDim*/
//                           1.0, /*cutoff*/
//                           1.12, /*particleSpacing*/
//                           1.0, /*epsilon*/
//                           1.0, /*sigma*/
//                           1.0, /*mass*/
//                           20, /*number of Iterations*/
//                           0.002, /*delta_t*/
//                           "P:2*2_Cutoff:1.0_Sigma:1.0_PS:1.12" /*filename*/ );
