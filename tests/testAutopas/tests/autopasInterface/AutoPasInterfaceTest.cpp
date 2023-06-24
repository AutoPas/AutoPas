/**
 * @file AutoPasInterfaceTest.cpp
 * @author seckler
 * @date 13.05.19
 */

#include "AutoPasInterfaceTest.h"

#include "autopas/AutoPasDecl.h"
#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/selectors/ContainerSelector.h"
#include "autopas/tuning/selectors/ContainerSelectorInfo.h"
#include "autopas/tuning/selectors/TraversalSelector.h"
#include "autopas/tuning/selectors/TraversalSelectorInfo.h"
#include "autopas/utils/StaticCellSelector.h"
#include "testingHelpers/NumThreadGuard.h"
#include "testingHelpers/commonTypedefs.h"

extern template class autopas::AutoPas<Molecule>;
using LJFunctorGlobals =
    autopas::LJFunctor<Molecule, /* shifting */ true, /*mixing*/ false, autopas::FunctorN3Modes::Both,
                       /*globals*/ true>;
extern template bool autopas::AutoPas<Molecule>::iteratePairwise(LJFunctorGlobals *);

constexpr double cutoff = 1.1;
constexpr double skinPerTimestep = 0.2;
constexpr unsigned int rebuildFrequency = 3;
constexpr std::array<double, 3> boxMin{0., 0., 0.};
constexpr std::array<double, 3> boxMax{10., 10., 10.};

constexpr std::array<double, 3> zeroArr = {0., 0., 0.};

template <typename AutoPasT>
void defaultInit(AutoPasT &autoPas) {
  autoPas.setBoxMin(boxMin);
  autoPas.setBoxMax(boxMax);
  autoPas.setCutoff(cutoff);
  autoPas.setVerletSkinPerTimestep(skinPerTimestep);
  autoPas.setVerletRebuildFrequency(rebuildFrequency);
  autoPas.setNumSamples(3);

  // init autopas
  autoPas.init();
}

template <typename AutoPasT>
void defaultInit(AutoPasT &autoPas1, AutoPasT &autoPas2, size_t direction) {
  autoPas1.setBoxMin(boxMin);
  autoPas2.setBoxMax(boxMax);

  auto midLow = boxMin, midHigh = boxMax;
  midLow[direction] = (boxMax[direction] + boxMin[direction]) / 2;
  midHigh[direction] = (boxMax[direction] + boxMin[direction]) / 2;
  autoPas1.setBoxMax(midHigh);
  autoPas2.setBoxMin(midLow);

  for (auto &aP : {&autoPas1, &autoPas2}) {
    aP->setCutoff(cutoff);
    aP->setVerletSkinPerTimestep(skinPerTimestep);
    aP->setVerletRebuildFrequency(2);
    aP->setNumSamples(2);
    // init autopas
    aP->init();
  }
}

/**
 * Convert the leaving particle to entering particles.
 * Hereby the periodic boundary position change is done.
 * @param leavingParticles
 * @return vector of particles that will enter the container.
 */
std::vector<Molecule> convertToEnteringParticles(const std::vector<Molecule> &leavingParticles) {
  std::vector<Molecule> enteringParticles{leavingParticles};
  for (auto &p : enteringParticles) {
    auto pos = p.getR();
    for (auto dim = 0; dim < 3; dim++) {
      if (pos[dim] < boxMin[dim]) {
        // has to be smaller than boxMax
        pos[dim] = std::min(std::nextafter(boxMax[dim], -1), pos[dim] + (boxMax[dim] - boxMin[dim]));
      } else if (pos[dim] >= boxMax[dim]) {
        // should at least be boxMin
        pos[dim] = std::max(boxMin[dim], pos[dim] - (boxMax[dim] - boxMin[dim]));
      }
    }
    p.setR(pos);
  }
  return enteringParticles;
}

/**
 * Identifies and sends particles that are in the halo of neighboring AutoPas instances or the same instance (periodic
 * boundaries).
 * @param autoPas
 * @return vector of particles that are already shifted for the next process.
 */
auto identifyAndSendHaloParticles(autopas::AutoPas<Molecule> &autoPas) {
  std::vector<Molecule> haloParticles;

  for (const int x : {-1, 0, 1}) {
    for (const int y : {-1, 0, 1}) {
      for (const int z : {-1, 0, 1}) {
        if (x == 0 and y == 0 and z == 0) continue;
        std::array<int, 3> direction{x, y, z};
        std::array<double, 3> min{}, max{}, shiftVec{};
        for (size_t dim = 0; dim < 3; ++dim) {
          // The search domain has to be enlarged as the position of the particles is not certain.
          bool needsShift = false;
          if (direction[dim] == -1) {
            min[dim] = autoPas.getBoxMin()[dim];
            max[dim] = autoPas.getBoxMin()[dim] + cutoff;
            if (autoPas.getBoxMin()[dim] == boxMin[dim]) {
              needsShift = true;
            }
          } else if (direction[dim] == 1) {
            min[dim] = autoPas.getBoxMax()[dim] - cutoff;
            max[dim] = autoPas.getBoxMax()[dim];
            if (autoPas.getBoxMax()[dim] == boxMax[dim]) {
              needsShift = true;
            }
          } else {  // 0
            min[dim] = autoPas.getBoxMin()[dim];
            max[dim] = autoPas.getBoxMax()[dim];
          }
          if (needsShift) {
            shiftVec[dim] = -(boxMax[dim] - boxMin[dim]) * direction[dim];
          } else {
            shiftVec[dim] = 0;
          }
        }
        // here it is important to only iterate over the owned particles!
        for (auto iter = autoPas.getRegionIterator(min, max, autopas::IteratorBehavior::owned); iter.isValid();
             ++iter) {
          auto particleCopy = *iter;
          particleCopy.addR(shiftVec);
          haloParticles.push_back(particleCopy);
        }
      }
    }
  }

  return haloParticles;
}

size_t addEnteringParticles(autopas::AutoPas<Molecule> &autoPas, const std::vector<Molecule> &enteringParticles) {
  size_t numAdded = 0;
  for (const auto &p : enteringParticles) {
    if (autopas::utils::inBox(p.getR(), autoPas.getBoxMin(), autoPas.getBoxMax())) {
      autoPas.addParticle(p);
      ++numAdded;
    }
  }
  return numAdded;
}

void addHaloParticles(autopas::AutoPas<Molecule> &autoPas, const std::vector<Molecule> &haloParticles) {
  for (const auto &p : haloParticles) {
    autoPas.addHaloParticle(p);
  }
}

template <typename Functor>
void doSimulationLoop(autopas::AutoPas<Molecule> &autoPas, Functor *functor) {
  // 1. update Container; return value is vector of invalid == leaving particles!
  auto invalidParticles = autoPas.updateContainer();

  // 2. leaving and entering particles
  const auto &sendLeavingParticles = invalidParticles;
  // 2b. get+add entering particles (addParticle)
  const auto &enteringParticles = convertToEnteringParticles(sendLeavingParticles);
  auto numAdded = addEnteringParticles(autoPas, enteringParticles);

  EXPECT_EQ(numAdded, enteringParticles.size());

  // 3. halo particles
  // 3a. identify and send inner particles that are in the halo of other autopas instances or itself.
  auto sendHaloParticles = identifyAndSendHaloParticles(autoPas);

  // 3b. get halo particles
  const auto &recvHaloParticles = sendHaloParticles;
  addHaloParticles(autoPas, recvHaloParticles);

  // 4. iteratePairwise
  autoPas.iteratePairwise(functor);
}

template <typename Functor>
void doSimulationLoop(autopas::AutoPas<Molecule> &autoPas1, autopas::AutoPas<Molecule> &autoPas2, Functor *functor1,
                      Functor *functor2) {
  // 1. update Container; return value is vector of invalid = leaving particles!
  auto invalidParticles1 = autoPas1.updateContainer();
  auto invalidParticles2 = autoPas2.updateContainer();

  // 2. leaving and entering particles
  const auto &sendLeavingParticles1 = invalidParticles1;
  const auto &sendLeavingParticles2 = invalidParticles2;
  // 2b. get+add entering particles (addParticle)
  const auto &enteringParticles2 = convertToEnteringParticles(sendLeavingParticles1);
  const auto &enteringParticles1 = convertToEnteringParticles(sendLeavingParticles2);

  // the particles may either still be in the same container (just going over periodic boundaries) or in the other.
  size_t numAdded = 0;
  numAdded += addEnteringParticles(autoPas1, enteringParticles1);
  numAdded += addEnteringParticles(autoPas1, enteringParticles2);
  numAdded += addEnteringParticles(autoPas2, enteringParticles1);
  numAdded += addEnteringParticles(autoPas2, enteringParticles2);

  ASSERT_EQ(numAdded, enteringParticles1.size() + enteringParticles2.size());

  // 3. halo particles
  // 3a. identify and send inner particles that are in the halo of other autopas instances or itself.
  auto sendHaloParticles1 = identifyAndSendHaloParticles(autoPas1);
  auto sendHaloParticles2 = identifyAndSendHaloParticles(autoPas2);

  // 3b. get halo particles
  const auto &recvHaloParticles2 = sendHaloParticles1;
  const auto &recvHaloParticles1 = sendHaloParticles2;
  addHaloParticles(autoPas1, recvHaloParticles1);
  addHaloParticles(autoPas2, recvHaloParticles2);

  // 4. iteratePairwise
  autoPas1.iteratePairwise(functor1);
  autoPas2.iteratePairwise(functor2);
}

template <typename Functor>
void doAssertions(autopas::AutoPas<Molecule> &autoPas, Functor *functor, unsigned long numParticlesExpected, int line) {
  std::vector<Molecule> molecules(numParticlesExpected);
  size_t numParticles = 0;
  for (auto iter = autoPas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    ASSERT_LT(numParticles, numParticlesExpected) << "Too many particles owned by this container." << std::endl
                                                  << "Called from line: " << line;
    molecules[numParticles++] = *iter;
  }
  ASSERT_EQ(numParticles, numParticlesExpected)
      << "The container should own exactly " << numParticlesExpected << " particles!" << std::endl
      << "Called from line: " << line;

  for (auto &mol : molecules) {
    EXPECT_NEAR(autopas::utils::ArrayMath::dot(mol.getF(), mol.getF()), 390144. * 390144., 1.)
        << "wrong force calculated for particle: " << mol.toString() << std::endl
        << "Called from line: " << line;
  }

  EXPECT_NEAR(functor->getUpot(), 16128.983372449373 * numParticles / 2., 1e-5) << "wrong upot calculated" << std::endl
                                                                                << "Called from line: " << line;
  EXPECT_NEAR(functor->getVirial(), 195072. * numParticles / 2., 1e-5) << "wrong virial calculated" << std::endl
                                                                       << "Called from line: " << line;
}

template <typename Functor>
void doAssertions(autopas::AutoPas<Molecule> &autoPas1, autopas::AutoPas<Molecule> &autoPas2, Functor *functor1,
                  Functor *functor2, int line) {
  std::array<Molecule, 2> molecules{};
  size_t numParticles = 0;
  for (auto iter = autoPas1.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    ASSERT_LT(numParticles, 2) << "Too many owned particles." << std::endl << "Called from line: " << line;
    molecules[numParticles++] = *iter;
  }
  for (auto iter = autoPas2.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    ASSERT_LT(numParticles, 2) << "Too many owned particles." << std::endl << "Called from line: " << line;
    molecules[numParticles++] = *iter;
  }
  ASSERT_EQ(numParticles, 2) << "There should be exactly two owned particles!" << std::endl
                             << "Called from line: " << line;

  for (const auto &mol : molecules) {
    EXPECT_DOUBLE_EQ(autopas::utils::ArrayMath::dot(mol.getF(), mol.getF()), 390144. * 390144)
        << "wrong force calculated.";
  }

  EXPECT_DOUBLE_EQ(functor1->getUpot() + functor2->getUpot(), 16128.983372449373)
      << "wrong upot calculated" << std::endl
      << "Called from line: " << line;
  EXPECT_DOUBLE_EQ(functor1->getVirial() + functor2->getVirial(), 195072.) << "wrong virial calculated" << std::endl
                                                                           << "Called from line: " << line;
}

void setFromOptions(const testingTuple &options, autopas::AutoPas<Molecule> &autoPas) {
  auto containerOption = std::get<0>(std::get<0>(options));
  auto traversalOption = std::get<1>(std::get<0>(options));
  auto loadEstimatorOption = std::get<2>(std::get<0>(options));
  auto dataLayoutOption = std::get<1>(options);
  auto newton3Option = std::get<2>(options);
  auto cellSizeOption = std::get<3>(options);

  autoPas.setAllowedContainers({containerOption});
  autoPas.setAllowedTraversals({traversalOption});
  autoPas.setAllowedLoadEstimators({loadEstimatorOption});
  autoPas.setAllowedDataLayouts({dataLayoutOption});
  autoPas.setAllowedNewton3Options({newton3Option});
  autoPas.setAllowedCellSizeFactors(autopas::NumberSetFinite<double>(std::set<double>({cellSizeOption})));
}

void testSimulationLoop(testingTuple options) {
  // create AutoPas object
  autopas::AutoPas<Molecule> autoPas;

  setFromOptions(options, autoPas);

  defaultInit(autoPas);

  int numParticles = 0;
  int maxID = 0;

  auto addParticlePair = [&autoPas, &numParticles, &maxID](std::array<double, 3> pos1) {
    using namespace autopas::utils::ArrayMath::literals;
    // create two particles with distance .5
    double distance = .5;
    std::array<double, 3> distVec{0., distance, 0.};
    std::array<double, 3> pos2 = pos1 + distVec;

    Molecule particle1(pos1, {0., 0., 0.}, maxID++, 0);
    Molecule particle2(pos2, {0., 0., 0.}, maxID++, 0);
    numParticles += 2;

    // add the two particles!
    autoPas.addParticle(particle1);
    autoPas.addParticle(particle2);
  };

  auto moveParticlesAndResetF = [&autoPas](std::array<double, 3> moveVec) {
    using namespace autopas::utils::ArrayMath::literals;
    for (auto iter = autoPas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
      iter->setR(iter->getR() + moveVec);
      iter->setF(zeroArr);
    }
  };

  auto deleteIDs = [&autoPas, &numParticles](std::set<unsigned long> ids) {
    for (auto iter = autoPas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
      if (ids.find(iter->getID()) != ids.end()) {
        autoPas.deleteParticle(iter);
        --numParticles;
      }
    }
  };

  addParticlePair({9.99, 5., 5.});

  LJFunctorGlobals functor(cutoff);
  functor.setParticleProperties(24.0, 1);
  // do first simulation loop
  doSimulationLoop(autoPas, &functor);

  doAssertions(autoPas, &functor, numParticles, __LINE__);

  moveParticlesAndResetF({autoPas.getVerletSkin() / 6, 0., 0.});
  addParticlePair({9.99, 1., 5.});

  // do second simulation loop
  doSimulationLoop(autoPas, &functor);

  doAssertions(autoPas, &functor, numParticles, __LINE__);

  moveParticlesAndResetF({-autoPas.getVerletSkin() / 6, 0., 0.});
  addParticlePair({9.99, 7., 5.});
  deleteIDs({2, 3});

  // do third simulation loop.
  doSimulationLoop(autoPas, &functor);

  doAssertions(autoPas, &functor, numParticles, __LINE__);

  // update positions a bit (outside of domain!) + reset F
  moveParticlesAndResetF({autoPas.getVerletSkin() / 6, 0., 0.});

  // do fourth simulation loop, tests rebuilding of container.
  doSimulationLoop(autoPas, &functor);

  doAssertions(autoPas, &functor, numParticles, __LINE__);
}

/**
 * This test checks the correct calculation of an AutoPas container, where:
 * The container owns 26 particles, that lie with a distance of 0.25 away from the centers of the faces and edges, as
 * well as from the corners of the box. They always lie on a straight line `a` from these points to the center of the
 * box. Additionally, the container owns 26 halo particles that lie outside of the container on the same line `a` with a
 * distance of 0.5 from the owned particle.
 * No additional halo exchange is simulated in this test.
 * @param options
 */
void testHaloCalculation(testingTuple options) {
  using namespace autopas::utils::ArrayMath::literals;

  // create AutoPas object
  autopas::AutoPas<Molecule> autoPas;

  setFromOptions(options, autoPas);

  defaultInit(autoPas);

  // create particle pairs with distance .5
  double distance = .5;
  unsigned long id = 0;
  for (int x_diff : {-1, 0, 1}) {
    for (int y_diff : {-1, 0, 1}) {
      //    for (int y_diff : {0}) {
      for (int z_diff : {-1, 0, 1}) {
        //      for (int z_diff : {0}) {
        if (x_diff == 0 and y_diff == 0 and z_diff == 0) {
          continue;
        }
        double mul = 5.;
        double mid = 5.;
        std::array<double, 3> edge{x_diff * mul + mid, y_diff * mul + mid, z_diff * mul + mid};

        std::array<double, 3> diff = {x_diff * 1., y_diff * 1., z_diff * 1.};
        diff = autopas::utils::ArrayMath::normalize(diff) * (distance / 2.);

        auto pos1 = edge - diff;
        auto pos2 = edge + diff;

        Molecule particle1(pos1, {0., 0., 0.}, id++);
        autoPas.addParticle(particle1);
        Molecule particle2(pos2, {0., 0., 0.}, id++);
        autoPas.addHaloParticle(particle2);
      }
    }
  }

  LJFunctorGlobals functor(cutoff);
  functor.setParticleProperties(24, 1);

  autoPas.iteratePairwise(&functor);

  doAssertions(autoPas, &functor, 26, __LINE__);
}

TEST_P(AutoPasInterfaceTest, SimulationLoopTest) {
  // this test checks the correct behavior of the autopas interface.
  auto options = GetParam();
  try {
    testSimulationLoop(options);
  } catch (autopas::utils::ExceptionHandler::AutoPasException &autoPasException) {
    std::string str = autoPasException.what();
    if (str.find("Rejected the only configuration in the search space!") != std::string::npos) {
      GTEST_SKIP() << "skipped with exception: " << autoPasException.what() << std::endl;
    } else {
      // rethrow
      throw;
    }
  }
}

/**
 * This test checks the correct behavior of the AutoPas interface with respect to halo calculations, see also the
 * comments of testHaloCalculation() for a more detailed description.
 */
TEST_P(AutoPasInterfaceTest, HaloCalculationTest) {
  auto options = GetParam();
  try {
    testHaloCalculation(options);
  } catch (autopas::utils::ExceptionHandler::AutoPasException &autoPasException) {
    std::string str = autoPasException.what();
    if (str.find("Rejected the only configuration in the search space!") != std::string::npos) {
      GTEST_SKIP() << "skipped with exception: " << autoPasException.what() << std::endl;
    } else {
      // rethrow
      throw;
    }
  }
}

using ::testing::Combine;
using ::testing::UnorderedElementsAreArray;
using ::testing::Values;
using ::testing::ValuesIn;

INSTANTIATE_TEST_SUITE_P(
    Generated, AutoPasInterfaceTest,
    // proper indent
    Combine(ValuesIn([]() -> std::vector<std::tuple<autopas::ContainerOption, autopas::TraversalOption,
                                                    autopas::LoadEstimatorOption>> {
              std::vector<std::tuple<autopas::ContainerOption, autopas::TraversalOption, autopas::LoadEstimatorOption>>
                  tupleVector;
              for (const auto &containerOption : autopas::ContainerOption::getAllOptions()) {
                auto compatibleTraversals = autopas::compatibleTraversals::allCompatibleTraversals(containerOption);
                for (const auto &traversalOption : compatibleTraversals) {
                  auto compatibleLoadEstimators = autopas::loadEstimators::getApplicableLoadEstimators(
                      containerOption, traversalOption, autopas::LoadEstimatorOption::getAllOptions());
                  for (const auto &loadEstimatorOption : compatibleLoadEstimators) {
                    tupleVector.emplace_back(containerOption, traversalOption, loadEstimatorOption);
                  }
                }
              }
              return tupleVector;
            }()),
            ValuesIn(autopas::DataLayoutOption::getAllOptions()), ValuesIn(autopas::Newton3Option::getAllOptions()),
            Values(0.5, 1., 1.5)),
    AutoPasInterfaceTest::PrintToStringParamName());

////////////////////////////////////////////// FOR EVERY SINGLE CONTAINER //////////////////////////////////////////////

TEST_P(AutoPasInterface1ContainersTest, testResize) {
  // init an autopas object
  autopas::AutoPas<Molecule> autoPas;
  autoPas.setBoxMin({0, 0, 0});
  autoPas.setBoxMax({10, 10, 10});
  autoPas.setCutoff(1);
  autoPas.setVerletSkinPerTimestep(0.1);

  const auto &containerOp = GetParam();
  autoPas.setAllowedContainers({containerOp});
  autoPas.setAllowedTraversals(autopas::compatibleTraversals::allCompatibleTraversals(containerOp));

  autoPas.init();

  ASSERT_EQ(autoPas.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 0)
      << "Container was not initialized empty!";

  // add three purposely placed particles
  auto expectedParticles = addParticlesMinMidMax(autoPas);
  ASSERT_EQ(autoPas.getNumberOfParticles(autopas::IteratorBehavior::owned), expectedParticles.size())
      << "Container did not receive all particles before resize()!";

  auto boxMinNew = autopas::utils::ArrayMath::add(autoPas.getBoxMin(), {.5, .5, .5});
  auto boxMaxNew = autopas::utils::ArrayMath::add(autoPas.getBoxMax(), {1, 1, 1});

  auto particlesOutside = autoPas.resizeBox(boxMinNew, boxMaxNew);

  // remove particles that are now outside from the expectation
  for (auto &p : particlesOutside) {
    auto pInExpected = std::find(expectedParticles.begin(), expectedParticles.end(), p);
    EXPECT_NE(pInExpected, expectedParticles.end())
        << "Particle that was returned as \"outside\" is not one of the initially expected particles!";
    if (pInExpected != expectedParticles.end()) {
      expectedParticles.erase(pInExpected);
    }
  }

  ASSERT_EQ(autoPas.getNumberOfParticles(), expectedParticles.size())
      << "Container does not contain all particles after resize!";

  std::vector<Molecule> particlesInsideAfterResize{};
  particlesInsideAfterResize.reserve(autoPas.getNumberOfParticles());
  for (auto &p : autoPas) {
    particlesInsideAfterResize.push_back(p);
  }

  EXPECT_THAT(particlesInsideAfterResize, ::testing::UnorderedElementsAreArray(expectedParticles));
}

INSTANTIATE_TEST_SUITE_P(Generated, AutoPasInterface1ContainersTest, ValuesIn(getTestableContainerOptions()),
                         AutoPasInterface1ContainersTest::PrintToStringParamName());

/////////////////////////////////////// FOR EVERY COMBINATION OF TWO CONTAINERS ////////////////////////////////////////

void testSimulationLoop(autopas::ContainerOption containerOption1, autopas::ContainerOption containerOption2,
                        size_t autoPasDirection) {
  using namespace autopas::utils::ArrayMath::literals;
  // create AutoPas object
  autopas::AutoPas<Molecule> autoPas1;
  autoPas1.setOutputSuffix("1_");
  autoPas1.setAllowedContainers(std::set<autopas::ContainerOption>{containerOption1});
  autoPas1.setAllowedTraversals(autopas::compatibleTraversals::allCompatibleTraversals(containerOption1));
  autopas::AutoPas<Molecule> autoPas2;
  autoPas2.setOutputSuffix("2_");
  autoPas2.setAllowedContainers(std::set<autopas::ContainerOption>{containerOption2});
  autoPas2.setAllowedTraversals(autopas::compatibleTraversals::allCompatibleTraversals(containerOption2));

  defaultInit(autoPas1, autoPas2, autoPasDirection);

  // create two particles with distance .5
  double distance = .5;
  std::array<double, 3> pos1{9.99, 5., 5.};
  std::array<double, 3> distVec{0., distance, 0.};
  std::array<double, 3> pos2 = pos1 + distVec;

  {
    Molecule particle1(pos1, {0., 0., 0.}, 0, 0);
    Molecule particle2(pos2, {0., 0., 0.}, 1, 0);

    // add the two particles!
    for (auto *p : {&particle1, &particle2}) {
      if (autopas::utils::inBox(p->getR(), autoPas1.getBoxMin(), autoPas1.getBoxMax())) {
        autoPas1.addParticle(*p);
      } else {
        autoPas2.addParticle(*p);
      }
    }
  }

  LJFunctorGlobals functor1(cutoff);
  functor1.setParticleProperties(24.0, 1);
  LJFunctorGlobals functor2(cutoff);
  functor2.setParticleProperties(24.0, 1);
  // do first simulation loop
  doSimulationLoop(autoPas1, autoPas2, &functor1, &functor2);

  doAssertions(autoPas1, autoPas2, &functor1, &functor2, __LINE__);

  // update positions a bit (outside of domain!) + reset F
  {
    std::array<double, 3> moveVec{skinPerTimestep * rebuildFrequency / 3., 0., 0.};
    for (auto *aP : {&autoPas1, &autoPas2}) {
      for (auto iter = aP->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
        iter->setR(iter->getR() + moveVec);
        iter->setF(zeroArr);
      }
    }
  }

  // do second simulation loop
  doSimulationLoop(autoPas1, autoPas2, &functor1, &functor2);

  doAssertions(autoPas1, autoPas2, &functor1, &functor2, __LINE__);

  // update positions a bit (outside of domain!) + reset F
  {
    std::array<double, 3> moveVec{-skinPerTimestep * rebuildFrequency / 3., 0., 0.};
    for (auto *aP : {&autoPas1, &autoPas2}) {
      for (auto iter = aP->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
        iter->setR(iter->getR() + moveVec);
        iter->setF(zeroArr);
      }
    }
  }

  // do third simulation loop
  doSimulationLoop(autoPas1, autoPas2, &functor1, &functor2);

  doAssertions(autoPas1, autoPas2, &functor1, &functor2, __LINE__);

  // update positions a bit (outside of domain!) + reset F
  {
    std::array<double, 3> moveVec{skinPerTimestep * rebuildFrequency / 3., 0., 0.};
    for (auto *aP : {&autoPas1, &autoPas2}) {
      for (auto iter = aP->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
        iter->setR(iter->getR() + moveVec);
        iter->setF(zeroArr);
      }
    }
  }

  // do fourth simulation loop
  doSimulationLoop(autoPas1, autoPas2, &functor1, &functor2);

  doAssertions(autoPas1, autoPas2, &functor1, &functor2, __LINE__);
}

TEST_P(AutoPasInterface2ContainersTest, SimulationLoopTest) {
  // this test checks the correct behavior of the autopas interface.
  auto containerOptionTuple = GetParam();
  testSimulationLoop(std::get<0>(containerOptionTuple), std::get<1>(containerOptionTuple), 0);
}

using ::testing::Combine;
using ::testing::UnorderedElementsAreArray;
using ::testing::ValuesIn;

INSTANTIATE_TEST_SUITE_P(Generated, AutoPasInterface2ContainersTest,
                         Combine(ValuesIn(getTestableContainerOptions()), ValuesIn(getTestableContainerOptions())),
                         AutoPasInterface2ContainersTest::PrintToStringParamName());
