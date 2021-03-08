/**
 * @file RegionParticleIteratorTest.cpp
 * @author seckler
 * @date 03.04.18
 */
#include "RegionParticleIteratorTest.h"

#include <autopas/AutoPas.h>

#include "IteratorTestHelper.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/WrapOpenMP.h"
#include "testingHelpers/EmptyFunctor.h"

using namespace autopas;

//////////////////////////////// NEW TESTS ///////////////////////////////////

template <typename AutoPasT>
auto RegionParticleIteratorTest::defaultInit(AutoPasT &autoPas, autopas::ContainerOption &containerOption,
                                             double cellSizeFactor) {
  autoPas.setBoxMin({0., 0., 0.});
  autoPas.setBoxMax({10., 10., 10.});
  autoPas.setCutoff(1);
  autoPas.setVerletSkin(0.2);
  autoPas.setVerletRebuildFrequency(2);
  autoPas.setNumSamples(2);
  autoPas.setAllowedContainers(std::set<autopas::ContainerOption>{containerOption});
  autoPas.setAllowedTraversals(autopas::compatibleTraversals::allCompatibleTraversals(containerOption));
  autoPas.setAllowedCellSizeFactors(autopas::NumberSetFinite<double>(std::set<double>({cellSizeFactor})));

#ifdef AUTOPAS_CUDA
  autoPas.setVerletClusterSize(32);
#endif

  autoPas.init();

  auto haloBoxMin =
      autopas::utils::ArrayMath::subScalar(autoPas.getBoxMin(), autoPas.getVerletSkin() + autoPas.getCutoff());
  auto haloBoxMax =
      autopas::utils::ArrayMath::addScalar(autoPas.getBoxMax(), autoPas.getVerletSkin() + autoPas.getCutoff());

  return std::make_tuple(haloBoxMin, haloBoxMax);
}

TEST_P(RegionParticleIteratorTest, testRegionAroundCorner) {
  auto [containerOption, cellSizeFactor, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);

  using ::autopas::utils::ArrayMath::add;
  using ::autopas::utils::ArrayMath::mulScalar;
  using ::autopas::utils::ArrayMath::sub;
  auto domainLength = sub(autoPas.getBoxMax(), autoPas.getBoxMin());
  // draw a box around the lower corner of the domain
  auto searchBoxLengthHalf = mulScalar(domainLength, 0.3);
  std::array<double, 3> searchBoxMin = sub(autoPas.getBoxMin(), searchBoxLengthHalf);
  std::array<double, 3> searchBoxMax = add(autoPas.getBoxMin(), searchBoxLengthHalf);

  auto [particleIDsOwned, particleIDsHalo, particleIDsInBoxOwned, particleIDsInBoxHalo] =
      IteratorTestHelper::fillContainerAroundBoundary(autoPas, searchBoxMin, searchBoxMax);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyFunctor<Molecule> eFunctor;
    autoPas.iteratePairwise(&eFunctor);
  }

  std::vector<size_t> expectedIDs;
  switch (behavior) {
    case autopas::IteratorBehavior::ownedOnly: {
      expectedIDs = particleIDsInBoxOwned;
      break;
    }
    case autopas::IteratorBehavior::haloOnly: {
      expectedIDs = particleIDsInBoxHalo;
      break;
    }
    case autopas::IteratorBehavior::haloAndOwned: {
      expectedIDs = particleIDsInBoxOwned;
      expectedIDs.insert(expectedIDs.end(), particleIDsInBoxHalo.begin(), particleIDsInBoxHalo.end());
      break;
    }
    case autopas::IteratorBehavior::haloOwnedAndDummy: {
      GTEST_FAIL() << "IteratorBehavior::haloOwnedAndDummy should not be tested through this test"
                      " as container behavior with dummy particles is not uniform.";
      break;
    }
  }

  // sanity check: there should be particles in the expected region
  ASSERT_THAT(expectedIDs, ::testing::Not(::testing::IsEmpty()));

  // actual test
  IteratorTestHelper::provideRegionIterator(
      useConstIterator, autoPas, behavior, searchBoxMin, searchBoxMax,
      [&](const auto &autopas, auto &iter) { IteratorTestHelper::findParticles(autoPas, iter, expectedIDs); });
}

using ::testing::Combine;
using ::testing::UnorderedElementsAreArray;
using ::testing::Values;
using ::testing::ValuesIn;

static inline auto getTestableContainerOptions() {
#ifdef AUTOPAS_CUDA
  return autopas::ContainerOption::getAllOptions();
#else
  auto containerOptions = autopas::ContainerOption::getAllOptions();
  containerOptions.erase(containerOptions.find(autopas::ContainerOption::verletClusterCells));
  return containerOptions;
#endif
}

static inline auto getIteratorBehaviorOptions() {
  auto allOptions = autopas::IteratorBehavior::getAllOptions();
  std::set<autopas::IteratorBehavior> retSet;
  // we ignore dummy particles in the general tests because they can behave differently depending on the container
  std::set<autopas::IteratorBehavior> ignoredOptions = {autopas::IteratorBehavior::haloOwnedAndDummy};
  std::set_difference(allOptions.begin(), allOptions.end(), ignoredOptions.begin(), ignoredOptions.end(),
                      std::inserter(retSet, retSet.begin()));
  return retSet;
}

INSTANTIATE_TEST_SUITE_P(Generated, RegionParticleIteratorTest,
                         Combine(ValuesIn(getTestableContainerOptions()), /*cell size factor*/ Values(0.5, 1., 1.5),
                                 /*use const*/ Values(true, false), /*prior force calc*/ Values(true, false),
                                 ValuesIn(getIteratorBehaviorOptions())),
                         RegionParticleIteratorTest::PrintToStringParamName());

//////////////////////////////// OLD TESTS ///////////////////////////////////
/********************************** Linked Cells Tests **********************************/
//
// TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIterator) {
//  LinkedCells<TouchableParticle> lcContainer(_boxMin, _boxMax, _cutoff, 0., 1.);
//
//  // add a number of particles
//  autopasTools::generators::RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0),
//                                                               lcContainer.getBoxMin(), lcContainer.getBoxMax(), 100);
//#ifdef AUTOPAS_OPENMP
//#pragma omp parallel default(none) shared(lcContainer)
//#endif
//  // touch them using the regionIterator
//  for (auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax); iterator.isValid(); ++iterator) {
//    iterator->touch();
//  }
//
//  checkTouches(lcContainer, _regionMin, _regionMax);
//}
//
// void RegionParticleIteratorTest::testLinkedCellsRegionParticleIteratorBehaviorOwned() {
//  LinkedCells<TouchableParticle> lcContainer(_boxMin, _boxMax, _cutoff, 0., 1.);
//
//  // add a number of particles
//  autopasTools::generators::RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0),
//                                                               lcContainer.getBoxMin(), lcContainer.getBoxMax(), 100);
//
//  TouchableParticle part(utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), 100);
//  lcContainer.addHaloParticle(part);
//
//  // touch them using the regionIterator
//  auto testRegionMin = utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5);
//#ifdef AUTOPAS_OPENMP
//#pragma omp parallel default(none) shared(lcContainer, testRegionMin)
//#endif
//  for (auto iterator = lcContainer.getRegionIterator(testRegionMin, _regionMax, autopas::IteratorBehavior::ownedOnly);
//       iterator.isValid(); ++iterator) {
//    iterator->touch();
//  }
//  // owned cells only start at [0, 0, 0]!
//  std::array<double, 3> realMin = {0, 0, 0};
//  checkTouches(lcContainer, realMin, _regionMax);
//}
//
// TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorOwned1Thread) {
//#ifdef AUTOPAS_OPENMP
//  int before = omp_get_max_threads();
//  omp_set_num_threads(1);
//#endif
//  testLinkedCellsRegionParticleIteratorBehaviorOwned();
//#ifdef AUTOPAS_OPENMP
//  omp_set_num_threads(before);
//#endif
//}
//#ifdef AUTOPAS_OPENMP
// TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorOwned16Thread) {
//  int before = omp_get_max_threads();
//  omp_set_num_threads(16);
//  testLinkedCellsRegionParticleIteratorBehaviorOwned();
//  omp_set_num_threads(before);
//}
// TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorOwned32Thread) {
//  int before = omp_get_max_threads();
//  omp_set_num_threads(32);
//  testLinkedCellsRegionParticleIteratorBehaviorOwned();
//  omp_set_num_threads(before);
//}
// TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorOwned50Thread) {
//  int before = omp_get_max_threads();
//  omp_set_num_threads(50);
//  testLinkedCellsRegionParticleIteratorBehaviorOwned();
//  omp_set_num_threads(before);
//}
// TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorOwned240Thread) {
//  int before = omp_get_max_threads();
//  omp_set_num_threads(240);
//  testLinkedCellsRegionParticleIteratorBehaviorOwned();
//  omp_set_num_threads(before);
//}
//#endif
//
// void RegionParticleIteratorTest::testLinkedCellsRegionParticleIteratorBehaviorHalo() {
//  LinkedCells<TouchableParticle> lcContainer(_boxMin, _boxMax, _cutoff, 0., 1.);
//
//  // add a number of particles
//  autopasTools::generators::RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0),
//                                                               lcContainer.getBoxMin(), lcContainer.getBoxMax(), 100);
//
//  TouchableParticle part(utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), 100);
//  lcContainer.addHaloParticle(part);
//
//  auto testRegionMin = utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5);
//  // touch them using the regionIterator
//#ifdef AUTOPAS_OPENMP
//#pragma omp parallel default(none) shared(lcContainer, testRegionMin)
//#endif
//  for (auto iterator = lcContainer.getRegionIterator(testRegionMin, _regionMax, autopas::IteratorBehavior::haloOnly);
//       iterator.isValid(); ++iterator) {
//    iterator->touch();
//    bool isInRegionOfInterest = utils::inBox(iterator->getR(), testRegionMin, _regionMax);
//    bool isInHalo = not utils::inBox(iterator->getR(), _boxMin, _boxMax);
//    EXPECT_TRUE(isInRegionOfInterest and isInHalo)
//        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]"
//        << " in thread: " << autopas_get_thread_num() << std::endl;
//  }
//
//  // check the touch using the normal iterator (needed to check whether no particle was forgotten)
//  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
//    // this is a test for halo only! so we first check whether it's within our region of interest and then whether
//    it's
//    // not in the halo
//    bool isInRegionOfInterest = utils::inBox(iterator->getR(), testRegionMin, _regionMax);
//    bool isInHalo = not utils::inBox(iterator->getR(), _boxMin, _boxMax);
//    EXPECT_EQ(isInRegionOfInterest and isInHalo ? 1 : 0, iterator->getNumTouched())
//        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
//  }
//}
//
// TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorHalo1Thread) {
//#ifdef AUTOPAS_OPENMP
//  int before = omp_get_max_threads();
//  omp_set_num_threads(1);
//#endif
//  testLinkedCellsRegionParticleIteratorBehaviorHalo();
//#ifdef AUTOPAS_OPENMP
//  omp_set_num_threads(before);
//#endif
//}
//
//#ifdef AUTOPAS_OPENMP
// TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorHalo16Threads) {
//  int before = omp_get_max_threads();
//  omp_set_num_threads(16);
//  testLinkedCellsRegionParticleIteratorBehaviorHalo();
//  omp_set_num_threads(before);
//}
//
// TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorHalo32Threads) {
//  int before = omp_get_max_threads();
//  omp_set_num_threads(32);
//  testLinkedCellsRegionParticleIteratorBehaviorHalo();
//  omp_set_num_threads(before);
//}
//
// TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorHalo50Threads) {
//  int before = omp_get_max_threads();
//  omp_set_num_threads(50);
//  testLinkedCellsRegionParticleIteratorBehaviorHalo();
//  omp_set_num_threads(before);
//}
// TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorHalo240Threads) {
//  int before = omp_get_max_threads();
//  omp_set_num_threads(240);
//  testLinkedCellsRegionParticleIteratorBehaviorHalo();
//  omp_set_num_threads(before);
//}
//#endif
//
// TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorCopyConstructor) {
//  LinkedCells<TouchableParticle> lcContainer(_boxMin, _boxMax, _cutoff, 0., 1.);
//
//  // add a number of particles
//  autopasTools::generators::RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0),
//                                                               lcContainer.getBoxMin(), lcContainer.getBoxMax(), 100);
//
//  {
//    auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax);
//    auto iterator2 = iterator;
//
//    // touch them using the regionIterator
//    for (; iterator2.isValid(); ++iterator2) {
//      iterator2->touch();
//    }
//  }
//
//  // check the touch using the normal iterator
//  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
//    //  std::cout << "id: " << iterator->getID() << " at [" <<
//    //  iterator->getR()[0]
//    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
//    //              << "] touched:" << iterator->getNumTouched() << std::endl;
//
//    EXPECT_EQ(utils::inBox(iterator->getR(), _regionMin, _regionMax) ? 1 : 0, iterator->getNumTouched())
//        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
//  }
//}
//
// TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorCopyAssignment) {
//  LinkedCells<TouchableParticle> lcContainer(_boxMin, _boxMax, _cutoff, 0., 1.);
//
//  // add a number of particles
//  autopasTools::generators::RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0),
//                                                               lcContainer.getBoxMin(), lcContainer.getBoxMax(), 100);
//
//  {
//    auto iterator2 = lcContainer.getRegionIterator(_regionMin, _regionMax);
//    // touch them using the regionIterator
//    for (; iterator2.isValid(); ++iterator2) {
//      iterator2->touch();
//    }
//
//    auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax);
//    // touch them using the regionIterator
//    for (; iterator.isValid(); ++iterator) {
//      iterator->touch();
//    }
//
//    iterator2 = lcContainer.getRegionIterator(_regionMin, _regionMax);
//    // touch them using the regionIterator
//    for (; iterator2.isValid(); ++iterator2) {
//      iterator2->touch();
//    }
//  }
//
//  // check the touch using the normal iterator
//  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
//    //  std::cout << "id: " << iterator->getID() << " at [" <<
//    //  iterator->getR()[0]
//    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
//    //              << "] touched:" << iterator->getNumTouched() << std::endl;
//
//    EXPECT_EQ(utils::inBox(iterator->getR(), _regionMin, _regionMax) ? 3 : 0, iterator->getNumTouched())
//        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
//  }
//}
//
//
///*********************************** Direct Sum Tests ***********************************/
//
//
// TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIterator) {
//  DirectSum<TouchableParticle> container(_boxMin, _boxMax, _cutoff, 0.);
//
//  // add a number of particles
//  autopasTools::generators::RandomGenerator::fillWithParticles(container, TouchableParticle({0., 0., 0.}, 0),
//                                                               container.getBoxMin(), container.getBoxMax(), 100);
//
//  // touch them using the regionIterator
//#ifdef AUTOPAS_OPENMP
//#pragma omp parallel default(none) shared(container)
//#endif
//  for (auto iterator = container.getRegionIterator(_regionMin, _regionMax); iterator.isValid(); ++iterator) {
//    iterator->touch();
//  }
//
//  // check the touch. Iterating over cells provides more debug info than normal iterator.
//  // check the touch using the normal iterator
//  for (auto iterator = container.begin(); iterator.isValid(); ++iterator) {
//    //  std::cout << "id: " << iterator->getID() << " at [" <<
//    //  iterator->getR()[0]
//    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
//    //              << "] touched:" << iterator->getNumTouched() << std::endl;
//
//    EXPECT_EQ(utils::inBox(iterator->getR(), _regionMin, _regionMax) ? 1 : 0, iterator->getNumTouched())
//        << "at: [" << utils::ArrayUtils::to_string(iterator->getR()) << "]";
//  }
//}
//
// TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorBehaviorOwned) {
//  DirectSum<TouchableParticle> container(_boxMin, _boxMax, _cutoff, 0.);
//
//  // add a number of particles
//  autopasTools::generators::RandomGenerator::fillWithParticles(container, TouchableParticle({0., 0., 0.}, 0),
//                                                               container.getBoxMin(), container.getBoxMax(), 100);
//
//  TouchableParticle part(utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), 100);
//  container.addHaloParticle(part);
//
//  // touch them using the regionIterator
//#ifdef AUTOPAS_OPENMP
//#pragma omp parallel default(none) shared(container)
//#endif
//  for (auto iterator = container.getRegionIterator(utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), _regionMax,
//                                                   autopas::IteratorBehavior::ownedOnly);
//       iterator.isValid(); ++iterator) {
//    iterator->touch();
//  }
//
//  // check the touch using the normal iterator
//  for (auto iterator = container.begin(); iterator.isValid(); ++iterator) {
//    //  std::cout << "id: " << iterator->getID() << " at [" <<
//    //  iterator->getR()[0]
//    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
//    //              << "] touched:" << iterator->getNumTouched() << std::endl;
//
//    EXPECT_EQ(utils::inBox(iterator->getR(), _boxMin, _regionMax) ? 1 : 0, iterator->getNumTouched())
//        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
//  }
//}
//
// TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorBehaviorHalo) {
//  DirectSum<TouchableParticle> container(_boxMin, _boxMax, _cutoff, 0.);
//
//  // add a number of particles
//  autopasTools::generators::RandomGenerator::fillWithParticles(container, TouchableParticle({0., 0., 0.}, 0),
//                                                               container.getBoxMin(), container.getBoxMax(), 100);
//
//  TouchableParticle part(utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), 100);
//  container.addHaloParticle(part);
//
//  // touch them using the regionIterator
//#ifdef AUTOPAS_OPENMP
//#pragma omp parallel default(none) shared(container)
//#endif
//  for (auto iterator = container.getRegionIterator(utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), _regionMax,
//                                                   autopas::IteratorBehavior::haloOnly);
//       iterator.isValid(); ++iterator) {
//    iterator->touch();
//  }
//
//  // check the touch using the normal iterator
//  for (auto iterator = container.begin(); iterator.isValid(); ++iterator) {
//    //  std::cout << "id: " << iterator->getID() << " at [" <<
//    //  iterator->getR()[0]
//    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
//    //              << "] touched:" << iterator->getNumTouched() << std::endl;
//
//    EXPECT_EQ(utils::inBox(iterator->getR(), utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), _regionMax)
//                  ? (utils::inBox(iterator->getR(), _boxMin, _regionMax) ? 0 : 1)
//                  : 0,
//              iterator->getNumTouched())
//        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
//  }
//}
//
// TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorCopyConstructor) {
//  DirectSum<TouchableParticle> container(_boxMin, _boxMax, _cutoff, 0.);
//
//  // add a number of particles
//  autopasTools::generators::RandomGenerator::fillWithParticles(container, TouchableParticle({0., 0., 0.}, 0),
//                                                               container.getBoxMin(), container.getBoxMax(), 100);
//
//  {
//    auto iterator = container.getRegionIterator(_regionMin, _regionMax);
//    auto iterator2 = iterator;
//
//    // touch them using the regionIterator
//    for (; iterator2.isValid(); ++iterator2) {
//      iterator2->touch();
//    }
//  }
//
//  // check the touch using the normal iterator
//  for (auto iterator = container.begin(); iterator.isValid(); ++iterator) {
//    //  std::cout << "id: " << iterator->getID() << " at [" <<
//    //  iterator->getR()[0]
//    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
//    //              << "] touched:" << iterator->getNumTouched() << std::endl;
//
//    EXPECT_EQ(utils::inBox(iterator->getR(), _regionMin, _regionMax) ? 1 : 0, iterator->getNumTouched())
//        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
//  }
//}
//
// TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorCopyAssignment) {
//  DirectSum<TouchableParticle> container(_boxMin, _boxMax, _cutoff, 0.);
//
//  // add a number of particles
//  autopasTools::generators::RandomGenerator::fillWithParticles(container, TouchableParticle({0., 0., 0.}, 0),
//                                                               container.getBoxMin(), container.getBoxMax(), 100);
//
//  {
//    auto iterator2 = container.getRegionIterator(_regionMin, _regionMax);
//    // touch them using the regionIterator
//    for (; iterator2.isValid(); ++iterator2) {
//      iterator2->touch();
//    }
//
//    auto iterator = container.getRegionIterator(_regionMin, _regionMax);
//    // touch them using the regionIterator
//    for (; iterator.isValid(); ++iterator) {
//      iterator->touch();
//    }
//
//    iterator2 = container.getRegionIterator(_regionMin, _regionMax);
//    // touch them using the regionIterator
//    for (; iterator2.isValid(); ++iterator2) {
//      iterator2->touch();
//    }
//  }
//
//  // check the touch using the normal iterator
//  for (auto iterator = container.begin(); iterator.isValid(); ++iterator) {
//    //  std::cout << "id: " << iterator->getID() << " at [" <<
//    //  iterator->getR()[0]
//    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
//    //              << "] touched:" << iterator->getNumTouched() << std::endl;
//
//    EXPECT_EQ(utils::inBox(iterator->getR(), _regionMin, _regionMax) ? 3 : 0, iterator->getNumTouched())
//        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
//  }
//}
//
///*********************************** VerletList Tests ***********************************/
//
// void RegionParticleIteratorTest::checkTouches(LCTouch &lcContainer, std::array<double, 3> &regionMin,
//                                              std::array<double, 3> &regionMax) {
//  int numTouches = 0;
//  // check the touch. Iterating over cells provides more debug info than normal iterator.
//  for (size_t cellId = 0; cellId < lcContainer.getCells().size(); ++cellId) {
//    const auto cellId3D = utils::ThreeDimensionalMapping::oneToThreeD(cellId, {7, 7, 7});
//    for (auto pIter = lcContainer.getCells()[cellId].begin(); pIter.isValid(); ++pIter) {
//      ++numTouches;
//      EXPECT_EQ(utils::inBox(pIter->getR(), regionMin, regionMax) ? 1 : 0, pIter->getNumTouched())
//          << "at: [" << utils::ArrayUtils::to_string(pIter->getR()) << "]" << std::endl
//          << "in cell: " << cellId << " [" << utils::ArrayUtils::to_string(cellId3D) << "]";
//    }
//  }
//  EXPECT_GE(numTouches, 0) << "No Particles were checked!";
//}
