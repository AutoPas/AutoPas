/**
 * @file RegionParticleIteratorTest.cpp
 * @author seckler
 * @date 03.04.18
 */
#include "RegionParticleIteratorTest.h"

#include "autopas/containers/directSum/DirectSum.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/WrapOpenMP.h"

using namespace autopas;

/********************************** Linked Cells Tests **********************************/

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIterator) {
  LinkedCells<TouchableParticle> lcContainer(_boxMin, _boxMax, _cutoff, 0., 1.);

  // add a number of particles
  autopasTools::generators::RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0),
                                                               lcContainer.getBoxMin(), lcContainer.getBoxMax(), 100);
#ifdef AUTOPAS_OPENMP
#pragma omp parallel default(none) shared(lcContainer)
#endif
  // touch them using the regionIterator
  for (auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax); iterator.isValid(); ++iterator) {
    iterator->touch();
  }

  checkTouches(lcContainer, _regionMin, _regionMax);
}

void RegionParticleIteratorTest::testLinkedCellsRegionParticleIteratorBehaviorOwned() {
  LinkedCells<TouchableParticle> lcContainer(_boxMin, _boxMax, _cutoff, 0., 1.);

  // add a number of particles
  autopasTools::generators::RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0),
                                                               lcContainer.getBoxMin(), lcContainer.getBoxMax(), 100);

  TouchableParticle part(utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), 100);
  lcContainer.addHaloParticle(part);

  // touch them using the regionIterator
  auto testRegionMin = utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5);
#ifdef AUTOPAS_OPENMP
#pragma omp parallel default(none) shared(lcContainer, testRegionMin)
#endif
  for (auto iterator = lcContainer.getRegionIterator(testRegionMin, _regionMax, autopas::IteratorBehavior::ownedOnly);
       iterator.isValid(); ++iterator) {
    iterator->touch();
  }
  // owned cells only start at [0, 0, 0]!
  std::array<double, 3> realMin = {0, 0, 0};
  checkTouches(lcContainer, realMin, _regionMax);
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorOwned1Thread) {
#ifdef AUTOPAS_OPENMP
  int before = omp_get_max_threads();
  omp_set_num_threads(1);
#endif
  testLinkedCellsRegionParticleIteratorBehaviorOwned();
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(before);
#endif
}
#ifdef AUTOPAS_OPENMP
TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorOwned16Thread) {
  int before = omp_get_max_threads();
  omp_set_num_threads(16);
  testLinkedCellsRegionParticleIteratorBehaviorOwned();
  omp_set_num_threads(before);
}
TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorOwned32Thread) {
  int before = omp_get_max_threads();
  omp_set_num_threads(32);
  testLinkedCellsRegionParticleIteratorBehaviorOwned();
  omp_set_num_threads(before);
}
TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorOwned50Thread) {
  int before = omp_get_max_threads();
  omp_set_num_threads(50);
  testLinkedCellsRegionParticleIteratorBehaviorOwned();
  omp_set_num_threads(before);
}
TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorOwned240Thread) {
  int before = omp_get_max_threads();
  omp_set_num_threads(240);
  testLinkedCellsRegionParticleIteratorBehaviorOwned();
  omp_set_num_threads(before);
}
#endif

void RegionParticleIteratorTest::testLinkedCellsRegionParticleIteratorBehaviorHalo() {
  LinkedCells<TouchableParticle> lcContainer(_boxMin, _boxMax, _cutoff, 0., 1.);

  // add a number of particles
  autopasTools::generators::RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0),
                                                               lcContainer.getBoxMin(), lcContainer.getBoxMax(), 100);

  TouchableParticle part(utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), 100);
  lcContainer.addHaloParticle(part);

  auto testRegionMin = utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5);
  // touch them using the regionIterator
#ifdef AUTOPAS_OPENMP
#pragma omp parallel default(none) shared(lcContainer, testRegionMin)
#endif
  for (auto iterator = lcContainer.getRegionIterator(testRegionMin, _regionMax, autopas::IteratorBehavior::haloOnly);
       iterator.isValid(); ++iterator) {
    iterator->touch();
    bool isInRegionOfInterest = utils::inBox(iterator->getR(), testRegionMin, _regionMax);
    bool isInHalo = not utils::inBox(iterator->getR(), _boxMin, _boxMax);
    EXPECT_TRUE(isInRegionOfInterest and isInHalo)
        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]"
        << " in thread: " << autopas_get_thread_num() << std::endl;
  }

  // check the touch using the normal iterator (needed to check whether no particle was forgotten)
  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    // this is a test for halo only! so we first check whether it's within our region of interest and then whether it's
    // not in the halo
    bool isInRegionOfInterest = utils::inBox(iterator->getR(), testRegionMin, _regionMax);
    bool isInHalo = not utils::inBox(iterator->getR(), _boxMin, _boxMax);
    EXPECT_EQ(isInRegionOfInterest and isInHalo ? 1 : 0, iterator->getNumTouched())
        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
  }
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorHalo1Thread) {
#ifdef AUTOPAS_OPENMP
  int before = omp_get_max_threads();
  omp_set_num_threads(1);
#endif
  testLinkedCellsRegionParticleIteratorBehaviorHalo();
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(before);
#endif
}

#ifdef AUTOPAS_OPENMP
TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorHalo16Threads) {
  int before = omp_get_max_threads();
  omp_set_num_threads(16);
  testLinkedCellsRegionParticleIteratorBehaviorHalo();
  omp_set_num_threads(before);
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorHalo32Threads) {
  int before = omp_get_max_threads();
  omp_set_num_threads(32);
  testLinkedCellsRegionParticleIteratorBehaviorHalo();
  omp_set_num_threads(before);
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorHalo50Threads) {
  int before = omp_get_max_threads();
  omp_set_num_threads(50);
  testLinkedCellsRegionParticleIteratorBehaviorHalo();
  omp_set_num_threads(before);
}
TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorHalo240Threads) {
  int before = omp_get_max_threads();
  omp_set_num_threads(240);
  testLinkedCellsRegionParticleIteratorBehaviorHalo();
  omp_set_num_threads(before);
}
#endif

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorEmpty) {
  LinkedCells<TouchableParticle> lcContainer(_boxMin, _boxMax, _cutoff, 0., 1.);

  // add no particles

  // touch them using the regionIterator
#ifdef AUTOPAS_OPENMP
#pragma omp parallel default(none) shared(lcContainer)
#endif
  for (auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax); iterator.isValid(); ++iterator) {
    iterator->touch();
  }

  // no particles, hence the iterator should be invalid from the start
  EXPECT_FALSE(lcContainer.begin().isValid());

  // valid area is empty since nothing should have been touched
  checkTouches(lcContainer, _regionMin, _regionMin);
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorCopyConstructor) {
  LinkedCells<TouchableParticle> lcContainer(_boxMin, _boxMax, _cutoff, 0., 1.);

  // add a number of particles
  autopasTools::generators::RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0),
                                                               lcContainer.getBoxMin(), lcContainer.getBoxMax(), 100);

  {
    auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax);
    auto iterator2 = iterator;

    // touch them using the regionIterator
    for (; iterator2.isValid(); ++iterator2) {
      iterator2->touch();
    }
  }

  // check the touch using the normal iterator
  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    EXPECT_EQ(utils::inBox(iterator->getR(), _regionMin, _regionMax) ? 1 : 0, iterator->getNumTouched())
        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
  }
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorCopyAssignment) {
  LinkedCells<TouchableParticle> lcContainer(_boxMin, _boxMax, _cutoff, 0., 1.);

  // add a number of particles
  autopasTools::generators::RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0),
                                                               lcContainer.getBoxMin(), lcContainer.getBoxMax(), 100);

  {
    auto iterator2 = lcContainer.getRegionIterator(_regionMin, _regionMax);
    // touch them using the regionIterator
    for (; iterator2.isValid(); ++iterator2) {
      iterator2->touch();
    }

    auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax);
    // touch them using the regionIterator
    for (; iterator.isValid(); ++iterator) {
      iterator->touch();
    }

    iterator2 = lcContainer.getRegionIterator(_regionMin, _regionMax);
    // touch them using the regionIterator
    for (; iterator2.isValid(); ++iterator2) {
      iterator2->touch();
    }
  }

  // check the touch using the normal iterator
  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    EXPECT_EQ(utils::inBox(iterator->getR(), _regionMin, _regionMax) ? 3 : 0, iterator->getNumTouched())
        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
  }
}

/**
 * Tests for correct iterator behavior when some of the cells in and outside of the region of interest are empty
 */
TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorSparseDomain) {
  // box goes from {0,0,0} to {5,5,5} + one halo layer
  LinkedCells<TouchableParticle> lcContainer(_boxMin, _boxMax, _cutoff, 0., 1.);

  size_t idShouldTouch = 0;
  TouchableParticle p({0, 0, 0}, idShouldTouch);
  p.setR({2, 2, 2});
  p.setID(idShouldTouch++);
  lcContainer.addParticle(p);
  p.setR({2, 3, 2});
  p.setID(idShouldTouch++);
  lcContainer.addParticle(p);
  p.setR({2, 3, 3});
  p.setID(idShouldTouch++);
  lcContainer.addParticle(p);

  size_t idOffset = 1000;
  size_t idShouldNotTouch = idOffset;
  p.setR({1, 1, 1});
  p.setID(idShouldNotTouch++);
  lcContainer.addParticle(p);
  p.setR({2, 4.5, 2});
  p.setID(idShouldNotTouch++);
  lcContainer.addParticle(p);
  p.setR({4, 4, 4});
  p.setID(idShouldNotTouch++);
  lcContainer.addParticle(p);

  std::array<double, 3> regionOfInterstMin = {2, 2, 2};
  std::array<double, 3> regionOfInterstMax = {3, 4, 4};

  int particlesTouched = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : particlesTouched) default(none) shared(lcContainer, regionOfInterstMin, regionOfInterstMax)
#endif
  for (auto iterator = lcContainer.getRegionIterator(regionOfInterstMin, regionOfInterstMax); iterator.isValid();
       ++iterator) {
    iterator->touch();
    ++particlesTouched;
  }
  EXPECT_EQ(particlesTouched, idShouldTouch);

  int particlesChecked = 0;
  // no openmp for check!
  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    EXPECT_EQ(utils::inBox(iterator->getR(), regionOfInterstMin, regionOfInterstMax) ? 1 : 0, iterator->getNumTouched())
        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
    ++particlesChecked;
  }
  EXPECT_EQ(particlesChecked, idShouldTouch + idShouldNotTouch - idOffset);
}

/*********************************** Direct Sum Tests ***********************************/

TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorSparseDomain) {
  // box goes from {0,0,0} to {5,5,5} + one halo layer
  DirectSum<TouchableParticle> dsContainer(_boxMin, _boxMax, _cutoff, 0.);

  size_t idShouldTouch = 0;
  TouchableParticle p({0, 0, 0}, idShouldTouch);
  p.setR({2, 2, 2});
  p.setID(idShouldTouch++);
  dsContainer.addParticle(p);
  p.setR({2, 3, 2});
  p.setID(idShouldTouch++);
  dsContainer.addParticle(p);
  p.setR({2, 3, 3});
  p.setID(idShouldTouch++);
  dsContainer.addParticle(p);

  size_t idOffset = 1000;
  size_t idShouldNotTouch = idOffset;
  p.setR({1, 1, 1});
  p.setID(idShouldNotTouch++);
  dsContainer.addParticle(p);
  p.setR({2, 4.5, 2});
  p.setID(idShouldNotTouch++);
  dsContainer.addParticle(p);
  p.setR({4, 4, 4});
  p.setID(idShouldNotTouch++);
  dsContainer.addParticle(p);

  std::array<double, 3> regionOfInterstMin = {2, 2, 2};
  std::array<double, 3> regionOfInterstMax = {3, 4, 4};

  int particlesTouched = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : particlesTouched) default(none) shared(dsContainer, regionOfInterstMin, regionOfInterstMax)
#endif
  for (auto iterator = dsContainer.getRegionIterator(regionOfInterstMin, regionOfInterstMax); iterator.isValid();
       ++iterator) {
    iterator->touch();
    ++particlesTouched;
  }
  EXPECT_EQ(particlesTouched, idShouldTouch);

  int particlesChecked = 0;
  for (auto iterator = dsContainer.begin(); iterator.isValid(); ++iterator) {
    EXPECT_EQ(utils::inBox(iterator->getR(), regionOfInterstMin, regionOfInterstMax) ? 1 : 0, iterator->getNumTouched())
        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
    ++particlesChecked;
  }
  EXPECT_EQ(particlesChecked, idShouldTouch + idShouldNotTouch - idOffset);
}

TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIterator) {
  DirectSum<TouchableParticle> container(_boxMin, _boxMax, _cutoff, 0.);

  // add a number of particles
  autopasTools::generators::RandomGenerator::fillWithParticles(container, TouchableParticle({0., 0., 0.}, 0),
                                                               container.getBoxMin(), container.getBoxMax(), 100);

  // touch them using the regionIterator
#ifdef AUTOPAS_OPENMP
#pragma omp parallel default(none) shared(container)
#endif
  for (auto iterator = container.getRegionIterator(_regionMin, _regionMax); iterator.isValid(); ++iterator) {
    iterator->touch();
  }

  // check the touch. Iterating over cells provides more debug info than normal iterator.
  // check the touch using the normal iterator
  for (auto iterator = container.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    EXPECT_EQ(utils::inBox(iterator->getR(), _regionMin, _regionMax) ? 1 : 0, iterator->getNumTouched())
        << "at: [" << utils::ArrayUtils::to_string(iterator->getR()) << "]";
  }
}

TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorBehaviorOwned) {
  DirectSum<TouchableParticle> container(_boxMin, _boxMax, _cutoff, 0.);

  // add a number of particles
  autopasTools::generators::RandomGenerator::fillWithParticles(container, TouchableParticle({0., 0., 0.}, 0),
                                                               container.getBoxMin(), container.getBoxMax(), 100);

  TouchableParticle part(utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), 100);
  container.addHaloParticle(part);

  // touch them using the regionIterator
#ifdef AUTOPAS_OPENMP
#pragma omp parallel default(none) shared(container)
#endif
  for (auto iterator = container.getRegionIterator(utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), _regionMax,
                                                   autopas::IteratorBehavior::ownedOnly);
       iterator.isValid(); ++iterator) {
    iterator->touch();
  }

  // check the touch using the normal iterator
  for (auto iterator = container.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    EXPECT_EQ(utils::inBox(iterator->getR(), _boxMin, _regionMax) ? 1 : 0, iterator->getNumTouched())
        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
  }
}

TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorBehaviorHalo) {
  DirectSum<TouchableParticle> container(_boxMin, _boxMax, _cutoff, 0.);

  // add a number of particles
  autopasTools::generators::RandomGenerator::fillWithParticles(container, TouchableParticle({0., 0., 0.}, 0),
                                                               container.getBoxMin(), container.getBoxMax(), 100);

  TouchableParticle part(utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), 100);
  container.addHaloParticle(part);

  // touch them using the regionIterator
#ifdef AUTOPAS_OPENMP
#pragma omp parallel default(none) shared(container)
#endif
  for (auto iterator = container.getRegionIterator(utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), _regionMax,
                                                   autopas::IteratorBehavior::haloOnly);
       iterator.isValid(); ++iterator) {
    iterator->touch();
  }

  // check the touch using the normal iterator
  for (auto iterator = container.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    EXPECT_EQ(utils::inBox(iterator->getR(), utils::ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), _regionMax)
                  ? (utils::inBox(iterator->getR(), _boxMin, _regionMax) ? 0 : 1)
                  : 0,
              iterator->getNumTouched())
        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
  }
}

TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorEmpty) {
  DirectSum<TouchableParticle> container(_boxMin, _boxMax, _cutoff, 0.);

  // add no particles

  int i = 0;
  // touch them using the regionIterator
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : i) default(none) shared(container)
#endif
  for (auto iterator = container.getRegionIterator(_regionMin, _regionMax); iterator.isValid(); ++iterator) {
    iterator->touch();
    i++;
  }

  // check the touch using the normal iterator
  for (auto iterator = container.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    EXPECT_EQ(utils::inBox(iterator->getR(), _regionMin, _regionMax) ? 1 : 0, iterator->getNumTouched())
        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
    i++;
  }
  EXPECT_EQ(i, 0);
}

TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorCopyConstructor) {
  DirectSum<TouchableParticle> container(_boxMin, _boxMax, _cutoff, 0.);

  // add a number of particles
  autopasTools::generators::RandomGenerator::fillWithParticles(container, TouchableParticle({0., 0., 0.}, 0),
                                                               container.getBoxMin(), container.getBoxMax(), 100);

  {
    auto iterator = container.getRegionIterator(_regionMin, _regionMax);
    auto iterator2 = iterator;

    // touch them using the regionIterator
    for (; iterator2.isValid(); ++iterator2) {
      iterator2->touch();
    }
  }

  // check the touch using the normal iterator
  for (auto iterator = container.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    EXPECT_EQ(utils::inBox(iterator->getR(), _regionMin, _regionMax) ? 1 : 0, iterator->getNumTouched())
        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
  }
}

TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorCopyAssignment) {
  DirectSum<TouchableParticle> container(_boxMin, _boxMax, _cutoff, 0.);

  // add a number of particles
  autopasTools::generators::RandomGenerator::fillWithParticles(container, TouchableParticle({0., 0., 0.}, 0),
                                                               container.getBoxMin(), container.getBoxMax(), 100);

  {
    auto iterator2 = container.getRegionIterator(_regionMin, _regionMax);
    // touch them using the regionIterator
    for (; iterator2.isValid(); ++iterator2) {
      iterator2->touch();
    }

    auto iterator = container.getRegionIterator(_regionMin, _regionMax);
    // touch them using the regionIterator
    for (; iterator.isValid(); ++iterator) {
      iterator->touch();
    }

    iterator2 = container.getRegionIterator(_regionMin, _regionMax);
    // touch them using the regionIterator
    for (; iterator2.isValid(); ++iterator2) {
      iterator2->touch();
    }
  }

  // check the touch using the normal iterator
  for (auto iterator = container.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    EXPECT_EQ(utils::inBox(iterator->getR(), _regionMin, _regionMax) ? 3 : 0, iterator->getNumTouched())
        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
  }
}

/*********************************** VerletList Tests ***********************************/

TEST_F(RegionParticleIteratorTest, testVerletRegionParticleIteratorSparseDomain) {
  // box goes from {0,0,0} to {5,5,5} + one halo layer
  VerletLists<TouchableParticle> vlContainer(_boxMin, _boxMax, _cutoff, _cutoff / 3);

  size_t idShouldTouch = 0;
  TouchableParticle p({0, 0, 0}, idShouldTouch);
  p.setR({2, 2, 1});
  p.setID(idShouldTouch++);
  vlContainer.addParticle(p);
  p.setR({2, 2, 2});
  p.setID(idShouldTouch++);
  vlContainer.addParticle(p);

  size_t idOffset = 1000;
  size_t idShouldNotTouch = idOffset;
  p.setR({2, 2, 2});
  p.setID(idShouldNotTouch++);
  vlContainer.addParticle(p);

  // move stuff
#ifdef AUTOPAS_OPENMP
#pragma omp parallel default(none) shared(vlContainer)
#endif
  for (auto pIter = vlContainer.begin(); pIter.isValid(); ++pIter) {
    switch (pIter->getID()) {
      case 0: {
        // move particle 0 into region
        pIter->setR({2, 2, 1.6});
        break;
      }
      case 1000: {
        // move particle 1000 out of region
        pIter->setR({2, 2, 3});
        break;
      }
    }
  }

  std::array<double, 3> regionOfInterstMin = {1.5, 1.5, 1.5};
  std::array<double, 3> regionOfInterstMax = {2.5, 2.5, 2.5};

  int particlesTouched = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : particlesTouched) default(none) shared(vlContainer, regionOfInterstMin, regionOfInterstMax)
#endif
  for (auto iterator = vlContainer.getRegionIterator(regionOfInterstMin, regionOfInterstMax); iterator.isValid();
       ++iterator) {
    iterator->touch();
    ++particlesTouched;
  }
  EXPECT_EQ(particlesTouched, idShouldTouch);

  int particlesChecked = 0;
  for (auto iterator = vlContainer.begin(); iterator.isValid(); ++iterator) {
    EXPECT_EQ(utils::inBox(iterator->getR(), regionOfInterstMin, regionOfInterstMax) ? 1 : 0, iterator->getNumTouched())
        << " particle at [" << utils::ArrayUtils::to_string(iterator->getR()) << "]" << std::endl;
    ++particlesChecked;
  }
  EXPECT_EQ(particlesChecked, idShouldTouch + idShouldNotTouch - idOffset);
}

void RegionParticleIteratorTest::checkTouches(LCTouch &lcContainer, std::array<double, 3> &regionMin,
                                              std::array<double, 3> &regionMax) {
  int numTouches = 0;
  // check the touch. Iterating over cells provides more debug info than normal iterator.
  for (size_t cellId = 0; cellId < lcContainer.getCells().size(); ++cellId) {
    const auto cellId3D = utils::ThreeDimensionalMapping::oneToThreeD(cellId, {7, 7, 7});
    for (auto pIter = lcContainer.getCells()[cellId].begin(); pIter.isValid(); ++pIter) {
      ++numTouches;
      EXPECT_EQ(utils::inBox(pIter->getR(), regionMin, regionMax) ? 1 : 0, pIter->getNumTouched())
          << "at: [" << utils::ArrayUtils::to_string(pIter->getR()) << "]" << std::endl
          << "in cell: " << cellId << " [" << utils::ArrayUtils::to_string(cellId3D) << "]";
    }
  }
  EXPECT_GE(numTouches, 0) << "No Particles were checked!";
}
