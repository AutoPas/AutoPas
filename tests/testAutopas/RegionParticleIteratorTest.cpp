/**
 * @file RegionParticleIteratorTest.cpp
 * @author seckler
 * @date 03.04.18
 */

#include "RegionParticleIteratorTest.h"

using namespace autopas;

/********************************** Linked Cells Tests **********************************/

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIterator) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>> lcContainer(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0));
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  // touch them using the regionIterator
  for (auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax); iterator.isValid(); ++iterator) {
    iterator->touch();
  }

  checkTouches(lcContainer, _regionMin, _regionMax);
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorOwned) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>> lcContainer(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0), 100);

  TouchableParticle part(ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), 100);
  lcContainer.addHaloParticle(part);

  // touch them using the regionIterator
  auto testRegionMin = ArrayMath::addScalar(_boxMin, -_cutoff * 0.5);
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iterator = lcContainer.getRegionIterator(testRegionMin, _regionMax, autopas::IteratorBehavior::ownedOnly);
       iterator.isValid(); ++iterator) {
    iterator->touch();
  }

  checkTouches(lcContainer, testRegionMin, _regionMax);
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorHalo) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>> lcContainer(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0), 100);

  TouchableParticle part(ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), 100);
  lcContainer.addHaloParticle(part);

  auto testRegionMin = ArrayMath::addScalar(_boxMin, -_cutoff * 0.5);
  // touch them using the regionIterator
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iterator = lcContainer.getRegionIterator(testRegionMin, _regionMax, autopas::IteratorBehavior::haloOnly);
       iterator.isValid(); ++iterator) {
    iterator->touch();
    EXPECT_TRUE(utils::inBox(iterator->getR(), testRegionMin, _regionMax)
              ? (utils::inBox(iterator->getR(), _boxMin, _boxMax) ? 0 : 1)
              : 0)
              << " particle at [" << iterator->getR()[0] << ", " << iterator->getR()[1] << ", " << iterator->getR()[2] << "]"
              << std::endl;
  }

  // check the touch using the normal iterator
  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    // this is a test for halo only! so we first check whether it's within our region of interest and then whether it's
    // not in the halo
    EXPECT_EQ(utils::inBox(iterator->getR(), testRegionMin, _regionMax)
                  ? (utils::inBox(iterator->getR(), _boxMin, _boxMax) ? 0 : 1)
                  : 0,
              iterator->getNumTouched())
        << " particle at [" << iterator->getR()[0] << ", " << iterator->getR()[1] << ", " << iterator->getR()[2] << "]"
        << std::endl;
  }
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorEmpty) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>> lcContainer(_boxMin, _boxMax, _cutoff);

  // add no particles

  // touch them using the regionIterator
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
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
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>> lcContainer(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0));

  auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax);
  auto iterator2 = iterator;

  // touch them using the regionIterator
  for (; iterator2.isValid(); ++iterator2) {
    iterator2->touch();
  }

  // check the touch using the normal iterator
  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    EXPECT_EQ(utils::inBox(iterator->getR(), _regionMin, _regionMax) ? 1 : 0, iterator->getNumTouched())
        << " particle at [" << iterator->getR()[0] << ", " << iterator->getR()[1] << ", " << iterator->getR()[2] << "]"
        << std::endl;
  }
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorCopyAssignment) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>> lcContainer(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0));

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

  // check the touch using the normal iterator
  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    EXPECT_EQ(utils::inBox(iterator->getR(), _regionMin, _regionMax) ? 3 : 0, iterator->getNumTouched())
        << " particle at [" << iterator->getR()[0] << ", " << iterator->getR()[1] << ", " << iterator->getR()[2] << "]"
        << std::endl;
  }
}

/**
 * Tests for correct iterator behavior when some of the cells in and outside of the region of interest are empty
 */
TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorSparseDomain) {
  // box goes from {0,0,0} to {5,5,5} + one halo layer
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>> lcContainer(_boxMin, _boxMax, _cutoff);

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
#pragma omp parallel reduction(+ : particlesTouched)
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
        << " particle at [" << iterator->getR()[0] << ", " << iterator->getR()[1] << ", " << iterator->getR()[2] << "]"
        << std::endl;
    ++particlesChecked;
  }
  EXPECT_EQ(particlesChecked, idShouldTouch + idShouldNotTouch - idOffset);
}

/*********************************** Direct Sum Tests ***********************************/

TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorSparseDomain) {
  // box goes from {0,0,0} to {5,5,5} + one halo layer
  DirectSum<TouchableParticle, FullParticleCell<TouchableParticle>> dsContainer(_boxMin, _boxMax, _cutoff);

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
#pragma omp parallel reduction(+ : particlesTouched)
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
        << " particle at [" << iterator->getR()[0] << ", " << iterator->getR()[1] << ", " << iterator->getR()[2] << "]"
        << std::endl;
    ++particlesChecked;
  }
  EXPECT_EQ(particlesChecked, idShouldTouch + idShouldNotTouch - idOffset);
}

TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIterator) {
  DirectSum<TouchableParticle, FullParticleCell<TouchableParticle>> container(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(container, TouchableParticle({0., 0., 0.}, 0));

  // touch them using the regionIterator
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
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
        << "at: [" << iterator->getR()[0] << ", " << iterator->getR()[1] << ", " << iterator->getR()[2] << "]";
  }
}

TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorBehaviorOwned) {
  DirectSum<TouchableParticle, FullParticleCell<TouchableParticle>> container(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(container, TouchableParticle({0., 0., 0.}, 0), 100);

  TouchableParticle part(ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), 100);
  container.addHaloParticle(part);

  // touch them using the regionIterator
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iterator = container.getRegionIterator(ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), _regionMax,
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
        << " particle at [" << iterator->getR()[0] << ", " << iterator->getR()[1] << ", " << iterator->getR()[2] << "]"
        << std::endl;
  }
}

TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorBehaviorHalo) {
  DirectSum<TouchableParticle, FullParticleCell<TouchableParticle>> container(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(container, TouchableParticle({0., 0., 0.}, 0), 100);

  TouchableParticle part(ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), 100);
  container.addHaloParticle(part);

  // touch them using the regionIterator
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iterator = container.getRegionIterator(ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), _regionMax,
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

    EXPECT_EQ(utils::inBox(iterator->getR(), ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), _regionMax)
                  ? (utils::inBox(iterator->getR(), _boxMin, _regionMax) ? 0 : 1)
                  : 0,
              iterator->getNumTouched())
        << " particle at [" << iterator->getR()[0] << ", " << iterator->getR()[1] << ", " << iterator->getR()[2] << "]"
        << std::endl;
  }
}

TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorEmpty) {
  DirectSum<TouchableParticle, FullParticleCell<TouchableParticle>> container(_boxMin, _boxMax, _cutoff);

  // add no particles

  int i = 0;
  // touch them using the regionIterator
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : i)
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
        << " particle at [" << iterator->getR()[0] << ", " << iterator->getR()[1] << ", " << iterator->getR()[2] << "]"
        << std::endl;
    i++;
  }
  EXPECT_EQ(i, 0);
}

TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorCopyConstructor) {
  DirectSum<TouchableParticle, FullParticleCell<TouchableParticle>> container(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(container, TouchableParticle({0., 0., 0.}, 0));

  auto iterator = container.getRegionIterator(_regionMin, _regionMax);
  auto iterator2 = iterator;

  // touch them using the regionIterator
  for (; iterator2.isValid(); ++iterator2) {
    iterator2->touch();
  }

  // check the touch using the normal iterator
  for (auto iterator = container.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    EXPECT_EQ(utils::inBox(iterator->getR(), _regionMin, _regionMax) ? 1 : 0, iterator->getNumTouched())
        << " particle at [" << iterator->getR()[0] << ", " << iterator->getR()[1] << ", " << iterator->getR()[2] << "]"
        << std::endl;
  }
}

TEST_F(RegionParticleIteratorTest, testDirectSumRegionParticleIteratorCopyAssignment) {
  DirectSum<TouchableParticle, FullParticleCell<TouchableParticle>> container(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(container, TouchableParticle({0., 0., 0.}, 0));

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

  // check the touch using the normal iterator
  for (auto iterator = container.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    EXPECT_EQ(utils::inBox(iterator->getR(), _regionMin, _regionMax) ? 3 : 0, iterator->getNumTouched())
        << " particle at [" << iterator->getR()[0] << ", " << iterator->getR()[1] << ", " << iterator->getR()[2] << "]"
        << std::endl;
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
#pragma omp parallel
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
#pragma omp parallel reduction(+ : particlesTouched)
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
        << " particle at [" << iterator->getR()[0] << ", " << iterator->getR()[1] << ", " << iterator->getR()[2] << "]"
        << std::endl;
    ++particlesChecked;
  }
  EXPECT_EQ(particlesChecked, idShouldTouch + idShouldNotTouch - idOffset);
}

void RegionParticleIteratorTest::checkTouches(LCTouch &lcContainer, std::array<double, 3> &regionMin,
                                              std::array<double, 3> &regionMax) {
  int numTouches = 0;
  // check the touch. Iterating over cells provides more debug info than normal iterator.
  for (size_t cellId = 0; cellId < lcContainer.getCells().size(); ++cellId) {
    auto cellId3D = utils::ThreeDimensionalMapping::oneToThreeD(cellId, {7, 7, 7});
    for (auto pIter = lcContainer.getCells()[cellId].begin(); pIter.isValid(); ++pIter) {
      ++numTouches;
      EXPECT_EQ(utils::inBox(pIter->getR(), regionMin, regionMax) ? 1 : 0, pIter->getNumTouched())
          << "at: [" << pIter->getR()[0] << ", " << pIter->getR()[1] << ", " << pIter->getR()[2] << "]" << std::endl
          << "in cell: " << cellId << " [" << cellId3D[0] << " | " << cellId3D[1] << " | " << cellId3D[2] << "]";
    }
  }
  EXPECT_GE(numTouches, 0) << "No Particles were checked!";
}
