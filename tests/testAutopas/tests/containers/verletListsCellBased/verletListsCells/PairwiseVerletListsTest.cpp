/**
 * @file PairwiseVerletListsTest.cpp
 * @author tirgendetwas
 * @date 22.11.20
 */
#include "PairwiseVerletListsTest.h"

#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"
using ::testing::_;
using ::testing::AtLeast;
using ::testing::Values;

/**
 * Tests if the correct number of interactions is taking place.
 * Sums up the lengths of all particles' neighbor lists
 * and compares the result (result additionally halved if useNewton3 is off) to the number of interactions.
 * This scenario is an interaction between two particles and the result should be 1.
 */
TEST_P(PairwiseVerletListsTest, testTwoParticles) {
  MockFunctor<Particle> emptyFunctor;
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  auto paramPair = GetParam();
  double cellSizeFactor = paramPair.first;
  bool useNewton3 = paramPair.second;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, useNewton3)).Times(AtLeast(1));
  autopas::VerletListsCells<Particle, autopas::VLCCellPairNeighborList<Particle>> verletLists(
      min, max, cutoff, autopas::TraversalOption::lc_c18, skin, cellSizeFactor);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  if (useNewton3 == true) {
    autopas::VLCC18Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true,
                             autopas::VLCCellPairNeighborList<Particle>, autopas::ContainerOption::pairwiseVerletLists>
        traversal(verletLists.getCellsPerDimension(), &emptyFunctor, verletLists.getInteractionLength(),
                  verletLists.getCellLength());

    verletLists.rebuildNeighborLists(&traversal);
    verletLists.iteratePairwise(&traversal);
  }

  else if (useNewton3 == false) {
    autopas::VLCC18Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false,
                             autopas::VLCCellPairNeighborList<Particle>, autopas::ContainerOption::pairwiseVerletLists>
        traversal(verletLists.getCellsPerDimension(), &emptyFunctor, verletLists.getInteractionLength(),
                  verletLists.getCellLength());

    verletLists.rebuildNeighborLists(&traversal);
    verletLists.iteratePairwise(&traversal);
  }

  std::vector<Particle *> list;
  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) list.push_back(&*iter);

  size_t partners = 0;
  for (auto &pl : list) {
    partners += verletLists.getNumberOfPartners(pl);
  }

  if (useNewton3 == false) {
    partners = partners / 2;
  }

  EXPECT_EQ(partners, 1);
}

/**
 * Tests if the correct number of interactions is taking place.
 * Sums up the lengths of all particles' neighbor lists
 * and compares the result (result additionally halved if useNewton3 is off) to the number of interactions.
 * This scenario includes three particles where two pairs are interacting,
 * while the third pair of particles are not within range of each other.
 * The result should be 2.
 */
TEST_P(PairwiseVerletListsTest, testThreeParticlesOneFar) {
  MockFunctor<Particle> emptyFunctorOther;
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {5, 5, 5};
  double cutoff = 1.;
  double skin = 0.2;
  auto paramPair = GetParam();
  double cellSizeFactor = paramPair.first;
  bool useNewton3 = paramPair.second;

  EXPECT_CALL(emptyFunctorOther, AoSFunctor(_, _, useNewton3)).Times(AtLeast(1));

  autopas::VerletListsCells<Particle, autopas::VLCCellPairNeighborList<Particle>> verletLists(
      min, max, cutoff, autopas::TraversalOption::lc_c18, skin, cellSizeFactor);

  std::array<double, 3> r = {2.5, 2.5, 2.5};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {2.5, 2.5, 3.5};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);
  std::array<double, 3> r3 = {3, 3, 3.5};
  Particle p3(r3, {0., 0., 0.}, 1);
  verletLists.addParticle(p3);

  if (useNewton3 == true) {
    autopas::VLCC18Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true,
                             autopas::VLCCellPairNeighborList<Particle>, autopas::ContainerOption::pairwiseVerletLists>
        traversal(verletLists.getCellsPerDimension(), &emptyFunctorOther, verletLists.getInteractionLength(),
                  verletLists.getCellLength());

    verletLists.rebuildNeighborLists(&traversal);
    verletLists.iteratePairwise(&traversal);
  }

  else if (useNewton3 == false) {
    autopas::VLCC18Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false,
                             autopas::VLCCellPairNeighborList<Particle>, autopas::ContainerOption::pairwiseVerletLists>
        traversal(verletLists.getCellsPerDimension(), &emptyFunctorOther, verletLists.getInteractionLength(),
                  verletLists.getCellLength());

    verletLists.rebuildNeighborLists(&traversal);
    verletLists.iteratePairwise(&traversal);
  }

  std::vector<Particle *> list;
  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) list.push_back(&*iter);

  size_t partners = 0;
  auto iter = verletLists.begin();
  while (iter.isValid()) {
    size_t current = verletLists.getNumberOfPartners(&*iter);
    partners += current;
    ++iter;
  }

  if (useNewton3 == false) {
    partners = partners / 2;
  }

  EXPECT_EQ(partners, 2);
}

/**
 * Tests if the correct number of interactions is taking place.
 * Sums up the lengths of all particles' neighbor lists
 * and compares the result (result additionally halved if useNewton3 is off) to the number of interactions.
 * This scenario includes three particles interacting with each other and the result should be 3.
 */
TEST_P(PairwiseVerletListsTest, testThreeParticlesClose) {
  MockFunctor<Particle> mock;
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {5, 5, 5};
  double cutoff = 1.;
  double skin = 0.2;
  auto paramPair = GetParam();
  double cellSizeFactor = paramPair.first;
  bool useNewton3 = paramPair.second;
  EXPECT_CALL(mock, AoSFunctor(_, _, useNewton3)).Times(AtLeast(1));
  autopas::VerletListsCells<Particle, autopas::VLCCellPairNeighborList<Particle>> verletLists(
      min, max, cutoff, autopas::TraversalOption::lc_c18, skin, cellSizeFactor);

  std::array<double, 3> r = {2.5, 2.5, 2.5};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {2.5, 2.5, 3.5};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);
  std::array<double, 3> r3 = {2, 2, 3};
  Particle p3(r3, {0., 0., 0.}, 1);
  verletLists.addParticle(p3);

  if (useNewton3 == true) {
    autopas::VLCC18Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true,
                             autopas::VLCCellPairNeighborList<Particle>, autopas::ContainerOption::pairwiseVerletLists>
        traversal(verletLists.getCellsPerDimension(), &mock, verletLists.getInteractionLength(),
                  verletLists.getCellLength());

    verletLists.rebuildNeighborLists(&traversal);
    verletLists.iteratePairwise(&traversal);
  }

  else if (useNewton3 == false) {
    autopas::VLCC18Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false,
                             autopas::VLCCellPairNeighborList<Particle>, autopas::ContainerOption::pairwiseVerletLists>
        traversal(verletLists.getCellsPerDimension(), &mock, verletLists.getInteractionLength(),
                  verletLists.getCellLength());

    verletLists.rebuildNeighborLists(&traversal);
    verletLists.iteratePairwise(&traversal);
  }

  std::vector<Particle *> list;
  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) list.push_back(&*iter);

  size_t partners = 0;
  auto iter = verletLists.begin();
  while (iter.isValid()) {
    size_t current = verletLists.getNumberOfPartners(&*iter);
    partners += current;
    ++iter;
  }

  if (useNewton3 == false) {
    partners = partners / 2;
  }

  EXPECT_EQ(partners, 3);
}

/**
 * Tests if the correct number of interactions is taking place.
 * Sums up the lengths of all particles' neighbor lists
 * and compares the result (result additionally halved if useNewton3 is off) to the number of interactions.
 * This scenario includes a single particle and the result should be 0.
 */
TEST_P(PairwiseVerletListsTest, testOneParticle) {
  MockFunctor<Particle> mock;
  // EXPECT_CALL(mock, AoSFunctor(_, _, true)); ?????
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {5, 5, 5};
  double cutoff = 1.;
  double skin = 0.2;
  auto paramPair = GetParam();
  double cellSizeFactor = paramPair.first;
  bool useNewton3 = paramPair.second;
  autopas::VerletListsCells<Particle, autopas::VLCCellPairNeighborList<Particle>> verletLists(
      min, max, cutoff, autopas::TraversalOption::lc_c18, skin, cellSizeFactor);

  std::array<double, 3> r = {2.5, 2.5, 2.5};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);

  if (useNewton3 == true) {
    autopas::VLCC18Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true,
                             autopas::VLCCellPairNeighborList<Particle>, autopas::ContainerOption::pairwiseVerletLists>
        traversal(verletLists.getCellsPerDimension(), &mock, verletLists.getInteractionLength(),
                  verletLists.getCellLength());

    verletLists.rebuildNeighborLists(&traversal);
    verletLists.iteratePairwise(&traversal);
  }

  else if (useNewton3 == false) {
    autopas::VLCC18Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false,
                             autopas::VLCCellPairNeighborList<Particle>, autopas::ContainerOption::pairwiseVerletLists>
        traversal(verletLists.getCellsPerDimension(), &mock, verletLists.getInteractionLength(),
                  verletLists.getCellLength());

    verletLists.rebuildNeighborLists(&traversal);
    verletLists.iteratePairwise(&traversal);
  }

  std::vector<Particle *> list;
  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) list.push_back(&*iter);

  size_t partners = 0;
  auto iter = verletLists.begin();
  while (iter.isValid()) {
    size_t current = verletLists.getNumberOfPartners(&*iter);
    partners += current;
    ++iter;
  }

  EXPECT_EQ(partners, 0);
}

INSTANTIATE_TEST_SUITE_P(Generated, PairwiseVerletListsTest,
                         Values(std::make_pair(1.0, true), std::make_pair(2.0, true), std::make_pair(1.0, false),
                                std::make_pair(2.0, false)),
                         PairwiseVerletListsTest::PrintToStringParamName());