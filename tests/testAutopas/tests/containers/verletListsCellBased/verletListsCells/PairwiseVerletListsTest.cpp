/**
 * @file PairwiseVerletListsTest.h
 * @author tirgendetwas
 * @date 22.11.20
 */
#include "PairwiseVerletListsTest.h"

#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"
using ::testing::_;
using ::testing::AtLeast;
using ::testing::Each;
using ::testing::Eq;
using ::testing::Invoke;
using ::testing::Values;

TEST_P(PairwiseVerletListsTest, testTwoParticles) {
  MockFunctor<Particle> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(1);
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  double cellSizeFactor = GetParam();
  autopas::VerletListsCells<Particle, autopas::PairwiseVerletNeighborList<Particle>> verletLists(
      min, max, cutoff, autopas::TraversalOption::lc_c18, skin, cellSizeFactor);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  autopas::VLCC18Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true,
                           typename autopas::VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType, 1>
      traversal(verletLists.getCellsPerDimension(), &emptyFunctor, verletLists.getInteractionLength(),
                verletLists.getCellLength());

  verletLists.rebuildNeighborLists(&traversal);
  verletLists.iteratePairwise(&traversal);

  std::vector<Particle *> list;
  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) list.push_back(&*iter);

  size_t partners = 0;
  for (auto &pl : list) {
    partners += verletLists.getNumberOfPartners(pl);
  }

  EXPECT_EQ(partners, 1);
}

TEST_P(PairwiseVerletListsTest, testThreeParticlesOneFar) {
  MockFunctor<Particle> emptyFunctorOther;
  EXPECT_CALL(emptyFunctorOther, AoSFunctor(_, _, true)).Times(2);
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {5, 5, 5};
  double cutoff = 1.;
  double skin = 0.2;
  double cellSizeFactor = GetParam();
  autopas::VerletListsCells<Particle, autopas::PairwiseVerletNeighborList<Particle>> verletLists(
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

  autopas::VLCC18Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true,
                           typename autopas::VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType, 1>
      traversal(verletLists.getCellsPerDimension(), &emptyFunctorOther, verletLists.getInteractionLength(),
                verletLists.getCellLength());

  verletLists.rebuildNeighborLists(&traversal);
  verletLists.iteratePairwise(&traversal);

  std::vector<Particle *> list;
  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) list.push_back(&*iter);

  size_t partners = 0;
  auto iter = verletLists.begin();
  while (iter.isValid()) {
    size_t current = verletLists.getNumberOfPartners(&*iter);
    std::cout << current << std::endl;
    partners += current;
    ++iter;
  }

  EXPECT_EQ(partners, 2);
}

TEST_P(PairwiseVerletListsTest, testThreeParticlesClose) {
  MockFunctor<Particle> mock;
  EXPECT_CALL(mock, AoSFunctor(_, _, true)).Times(3);
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {5, 5, 5};
  double cutoff = 1.;
  double skin = 0.2;
  double cellSizeFactor = GetParam();
  autopas::VerletListsCells<Particle, autopas::PairwiseVerletNeighborList<Particle>> verletLists(
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

  autopas::VLCC18Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true,
                           typename autopas::VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType, 1>
      traversal(verletLists.getCellsPerDimension(), &mock, verletLists.getInteractionLength(),
                verletLists.getCellLength());

  verletLists.rebuildNeighborLists(&traversal);
  verletLists.iteratePairwise(&traversal);

  std::vector<Particle *> list;
  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) list.push_back(&*iter);

  size_t partners = 0;
  auto iter = verletLists.begin();
  while (iter.isValid()) {
    size_t current = verletLists.getNumberOfPartners(&*iter);
    partners += current;
    ++iter;
  }

  EXPECT_EQ(partners, 3);
}

TEST_P(PairwiseVerletListsTest, testOneParticle) {
  MockFunctor<Particle> mock;
  // EXPECT_CALL(mock, AoSFunctor(_, _, true)).Times(3);
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {5, 5, 5};
  double cutoff = 1.;
  double skin = 0.2;
  double cellSizeFactor = GetParam();
  autopas::VerletListsCells<Particle, autopas::PairwiseVerletNeighborList<Particle>> verletLists(
      min, max, cutoff, autopas::TraversalOption::lc_c18, skin, cellSizeFactor);

  std::array<double, 3> r = {2.5, 2.5, 2.5};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);

  autopas::VLCC18Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true,
                           typename autopas::VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType, 1>
      traversal(verletLists.getCellsPerDimension(), &mock, verletLists.getInteractionLength(),
                verletLists.getCellLength());

  verletLists.rebuildNeighborLists(&traversal);
  verletLists.iteratePairwise(&traversal);

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

INSTANTIATE_TEST_SUITE_P(Generated, PairwiseVerletListsTest, Values(1.0, 2.0),
                         VerletListsTest::PrintToStringParamName());