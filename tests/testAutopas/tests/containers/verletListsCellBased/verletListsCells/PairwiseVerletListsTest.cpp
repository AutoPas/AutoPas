/**
 * @file PairwiseVerletListsTest.cpp
 * @author tirgendetwas
 * @date 22.11.20
 */
#include "PairwiseVerletListsTest.h"

#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCCellPairNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCC18Traversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"

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
  MockPairwiseFunctor<ParticleFP64> emptyFunctor;
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  auto params = GetParam();
  double cellSizeFactor = std::get<0>(params);
  bool useNewton3 = std::get<1>(params);
  auto buildType = std::get<2>(params);
  const autopas::LoadEstimatorOption loadEstimator = autopas::LoadEstimatorOption::none;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, useNewton3)).Times(AtLeast(1));
  autopas::VerletListsCells<ParticleFP64, autopas::VLCCellPairNeighborList<ParticleFP64>> verletLists(
      min, max, cutoff, skin, cellSizeFactor, loadEstimator, buildType);

  std::array<double, 3> r = {2, 2, 2};
  ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  autopas::VLCC18Traversal<FPCell, MPairwiseFunctor, autopas::VLCCellPairNeighborList<ParticleFP64>> traversal(
      verletLists.getCellsPerDimension(), &emptyFunctor, verletLists.getInteractionLength(),
      verletLists.getCellLength(), autopas::DataLayoutOption::aos, useNewton3,
      autopas::ContainerOption::pairwiseVerletLists);

  verletLists.rebuildNeighborLists(&traversal);
  verletLists.computeInteractions(&traversal);

  std::vector<ParticleFP64 *> list;
  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) list.push_back(&*iter);

  size_t partners = 0;
  for (auto &pl : list) {
    partners += verletLists.getNumberOfPartners(pl);
  }

  if (not useNewton3) {
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
  MockPairwiseFunctor<ParticleFP64> emptyFunctorOther;
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {5, 5, 5};
  double cutoff = 1.;
  double skin = 0.2;
  auto params = GetParam();
  double cellSizeFactor = std::get<0>(params);
  bool useNewton3 = std::get<1>(params);
  auto buildType = std::get<2>(params);
  const autopas::LoadEstimatorOption loadEstimator = autopas::LoadEstimatorOption::none;

  EXPECT_CALL(emptyFunctorOther, AoSFunctor(_, _, useNewton3)).Times(AtLeast(1));

  autopas::VerletListsCells<ParticleFP64, autopas::VLCCellPairNeighborList<ParticleFP64>> verletLists(
      min, max, cutoff, skin, cellSizeFactor, loadEstimator, buildType);

  std::array<double, 3> r = {2.5, 2.5, 2.5};
  ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {2.5, 2.5, 3.5};
  ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);
  std::array<double, 3> r3 = {3, 3, 3.5};
  ParticleFP64 p3(r3, {0., 0., 0.}, 1);
  verletLists.addParticle(p3);

  autopas::VLCC18Traversal<FPCell, MPairwiseFunctor, autopas::VLCCellPairNeighborList<ParticleFP64>> traversal(
      verletLists.getCellsPerDimension(), &emptyFunctorOther, verletLists.getInteractionLength(),
      verletLists.getCellLength(), autopas::DataLayoutOption::aos, useNewton3,
      autopas::ContainerOption::pairwiseVerletLists);

  verletLists.rebuildNeighborLists(&traversal);
  verletLists.computeInteractions(&traversal);

  std::vector<ParticleFP64 *> list;
  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) list.push_back(&*iter);

  size_t partners = 0;
  auto iter = verletLists.begin();
  while (iter.isValid()) {
    size_t current = verletLists.getNumberOfPartners(&*iter);
    partners += current;
    ++iter;
  }

  if (not useNewton3) {
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
  MockPairwiseFunctor<ParticleFP64> mock;
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {5, 5, 5};
  double cutoff = 1.;
  double skin = 0.2;
  auto params = GetParam();
  double cellSizeFactor = std::get<0>(params);
  bool useNewton3 = std::get<1>(params);
  auto buildType = std::get<2>(params);
  const autopas::LoadEstimatorOption loadEstimator = autopas::LoadEstimatorOption::none;

  EXPECT_CALL(mock, AoSFunctor(_, _, useNewton3)).Times(AtLeast(1));
  autopas::VerletListsCells<ParticleFP64, autopas::VLCCellPairNeighborList<ParticleFP64>> verletLists(
      min, max, cutoff, skin, cellSizeFactor, loadEstimator, buildType);

  std::array<double, 3> r = {2.5, 2.5, 2.5};
  ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {2.5, 2.5, 3.5};
  ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);
  std::array<double, 3> r3 = {2, 2, 3};
  ParticleFP64 p3(r3, {0., 0., 0.}, 1);
  verletLists.addParticle(p3);

  autopas::VLCC18Traversal<FPCell, MPairwiseFunctor, autopas::VLCCellPairNeighborList<ParticleFP64>> traversal(
      verletLists.getCellsPerDimension(), &mock, verletLists.getInteractionLength(), verletLists.getCellLength(),
      autopas::DataLayoutOption::aos, useNewton3, autopas::ContainerOption::pairwiseVerletLists);

  verletLists.rebuildNeighborLists(&traversal);
  verletLists.computeInteractions(&traversal);

  std::vector<ParticleFP64 *> list;
  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) list.push_back(&*iter);

  size_t partners = 0;
  auto iter = verletLists.begin();
  while (iter.isValid()) {
    size_t current = verletLists.getNumberOfPartners(&*iter);
    partners += current;
    ++iter;
  }

  if (not useNewton3) {
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
  MockPairwiseFunctor<ParticleFP64> mock;
  // EXPECT_CALL(mock, AoSFunctor(_, _, true)); ?????
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {5, 5, 5};
  double cutoff = 1.;
  double skin = 0.2;
  auto params = GetParam();
  double cellSizeFactor = std::get<0>(params);
  bool useNewton3 = std::get<1>(params);
  auto buildType = std::get<2>(params);
  const autopas::LoadEstimatorOption loadEstimator = autopas::LoadEstimatorOption::none;
  autopas::VerletListsCells<ParticleFP64, autopas::VLCCellPairNeighborList<ParticleFP64>> verletLists(
      min, max, cutoff, skin, cellSizeFactor, loadEstimator, buildType);

  std::array<double, 3> r = {2.5, 2.5, 2.5};
  ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);

  autopas::VLCC18Traversal<FPCell, MPairwiseFunctor, autopas::VLCCellPairNeighborList<ParticleFP64>> traversal(
      verletLists.getCellsPerDimension(), &mock, verletLists.getInteractionLength(), verletLists.getCellLength(),
      autopas::DataLayoutOption::aos, useNewton3, autopas::ContainerOption::pairwiseVerletLists);

  verletLists.rebuildNeighborLists(&traversal);
  verletLists.computeInteractions(&traversal);

  std::vector<ParticleFP64 *> list;
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

/**
 * Fills an AoS and an SoA list with particles, executes one iteration
 * and compares the positions of each pair of corresponding particles.
 * */
TEST_P(PairwiseVerletListsTest, SoAvsAoSLJ) {
  const double cutoff = 2.;
  auto params = GetParam();
  double cellSizeFactor = std::get<0>(params);
  bool useNewton3 = std::get<1>(params);

  // changing type from Particle to Molecule
  auto oldBuildType = std::get<2>(params);
  auto buildType = autopas::VerletListsCellsHelpers::VLCBuildType::aosBuild;
  if (oldBuildType == autopas::VerletListsCellsHelpers::VLCBuildType::soaBuild) {
    buildType = autopas::VerletListsCellsHelpers::VLCBuildType::soaBuild;
  }

  const autopas::LoadEstimatorOption loadEstimator = autopas::LoadEstimatorOption::none;
  std::array<double, 3> min = {0, 0, 0};
  std::array<double, 3> max = {10, 10, 10};
  autopas::VerletListsCells<Molecule, autopas::VLCCellPairNeighborList<Molecule>> verletLists1(
      min, max, cutoff, 0.01, cellSizeFactor, loadEstimator, buildType);
  autopas::VerletListsCells<Molecule, autopas::VLCCellPairNeighborList<Molecule>> verletLists2(
      min, max, cutoff, 0.01, cellSizeFactor, loadEstimator, buildType);

  Molecule defaultParticle({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(verletLists1, defaultParticle, verletLists1.getBoxMin(),
                                                                verletLists1.getBoxMax(), 100);
  autopasTools::generators::UniformGenerator::fillWithParticles(verletLists2, defaultParticle, verletLists2.getBoxMin(),
                                                                verletLists2.getBoxMax(), 100);
  LJFunctorType<> ljFunctor(cutoff);
  ljFunctor.setParticleProperties(1., 1.);

  autopas::VLCC18Traversal<FMCell, LJFunctorType<>, autopas::VLCCellPairNeighborList<Molecule>> verletTraversal1(
      verletLists1.getCellsPerDimension(), &ljFunctor, verletLists1.getInteractionLength(),
      verletLists1.getCellLength(), autopas::DataLayoutOption::aos, useNewton3,
      autopas::ContainerOption::Value::pairwiseVerletLists);
  autopas::VLCC18Traversal<FMCell, LJFunctorType<>, autopas::VLCCellPairNeighborList<Molecule>> soaTraversal(
      verletLists2.getCellsPerDimension(), &ljFunctor, verletLists2.getInteractionLength(),
      verletLists2.getCellLength(), autopas::DataLayoutOption::soa, useNewton3,
      autopas::ContainerOption::Value::pairwiseVerletLists);

  verletLists1.rebuildNeighborLists(&verletTraversal1);
  verletLists2.rebuildNeighborLists(&soaTraversal);
  verletLists1.computeInteractions(&verletTraversal1);
  verletLists2.computeInteractions(&soaTraversal);

  auto iter1 = verletLists1.begin();
  auto iter2 = verletLists2.begin();

  for (; iter1.isValid() && iter2.isValid(); ++iter1, ++iter2) {
    for (unsigned int dim = 0; dim < 3; dim++) {
      ASSERT_NEAR(iter1->getR()[dim], iter2->getR()[dim], fabs(iter1->getR()[dim] * 1e-7));
      EXPECT_NEAR(iter1->getF()[dim], iter2->getF()[dim], fabs(iter1->getF()[dim] * 1e-7));
    }
  }
  EXPECT_FALSE(iter1.isValid());
  EXPECT_FALSE(iter2.isValid());
}

INSTANTIATE_TEST_SUITE_P(Generated, PairwiseVerletListsTest,
                         Values(std::make_tuple(1.0, true, autopas::VerletListsCellsHelpers::VLCBuildType::aosBuild),
                                std::make_tuple(2.0, true, autopas::VerletListsCellsHelpers::VLCBuildType::aosBuild),
                                std::make_tuple(1.0, false, autopas::VerletListsCellsHelpers::VLCBuildType::aosBuild),
                                std::make_tuple(2.0, false, autopas::VerletListsCellsHelpers::VLCBuildType::aosBuild),
                                std::make_tuple(1.0, true, autopas::VerletListsCellsHelpers::VLCBuildType::soaBuild),
                                std::make_tuple(2.0, true, autopas::VerletListsCellsHelpers::VLCBuildType::soaBuild),
                                std::make_tuple(1.0, false, autopas::VerletListsCellsHelpers::VLCBuildType::soaBuild),
                                std::make_tuple(2.0, false, autopas::VerletListsCellsHelpers::VLCBuildType::soaBuild)),
                         PairwiseVerletListsTest::PrintToStringParamName());