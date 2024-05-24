/**
 * @file VerletListsCellsTest.cpp
 * @author nguyen
 * @date 02.09.18
 */
#include "VerletListsCellsTest.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCAllCellsNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCC18Traversal.h"
#include "autopas/particles/Particle.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "mocks/MockFunctor.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::_;
using ::testing::AtLeast;

void applyFunctor(MockFunctor<Particle> &functor, const double cellSizefactor,
                  autopas::VerletListsCellsHelpers<Particle>::VLCBuildType buildType) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skinPerTimestep = 0.01;
  unsigned int rebuildFrequency = 20;
  const autopas::LoadEstimatorOption loadEstimator = autopas::LoadEstimatorOption::none;
  autopas::VerletListsCells<Particle, autopas::VLCAllCellsNeighborList<Particle>> verletLists(
      min, max, cutoff, skinPerTimestep, rebuildFrequency, cellSizefactor, loadEstimator, buildType);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  autopas::VLCC18Traversal<FPCell, MFunctor, autopas::VLCAllCellsNeighborList<Particle>> traversal(
      verletLists.getCellsPerDimension(), &functor, verletLists.getInteractionLength(), verletLists.getCellLength(),
      autopas::DataLayoutOption::aos, true, autopas::ContainerOption::verletListsCells);

  verletLists.rebuildNeighborLists(&traversal);
  verletLists.iteratePairwise(&traversal);

  std::vector<Particle *> list;
  for (auto iter = verletLists.begin(autopas::IteratorBehavior::ownedOrHalo); iter.isValid(); ++iter) {
    list.push_back(&*iter);
  }

  EXPECT_EQ(list.size(), 2);
  size_t partners = 0;

  for (auto &pl : list) {
    partners += verletLists.getNumberOfPartners(pl);
  }
  EXPECT_EQ(partners, 1);
}

/**
 * Fills an AoS and an SoA list with particles, executes one iteration
 * and compares the positions of each pair of corresponding particles.
 * */
void soaTest(const double cellSizeFactor, autopas::VerletListsCellsHelpers<Particle>::VLCBuildType oldBuildType) {
  const double cutoff = 2.;
  const autopas::LoadEstimatorOption loadEstimator = autopas::LoadEstimatorOption::none;
  std::array<double, 3> min = {0, 0, 0};
  std::array<double, 3> max = {10, 10, 10};

  auto buildType = autopas::VerletListsCellsHelpers<Molecule>::VLCBuildType::aosBuild;
  if (oldBuildType == autopas::VerletListsCellsHelpers<Particle>::VLCBuildType::soaBuild) {
    buildType = autopas::VerletListsCellsHelpers<Molecule>::VLCBuildType::soaBuild;
  }

  autopas::VerletListsCells<Molecule, autopas::VLCCellPairNeighborList<Molecule>> verletLists1(
      min, max, cutoff, 0.01, 30, cellSizeFactor, loadEstimator, buildType);
  autopas::VerletListsCells<Molecule, autopas::VLCCellPairNeighborList<Molecule>> verletLists2(
      min, max, cutoff, 0.01, 30, cellSizeFactor, loadEstimator, buildType);

  Molecule defaultParticle({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(verletLists1, defaultParticle, verletLists1.getBoxMin(),
                                                               verletLists1.getBoxMax(), 100);
  autopasTools::generators::RandomGenerator::fillWithParticles(verletLists2, defaultParticle, verletLists2.getBoxMin(),
                                                               verletLists2.getBoxMax(), 100);
  LJFunctorType<> ljFunctor(cutoff);
  ljFunctor.setParticleProperties(1., 1.);

  autopas::VLCC18Traversal<FMCell, LJFunctorType<>, autopas::VLCCellPairNeighborList<Molecule>>
      verletTraversal1(verletLists1.getCellsPerDimension(), &ljFunctor, verletLists1.getInteractionLength(),
                       verletLists1.getCellLength(), autopas::DataLayoutOption::aos, true,
                       autopas::ContainerOption::verletListsCells);
  autopas::VLCC18Traversal<FMCell, LJFunctorType<>, autopas::VLCCellPairNeighborList<Molecule>> soaTraversal(
      verletLists2.getCellsPerDimension(), &ljFunctor, verletLists2.getInteractionLength(),
      verletLists2.getCellLength(), autopas::DataLayoutOption::soa, true, autopas::ContainerOption::verletListsCells);

  verletLists1.rebuildNeighborLists(&verletTraversal1);
  verletLists2.rebuildNeighborLists(&soaTraversal);
  verletLists1.iteratePairwise(&verletTraversal1);
  verletLists2.iteratePairwise(&soaTraversal);

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

TEST_F(VerletListsCellsTest, testVerletListBuild) {
  MockFunctor<Particle> emptyFunctorAoSBuild;
  EXPECT_CALL(emptyFunctorAoSBuild, AoSFunctor(_, _, true)).Times(1);

  applyFunctor(emptyFunctorAoSBuild, 1.0, autopas::VerletListsCellsHelpers<Particle>::VLCBuildType::aosBuild);

  MockFunctor<Particle> emptyFunctorSoABuild;
  EXPECT_CALL(emptyFunctorSoABuild, AoSFunctor(_, _, true)).Times(1);

  applyFunctor(emptyFunctorSoABuild, 1.0, autopas::VerletListsCellsHelpers<Particle>::VLCBuildType::soaBuild);

  MockFunctor<Particle> emptyFunctorAoSBuild_cs2;
  EXPECT_CALL(emptyFunctorAoSBuild_cs2, AoSFunctor(_, _, true)).Times(1);

  applyFunctor(emptyFunctorAoSBuild_cs2, 2.0, autopas::VerletListsCellsHelpers<Particle>::VLCBuildType::aosBuild);

  MockFunctor<Particle> emptyFunctorSoABuild_cs2;
  EXPECT_CALL(emptyFunctorSoABuild_cs2, AoSFunctor(_, _, true)).Times(1);

  applyFunctor(emptyFunctorSoABuild_cs2, 2.0, autopas::VerletListsCellsHelpers<Particle>::VLCBuildType::soaBuild);

  soaTest(1.0, autopas::VerletListsCellsHelpers<Particle>::VLCBuildType::soaBuild);
  soaTest(2.0, autopas::VerletListsCellsHelpers<Particle>::VLCBuildType::soaBuild);
  soaTest(1.0, autopas::VerletListsCellsHelpers<Particle>::VLCBuildType::aosBuild);
  soaTest(2.0, autopas::VerletListsCellsHelpers<Particle>::VLCBuildType::aosBuild);
}
