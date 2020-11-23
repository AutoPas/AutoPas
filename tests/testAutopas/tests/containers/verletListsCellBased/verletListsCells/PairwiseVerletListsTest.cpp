#include "PairwiseVerletListsTest.h"

#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"
using ::testing::_;

TEST_F(PairwiseVerletListsTest, testTwoParticles)
{
  MockFunctor<Particle> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(1);
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletListsCells<Particle, autopas::PairwiseVerletNeighborList<Particle>> verletLists(min, max, cutoff, autopas::TraversalOption::lc_c18, skin,
                                                                                                   1.0);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  autopas::VLCC18Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true, typename autopas::VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType, false>
      traversal(verletLists.getCellsPerDimension(), &emptyFunctor, verletLists.getInteractionLength(), verletLists.getCellLength());

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