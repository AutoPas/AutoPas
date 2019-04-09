/**
 * @file VerletListsCellsTest.cpp
 * @author nguyen
 * @date 02.09.18
 */

#include "VerletListsCellsTest.h"

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Invoke;

TEST_F(VerletListsCellsTest, VerletListConstructor) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletListsCells<Particle> verletLists(min, max, cutoff, autopas::TraversalOption::c18, skin);
}

TEST_F(VerletListsCellsTest, testVerletListBuild) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletListsCells<Particle> verletLists(min, max, cutoff, autopas::TraversalOption::c18, skin);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  autopas::C18TraversalVerlet<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> traversal(
      verletLists.getCellsPerDimension(), &emptyFunctor);
  verletLists.iteratePairwise(&emptyFunctor, &traversal, true);

  std::vector<Particle*> list;
  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) list.push_back(&*iter);

  EXPECT_EQ(list.size(), 2);
  int partners = 0;
  for (auto p : list) {
    partners += verletLists.getVerletList(p).size();
  }
  EXPECT_EQ(partners, 1);
}

TEST_F(VerletListsCellsTest, testVerletList) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletListsCells<Particle> verletLists(min, max, cutoff, autopas::TraversalOption::c18, skin);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> mockFunctor;
  using ::testing::_;  // anything is ok
  EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true));

  autopas::C18TraversalVerlet<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> traversal(
      verletLists.getCellsPerDimension(), &mockFunctor);
  verletLists.iteratePairwise(&mockFunctor, &traversal, true);

  std::vector<Particle*> list;
  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) list.push_back(&*iter);

  EXPECT_EQ(list.size(), 2);
  int partners = 0;
  for (auto p : list) {
    partners += verletLists.getVerletList(p).size();
  }
  EXPECT_EQ(partners, 1);
}
