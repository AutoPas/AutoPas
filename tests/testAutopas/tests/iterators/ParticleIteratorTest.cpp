/**
 * @file ParticleIteratorTest.cpp
 * @author tchipev
 * @date 19.01.18
 */

#include "ParticleIteratorTest.h"

#include "autopas/containers/directSum/DirectSum.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/iterators/ParticleIterator.h"
#include "testingHelpers/commonTypedefs.h"

namespace ParticleIteratorTest {

using namespace autopas;
using namespace autopas::internal;

void ParticleIteratorTest::SetUp() {
  for (int i = 0; i < 20; ++i) {
    std::array<double, 3> arr{};
    for (auto &a : arr) {
      a = static_cast<double>(i);
    }
    Molecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), 0);
    _vecOfMolecules.push_back(m);
  }

  _currentIndex = 0;
}

void ParticleIteratorTest::TearDown() {}

TEST_F(ParticleIteratorTest, testFullIterator_EFEFFEEFEF) {
  // Empty Full Empty Full Full Empty Empty Full Empty Full
  std::vector<FMCell> data(10);

  for (auto i : {1ul, 3ul, 4ul, 7ul, 9ul}) {
    fillWithParticles(&data.at(i));
  }

  std::vector<size_t> foundParticles;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(vecMerge : foundParticles)
#endif
  {
    for (auto iter = ParticleIterator<Molecule, FMCell, true>(&data); iter.isValid(); ++iter) {
      auto particleID = iter->getID();
      foundParticles.push_back(particleID);
      for (int d = 0; d < 3; ++d) {
        EXPECT_DOUBLE_EQ(iter->getR()[d], particleID);
      }
    }
  }

  ASSERT_EQ(foundParticles.size(), 20);
  std::sort(foundParticles.begin(), foundParticles.end());
  for (size_t i = 0; i < 20; ++i) {
    ASSERT_EQ(i, foundParticles[i]);
  }
}

TEST_F(ParticleIteratorTest, testFullIterator_FEFEEFFEFE) {
  // Full Empty Full Empty Empty Full Full Empty Full Empty
  std::vector<FMCell> data(10);

  for (auto i : {0ul, 2ul, 5ul, 6ul, 8ul}) {
    fillWithParticles(&data.at(i));
  }

  std::vector<size_t> foundParticles;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(vecMerge : foundParticles)
#endif
  {
    for (auto iter = ParticleIterator<Molecule, FMCell, true>(&data); iter.isValid(); ++iter) {
      auto particleID = iter->getID();
      foundParticles.push_back(particleID);
      for (int d = 0; d < 3; ++d) {
        EXPECT_DOUBLE_EQ(iter->getR()[d], particleID);
      }
    }
  }

  ASSERT_EQ(foundParticles.size(), 20);
  std::sort(foundParticles.begin(), foundParticles.end());
  for (size_t i = 0; i < 20; ++i) {
    ASSERT_EQ(i, foundParticles[i]);
  }
}

TEST_F(ParticleIteratorTest, testFullIterator_deletion) {
  // Full Empty Full Empty Empty Full Full Empty Full Empty
  std::vector<FMCell> data(10);

  for (auto i : {0ul, 2ul, 5ul, 6ul, 8ul}) {
    fillWithParticles(&data.at(i));
  }

  int numFoundParticles = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : numFoundParticles)
#endif
  {
    ParticleIterator<Molecule, FMCell, true> iter(&data);
    for (; iter.isValid(); ++iter, ++numFoundParticles) {
      autopas::internal::deleteParticle(iter);
    }
  }
  ASSERT_EQ(numFoundParticles, 20);

  numFoundParticles = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : numFoundParticles)
#endif
  {
    ParticleIterator<Molecule, FMCell, true> iter(&data);
    for (; iter.isValid(); ++iter) {
      ++numFoundParticles;
    }
  }
  ASSERT_EQ(numFoundParticles, 0);
}

TEST_F(ParticleIteratorTest, testFullIterator_mutable) {
  // Full Empty Full Empty Empty Full Full Empty Full Empty
  std::vector<FMCell> data(10);

  for (auto i : {0ul, 2ul, 5ul, 6ul, 8ul}) {
    fillWithParticles(&data.at(i));
  }

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  {
    ParticleIterator<Molecule, FMCell, true> iter(&data);
    for (; iter.isValid(); ++iter) {
      double newVel = iter->getID() + 1;
      std::array<double, 3> newVelArr = {newVel, newVel, newVel};
      iter->setV(newVelArr);
    }
  }

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  {
    ParticleIterator<Molecule, FMCell, true> iter(&data);
    for (; iter.isValid(); ++iter) {
      double expectedVel = iter->getID() + 1;
      auto vel = iter->getV();
      EXPECT_EQ(vel[0], expectedVel);
      EXPECT_EQ(vel[1], expectedVel);
      EXPECT_EQ(vel[2], expectedVel);
    }
  }
}

/**
 * Test the iterator behavior for owning and halo molecules.
 * this function expects that two molecules are already added to container:
 * mol should be added as normal particle.
 * haloMol as haloParticle.
 * @tparam Container
 * @tparam Molecule
 * @param container should have already an added mol (as owning molecule) and
 * haloMol (als halo molecule)
 * @param mol
 * @param haloMol
 */
template <class Container, class Molecule>
void testContainerIteratorBehavior(Container &container, Molecule &mol, Molecule &haloMol) {
  // default
  int count = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : count)
#endif
  for (auto iter = container.begin(); iter.isValid(); ++iter) {
    count++;
  }
  EXPECT_EQ(count, 2);

  // haloAndOwned (same as default)
  count = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : count)
#endif
  for (auto iter = container.begin(IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
    count++;
  }
  EXPECT_EQ(count, 2);

  // owned only
  count = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : count)
#endif
  for (auto iter = container.begin(IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    count++;
    EXPECT_EQ(iter->getID(), mol.getID());
  }
  EXPECT_EQ(count, 1);

  // halo only
  count = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : count)
#endif
  for (auto iter = container.begin(IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
    count++;
    EXPECT_EQ(iter->getID(), haloMol.getID());
  }
  EXPECT_EQ(count, 1);
}

TEST_F(ParticleIteratorTest, testIteratorBehaviorDirectSum) {
  DirectSum<Molecule> ds({0., 0., 0.}, {10., 10., 10.}, 3, 0.);
  Molecule mol({1., 1., 1.}, {0., 0., 0.}, 1);
  ds.addParticle(mol);
  Molecule haloMol({-1., 1., 1.}, {0., 0., 0.}, 2, 0);
  ds.addHaloParticle(haloMol);

  testContainerIteratorBehavior(ds, mol, haloMol);
}

TEST_F(ParticleIteratorTest, testIteratorBehaviorLinkedCells) {
  LinkedCells<Molecule> linkedCells({0., 0., 0.}, {10., 10., 10.}, 3, 0., 1.);
  Molecule mol({1., 1., 1.}, {0., 0., 0.}, 1);
  linkedCells.addParticle(mol);
  Molecule haloMol({-1., 1., 1.}, {0., 0., 0.}, 2, 0);
  linkedCells.addHaloParticle(haloMol);

  testContainerIteratorBehavior(linkedCells, mol, haloMol);
}

TEST_F(ParticleIteratorTest, testIteratorBehaviorVerletLists) {
  VerletLists<Molecule> verletLists({0., 0., 0.}, {10., 10., 10.}, 3, 0.);
  Molecule mol({1., 1., 1.}, {0., 0., 0.}, 1, 0);
  verletLists.addParticle(mol);
  Molecule haloMol({-1., 1., 1.}, {0., 0., 0.}, 2, 0);
  verletLists.addHaloParticle(haloMol);

  // test normally
  testContainerIteratorBehavior(verletLists, mol, haloMol);

  // swap everything around, test if it still valid :)
  haloMol.setR({1., 1., 1.});
  mol.setR({-1., 1., 1.});
  verletLists.begin(IteratorBehavior::ownedOnly)->setR({-1., 1., 1.});
  verletLists.begin(IteratorBehavior::haloOnly)->setR({1., 1., 1.});

  testContainerIteratorBehavior(verletLists, mol, haloMol);
}
}  // end namespace ParticleIteratorTest
