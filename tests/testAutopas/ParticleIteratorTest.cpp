/*
 * ParticleIteratorTest.cpp
 *
 *  Created on: 19 Jan 2018
 *      Author: tchipevn
 */

#include "ParticleIteratorTest.h"
#include "autopasIncludes.h"
#include "gtest/gtest.h"

using namespace autopas;
using namespace autopas::internal;

void ParticleIteratorTest::SetUp() {
  for (int i = 0; i < 20; ++i) {
    std::array<double, 3> arr{};
    for (auto& a : arr) {
      a = static_cast<double>(i);
    }
    MoleculeLJ m(arr, {0., 0., 0.}, static_cast<unsigned long>(i));
    _vecOfMolecules.push_back(m);
  }

  _currentIndex = 0;
}

void ParticleIteratorTest::TearDown() {}

TEST_F(ParticleIteratorTest, testFullIterator_EFEFFEEFEF) {
  // Empty Full Empty Full Full Empty Empty Full Empty Full
  std::vector<FullParticleCell<MoleculeLJ>> data(10);

  for (auto i : {1ul, 3ul, 4ul, 7ul, 9ul}) {
    fillWithParticles(&data.at(i));
  }

  ParticleIterator<MoleculeLJ, FullParticleCell<MoleculeLJ>> iter(&data);
  int i = 0;
  for (; iter.isValid(); ++iter, ++i) {
    for (int d = 0; d < 3; ++d) {
      ASSERT_DOUBLE_EQ(iter->getR()[d], _vecOfMolecules[i].getR()[d]);
    }
    ASSERT_EQ(iter->getID(), _vecOfMolecules[i].getID());
  }
}

TEST_F(ParticleIteratorTest, testFullIterator_FEFEEFFEFE) {
  // Full Empty Full Empty Empty Full Full Empty Full Empty
  std::vector<FullParticleCell<MoleculeLJ>> data(10);

  for (auto i : {0ul, 2ul, 5ul, 6ul, 8ul}) {
    fillWithParticles(&data.at(i));
  }

  ParticleIterator<MoleculeLJ, FullParticleCell<MoleculeLJ>> iter(&data);
  int i = 0;
  for (; iter.isValid(); ++iter, ++i) {
    for (int d = 0; d < 3; ++d) {
      ASSERT_DOUBLE_EQ(iter->getR()[d], _vecOfMolecules[i].getR()[d]);
    }
    ASSERT_EQ(iter->getID(), _vecOfMolecules[i].getID());
  }
}

TEST_F(ParticleIteratorTest, testFullIterator_deletion) {
  // Full Empty Full Empty Empty Full Full Empty Full Empty
  std::vector<FullParticleCell<MoleculeLJ>> data(10);

  for (auto i : {0ul, 2ul, 5ul, 6ul, 8ul}) {
    fillWithParticles(&data.at(i));
  }
  {
    ParticleIterator<MoleculeLJ, FullParticleCell<MoleculeLJ>> iter(&data);
    int i = 0;
    for (; iter.isValid(); ++iter, ++i) {
      iter.deleteCurrentParticle();
    }
    ASSERT_EQ(i, 20);
  }

  {
    ParticleIterator<MoleculeLJ, FullParticleCell<MoleculeLJ>> iter(&data);
    int i = 0;
    for (; iter.isValid(); ++iter) {
      ++i;
    }
    ASSERT_EQ(i, 0);
  }
}

TEST_F(ParticleIteratorTest, testFullIterator_mutable) {
  // Full Empty Full Empty Empty Full Full Empty Full Empty
  std::vector<FullParticleCell<MoleculeLJ>> data(10);

  for (auto i : {0ul, 2ul, 5ul, 6ul, 8ul}) {
    fillWithParticles(&data.at(i));
  }
  {
    ParticleIterator<MoleculeLJ, FullParticleCell<MoleculeLJ>> iter(&data);
    double i = 0.;
    for (; iter.isValid(); ++iter, ++i) {
      iter->setV({i, i, i});
    }
  }

  {
    ParticleIterator<MoleculeLJ, FullParticleCell<MoleculeLJ>> iter(&data);
    double i = 0.;
    for (; iter.isValid(); ++iter, ++i) {
      auto vel = iter->getV();
      ASSERT_EQ(vel[0], i);
      ASSERT_EQ(vel[1], i);
      ASSERT_EQ(vel[2], i);
    }
  }
}

TEST_F(ParticleIteratorTest, testRMMIterator_EFEFFEEFEF) {
  // Empty Full Empty Full Full Empty Empty Full Empty Full
  std::vector<RMMParticleCell<MoleculeLJ>> data(10);

  for (auto i : {1u, 3u, 4u, 7u, 9u}) {
    fillWithParticles(&data.at(i));
  }

  ParticleIterator<MoleculeLJ, RMMParticleCell<MoleculeLJ>> iter(&data);
  int i = 0;
  for (; iter.isValid(); ++iter, ++i) {
    for (int d = 0; d < 3; ++d) {
      ASSERT_DOUBLE_EQ(iter->getR()[d], _vecOfMolecules[i].getR()[d]);
    }
    //		ASSERT_EQ(iter->getID(), _vecOfMolecules[i].getID());
  }
}

TEST_F(ParticleIteratorTest, testRMMIterator_FEFEEFFEFE) {
  // Full Empty Full Empty Empty Full Full Empty Full Empty
  std::vector<RMMParticleCell<MoleculeLJ>> data(10);

  for (auto i : {0u, 2u, 5u, 6u, 8u}) {
    fillWithParticles(&data.at(i));
  }

  ParticleIterator<MoleculeLJ, RMMParticleCell<MoleculeLJ>> iter(&data);
  int i = 0;
  for (; iter.isValid(); ++iter, ++i) {
    for (int d = 0; d < 3; ++d) {
      ASSERT_DOUBLE_EQ(iter->getR()[d], _vecOfMolecules[i].getR()[d]);
    }
    //		ASSERT_EQ(iter->getID(), _vecOfMolecules[i].getID());
  }
}

TEST_F(ParticleIteratorTest, testRMMIterator_deletion) {
  // Full Empty Full Empty Empty Full Full Empty Full Empty
  std::vector<RMMParticleCell<MoleculeLJ>> data(10);

  for (auto i : {0u, 2u, 5u, 6u, 8u}) {
    fillWithParticles(&data.at(i));
  }
  {
    ParticleIterator<MoleculeLJ, RMMParticleCell<MoleculeLJ>> iter(&data);
    int i = 0;
    for (; iter.isValid(); ++iter, ++i) {
      iter.deleteCurrentParticle();
    }
    ASSERT_EQ(i, 20);
  }
  {
    ParticleIterator<MoleculeLJ, RMMParticleCell<MoleculeLJ>> iter(&data);
    int i = 0;
    for (; iter.isValid(); ++iter) {
      ++i;
    }
    ASSERT_EQ(i, 0);
  }
}

TEST_F(ParticleIteratorTest, testRMMIterator_mutable) {
  // Full Empty Full Empty Empty Full Full Empty Full Empty
  std::vector<RMMParticleCell<MoleculeLJ>> data(10);

  for (auto i : {0ul, 2ul, 5ul, 6ul, 8ul}) {
    fillWithParticles(&data.at(i));
  }
  {
    ParticleIterator<MoleculeLJ, RMMParticleCell<MoleculeLJ>> iter(&data);
    double i = 0.;
    for (; iter.isValid(); ++iter, ++i) {
      iter->setF({i, i, i});
    }
  }

  {
    ParticleIterator<MoleculeLJ, RMMParticleCell<MoleculeLJ>> iter(&data);
    double i = 0.;
    for (; iter.isValid(); ++iter, ++i) {
      auto force = iter->getF();
      ASSERT_EQ(force[0], i);
      ASSERT_EQ(force[1], i);
      ASSERT_EQ(force[2], i);
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
void testContainerIteratorBehavior(Container& container, Molecule& mol, Molecule& haloMol) {
  // default
  int count = 0;
  for (auto iter = container.begin(); iter.isValid(); ++iter) {
    count++;
  }
  EXPECT_EQ(count, 2);

  // haloAndOwned (same as default)
  count = 0;
  for (auto iter = container.begin(IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
    count++;
  }
  EXPECT_EQ(count, 2);

  // owned only
  count = 0;
  for (auto iter = container.begin(IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    count++;
    EXPECT_EQ(iter->getID(), mol.getID());
  }
  EXPECT_EQ(count, 1);

  // halo only
  count = 0;
  for (auto iter = container.begin(IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
    count++;
    EXPECT_EQ(iter->getID(), haloMol.getID());
  }
  EXPECT_EQ(count, 1);
}

TEST_F(ParticleIteratorTest, testIteratorBehaviorDirectSum) {
  DirectSum<MoleculeLJ, FullParticleCell<MoleculeLJ>> ds({0., 0., 0.}, {10., 10., 10.}, 3);
  MoleculeLJ mol({1., 1., 1.}, {0., 0., 0.}, 1);
  ds.addParticle(mol);
  MoleculeLJ haloMol({-1., 1., 1.}, {0., 0., 0.}, 2);
  ds.addHaloParticle(haloMol);

  testContainerIteratorBehavior(ds, mol, haloMol);
}

TEST_F(ParticleIteratorTest, testIteratorBehaviorLinkedCells) {
  LinkedCells<MoleculeLJ, FullParticleCell<MoleculeLJ>> linkedCells({0., 0., 0.}, {10., 10., 10.}, 3);
  MoleculeLJ mol({1., 1., 1.}, {0., 0., 0.}, 1);
  linkedCells.addParticle(mol);
  MoleculeLJ haloMol({-1., 1., 1.}, {0., 0., 0.}, 2);
  linkedCells.addHaloParticle(haloMol);

  testContainerIteratorBehavior(linkedCells, mol, haloMol);
}

TEST_F(ParticleIteratorTest, testIteratorBehaviorVerletLists) {
  VerletLists<MoleculeLJ> verletLists({0., 0., 0.}, {10., 10., 10.}, 3, 0., 1);
  MoleculeLJ mol({1., 1., 1.}, {0., 0., 0.}, 1);
  verletLists.addParticle(mol);
  MoleculeLJ haloMol({-1., 1., 1.}, {0., 0., 0.}, 2);
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