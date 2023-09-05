/**
 * @file ParticleCounterTest.cpp
 * @author D. Martin
 * @date 01.09.23
 */

#include "ParticleCounterTest.h"

TEST_F(ParticleCounterTest, testParticleCountersFullParticleCell) {
  autopas::FullParticleCell<Particle> cell({1., 1., 1.});

  // expect ownershipstate ownedOrHalo
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == (autopas::OwnershipState::owned | autopas::OwnershipState::halo));
// insert 16 owned particles
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < 16; i++) {
    Particle p({(double)i / 32.0, (double)i / 32.0, (double)i / 32.0}, {0., 0., 0.}, i);
    p.setOwnershipState(autopas::OwnershipState::owned);
    cell.addParticle(p);
  }

  // expect ownershipstate owned
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == autopas::OwnershipState::owned);

// insert 16 halo particles
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for
#endif
  for (int i = 16; i < 32; i++) {
    Particle p({(double)i / 32.0, (double)i / 32.0, (double)i / 32.0}, {0., 0., 0.}, i);
    p.setOwnershipState(autopas::OwnershipState::halo);
    cell.addParticle(p);
  }

  EXPECT_EQ(cell.getNumberOfOwnedParticles(), 16);
  EXPECT_EQ(cell.getNumberOfHaloParticles(), 16);

  cell.deleteDummyParticles();

  EXPECT_EQ(cell.getNumberOfOwnedParticles(), 16);
  EXPECT_EQ(cell.getNumberOfHaloParticles(), 16);

  // expect ownershipstate ownedOrHalo
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == (autopas::OwnershipState::owned | autopas::OwnershipState::halo));

// remove 16 owned particles
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < 16; i++) {
    cell.deleteByIndex(i);
  }

  // expect ownershipstate halo
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == autopas::OwnershipState::halo);

// remove 8 halo particles
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < 8; i++) {
    cell.deleteByIndex(i);
  }

  EXPECT_EQ(cell.getNumberOfOwnedParticles(), 0);
  EXPECT_EQ(cell.getNumberOfHaloParticles(), 8);

  cell.clear();

  EXPECT_EQ(cell.getNumberOfOwnedParticles(), 0);
  EXPECT_EQ(cell.getNumberOfHaloParticles(), 0);

  // expect ownershipstate ownedOrHalo
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == (autopas::OwnershipState::owned | autopas::OwnershipState::halo));
}

TEST_F(ParticleCounterTest, testParticleCountersSortedCellView) {
  autopas::FullParticleCell<Particle> cell({1., 1., 1.});

  // expect ownershipstate ownedOrHalo
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == (autopas::OwnershipState::owned | autopas::OwnershipState::halo));
// insert 16 owned particles
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < 16; i++) {
    Particle p({(double)i / 32.0, (double)i / 32.0, (double)i / 32.0}, {0., 0., 0.}, i);
    p.setOwnershipState(autopas::OwnershipState::owned);
    cell.addParticle(p);
  }

  // expect ownershipstate owned
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == autopas::OwnershipState::owned);

// insert 16 halo particles
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for
#endif
  for (int i = 16; i < 32; i++) {
    Particle p({(double)i / 32.0, (double)i / 32.0, (double)i / 32.0}, {0., 0., 0.}, i);
    p.setOwnershipState(autopas::OwnershipState::halo);
    cell.addParticle(p);
  }

  auto scv = autopas::SortedCellView<Particle, autopas::FullParticleCell<Particle>>(cell, {1., 0., 0.});

  EXPECT_EQ(scv.getNumberOfOwnedParticles(), 16);
  EXPECT_EQ(scv.getNumberOfHaloParticles(), 16);

  scv.deleteDummyParticles();

  EXPECT_EQ(scv.getNumberOfOwnedParticles(), 16);
  EXPECT_EQ(scv.getNumberOfHaloParticles(), 16);

  // expect ownershipstate ownedOrHalo
  EXPECT_TRUE(scv.getPossibleParticleOwnerships() == (autopas::OwnershipState::owned | autopas::OwnershipState::halo));

  // remove 16 owned particles
  for (int i = 0; i < 16; i++) {
    scv.deleteByIndex(i);
  }

  // expect ownershipstate halo
  EXPECT_TRUE(scv.getPossibleParticleOwnerships() == autopas::OwnershipState::halo);

  // remove 8 halo particles
  for (int i = 0; i < 8; i++) {
    scv.deleteByIndex(i);
  }

  EXPECT_EQ(scv.getNumberOfOwnedParticles(), 0);
  EXPECT_EQ(scv.getNumberOfHaloParticles(), 8);

  scv.clear();

  EXPECT_EQ(scv.getNumberOfOwnedParticles(), 0);
  EXPECT_EQ(scv.getNumberOfHaloParticles(), 0);

  // expect ownershipstate ownedOrHalo
  EXPECT_TRUE(scv.getPossibleParticleOwnerships() == (autopas::OwnershipState::owned | autopas::OwnershipState::halo));
}

TEST_F(ParticleCounterTest, testParticleCountersReferenceParticleCell) {
  autopas::ReferenceParticleCell<Particle> cell({1., 1., 1.});

  std::vector<Particle> particles;

  // create Particles
  for (int i = 0; i < 32; i++) {
    Particle p({(double)i / 32.0, (double)i / 32.0, (double)i / 32.0}, {0., 0., 0.}, i);
    p.setOwnershipState(i < 16 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
    particles.push_back(p);
  }

  // expect ownershipstate ownedOrHalo
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == (autopas::OwnershipState::owned | autopas::OwnershipState::halo));
// insert 16 owned particles
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < 16; i++) {
    cell.addParticleReference(&particles[i]);
  }

  // expect ownershipstate owned
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == autopas::OwnershipState::owned);

// insert 16 halo particles
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for
#endif
  for (int i = 16; i < 32; i++) {
    cell.addParticleReference(&particles[i]);
  }

  EXPECT_EQ(cell.getNumberOfOwnedParticles(), 16);
  EXPECT_EQ(cell.getNumberOfHaloParticles(), 16);

  cell.deleteDummyParticles();

  EXPECT_EQ(cell.getNumberOfOwnedParticles(), 16);
  EXPECT_EQ(cell.getNumberOfHaloParticles(), 16);

  // expect ownershipstate ownedOrHalo
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == (autopas::OwnershipState::owned | autopas::OwnershipState::halo));

  // remove 16 owned particles
  for (int i = 0; i < 16; i++) {
    cell.deleteByIndex(i);
  }

  // expect ownershipstate halo
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == autopas::OwnershipState::halo);

  // remove 8 halo particles
  for (int i = 0; i < 8; i++) {
    cell.deleteByIndex(i);
  }

  EXPECT_EQ(cell.getNumberOfOwnedParticles(), 0);
  EXPECT_EQ(cell.getNumberOfHaloParticles(), 8);

  cell.clear();

  EXPECT_EQ(cell.getNumberOfOwnedParticles(), 0);
  EXPECT_EQ(cell.getNumberOfHaloParticles(), 0);

  // expect ownershipstate ownedOrHalo
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == (autopas::OwnershipState::owned | autopas::OwnershipState::halo));
}

TEST_F(ParticleCounterTest, testParticleCountersClusterTower) {
  autopas::internal::ClusterTower<Particle> cell(1);

  // expect ownershipstate ownedOrHalo
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == (autopas::OwnershipState::owned | autopas::OwnershipState::halo));
  // insert 16 owned particles
  for (int i = 0; i < 16; i++) {
    Particle p({(double)i / 32.0, (double)i / 32.0, (double)i / 32.0}, {0., 0., 0.}, i);
    p.setOwnershipState(autopas::OwnershipState::owned);
    cell.addParticle(p);
  }

  // expect ownershipstate owned
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == autopas::OwnershipState::owned);

  // insert 16 halo particles
  for (int i = 16; i < 32; i++) {
    Particle p({(double)i / 32.0, (double)i / 32.0, (double)i / 32.0}, {0., 0., 0.}, i);
    p.setOwnershipState(autopas::OwnershipState::halo);
    cell.addParticle(p);
  }

  EXPECT_EQ(cell.getNumberOfOwnedParticles(), 16);
  EXPECT_EQ(cell.getNumberOfHaloParticles(), 16);

  cell.deleteDummyParticles();

  EXPECT_EQ(cell.getNumberOfOwnedParticles(), 16);
  EXPECT_EQ(cell.getNumberOfHaloParticles(), 16);

  // expect ownershipstate ownedOrHalo
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == (autopas::OwnershipState::owned | autopas::OwnershipState::halo));

  // remove 16 owned particles
  for (int i = 0; i < 16; i++) {
    cell.deleteByIndex(i);
  }

  // expect ownershipstate halo
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == autopas::OwnershipState::halo);

  // remove 8 halo particles
  for (int i = 0; i < 8; i++) {
    cell.deleteByIndex(i);
  }

  EXPECT_EQ(cell.getNumberOfOwnedParticles(), 0);
  EXPECT_EQ(cell.getNumberOfHaloParticles(), 8);

  cell.clear();

  EXPECT_EQ(cell.getNumberOfOwnedParticles(), 0);
  EXPECT_EQ(cell.getNumberOfHaloParticles(), 0);

  // expect ownershipstate ownedOrHalo
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == (autopas::OwnershipState::owned | autopas::OwnershipState::halo));
}