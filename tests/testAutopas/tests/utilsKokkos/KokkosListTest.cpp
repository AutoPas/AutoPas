/**
 * @file KokkosListTest.cpp
 * @author Luis Gall
 * @date 30.10.2026
 *
 * Currently covered methods of KokkosList: operator<useDeviceView=false>()(index), addParticle(Particle), resize(size)
 */

#ifdef AUTOPAS_ENABLE_KOKKOS

#include "KokkosListTest.h"

TEST_F(KokkosListTest, testInitialization) {
  autopas::utilsKokkos::KokkosList<autopas::KokkosParticleBaseFP32> kokkosList{};
}

TEST_F(KokkosListTest, testConstructor) {
  autopas::utilsKokkos::KokkosList<autopas::KokkosParticleBaseFP32> kokkosList{autopas::DataLayoutOption::soa, 42};

  auto activeLayout = kokkosList.getActiveDataLayout();
  auto size = kokkosList.size();
  auto capacity = kokkosList.getCapacity();

  EXPECT_EQ(activeLayout, autopas::DataLayoutOption::soa);
  EXPECT_EQ(size, 0);
  EXPECT_EQ(capacity, 42);

  auto soa = kokkosList.getSoA();
  auto aos = kokkosList.getAoS();

  EXPECT_EQ(soa.getView<0>().extent(0), 42);
  EXPECT_EQ(aos.getView().extent(0), 42);
}

TEST_F(KokkosListTest, testResizeEmpty) {
  autopas::utilsKokkos::KokkosList<autopas::KokkosParticleBaseFP32> kokkosList{autopas::DataLayoutOption::soa, 42};

  // Case 1: downsize the list
  kokkosList.resize(16);

  auto sizeA = kokkosList.size();
  auto capacityA = kokkosList.getCapacity();

  EXPECT_EQ(sizeA, 0);
  EXPECT_EQ(capacityA, 16);
  EXPECT_EQ(kokkosList.getAoS().getView().extent(0), 16);
  EXPECT_EQ(kokkosList.getSoA().getView<0>().extent(0), 16);

  // todo: expect that warning is issued

  // Case 2: increase size of the list
  kokkosList.resize(32);

  auto sizeB = kokkosList.size();
  auto capacityB = kokkosList.getCapacity();

  EXPECT_EQ(sizeB, 0);
  EXPECT_EQ(capacityB, 32);
  EXPECT_EQ(kokkosList.getAoS().getView().extent(0), 32);
  EXPECT_EQ(kokkosList.getSoA().getView<0>().extent(0), 32);
}

TEST_F(KokkosListTest, testResizeFull) {
  autopas::utilsKokkos::KokkosList<autopas::KokkosParticleBaseFP32> kokkosList{autopas::DataLayoutOption::soa, 8};

  // Fill list with dummy particles
  autopas::KokkosParticleBaseFP32 particle{};
  for (int i = 0; i < 8; ++i) {
    kokkosList.addParticle(particle);
  }

  size_t sizeA = kokkosList.size();
  size_t capacityA = kokkosList.getCapacity();
  EXPECT_EQ(sizeA, 8);
  EXPECT_EQ(capacityA, 8);
  EXPECT_EQ(kokkosList.getAoS().getView().extent(0), 8);
  EXPECT_EQ(kokkosList.getSoA().getView<0>().extent(0), 8);

  kokkosList.resize(4);

  size_t sizeB = kokkosList.size();
  size_t capacityB = kokkosList.getCapacity();
  EXPECT_EQ(sizeB, 4);
  EXPECT_EQ(capacityB, 4);
  EXPECT_EQ(kokkosList.getAoS().getView().extent(0), 4);
  EXPECT_EQ(kokkosList.getSoA().getView<0>().extent(0), 4);
}

TEST_F(KokkosListTest, testAddParticleSoA) {

  autopas::utilsKokkos::KokkosList<autopas::KokkosParticleBaseFP32> kokkosList{autopas::DataLayoutOption::soa, 0};

  autopas::KokkosParticleBaseFP32 particle{};

  signed int id = 42;
  const std::array<float, 3> position{1.f, 2.f, 3.f};
  const std::array<float, 3> force{0.5f, 0.4f, 0.3f};
  const std::array<float, 3> velocity{0.123f, 0.456f, 0.789f};
  autopas::OwnershipState ownershipState{autopas::OwnershipState::halo};

  particle.setID(id);
  particle.setR(position);
  particle.setF(force);
  particle.setV(velocity);
  particle.setOwnershipState(ownershipState);

  // Adding a particle into an empty KokkosStorage object should increase the size
  kokkosList.addParticle(particle);

  size_t sizeA = kokkosList.size();
  size_t capacityA = kokkosList.getCapacity();

  EXPECT_EQ(sizeA, 1);
  EXPECT_EQ(capacityA, 2);

  // Test if the stored particle is correctly added to the SoA views of the storage object
  auto soa = kokkosList.getSoA();

  auto storedIdSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::id-1, false>(0);
  auto storedXSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posX-1, false>(0);
  auto storedYSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posY-1, false>(0);
  auto storedZSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posZ-1, false>(0);
  auto storedVxSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velX-1, false>(0);
  auto storedVySoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velY-1, false>(0);
  auto storedVzSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velZ-1, false>(0);
  auto storedFxSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceX-1, false>(0);
  auto storedFySoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceY-1, false>(0);
  auto storedFzSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceZ-1, false>(0);
  auto storedStateSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::ownershipState-1, false>(0);

  EXPECT_EQ(storedIdSoA, id);
  EXPECT_EQ(storedXSoA, position[0]);
  EXPECT_EQ(storedYSoA, position[1]);
  EXPECT_EQ(storedZSoA, position[2]);
  EXPECT_EQ(storedVxSoA, velocity[0]);
  EXPECT_EQ(storedVySoA, velocity[1]);
  EXPECT_EQ(storedVzSoA, velocity[2]);
  EXPECT_EQ(storedFxSoA, force[0]);
  EXPECT_EQ(storedFySoA, force[1]);
  EXPECT_EQ(storedFzSoA, force[2]);
  EXPECT_EQ(storedStateSoA, ownershipState);

  // Test if the view's operator() correctly accesses the right particle data
  auto storedIdList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::id, false>(0);
  auto storedXList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posX, false>(0);
  auto storedYList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posY, false>(0);
  auto storedZList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posZ, false>(0);
  auto storedVxList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velX, false>(0);
  auto storedVyList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velY, false>(0);
  auto storedVzList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velZ, false>(0);
  auto storedFxList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceX, false>(0);
  auto storedFyList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceY, false>(0);
  auto storedFzList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceZ, false>(0);
  auto storedStateList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::ownershipState, false>(0);

  EXPECT_EQ(storedIdList, id);
  EXPECT_EQ(storedXList, position[0]);
  EXPECT_EQ(storedYList, position[1]);
  EXPECT_EQ(storedZList, position[2]);
  EXPECT_EQ(storedVxList, velocity[0]);
  EXPECT_EQ(storedVyList, velocity[1]);
  EXPECT_EQ(storedVzList, velocity[2]);
  EXPECT_EQ(storedFxList, force[0]);
  EXPECT_EQ(storedFyList, force[1]);
  EXPECT_EQ(storedFzList, force[2]);
  EXPECT_EQ(storedStateList, ownershipState);
}

TEST_F(KokkosListTest, testAddParticleAoS) {

  autopas::utilsKokkos::KokkosList<autopas::KokkosParticleBaseFP32> kokkosList{autopas::DataLayoutOption::aos, 0};

  autopas::KokkosParticleBaseFP32 particle{};

  signed int id = 42;
  const std::array<float, 3> position{1.f, 2.f, 3.f};
  const std::array<float, 3> force{0.5f, 0.4f, 0.3f};
  const std::array<float, 3> velocity{0.123f, 0.456f, 0.789f};
  autopas::OwnershipState ownershipState{autopas::OwnershipState::halo};

  particle.setID(id);
  particle.setR(position);
  particle.setF(force);
  particle.setV(velocity);
  particle.setOwnershipState(ownershipState);

  // Adding a particle into an empty KokkosStorage object should increase the size
  kokkosList.addParticle(particle);

  size_t sizeA = kokkosList.size();
  size_t capacityA = kokkosList.getCapacity();

  EXPECT_EQ(sizeA, 1);
  EXPECT_EQ(capacityA, 2);

  // Test if the stored particle is correctly added to the SoA views of the storage object
  auto aos = kokkosList.getAoS();

  auto storedIdAoS = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::id, false>(0);
  auto storedXAoS = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posX, false>(0);
  auto storedYAoS = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posY, false>(0);
  auto storedZAoS = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posZ, false>(0);
  auto storedVxAoS = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velX, false>(0);
  auto storedVyAoS = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velY, false>(0);
  auto storedVzAoS = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velZ, false>(0);
  auto storedFxAoS = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceX, false>(0);
  auto storedFyAoS = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceY, false>(0);
  auto storedFzAoS = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceZ, false>(0);
  auto storedStateAoS = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::ownershipState, false>(0);

  EXPECT_EQ(storedIdAoS, id);
  EXPECT_EQ(storedXAoS, position[0]);
  EXPECT_EQ(storedYAoS, position[1]);
  EXPECT_EQ(storedZAoS, position[2]);
  EXPECT_EQ(storedVxAoS, velocity[0]);
  EXPECT_EQ(storedVyAoS, velocity[1]);
  EXPECT_EQ(storedVzAoS, velocity[2]);
  EXPECT_EQ(storedFxAoS, force[0]);
  EXPECT_EQ(storedFyAoS, force[1]);
  EXPECT_EQ(storedFzAoS, force[2]);
  EXPECT_EQ(storedStateAoS, ownershipState);

  // Test if the view's operator() correctly accesses the right particle data
  auto storedIdList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::id, false>(0);
  auto storedXList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posX, false>(0);
  auto storedYList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posY, false>(0);
  auto storedZList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posZ, false>(0);
  auto storedVxList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velX, false>(0);
  auto storedVyList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velY, false>(0);
  auto storedVzList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velZ, false>(0);
  auto storedFxList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceX, false>(0);
  auto storedFyList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceY, false>(0);
  auto storedFzList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceZ, false>(0);
  auto storedStateList = kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::ownershipState, false>(0);

  EXPECT_EQ(storedIdList, id);
  EXPECT_EQ(storedXList, position[0]);
  EXPECT_EQ(storedYList, position[1]);
  EXPECT_EQ(storedZList, position[2]);
  EXPECT_EQ(storedVxList, velocity[0]);
  EXPECT_EQ(storedVyList, velocity[1]);
  EXPECT_EQ(storedVzList, velocity[2]);
  EXPECT_EQ(storedFxList, force[0]);
  EXPECT_EQ(storedFyList, force[1]);
  EXPECT_EQ(storedFzList, force[2]);
  EXPECT_EQ(storedStateList, ownershipState);
}

#endif