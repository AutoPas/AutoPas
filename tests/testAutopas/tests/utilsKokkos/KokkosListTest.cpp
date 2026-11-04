/**
 * @file KokkosListTest.cpp
 * @author Luis Gall
 * @date 30.10.2026
 *
 * Currently covered methods of KokkosList: operator<useDeviceView>()(index), addParticle(Particle), resize(size),
 * syncAll<Space>(), modifyAll<Space>()
 */

#ifdef AUTOPAS_ENABLE_KOKKOS

#include "KokkosListTest.h"

#include "KokkosTestUtils.h"

using namespace kokkosTestUtils;

namespace {

using List = autopas::utilsKokkos::KokkosList<Particle>;

/**
 * Adds the particles with the indices [first, last) to the list.
 * @param list
 * @param first
 * @param last
 */
void addParticles(List &list, size_t first, size_t last) {
  for (size_t i = first; i < last; ++i) {
    list.addParticle(makeParticle(i));
  }
}

/**
 * Checks the raw device storage of the active data layout for the particles [0, numParticles).
 * The inactive data layout is not checked since it holds uninitialized data.
 * @param list
 * @param numParticles
 */
void expectDeviceStorage(const List &list, size_t numParticles) {
  if (list.getActiveDataLayout() == autopas::DataLayoutOption::soa) {
    kokkosTestUtils::expectDeviceStorage(list.getSoA(), numParticles);
  } else {
    kokkosTestUtils::expectDeviceStorage(list.getAoS(), numParticles);
  }
}

/**
 * Checks that the device views of both data layouts have the given extent.
 * @param list
 * @param extent
 */
void expectDeviceExtents(const List &list, size_t extent) {
  EXPECT_EQ(list.getAoS().getView().view_device().extent(0), extent);
  EXPECT_EQ(list.getSoA().template getView<0>().view_device().extent(0), extent);
  EXPECT_EQ(list.getSoA().template getView<Attr::ownershipState - 1>().view_device().extent(0), extent);
}

/**
 * Adds a single particle to an empty list and checks it on the device.
 * @param layout
 */
void testAddParticleDevice(autopas::DataLayoutOption layout) {
  List list{layout, 0};

  addParticles(list, 0, 1);
  EXPECT_EQ(list.size(), 1u);
  EXPECT_EQ(list.getCapacity(), 2u);

  // addParticle() already handles the host side, so only the sync to the device is left to do.
  list.syncAll<DeviceExecSpace>();

  expectDeviceStorage(list, 1);
  EXPECT_EQ(countDeviceMismatches(list, 0, 1), 0);
}

/**
 * Adds more particles than fit into the initial capacity, which reallocates the views several times.
 * @param layout
 */
void testAddParticleGrowthDevice(autopas::DataLayoutOption layout) {
  List list{layout, 0};

  addParticles(list, 0, 5);
  EXPECT_EQ(list.size(), 5u);
  EXPECT_EQ(list.getCapacity(), 8u);

  list.syncAll<DeviceExecSpace>();

  expectDeviceExtents(list, 8);
  expectDeviceStorage(list, 5);
  EXPECT_EQ(countDeviceMismatches(list, 0, 5), 0);
}

/**
 * Shrinks a full list and checks that the remaining particles are still correct on the device.
 * DualView::resize() takes a different code path depending on which side holds the newest data, so both are tested.
 * @param layout
 * @param syncBeforeResize If true the device is up to date before the resize, otherwise the host holds the newest data.
 */
void testResizeFullDevice(autopas::DataLayoutOption layout, bool syncBeforeResize) {
  List list{layout, 8};
  addParticles(list, 0, 8);

  if (syncBeforeResize) {
    list.syncAll<DeviceExecSpace>();
  }

  list.resize(4);
  list.syncAll<DeviceExecSpace>();

  EXPECT_EQ(list.size(), 4u);
  EXPECT_EQ(list.getCapacity(), 4u);
  expectDeviceExtents(list, 4);
  expectDeviceStorage(list, 4);
  EXPECT_EQ(countDeviceMismatches(list, 0, 4), 0);
}

/**
 * Modifies particles on the device, then adds a particle on the host which forces the views to grow.
 * The modifications and the new particle have to be available on both sides afterwards.
 * @param layout
 */
void testDeviceRoundTrip(autopas::DataLayoutOption layout) {
  const float offset = 1000.f;

  List list{layout, 0};
  addParticles(list, 0, 2);
  EXPECT_EQ(list.getCapacity(), 2u);

  // The device has to be up to date before it may be marked as modified. Otherwise, DualView::modify() aborts.
  list.syncAll<DeviceExecSpace>();
  overwriteOnDevice(list, 2, offset);
  list.modifyAll<DeviceExecSpace>();

  // addParticle() syncs the modifications back to the host before writing to it.
  addParticles(list, 2, 3);
  EXPECT_EQ(list.size(), 3u);
  EXPECT_EQ(list.getCapacity(), 4u);

  expectOnHost(list, 0, offset);
  expectOnHost(list, 1, offset);
  expectOnHost(list, 2);

  list.syncAll<DeviceExecSpace>();

  EXPECT_EQ(countDeviceMismatches(list, 0, 2, offset), 0);
  EXPECT_EQ(countDeviceMismatches(list, 2, 3), 0);
}

}  // namespace

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

TEST_F(KokkosListTest, testAddParticleSoAHost) {

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

TEST_F(KokkosListTest, testAddParticleAoSHost) {

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

TEST_F(KokkosListTest, testConstructorDevice) {
  const List list{autopas::DataLayoutOption::soa, 42};

  expectDeviceExtents(list, 42);
}

TEST_F(KokkosListTest, testResizeEmptyDevice) {
  List list{autopas::DataLayoutOption::soa, 42};

  // Case 1: downsize the list
  list.resize(16);
  expectDeviceExtents(list, 16);

  // Case 2: increase size of the list
  list.resize(32);
  expectDeviceExtents(list, 32);
}

TEST_F(KokkosListTest, testResizeFullSoADevice) {
  testResizeFullDevice(autopas::DataLayoutOption::soa, true);
  testResizeFullDevice(autopas::DataLayoutOption::soa, false);
}

TEST_F(KokkosListTest, testResizeFullAoSDevice) {
  testResizeFullDevice(autopas::DataLayoutOption::aos, true);
  testResizeFullDevice(autopas::DataLayoutOption::aos, false);
}

TEST_F(KokkosListTest, testAddParticleSoADevice) { testAddParticleDevice(autopas::DataLayoutOption::soa); }

TEST_F(KokkosListTest, testAddParticleAoSDevice) { testAddParticleDevice(autopas::DataLayoutOption::aos); }

TEST_F(KokkosListTest, testAddParticleGrowthSoADevice) { testAddParticleGrowthDevice(autopas::DataLayoutOption::soa); }

TEST_F(KokkosListTest, testAddParticleGrowthAoSDevice) { testAddParticleGrowthDevice(autopas::DataLayoutOption::aos); }

TEST_F(KokkosListTest, testDeviceRoundTripSoA) { testDeviceRoundTrip(autopas::DataLayoutOption::soa); }

TEST_F(KokkosListTest, testDeviceRoundTripAoS) { testDeviceRoundTrip(autopas::DataLayoutOption::aos); }

#endif