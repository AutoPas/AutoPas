/**
 * @file KokkosListTest.cpp
 * @author Luis Gall
 * @date 30.10.2026
 *
 * Currently covered methods of KokkosList: operator<useDeviceView>()(index), addParticle(Particle), resize(size),
 * syncAll<Space>(), modifyAll<Space>(), sync<Space, Attribute>(), modify<Space, Attribute>()
 */

#ifdef AUTOPAS_ENABLE_KOKKOS

#include "KokkosListTest.h"

#include "KokkosTestUtils.h"

using namespace kokkosTestUtils;

namespace {

using List = autopas::utils::kokkos::KokkosList<Particle>;

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
 * The inactive data layout is not checked since only the active one is allocated.
 * @param list
 * @param numParticles
 */
void expectDeviceStorage(const List &list, size_t numParticles, float offset = 0.f) {
  if (list.getActiveDataLayout() == autopas::DataLayoutOption::soa) {
    kokkosTestUtils::expectDeviceStorage(list.getSoA(), numParticles, offset);
  } else {
    kokkosTestUtils::expectDeviceStorage(list.getAoS(), numParticles, offset);
  }
}

/**
 * Checks that the device views of the active data layout have the given extent.
 * The inactive data layout is not checked since only the active one is allocated.
 * @param list
 * @param extent
 */
void expectDeviceExtents(const List &list, size_t extent) {
  if (list.getActiveDataLayout() == autopas::DataLayoutOption::soa) {
    EXPECT_EQ(list.getSoA().template getView<Attr::id>().view_device().extent(0), extent);
    EXPECT_EQ(list.getSoA().template getView<Attr::ownershipState>().view_device().extent(0), extent);
  } else {
    EXPECT_EQ(list.getAoS().getView().view_device().extent(0), extent);
  }
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
  expectDeviceStorage(list, 2, offset);
  EXPECT_EQ(countDeviceMismatches(list, 0, 2., offset), 0);

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

/**
 * Tests if the initialization of an empty KokkosList works with the default constructor
 */
TEST_F(KokkosListTest, testInitialization) { List kokkosList{}; }

/**
 * Tests if the initialization of a allocated KokkosList works with the custom constructor
 */
TEST_F(KokkosListTest, testConstructor) {
  List kokkosList{autopas::DataLayoutOption::soa, 42};

  auto activeLayout = kokkosList.getActiveDataLayout();
  auto size = kokkosList.size();
  auto capacity = kokkosList.getCapacity();

  EXPECT_EQ(activeLayout, autopas::DataLayoutOption::soa);
  EXPECT_EQ(size, 0);
  EXPECT_EQ(capacity, 42);
  expectDeviceExtents(kokkosList, 42);

  // Only the active data layout is allocated.
  auto soa = kokkosList.getSoA();
  auto aos = kokkosList.getAoS();
  EXPECT_EQ(soa.getView<Attr::id>().extent(0), 42);
  EXPECT_EQ(aos.getView().extent(0), 0);
}

/**
 * A default-constructed DataLayoutOption is invalid. The constructor's resize() has to reject it instead of silently
 * allocating nothing.
 */
TEST_F(KokkosListTest, testConstructorInvalidLayout) {
  // Other tests change the global behavior, so make sure that exceptions are thrown.
  autopas::utils::ExceptionHandler::setBehavior(autopas::utils::ExceptionBehavior::throwException);

  EXPECT_THROW((List{autopas::DataLayoutOption{}, 4}), autopas::utils::ExceptionHandler::AutoPasException);
}

/**
 * Tests multiple resize cases, e.g. increasing or decreasing size for an empty KokkosList
 */
TEST_F(KokkosListTest, testResizeEmpty) {
  List kokkosList{autopas::DataLayoutOption::soa, 42};

  // Case 1: downsize the list
  kokkosList.resize(16);  // todo: expect that warning is issued

  auto sizeA = kokkosList.size();
  auto capacityA = kokkosList.getCapacity();

  EXPECT_EQ(sizeA, 0);
  EXPECT_EQ(capacityA, 16);
  expectDeviceExtents(kokkosList, 16);

  // Case 2: increase size of the list
  kokkosList.resize(32);

  auto sizeB = kokkosList.size();
  auto capacityB = kokkosList.getCapacity();

  EXPECT_EQ(sizeB, 0);
  EXPECT_EQ(capacityB, 32);
  expectDeviceExtents(kokkosList, 32);
}

/**
 * Tests multiple resize cases, e.g. increasing or decreasing size for an full KokkosList and if the stored data
 * survives the resize
 */
TEST_F(KokkosListTest, testResizeFull) {
  List kokkosList{autopas::DataLayoutOption::soa, 8};

  addParticles(kokkosList, 0, 8);

  size_t sizeA = kokkosList.size();
  size_t capacityA = kokkosList.getCapacity();
  EXPECT_EQ(sizeA, 8);
  EXPECT_EQ(capacityA, 8);
  expectDeviceExtents(kokkosList, 8);

  kokkosList.resize(4);
  kokkosList.syncAll<DeviceExecSpace>();
  expectDeviceStorage(kokkosList, 4);
  EXPECT_EQ(countDeviceMismatches(kokkosList, 0, 4), 0);

  size_t sizeB = kokkosList.size();
  size_t capacityB = kokkosList.getCapacity();
  EXPECT_EQ(sizeB, 4);
  EXPECT_EQ(capacityB, 4);
  expectDeviceExtents(kokkosList, 4);
}

/**
 * Tests if a particle is correctly added to the KokkosList on the Host with SoA as the active storage format
 */
TEST_F(KokkosListTest, testAddParticleSoAHost) {
  List kokkosList{autopas::DataLayoutOption::soa, 0};

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

  auto storedIdSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::id, false>(0);
  auto storedXSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posX, false>(0);
  auto storedYSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posY, false>(0);
  auto storedZSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posZ, false>(0);
  auto storedVxSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velX, false>(0);
  auto storedVySoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velY, false>(0);
  auto storedVzSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velZ, false>(0);
  auto storedFxSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceX, false>(0);
  auto storedFySoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceY, false>(0);
  auto storedFzSoA = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceZ, false>(0);
  auto storedStateSoA =
      soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::ownershipState, false>(0);

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
  auto storedStateList =
      kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::ownershipState, false>(0);

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

/**
 * Tests if a particle is correctly added to the KokkosList on the Host with AoS as the active storage format
 */
TEST_F(KokkosListTest, testAddParticleAoSHost) {
  List kokkosList{autopas::DataLayoutOption::aos, 0};

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
  auto storedStateAoS =
      aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::ownershipState, false>(0);

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
  auto storedStateList =
      kokkosList.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::ownershipState, false>(0);

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

/**
 * Tests if the sync<Attribute>() function only syncs the respective attribute and no other Views
 */
TEST_F(KokkosListTest, testPerAttributeSyncSoA) {
  if constexpr (hostAndDeviceShareMemory) {
    GTEST_SKIP() << "Host and device share memory, so DualView sync flags are no-ops.";
  }

  const float newValue = -1.f;
  List list{autopas::DataLayoutOption::soa, 0};
  addParticles(list, 0, 2);

  list.syncAll<DeviceExecSpace>();

  // Both attributes change on the device, but only posX is marked as modified.
  writePositionsXYOnDevice(list, 2, newValue);
  list.modify<DeviceExecSpace, Attr::posX>();

  EXPECT_TRUE(list.getSoA().getView<Attr::posX>().need_sync_host());
  EXPECT_FALSE(list.getSoA().getView<Attr::posY>().need_sync_host());

  list.sync<HostExecSpace, Attr::posX>();

  EXPECT_FALSE(list.getSoA().getView<Attr::posX>().need_sync_host());

  for (size_t i = 0; i < 2; ++i) {
    SCOPED_TRACE("Particle " + std::to_string(i));
    EXPECT_EQ((list.operator()<Attr::posX, false>(i)), newValue);
    EXPECT_EQ((list.operator()<Attr::posY, false>(i)), expectedAt(i).y);
  }
}

/**
 * Adds particles to a default-constructed list. Its DualViews have no modified flags until the first resize, so the
 * first addParticle() has to mark the host as modified only after growing the views. Otherwise, the device is marked as
 * modified, and the next addParticle() overwrites the first particle on the host with uninitialized device memory.
 */
TEST_F(KokkosListTest, testAddParticleDefaultConstructedDevice) {
  List list{};

  addParticles(list, 0, 3);
  EXPECT_EQ(list.size(), 3u);
  EXPECT_EQ(list.getCapacity(), 4u);

  for (size_t i = 0; i < 3; ++i) {
    expectOnHost(list, i);
  }

  list.syncAll<DeviceExecSpace>();

  expectDeviceStorage(list, 3);
  EXPECT_EQ(countDeviceMismatches(list, 0, 3), 0);
}

/* The following tests are described in the respective test functions in separate comments */

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