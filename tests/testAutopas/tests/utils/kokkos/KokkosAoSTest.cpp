/**
 * @file KokkosAoSTest.cpp
 * @author Luis Gall
 * @date 28.10.2026
 *
 * Currently covered methods of KokkosAoS: operator<useDeviceView>()(index), addParticle(index, Particle), resize(size),
 * sync<Space>(), modify<Space>()
 */

#ifdef AUTOPAS_ENABLE_KOKKOS

#include "KokkosAoSTest.h"

#include "KokkosTestUtils.h"

using namespace kokkosTestUtils;

namespace {

/**
 * Creates an AoS with the given size and fills it with the particles [0, size) on the host.
 * In contrast to KokkosList, KokkosAoS::addParticle() does not handle the sync state on its own.
 * @param size
 * @return The filled AoS, in which the host holds the newest data.
 */
AoS makeFilledAoS(size_t size) {
  AoS aos{};
  aos.resize(size);

  // resize() marks the device as modified. Without this sync, the modify() below aborts.
  aos.sync<HostExecSpace>();
  for (size_t i = 0; i < size; ++i) {
    aos.addParticle(i, makeParticle(i));
  }
  aos.modify<HostExecSpace>();

  return aos;
}

}  // namespace

/**
 * Tests if the initialization of an empty KokkosAoS works with the default constructor
 */
TEST_F(KokkosAoSTest, testInitialization) { AoS aos; }

/**
 * Tests if three particles are correctly added to the AoS on the Host, synced to the Device, and remain as expected
 */
TEST_F(KokkosAoSTest, testAddParticleHostAndDevice) {
  AoS aos = makeFilledAoS(3);

  EXPECT_EQ(aos.getView().extent(0), 3u);
  for (size_t i = 0; i < 3; ++i) {
    expectOnHost(aos, i);
  }

  aos.sync<DeviceExecSpace>();

  expectDeviceStorage(aos, 3);
  EXPECT_EQ(countDeviceMismatches(aos, 0, 3), 0);
}

/**
 * Tests if three particles that are modified on the Device are as expected in the Device storage, and are correctly
 * synced back to the Host
 */
TEST_F(KokkosAoSTest, testDeviceRoundTrip) {
  const float offset = 1000.f;
  AoS aos = makeFilledAoS(3);

  // The device has to be up to date before it may be marked as modified. Otherwise, DualView::modify() aborts.
  aos.sync<DeviceExecSpace>();
  overwriteOnDevice(aos, 3, offset);
  aos.modify<DeviceExecSpace>();

  expectDeviceStorage(aos, 3, offset);
  EXPECT_EQ(countDeviceMismatches(aos, 0, 3, offset), 0);

  aos.sync<HostExecSpace>();

  for (size_t i = 0; i < 3; ++i) {
    expectOnHost(aos, i, offset);
  }
}

/**
 * Tests if the Host particles survive a resize in case the Host View holds the most up-to-date data
 */
TEST_F(KokkosAoSTest, testResizeHostNewest) {
  AoS aos = makeFilledAoS(2);

  aos.resize(8);
  EXPECT_EQ(aos.getView().view_host().extent(0), 8u);
  EXPECT_EQ(aos.getView().view_device().extent(0), 8u);

  for (size_t i = 0; i < 2; ++i) {
    expectOnHost(aos, i);
  }

  aos.sync<DeviceExecSpace>();

  expectDeviceStorage(aos, 2);
  EXPECT_EQ(countDeviceMismatches(aos, 0, 2), 0);
}

/**
 * Tests if the Device particles survive a resize in case the Device View holds the most up-to-date data
 */
TEST_F(KokkosAoSTest, testResizeDeviceNewest) {
  const float offset = 1000.f;
  AoS aos = makeFilledAoS(2);

  aos.sync<DeviceExecSpace>();
  overwriteOnDevice(aos, 2, offset);
  aos.modify<DeviceExecSpace>();

  aos.resize(8);
  EXPECT_EQ(aos.getView().view_host().extent(0), 8u);
  EXPECT_EQ(aos.getView().view_device().extent(0), 8u);

  expectDeviceStorage(aos, 2, offset);
  EXPECT_EQ(countDeviceMismatches(aos, 0, 2, offset), 0);

  aos.sync<HostExecSpace>();

  for (size_t i = 0; i < 2; ++i) {
    expectOnHost(aos, i, offset);
  }
}

/**
 * Tests if after a resize that downsizes a full list, the particles are still existing on Host and Device
 */
TEST_F(KokkosAoSTest, testShrink) {
  AoS aos = makeFilledAoS(8);

  aos.resize(4);
  for (size_t i = 0; i < 4; ++i) {
    expectOnHost(aos, i);
  }
  aos.sync<DeviceExecSpace>();

  EXPECT_EQ(aos.getView().view_device().extent(0), 4u);
  expectDeviceStorage(aos, 4);
  EXPECT_EQ(countDeviceMismatches(aos, 0, 4), 0);
}

#endif
