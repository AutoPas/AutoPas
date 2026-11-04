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

TEST_F(KokkosAoSTest, testInitialization) { AoS aos; }

TEST_F(KokkosAoSTest, testAddParticleHost) {
  const AoS aos = makeFilledAoS(3);

  EXPECT_EQ(aos.getView().extent(0), 3u);
  for (size_t i = 0; i < 3; ++i) {
    expectOnHost(aos, i);
  }
}

TEST_F(KokkosAoSTest, testAddParticleDevice) {
  AoS aos = makeFilledAoS(3);

  aos.sync<DeviceExecSpace>();

  expectDeviceStorage(aos, 3);
  EXPECT_EQ(countDeviceMismatches(aos, 0, 3), 0);
}

TEST_F(KokkosAoSTest, testDeviceRoundTrip) {
  const float offset = 1000.f;
  AoS aos = makeFilledAoS(3);

  // The device has to be up to date before it may be marked as modified. Otherwise, DualView::modify() aborts.
  aos.sync<DeviceExecSpace>();
  overwriteOnDevice(aos, 3, offset);
  aos.modify<DeviceExecSpace>();

  aos.sync<HostExecSpace>();

  for (size_t i = 0; i < 3; ++i) {
    expectOnHost(aos, i, offset);
  }
}

// DualView::resize() takes a different code path depending on which side holds the newest data.
TEST_F(KokkosAoSTest, testResizeHostNewest) {
  AoS aos = makeFilledAoS(2);

  aos.resize(8);
  EXPECT_EQ(aos.getView().view_host().extent(0), 8u);
  EXPECT_EQ(aos.getView().view_device().extent(0), 8u);

  aos.sync<DeviceExecSpace>();

  expectDeviceStorage(aos, 2);
  EXPECT_EQ(countDeviceMismatches(aos, 0, 2), 0);
}

TEST_F(KokkosAoSTest, testResizeDeviceNewest) {
  const float offset = 1000.f;
  AoS aos = makeFilledAoS(2);

  aos.sync<DeviceExecSpace>();
  overwriteOnDevice(aos, 2, offset);
  aos.modify<DeviceExecSpace>();

  aos.resize(8);
  EXPECT_EQ(aos.getView().view_host().extent(0), 8u);
  EXPECT_EQ(aos.getView().view_device().extent(0), 8u);

  aos.sync<HostExecSpace>();

  for (size_t i = 0; i < 2; ++i) {
    expectOnHost(aos, i, offset);
  }
}

TEST_F(KokkosAoSTest, testShrink) {
  AoS aos = makeFilledAoS(8);

  aos.resize(4);
  aos.sync<DeviceExecSpace>();

  EXPECT_EQ(aos.getView().view_device().extent(0), 4u);
  expectDeviceStorage(aos, 4);
  EXPECT_EQ(countDeviceMismatches(aos, 0, 4), 0);
}

#endif
