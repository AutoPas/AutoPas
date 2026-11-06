/**
 * @file KokkosSoATest.cpp
 * @author Luis Gall
 * @date 28.10.2026
 *
 * Currently covered methods of KokkosSoA: operator<useDeviceView>()(index), addParticle(index, Particle), resize(size),
 * getView<Attribute>(), sync<Space, Attribute>(), syncAll<Space>(), modify<Space, Attribute>(), modifyAll<Space>()
 */

#ifdef AUTOPAS_ENABLE_KOKKOS

#include "KokkosSoATest.h"

#include "KokkosTestUtils.h"

using namespace kokkosTestUtils;

namespace {

/**
 * Column names for the generic SoAs below. Nothing is skipped, so the first name is stored in the first column.
 */
enum class TestColumns { first, second, third, fourth };

using TestSoA =
    autopas::utilsKokkos::KokkosSoA<TestColumns, TestColumns::first, unsigned int *, float *, float *, float *>;

static_assert(
    std::is_same_v<TestSoA::Type,
                   std::tuple<Kokkos::DualView<unsigned int *, autopas::utilsKokkos::DeviceMemSpace::device_type>,
                              Kokkos::DualView<float *, autopas::utilsKokkos::DeviceMemSpace::device_type>,
                              Kokkos::DualView<float *, autopas::utilsKokkos::DeviceMemSpace::device_type>,
                              Kokkos::DualView<float *, autopas::utilsKokkos::DeviceMemSpace::device_type>>>,
    "KokkosSoA must store one DualView per attribute.");

static_assert(std::is_same_v<decltype(std::get<0>(std::declval<TestSoA::Type &>()).view_host().data()), unsigned int *>,
              "The host view of the first column must store unsigned int.");
static_assert(std::is_same_v<decltype(std::get<1>(std::declval<TestSoA::Type &>()).view_host().data()), float *>,
              "The host view of the second column must store float.");

/**
 * Creates an SoA with the given size and fills it with the particles [0, size) on the host.
 * In contrast to KokkosList, KokkosSoA::addParticle() does not handle the sync state on its own.
 * @param size
 * @return The filled SoA, in which the host holds the newest data.
 */
SoA makeFilledSoA(size_t size) {
  SoA soa{};
  soa.resize(size);

  // resize() marks the device as modified. Without this sync, the modifyAll() below aborts.
  soa.syncAll<HostExecSpace>();
  for (size_t i = 0; i < size; ++i) {
    soa.addParticle(i, makeParticle(i));
  }
  soa.modifyAll<HostExecSpace>();

  return soa;
}

}  // namespace

TEST_F(KokkosSoATest, testInitialization) { SoA soa; }

TEST_F(KokkosSoATest, testGetView) {
  autopas::utilsKokkos::KokkosSoA<TestColumns, TestColumns::first, unsigned int *, float *, double *> soa{};

  soa.resize(1);

  const auto &storedView1 = soa.getView<TestColumns::first>();
  const auto &storedView2 = soa.getView<TestColumns::second>();
  const auto &storedView3 = soa.getView<TestColumns::third>();

  EXPECT_EQ(storedView1.extent(0), 1u);
  EXPECT_EQ(storedView2.extent(0), 1u);
  EXPECT_EQ(storedView3.extent(0), 1u);

  static_assert(std::is_same_v<std::remove_cvref_t<decltype(storedView1)>::value_type, unsigned int>,
                "The first attribute must be unsigned int.");
  static_assert(std::is_same_v<std::remove_cvref_t<decltype(storedView2)>::value_type, float>,
                "The second attribute must be float.");
  static_assert(std::is_same_v<std::remove_cvref_t<decltype(storedView3)>::value_type, double>,
                "The third attribute must be double.");
}

TEST_F(KokkosSoATest, testAddParticleHost) {
  const SoA soa = makeFilledSoA(3);

  EXPECT_EQ(soa.getView<Attr::id>().extent(0), 3u);
  for (size_t i = 0; i < 3; ++i) {
    expectOnHost(soa, i);
  }
}

TEST_F(KokkosSoATest, testAddParticleDevice) {
  SoA soa = makeFilledSoA(3);

  soa.syncAll<DeviceExecSpace>();

  expectDeviceStorage(soa, 3);
  EXPECT_EQ(countDeviceMismatches(soa, 0, 3), 0);
}

TEST_F(KokkosSoATest, testDeviceRoundTrip) {
  const float offset = 1000.f;
  SoA soa = makeFilledSoA(3);

  // The device has to be up to date before it may be marked as modified. Otherwise, DualView::modify() aborts.
  soa.syncAll<DeviceExecSpace>();
  overwriteOnDevice(soa, 3, offset);
  soa.modifyAll<DeviceExecSpace>();

  soa.syncAll<HostExecSpace>();

  for (size_t i = 0; i < 3; ++i) {
    expectOnHost(soa, i, offset);
  }
}

// DualView::resize() takes a different code path depending on which side holds the newest data.
TEST_F(KokkosSoATest, testResizeHostNewest) {
  SoA soa = makeFilledSoA(2);

  soa.resize(8);
  EXPECT_EQ(soa.getView<Attr::id>().view_host().extent(0), 8u);
  EXPECT_EQ(soa.getView<Attr::id>().view_device().extent(0), 8u);

  soa.syncAll<DeviceExecSpace>();

  expectDeviceStorage(soa, 2);
  EXPECT_EQ(countDeviceMismatches(soa, 0, 2), 0);
}

TEST_F(KokkosSoATest, testResizeDeviceNewest) {
  const float offset = 1000.f;
  SoA soa = makeFilledSoA(2);

  soa.syncAll<DeviceExecSpace>();
  overwriteOnDevice(soa, 2, offset);
  soa.modifyAll<DeviceExecSpace>();

  soa.resize(8);
  EXPECT_EQ(soa.getView<Attr::id>().view_host().extent(0), 8u);
  EXPECT_EQ(soa.getView<Attr::id>().view_device().extent(0), 8u);

  soa.syncAll<HostExecSpace>();

  for (size_t i = 0; i < 2; ++i) {
    expectOnHost(soa, i, offset);
  }
}

TEST_F(KokkosSoATest, testShrink) {
  SoA soa = makeFilledSoA(8);

  soa.resize(4);
  soa.syncAll<DeviceExecSpace>();

  EXPECT_EQ(soa.getView<Attr::id>().view_device().extent(0), 4u);
  expectDeviceStorage(soa, 4);
  EXPECT_EQ(countDeviceMismatches(soa, 0, 4), 0);
}

// Columns have to sync independently of each other.
TEST_F(KokkosSoATest, testPerColumnSync) {
  if constexpr (hostAndDeviceShareMemory) {
    GTEST_SKIP() << "Host and device share memory, so DualView sync flags are no-ops.";
  }

  const float newValue = -1.f;
  SoA soa = makeFilledSoA(2);

  soa.syncAll<DeviceExecSpace>();

  // Both columns change on the device, but only posX is marked as modified.
  writePositionsXYOnDevice(soa, 2, newValue);
  soa.modify<DeviceExecSpace, Attr::posX>();

  EXPECT_TRUE(soa.getView<Attr::posX>().need_sync_host());
  EXPECT_FALSE(soa.getView<Attr::posY>().need_sync_host());

  soa.sync<HostExecSpace, Attr::posX>();

  EXPECT_FALSE(soa.getView<Attr::posX>().need_sync_host());

  for (size_t i = 0; i < 2; ++i) {
    SCOPED_TRACE("Particle " + std::to_string(i));
    EXPECT_EQ((soa.operator()<Attr::posX, false>(i)), newValue);
    EXPECT_EQ((soa.operator()<Attr::posY, false>(i)), expectedAt(i).y);
  }
}

#endif
