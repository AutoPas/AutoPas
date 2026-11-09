/**
 * @file KokkosTestUtils.h
 * @author Luis Gall
 * @date 04.11.2026
 *
 * Helpers shared by the tests of KokkosSoA, KokkosAoS, and KokkosList.
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <string>
#include <type_traits>

#include "autopas/particles/ParticleDefinitions.h"
#include "autopas/utilsKokkos/KokkosAoS.h"
#include "autopas/utilsKokkos/KokkosSoA.h"
#include "autopas/utilsKokkos/KokkosSpaces.h"

namespace kokkosTestUtils {

/**
 * Particle type used in all Kokkos storage tests.
 */
using Particle = autopas::KokkosParticleBaseFP32;
/**
 * Attribute names of the particle type.
 */
using Attr = Particle::AttributeNames;
/**
 * SoA storage of the particle type.
 */
using SoA = Particle::KokkosSoAArraysType;
/**
 * AoS storage of the particle type.
 */
using AoS = autopas::utilsKokkos::KokkosAoS<Particle>;
using autopas::utilsKokkos::DeviceExecSpace;
using autopas::utilsKokkos::HostExecSpace;

/**
 * True if host and device share one memory space (e.g. OpenMP-only builds).
 * DualViews are then single-device: modify() and sync() are no-ops and need_sync_host() / need_sync_device() are
 * always false. Tests that rely on stale data or sync flags have to be skipped in this case.
 */
inline constexpr bool hostAndDeviceShareMemory =
    std::is_same_v<autopas::utilsKokkos::DeviceMemSpace, autopas::utilsKokkos::HostMemSpace>;

/**
 * Attribute values a particle is expected to hold.
 * This only consists of scalars, so it can be used inside of device kernels.
 */
struct ExpectedParticle {
  /**
   * Particle id.
   */
  unsigned int id;
  /**
   * Position, velocity, and force components.
   */
  float x, y, z, vx, vy, vz, fx, fy, fz;
  /**
   * Ownership state.
   */
  autopas::OwnershipState state;
};

/**
 * Generates values that are distinct for every index, so that lost, shifted, or swapped particles are detected.
 * All values are exactly representable as float, hence they can be compared exactly on host and device.
 * @param i Index of the particle.
 * @param offset Added to all floating point attributes, used to distinguish data that was overwritten on the device.
 * @return The values of the particle.
 */
KOKKOS_INLINE_FUNCTION ExpectedParticle expectedAt(size_t i, float offset = 0.f) {
  const float base = static_cast<float>(i) * 100.f + offset;
  return {static_cast<unsigned int>(42 + i),
          base + 1.f,
          base + 2.f,
          base + 3.f,
          base + 4.f,
          base + 5.f,
          base + 6.f,
          base + 7.f,
          base + 8.f,
          base + 9.f,
          i % 2 == 0 ? autopas::OwnershipState::halo : autopas::OwnershipState::owned};
}

/**
 * Creates a particle holding the values of expectedAt().
 * @param i Index of the particle.
 * @return The particle.
 */
inline Particle makeParticle(size_t i) {
  const auto expected = expectedAt(i);
  Particle particle{};
  particle.setID(expected.id);
  particle.setR({expected.x, expected.y, expected.z});
  particle.setV({expected.vx, expected.vy, expected.vz});
  particle.setF({expected.fx, expected.fy, expected.fz});
  particle.setOwnershipState(expected.state);
  return particle;
}

/**
 * Reads all attributes of the particles [first, last) through the storage's operator() in a device kernel.
 *
 * This is a free function because the kernel must not be defined in TEST_F's private TestBody(), which nvcc rejects.
 *
 * @tparam Storage SoA, AoS, or KokkosList.
 * @param storage
 * @param first
 * @param last
 * @param offset Offset that was used to generate the expected values.
 * @return Number of attributes that do not have the expected value.
 */
template <typename Storage>
int countDeviceMismatches(const Storage &storage, size_t first, size_t last, float offset = 0.f) {
  int errors = 0;
  Kokkos::parallel_reduce(
      "kokkosTestUtils::countDeviceMismatches", Kokkos::RangePolicy<DeviceExecSpace>(first, last),
      KOKKOS_LAMBDA(const size_t i, int &err) {
        const ExpectedParticle expected = expectedAt(i, offset);
        if (storage.template operator()<Attr::id, true>(i) != expected.id) ++err;
        if (storage.template operator()<Attr::posX, true>(i) != expected.x) ++err;
        if (storage.template operator()<Attr::posY, true>(i) != expected.y) ++err;
        if (storage.template operator()<Attr::posZ, true>(i) != expected.z) ++err;
        if (storage.template operator()<Attr::velX, true>(i) != expected.vx) ++err;
        if (storage.template operator()<Attr::velY, true>(i) != expected.vy) ++err;
        if (storage.template operator()<Attr::velZ, true>(i) != expected.vz) ++err;
        if (storage.template operator()<Attr::forceX, true>(i) != expected.fx) ++err;
        if (storage.template operator()<Attr::forceY, true>(i) != expected.fy) ++err;
        if (storage.template operator()<Attr::forceZ, true>(i) != expected.fz) ++err;
        if (storage.template operator()<Attr::ownershipState, true>(i) != expected.state) ++err;
      },
      errors);
  return errors;
}

/**
 * Overwrites all floating point attributes of the particles [0, numParticles) in a device kernel.
 * The id and the ownership state remain untouched.
 * @tparam Storage SoA, AoS, or KokkosList.
 * @param storage
 * @param numParticles
 * @param offset Offset used to generate the new values.
 */
template <typename Storage>
void overwriteOnDevice(const Storage &storage, size_t numParticles, float offset) {
  Kokkos::parallel_for(
      "kokkosTestUtils::overwriteOnDevice", Kokkos::RangePolicy<DeviceExecSpace>(0, numParticles),
      KOKKOS_LAMBDA(const size_t i) {
        const ExpectedParticle expected = expectedAt(i, offset);
        storage.template operator()<Attr::posX, true>(i) = expected.x;
        storage.template operator()<Attr::posY, true>(i) = expected.y;
        storage.template operator()<Attr::posZ, true>(i) = expected.z;
        storage.template operator()<Attr::velX, true>(i) = expected.vx;
        storage.template operator()<Attr::velY, true>(i) = expected.vy;
        storage.template operator()<Attr::velZ, true>(i) = expected.vz;
        storage.template operator()<Attr::forceX, true>(i) = expected.fx;
        storage.template operator()<Attr::forceY, true>(i) = expected.fy;
        storage.template operator()<Attr::forceZ, true>(i) = expected.fz;
      });
}

/**
 * Writes to the posX and posY attributes of the particles [0, numParticles) in a device kernel, without touching the
 * sync state.
 * @tparam Storage SoA, AoS, or KokkosList.
 * @param storage
 * @param numParticles
 * @param value
 */
template <typename Storage>
void writePositionsXYOnDevice(const Storage &storage, size_t numParticles, float value) {
  Kokkos::parallel_for(
      "kokkosTestUtils::writePositionsXYOnDevice", Kokkos::RangePolicy<DeviceExecSpace>(0, numParticles),
      KOKKOS_LAMBDA(const size_t i) {
        storage.template operator()<Attr::posX, true>(i) = value;
        storage.template operator()<Attr::posY, true>(i) = value;
      });
}

/**
 * Checks all attributes of one particle through the storage's operator() on the host.
 * @tparam Storage SoA, AoS, or KokkosList.
 * @param storage
 * @param i Index of the particle.
 * @param offset Offset that was used to generate the expected values.
 */
template <typename Storage>
void expectOnHost(const Storage &storage, size_t i, float offset = 0.f) {
  SCOPED_TRACE("Particle " + std::to_string(i));
  const auto expected = expectedAt(i, offset);
  EXPECT_EQ((storage.template operator()<Attr::id, false>(i)), expected.id);
  EXPECT_EQ((storage.template operator()<Attr::posX, false>(i)), expected.x);
  EXPECT_EQ((storage.template operator()<Attr::posY, false>(i)), expected.y);
  EXPECT_EQ((storage.template operator()<Attr::posZ, false>(i)), expected.z);
  EXPECT_EQ((storage.template operator()<Attr::velX, false>(i)), expected.vx);
  EXPECT_EQ((storage.template operator()<Attr::velY, false>(i)), expected.vy);
  EXPECT_EQ((storage.template operator()<Attr::velZ, false>(i)), expected.vz);
  EXPECT_EQ((storage.template operator()<Attr::forceX, false>(i)), expected.fx);
  EXPECT_EQ((storage.template operator()<Attr::forceY, false>(i)), expected.fy);
  EXPECT_EQ((storage.template operator()<Attr::forceZ, false>(i)), expected.fz);
  EXPECT_EQ((storage.template operator()<Attr::ownershipState, false>(i)), expected.state);
}

/**
 * Copies one SoA column from the device to the host.
 * @tparam Attribute Attribute name of the column.
 * @param soa
 * @return Host copy of the device view.
 */
template <Attr Attribute>
auto copyDeviceColumnToHost(const SoA &soa) {
  return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, soa.template getView<Attribute>().view_device());
}

/**
 * Checks the raw device storage of the particles [0, numParticles), independent of any operator().
 * @param soa
 * @param numParticles
 */
inline void expectDeviceStorage(const SoA &soa, size_t numParticles, float offset = 0.f) {
  const auto ids = copyDeviceColumnToHost<Attr::id>(soa);
  const auto xs = copyDeviceColumnToHost<Attr::posX>(soa);
  const auto ys = copyDeviceColumnToHost<Attr::posY>(soa);
  const auto zs = copyDeviceColumnToHost<Attr::posZ>(soa);
  const auto vxs = copyDeviceColumnToHost<Attr::velX>(soa);
  const auto vys = copyDeviceColumnToHost<Attr::velY>(soa);
  const auto vzs = copyDeviceColumnToHost<Attr::velZ>(soa);
  const auto fxs = copyDeviceColumnToHost<Attr::forceX>(soa);
  const auto fys = copyDeviceColumnToHost<Attr::forceY>(soa);
  const auto fzs = copyDeviceColumnToHost<Attr::forceZ>(soa);
  const auto states = copyDeviceColumnToHost<Attr::ownershipState>(soa);

  for (size_t i = 0; i < numParticles; ++i) {
    SCOPED_TRACE("Particle " + std::to_string(i));
    const auto expected = expectedAt(i, offset);
    EXPECT_EQ(ids(i), expected.id);
    EXPECT_EQ(xs(i), expected.x);
    EXPECT_EQ(ys(i), expected.y);
    EXPECT_EQ(zs(i), expected.z);
    EXPECT_EQ(vxs(i), expected.vx);
    EXPECT_EQ(vys(i), expected.vy);
    EXPECT_EQ(vzs(i), expected.vz);
    EXPECT_EQ(fxs(i), expected.fx);
    EXPECT_EQ(fys(i), expected.fy);
    EXPECT_EQ(fzs(i), expected.fz);
    EXPECT_EQ(states(i), expected.state);
  }
}

/**
 * Checks the raw device storage of the particles [0, numParticles), independent of any operator().
 * @param aos
 * @param numParticles
 */
inline void expectDeviceStorage(const AoS &aos, size_t numParticles, float offset = 0.f) {
  const auto particles = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, aos.getView().view_device());

  for (size_t i = 0; i < numParticles; ++i) {
    SCOPED_TRACE("Particle " + std::to_string(i));
    const auto expected = expectedAt(i, offset);
    EXPECT_EQ(particles(i).get<Attr::id>(), expected.id);
    EXPECT_EQ(particles(i).get<Attr::posX>(), expected.x);
    EXPECT_EQ(particles(i).get<Attr::posY>(), expected.y);
    EXPECT_EQ(particles(i).get<Attr::posZ>(), expected.z);
    EXPECT_EQ(particles(i).get<Attr::velX>(), expected.vx);
    EXPECT_EQ(particles(i).get<Attr::velY>(), expected.vy);
    EXPECT_EQ(particles(i).get<Attr::velZ>(), expected.vz);
    EXPECT_EQ(particles(i).get<Attr::forceX>(), expected.fx);
    EXPECT_EQ(particles(i).get<Attr::forceY>(), expected.fy);
    EXPECT_EQ(particles(i).get<Attr::forceZ>(), expected.fz);
    EXPECT_EQ(particles(i).get<Attr::ownershipState>(), expected.state);
  }
}

}  // namespace kokkosTestUtils

#endif
