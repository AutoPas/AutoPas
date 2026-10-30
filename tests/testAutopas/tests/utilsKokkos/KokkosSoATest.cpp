/**
 * @file KokkosSoATest.cpp
 * @author Luis Gall
 * @date 28.10.2026
 *
 * Currently covered methods of KokkosSoA: operator<useDeviceView=false>()(index), addParticle(index, Particle)
 */

#ifdef AUTOPAS_ENABLE_KOKKOS

#include "KokkosSoATest.h"

TEST_F(KokkosSoATest, testInitialization) {
  autopas::utilsKokkos::KokkosSoA<autopas::KokkosParticleBaseFP32::KokkosSoAArraysType> soa;
}

TEST_F(KokkosSoATest, SoATypeTest) {
  static_assert(
      std::is_same<autopas::utilsKokkos::KokkosSoA<unsigned int *, float *, float *, float *>::Type,
                   std::tuple<Kokkos::DualView<unsigned int *, autopas::utilsKokkos::DeviceMemSpace::device_type>,
                              Kokkos::DualView<float *, autopas::utilsKokkos::DeviceMemSpace::device_type>,
                              Kokkos::DualView<float *, autopas::utilsKokkos::DeviceMemSpace::device_type>,
                              Kokkos::DualView<float *, autopas::utilsKokkos::DeviceMemSpace::device_type>>>::value,
      "type must be same");
}

using SoAType = autopas::utilsKokkos::KokkosSoA<unsigned int *, float *, float *, float *>::Type;

TEST_F(KokkosSoATest, SoAStorageTestGet) {
  SoAType soa{};

  EXPECT_EQ(std::get<0>(soa).extent(0), 0);
  EXPECT_EQ(std::get<1>(soa).extent(0), 0);
  EXPECT_EQ(std::get<2>(soa).extent(0), 0);
  EXPECT_EQ(std::get<3>(soa).extent(0), 0);

  size_t size1 = 4;
  size_t size2 = 5;
  size_t size3 = 6;
  size_t size4 = 7;

  std::get<0>(soa).resize(size1);
  std::get<1>(soa).resize(size2);
  std::get<2>(soa).resize(size3);
  std::get<3>(soa).resize(size4);

  EXPECT_EQ(std::get<0>(soa).extent(0), size1);
  EXPECT_EQ(std::get<1>(soa).extent(0), size2);
  EXPECT_EQ(std::get<2>(soa).extent(0), size3);
  EXPECT_EQ(std::get<3>(soa).extent(0), size4);

  static_assert(std::is_same<decltype(std::get<0>(soa).view_host().data()), unsigned int *>::value,
                "id type must be proper(unsigned int)");
  static_assert(std::is_same<decltype(std::get<1>(soa).view_host().data()), float *>::value,
                "position type must be proper(float)");
  static_assert(std::is_same<decltype(std::get<2>(soa).view_host().data()), float *>::value,
                "position type must be proper(float)");
  static_assert(std::is_same<decltype(std::get<3>(soa).view_host().data()), float *>::value,
                "position type must be proper(float)");
}

TEST_F(KokkosSoATest, SoATestGet) {
  autopas::utilsKokkos::KokkosSoA<unsigned int *, float *, double *> soa {};

  soa.resize(1);

  auto storedView1 = soa.getView<0>();
  auto storedView2 = soa.getView<1>();
  auto storedView3 = soa.getView<2>();

  EXPECT_EQ(storedView1.extent(0), 1);
  EXPECT_EQ(storedView2.extent(0), 1);
  EXPECT_EQ(storedView3.extent(0), 1);

  static_assert(std::is_same<decltype(storedView1)::value_type, unsigned int>::value,
                "first attribute type must be proper(unsigned int)");
  static_assert(std::is_same<decltype(storedView2)::value_type, float>::value,
                "second attribute type must be proper(float)");
  static_assert(std::is_same<decltype(storedView3)::value_type, double>::value,
                "third attribute type must be proper(double)");
}

TEST_F(KokkosSoATest, SoATestAddParticle) {
  autopas::KokkosParticleBaseFP32::KokkosSoAArraysType soa;
  soa.resize(1);

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

  soa.template addParticle<autopas::KokkosParticleBaseFP32>(0, particle);

  auto storedId = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::id - 1, false>(0);
  auto storedX = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posX - 1, false>(0);
  auto storedY = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posY - 1, false>(0);
  auto storedZ = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posZ - 1, false>(0);
  auto storedVx = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velX - 1, false>(0);
  auto storedVy = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velY - 1, false>(0);
  auto storedVz = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velZ - 1, false>(0);
  auto storedFx = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceX - 1, false>(0);
  auto storedFy = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceY - 1, false>(0);
  auto storedFz = soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceZ - 1, false>(0);
  auto storedState =
      soa.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::ownershipState - 1, false>(0);

  EXPECT_EQ(storedId, id);
  EXPECT_EQ(storedX, position[0]);
  EXPECT_EQ(storedY, position[1]);
  EXPECT_EQ(storedZ, position[2]);
  EXPECT_EQ(storedVx, velocity[0]);
  EXPECT_EQ(storedVy, velocity[1]);
  EXPECT_EQ(storedVz, velocity[2]);
  EXPECT_EQ(storedFx, force[0]);
  EXPECT_EQ(storedFy, force[1]);
  EXPECT_EQ(storedFz, force[2]);
  EXPECT_EQ(storedState, ownershipState);
}

#endif
