/**
 * @file KokkosSoATest.cpp
 * @author Luis Gall
 * @date 28.10.2026
 *
 */

#ifdef AUTOPAS_ENABLE_KOKKOS

#include "KokkosSoATest.h"

TEST_F(KokkosSoATest, testInitialization) { autopas::utilsKokkos::KokkosSoA<autopas::ParticleBaseFP32::SoAArraysType> soa; }

TEST_F(KokkosSoATest, SoATypeTest) {
  static_assert(std::is_same<autopas::utilsKokkos::KokkosSoA<unsigned int*, float*, float*, float*>::Type,
                             std::tuple<Kokkos::DualView<unsigned int*, autopas::utilsKokkos::DeviceMemSpace::device_type>,
                                        Kokkos::DualView<float*, autopas::utilsKokkos::DeviceMemSpace::device_type>,
                                        Kokkos::DualView<float*, autopas::utilsKokkos::DeviceMemSpace::device_type>,
                                        Kokkos::DualView<float*, autopas::utilsKokkos::DeviceMemSpace::device_type>>>::value,
                "type must be same");
}

using SoAType = autopas::utilsKokkos::KokkosSoA<unsigned int*, float*, float*, float*>::Type;

TEST_F(KokkosSoATest, SoAStorageTestGet) {

  Kokkos::initialize();
  {
    SoAType soa {};

    EXPECT_EQ(std::get<0>(soa).extent(0), 0);
    EXPECT_EQ(std::get<1>(soa).extent(0), 0);
    EXPECT_EQ(std::get<2>(soa).extent(0), 0);
    EXPECT_EQ(std::get<3>(soa).extent(0), 0);

    std::get<0>(soa).resize(4);
    std::get<1>(soa).resize(5);
    std::get<2>(soa).resize(6);
    std::get<3>(soa).resize(7);

    EXPECT_EQ(std::get<0>(soa).extent(0), 4);
    EXPECT_EQ(std::get<1>(soa).extent(0), 5);
    EXPECT_EQ(std::get<2>(soa).extent(0), 6);
    EXPECT_EQ(std::get<3>(soa).extent(0), 7);

    // TODO: check for correct type similar to SoATest.cpp SoAStorageTestGet
  }
  Kokkos::finalize();
}

TEST_F(KokkosSoATest, SoATestAddParticle) {
  Kokkos::initialize();
  {
    autopas::KokkosParticleBaseFP32::KokkosSoAArraysType soa;
    soa.resize(1);

    autopas::KokkosParticleBaseFP32 particle {};

    soa.template addParticle<autopas::KokkosParticleBaseFP32, false>(0, particle);l

  }
  Kokkos::finalize();
}

#endif
