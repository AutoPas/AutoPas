/**
 * @file KokkosListTest.cpp
 * @author Luis Gall
 * @date 30.10.2026
 *
 * Currently covered methods of KokkosList: operator()(index), addParticle(index, Particle), resize(size)
 */

#ifdef AUTOPAS_ENABLE_KOKKOS

#include "KokkosListTest.h"

TEST_F(KokkosListTest, testInitialization) {
  autopas::utilsKokkos::KokkosList<autopas::KokkosParticleBaseFP32> kokkosList{};
}

TEST_F(KokkosListTest, testConstructor) {
  autopas::utilsKokkos::KokkosList<autopas::KokkosParticleBaseFP32> kokkosList{autopas::DataLayoutOption::soa, 42};
}

#endif