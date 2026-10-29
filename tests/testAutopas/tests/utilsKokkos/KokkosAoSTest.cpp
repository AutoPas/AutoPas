/**
* @file KokkosAoSTest.cpp
 * @author Luis Gall
 * @date 28.10.2026
 *
 * Currently covered methods of KokkosAoS: operator()(index), addParticle(index, Particle), resize(size)
 */

#ifdef AUTOPAS_ENABLE_KOKKOS

#include "KokkosAoSTest.h"

TEST_F(KokkosAoSTest, testInitialization) { autopas::utilsKokkos::KokkosAoS<autopas::KokkosParticleBaseFP32> aos; }

TEST_F(KokkosAoSTest, AoSTestAddParticle) {
    autopas::utilsKokkos::KokkosAoS<autopas::KokkosParticleBaseFP32> aos {};
    aos.resize(1);

    autopas::KokkosParticleBaseFP32 particle {};

    signed int id = 42;
    const std::array<float, 3> position {1.f,2.f,3.f};
    const std::array<float, 3> force {0.5f,0.4f,0.3f};
    const std::array<float, 3> velocity {0.123f,0.456f,0.789f};
    autopas::OwnershipState ownershipState {autopas::OwnershipState::halo};

    particle.setID(id);
    particle.setR(position);
    particle.setF(force);
    particle.setV(velocity);
    particle.setOwnershipState(ownershipState);

    aos.addParticle(0, particle);

    auto storedId = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::id, false>(0);
    auto storedX = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posX, false>(0);
    auto storedY = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posY, false>(0);
    auto storedZ = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::posZ, false>(0);
    auto storedVx = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velX, false>(0);
    auto storedVy = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velY, false>(0);
    auto storedVz = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::velZ, false>(0);
    auto storedFx = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceX, false>(0);
    auto storedFy = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceY, false>(0);
    auto storedFz = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::forceZ, false>(0);
    auto storedState = aos.template operator()<autopas::KokkosParticleBaseFP32::AttributeNames::ownershipState, false>(0);

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
