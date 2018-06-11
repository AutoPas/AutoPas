/**
 * @file SoATest.cpp
 * @author seckler
 * @date 08.06.18
 */

#include "SoATest.h"
#include "particles/Particle.h"
#include "utils/SoA.h"

TEST_F(SoATest, testInitialization) { autopas::SoA<autopas::Particle> soa; }