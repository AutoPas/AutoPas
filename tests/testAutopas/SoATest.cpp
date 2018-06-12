/**
 * @file SoATest.cpp
 * @author seckler
 * @date 08.06.18
 */

#include "SoATest.h"
#include "particles/Particle.h"
#include "utils/SoA.h"
#include "utils/SoAType.h"

TEST_F(SoATest, testInitialization) { autopas::SoA<autopas::Particle> soa; }

TEST_F(SoATest, SoATypeTest) {
  static_assert(std::is_same<autopas::utils::SoAType<size_t, double, double, double>::Type,
                             std::tuple<std::vector<size_t, autopas::AlignedAllocator<size_t>>,
                                        std::vector<double, autopas::AlignedAllocator<double>>,
                                        std::vector<double, autopas::AlignedAllocator<double>>,
                                        std::vector<double, autopas::AlignedAllocator<double>>>>::value,
                "type must be same");
}