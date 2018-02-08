//
// Created by ga68cat on 2/8/18.
//

#ifndef AUTOPAS_AOSVSSOATest_H
#define AUTOPAS_AOSVSSOATest_H

#include <autopas.h>
#include <gtest/gtest.h>

class AoSvsSoATest : public ::testing::Test {
 public:
  AoSvsSoATest() = default;

  ~AoSvsSoATest() override = default;

  void generateParticles(std::vector<autopas::Particle> *particles);
};

#endif  // AUTOPAS_AOSVSSOATest_H
