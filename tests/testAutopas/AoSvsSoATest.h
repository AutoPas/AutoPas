//
// Created by ga68cat on 2/8/18.
//

#ifndef AUTOPAS_AOSVSSOATest_H
#define AUTOPAS_AOSVSSOATest_H

#include <autopasIncludes.h>
#include <gtest/gtest.h>
#include "AutoPasTest.h"

class AoSvsSoATest : public AutoPasTest {
 public:

  void generateParticles(std::vector<autopas::Particle> *particles);
};

#endif  // AUTOPAS_AOSVSSOATest_H
