//
// Created by ga68cat on 2/8/18.
//

#ifndef AUTOPAS_AOSVSSOATest_H
#define AUTOPAS_AOSVSSOATest_H

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopasIncludes.h"

class AoSvsSoATest : public AutoPasTestBase {
 public:
  void generateParticles(std::vector<autopas::Particle> *particles);
};

#endif  // AUTOPAS_AOSVSSOATest_H
