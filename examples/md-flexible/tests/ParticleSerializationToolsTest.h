/**
 * @file ParticleSerializationToolsTest.h
 * @author J. Körner
 * @date 19/05/21
 */
#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/molecularDynamics/MoleculeLJ.h"
#include "src/ParticleAttributes.h"
#include "src/ParticleSerializationTools.h"

class ParticleSerializationToolsTest : public AutoPasTestBase {
 public:
  ParticleSerializationToolsTest();

 protected:
  ParticleAttributes _particle;

  autopas::MoleculeLJ<double> _molecule;
};
