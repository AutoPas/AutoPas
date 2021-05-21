/**
 * @file ParticleSerializationToolsTest.h
 * @author J. Körner
 * @date 19/05/21
 */
#pragma once

#include <gtest/gtest.h>

#include "autopas/molecularDynamics/MoleculeLJ.h"
#include "AutoPasTestBase.h"
#include "src/ParticleSerializationTools.h"

#include "src/ParticleAttributes.h"

class ParticleSerializationToolsTest : public AutoPasTestBase {
 public:
  ParticleSerializationToolsTest();

 protected:
	ParticleAttributes _particle;

	autopas::MoleculeLJ<double> _molecule;
};
