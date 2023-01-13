/**
 * @file TimeDiscretizationTest.h
 * @author N. Fottner
 * @date 05/22/19.
 */

#pragma once

#include <gtest/gtest.h>

#include <vector>

#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "autopas/molecularDynamics/MultisiteMoleculeLJ.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopas/utils/ArrayMath.h"
#include "autopasTools/generators/GridGenerator.h"
#include "src/TimeDiscretization.h"
#include "src/configuration/MDFlexConfig.h"

using Molecule = autopas::MoleculeLJ;
using MultisiteMolecule = autopas::MultisiteMoleculeLJ;

namespace {
/**
 * Initializes an AutoPas container for a simple 5x5x5 domain, and fills it with particles.
 * @param autopasContainer
 */
void fillWithParticlesAndInit(autopas::AutoPas<ParticleType> &autopasContainer);

/**
 * Initialise particle properties library with a single site type (and multi-site molecule type, if needed)
 * @param PPL
 */
void initPPL(ParticlePropertiesLibrary<> &PPL);

/**
 * Shared implementation of testCalculateVelocities for both molecules types to reduce code duplication.
 * @tparam MoleculeType Either Molecule or MultisiteMolecule.
 */
template <class MoleculeType> void testCalculateVelocitiesImpl();
}

class TimeDiscretizationTest : public AutoPasTestBase {
 public:
  /**
   * Constructor.
   */
  TimeDiscretizationTest() = default;
};