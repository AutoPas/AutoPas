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
#include "autopas/utils/ArrayMath.h"
#include "autopasTools/generators/GridGenerator.h"
#include "molecularDynamicsLibrary/MultisiteMoleculeLJ.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "src/TimeDiscretization.h"
#include "src/configuration/MDFlexConfig.h"

#if MD_FLEXIBLE_MODE!=MULTISITE or not defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH)
#include "src/MoleculeContainer.h"
#endif

namespace {
/**
 * Initializes an AutoPas container for a simple 5x5x5 domain, and fills it with particles.
 * @param autopasContainer
 */
#if MD_FLEXIBLE_MODE!=MULTISITE or not defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH)
void fillWithParticlesAndInit(autopas::AutoPas<ParticleType> &autopasContainer);
#else
void fillWithParticlesAndInit(autopas::AutoPas<ParticleType> &autopasContainer, MoleculeContainer& moleculeContainer, ParticlePropertiesLibrary<> &PPL);
#endif

/**
 * Initialise particle properties library with a single site type (and multi-site molecule type, if needed)
 * @param PPL
 */
void initPPL(ParticlePropertiesLibrary<> &PPL);

}  // namespace

class TimeDiscretizationTest : public AutoPasTestBase {
 public:
  /**
   * Constructor.
   */
  TimeDiscretizationTest() = default;
};