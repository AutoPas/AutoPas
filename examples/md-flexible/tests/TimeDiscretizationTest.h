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
#include "autopas/molecularDynamics/MulticenteredMoleculeLJ.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopas/utils/ArrayMath.h"

using Molecule = autopas::MoleculeLJ;
using MultisiteMolecule = autopas::MulticenteredMoleculeLJ;

namespace {
template <class MoleculeType>
void fillWithParticlesAndInit(autopas::AutoPas<MoleculeType> &autopasContainer);

template<> void fillWithParticlesAndInit<MultisiteMolecule>(autopas::AutoPas<MultisiteMolecule> &autopasContainer);

/**
 * Initialise particle properties library.
 * This function should have a valid molecule type.
 * @tparam MoleculeType
 * @param PPL
 */
template <class MoleculeType>
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