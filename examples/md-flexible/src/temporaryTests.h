//
// Created by johnny on 01.12.23.
//
#pragma once
#include "TypeDefinitions.h"
#include "autopas/AutoPas.h"
#include "MoleculeContainer.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Quaternion.h"
#include "TimeDiscretization.h"

#if defined(I_JUST_WANNA_RUN_UNIT_TESTS)
namespace temporaryTests {

void fillWithParticlesAndInit(autopas::AutoPas<ParticleType> &autopasContainer, MoleculeContainer &moleculeContainer,
                              ParticlePropertiesLibrary<> &PPL);

void initPPL(ParticlePropertiesLibrary<> &PPL);

void testCalculateVelocities();

void testCalculateQuaternion();

void testCalculateAngularVelocities();

void testCalculatePositions();

}
#endif