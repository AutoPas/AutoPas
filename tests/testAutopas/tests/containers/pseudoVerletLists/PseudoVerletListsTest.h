/**
 * @file PseudoVerletListsTest.h
 * @author Lars Doll
 * @date 20.12.2025
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/containers/pseudoVerletLists/PseudoVerletLists.h"
#include "autopas/particles/ParticleDefinitions.h"

class PseudoVerletListsTest : public AutoPasTestBase, public ::testing::WithParamInterface<double> {};
