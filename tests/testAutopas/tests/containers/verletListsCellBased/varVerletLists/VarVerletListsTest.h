/**
 * @file VarVerletListsTest.h
 * @author humig
 * @date 21.05.19
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/utils/generators/UniformGenerator.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "mocks/MockPairwiseFunctor.h"
#include "testingHelpers/commonTypedefs.h"

class VarVerletListsTest : public AutoPasTestBase {};