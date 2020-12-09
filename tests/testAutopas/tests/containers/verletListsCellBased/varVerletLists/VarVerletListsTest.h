/**
 * @file VarVerletListsTest.h
 * @author humig
 * @date 21.05.19
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/Particle.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/commonTypedefs.h"

namespace VarVerletListsTest {

class VarVerletListsTest : public AutoPasTestBase {};
}  // end namespace VarVerletListsTest
