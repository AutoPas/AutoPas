/**
 * @file VerletListsCellsTest.h
 * @author nguyen
 * @date 02.09.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/Particle.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/commonTypedefs.h"

namespace VerletListsCellsTest {

class VerletListsCellsTest : public AutoPasTestBase {};

} // end namespace VerletListsCellsTest
