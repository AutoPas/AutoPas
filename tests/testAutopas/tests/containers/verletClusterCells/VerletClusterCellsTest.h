/**
 * @file VerletClusterCellsTest.h
 * @author jspahl
 * @date 6.4.19
 */
#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/Particle.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/commonTypedefs.h"

namespace VerletClusterCellsTest {

class VerletClusterCellsTest : public AutoPasTestBase {};

}  // end namespace VerletClusterCellsTest
