/**
 * @file PairwiseVerletListsTest.h
 * @author tirgendetwas
 * @date 22.11.20
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/Particle.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/commonTypedefs.h"

class PairwiseVerletListsTest : public AutoPasTestBase {};