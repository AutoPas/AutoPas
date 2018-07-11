/**
 * @file VerletListsTest.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/Particle.h"
#include "mocks/MockFunctor.h"
#include "mocks/MockVerletLists.h"
#include "testingHelpers/RandomGenerator.h"

class VerletListsTest : public AutoPasTestBase {};
