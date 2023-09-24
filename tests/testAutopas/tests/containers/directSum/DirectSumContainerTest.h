/**
 * @file DirectSumTest.h
 * @author seckler
 * @date 27.04.18
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/directSum/DirectSum.h"
#include "autopas/particles/Particle.h"
#include "testingHelpers/commonTypedefs.h"

using ParamType = bool;

class DirectSumContainerTest : public AutoPasTestBase, public ::testing::WithParamInterface<ParamType> {};
