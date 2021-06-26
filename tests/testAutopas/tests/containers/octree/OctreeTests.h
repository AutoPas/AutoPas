/**
 * @file OctreeTests.h
 * @author Johannes Spies
 * @date 15.04.2021
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
//#include "autopas/containers/CellBlock3D.h"
//#include "testingHelpers/commonTypedefs.h"

using GeneratorSpec = std::tuple<int unsigned, int unsigned>;

class OctreeTest : public AutoPasTestBase, public ::testing::WithParamInterface<GeneratorSpec> {};
