/*
 * SlicedTraversalTest.h
 *
 *  Created on: 22 Jan 2018
 *      Author: gratl
 */

#pragma once

#include <AutoPas.h>
#include <gtest/gtest.h>
#include <mocks/MockFunctor.h>
#include <testingHelpers/commonTypedefs.h>

class SlicedTraversalTest : public testing::Test {
 public:
  SlicedTraversalTest() = default;

  ~SlicedTraversalTest() override = default;

  void fillWithParticles(std::vector<FPCell> &cells, std::array<size_t, 3> particlesPerDim);
};
