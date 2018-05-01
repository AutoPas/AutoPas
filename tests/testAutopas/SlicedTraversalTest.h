/*
 * SlicedTraversalTest.h
 *
 *  Created on: 22 Jan 2018
 *      Author: gratl
 */

#pragma once

#include <gtest/gtest.h>
#include <autopasIncludes.h>
#include <mocks/MockFunctor.h>
//#include "../../examples/md/mdutils.h"

class MyParticle : public autopas::Particle {
 public:
  MyParticle(std::array<double, 3> pos, unsigned long id)
      : autopas::Particle(pos, {0, 0, 0}, id){};

 private:
};

class SlicedTraversalTest : public testing::Test {
 public:
  SlicedTraversalTest() = default;

  ~SlicedTraversalTest() override = default;

};
