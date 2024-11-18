/**
 * @file DSSequentialTraversalTest.h
 * @author F. Gratl
 * @date 11/23/18
 */

#pragma once

#include "AutoPasTestBase.h"
#include "testingHelpers/commonTypedefs.h"

class DSSequentialTraversalTest : public AutoPasTestBase {
 public:
  void testTraversalPairwise(bool useSoA);
  void testTraversalTriwise(bool useSoA);

 private:
  std::vector<FPCell> fillParticleCells(size_t numParticles, size_t numHaloParticlesPerCell);
};
