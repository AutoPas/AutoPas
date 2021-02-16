/**
 * @file ParticleIteratorTest.h
 * @author tchipev
 * @date 19.01.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/utils/WrapOpenMP.h"
#include "testingHelpers/commonTypedefs.h"

class ParticleIteratorTest : public AutoPasTestBase {
 protected:
  /**
   * Generates a given amount of cells where only indicated cells contain a given amount of particles.
   * @param numCells
   * @param cellsToFill
   * @param particlesPerCell
   * @return Vector of generated and filled cells.
   */
  std::vector<FMCell> generateCellsWithPattern(size_t numCells, const std::vector<size_t> &cellsToFill,
                                               size_t particlesPerCell);

  /**
   * Checks that all particles in a vector of cells are accessed by the iterator.
   * @param cellsWithParticles Indices of cells that shall contain particles.
   * @param testWithAdditionalParticles Whether the iterator shall also consider a vector of additional particles
   * besides the cells
   */
  void testAllParticlesFoundPattern(const std::vector<size_t> &cellsWithParticles, bool testWithAdditionalParticles);
};
