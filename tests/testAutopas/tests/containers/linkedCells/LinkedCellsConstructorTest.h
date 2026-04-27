/**
 * @file LinkedCellsConstructorTest.h
 * @author Alexander Glas
 * @date 27.04.26
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/CellBlock3D.h"
#include "testingHelpers/commonTypedefs.h"

class LinkedCellsConstructorTest : public AutoPasTestBase {
 public:
  static std::array<unsigned long, 3> cellsPerDimFromCellSizeFactor(const std::array<double, 3> &boxMin,
                                                                     const std::array<double, 3> &boxMax,
                                                                     double interactionLength,
                                                                     double cellSizeFactor);

 protected:
  static std::vector<std::array<double, 3>> representativePositions(const std::array<double, 3> &boxMin,
                                                                    const std::array<double, 3> &boxMax);

  static std::vector<ParticleFP64> representativeOwnedParticles(const std::array<double, 3> &boxMin,
                                                                const std::array<double, 3> &boxMax);

  static std::vector<ParticleFP64> representativeHaloParticles(const std::array<double, 3> &boxMin,
                                                               const std::array<double, 3> &boxMax);
};
