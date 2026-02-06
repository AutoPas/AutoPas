/**
 * @file ContainerSelectorTest.h
 * @author F. Gratl
 * @date 22.06.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"

class ContainerSelectorTest : public AutoPasTestBase {
 protected:
  const std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 10};
  const double cutoff = 1;
  const double cellSizeFactor = 1;
  const double verletSkin = 0.1;
  const unsigned int verletRebuildFrequency = 2;
  const bool orderCellsByMortonIndex = true;
  const bool useOptimizedLJFunctor = true;
  const bool useCompactAoS = true;
  const bool reserveVLSizes = true;
  const bool bucketSortParticles = true;
  const bool sortVerletLists = true;
  const size_t sortingFrequency = 1;
};
