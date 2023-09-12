/**
 * @file CellOwnershipTest.h
 * @author D. Martin
 * @date 10.09.23
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/cells/SortedCellView.h"
#include "autopas/containers/octree/OctreeLeafNode.h"
#include "autopas/containers/verletClusterLists/ClusterTower.h"
#include "testingHelpers/commonTypedefs.h"

template <typename T>
class CellOwnershipTestTyped : public AutoPasTestBase {
 public:
  CellOwnershipTestTyped() = default;

  ~CellOwnershipTestTyped() override = default;
};