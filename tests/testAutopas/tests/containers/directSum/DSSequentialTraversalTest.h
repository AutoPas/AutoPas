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
  void testTraversal(bool useSoA);
};
