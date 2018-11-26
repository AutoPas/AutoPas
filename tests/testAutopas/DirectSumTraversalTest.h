/**
 * @file DirectSumTraversalTest.h
 * @author F. Gratl
 * @date 11/23/18
 */

#pragma once

#include "AutoPasTestBase.h"
#include "testingHelpers/commonTypedefs.h"

class DirectSumTraversalTest : public AutoPasTestBase {
 public:
  void testTraversal(bool useSoA);
};
