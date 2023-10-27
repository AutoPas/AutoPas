/**
 * @file DSSequentialTraversalTest3B.h
 * @author muehlhaeusser
 * @date 26.10.23
 */

#pragma once

#include "AutoPasTestBase.h"
#include "testingHelpers/commonTypedefs.h"

class DSSequentialTraversalTest3B : public AutoPasTestBase {
 public:
  void testTraversal(bool useSoA);
};
