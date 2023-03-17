/**
 * @file ForceSequentialTest.cpp
 * @author F. Gratl
 * @date 15.02.2023
 */

#include "ForceSequentialTest.h"

#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"

TEST_P(ForceSequentialTest, testForceSequential) {
  const auto [containerOp] = GetParam();

  autopas::AutoPas<Molecule> autoPas;
  autoPas.setAllowedContainers({containerOp});
  autoPas.setAllowedTraversals(autopas::options::TraversalOption::getAllOptions());
  autoPas.init();

  // fill AP with
}