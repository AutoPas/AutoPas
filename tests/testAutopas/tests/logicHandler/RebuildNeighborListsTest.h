/**
* @file RebuildNeighborListsTest.h
* @author muehlhaeusser
* @date 29.11.2024
*/

#pragma once

#include <gtest/gtest.h>

#include "../../AutoPasTestBase.h"
#include "autopas/tuning/Configuration.h"

class RebuildNeighborListsTest : public AutoPasTestBase,
                                 public ::testing::WithParamInterface<std::tuple<autopas::Configuration, autopas::Configuration>> {

 public:

  static std::set<autopas::Configuration> getPairwiseConfigs();
  static std::set<autopas::Configuration> getTriwiseConfigs();

};
