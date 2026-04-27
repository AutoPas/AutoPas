/**
 * @file HGC08CellHandlerTest.h
 * @author Alexander Glas
 * @date 27.04.26
 */

#pragma once

#include <array>

#include "AutoPasTestBase.h"

class HGC08CellHandlerTest : public AutoPasTestBase {
 public:
  void runDecomposeScenario(const std::array<double, 3> &boxMin);
};