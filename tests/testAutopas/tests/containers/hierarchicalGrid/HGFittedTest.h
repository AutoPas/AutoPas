/**
 * @file HGFittedTest.h
 * @author Alexander Glas
 * @date 27.04.26
 */

#pragma once

#include <array>
#include <string>

#include "AutoPasTestBase.h"

class HGFittedTest : public AutoPasTestBase {
 public:
  static void expectFittedHierarchy(const std::string &hierarchyDump, const std::array<double, 3> &boxMin,
                                    const std::array<double, 3> &boxMax);
};