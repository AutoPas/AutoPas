/**
 * @file GeneratorsTest.h
 * @author N. Fottner
 * @date 02/08/19
 */
#pragma once
#include <gtest/gtest.h>

#include <array>

#include "AutoPasTestBase.h"

class GeneratorsTest : public AutoPasTestBase {
 public:
  GeneratorsTest() = default;

 protected:
  double epsilon{1.0};
  double sigma{1.0};
  double cutoff{1.};
  std::array<double, 3> boxmin{{0., 0., 0.}};
  std::array<double, 3> boxmax{{5., 5., 5.}};
};