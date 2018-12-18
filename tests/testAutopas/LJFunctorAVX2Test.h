/**
 * @file LJFunctorAVX2Test.h
 * @author F. Gratl
 * @date 12/17/18
 */

#pragma once

#include "AutoPasTestBase.h"

class LJFunctorAVX2Test : public AutoPasTestBase {
 public:
  LJFunctorAVX2Test()
      : AutoPasTestBase(), _cutoff{1.}, _epsilon{1.}, _sigma{1.}, _lowCorner{0, 0, 0}, _highCorner{2, 1, 1} {}

  const double _cutoff;
  const double _epsilon;
  const double _sigma;
  const std::array<double, 3> _lowCorner;
  const std::array<double, 3> _highCorner;
};
