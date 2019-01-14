/**
 * @file CellBlock3DTest.h
 * @author tchipev
 * @date 19.01.18
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"

class CellBlock3DTest : public AutoPasTestBase {
 public:
  CellBlock3DTest()
      : _cells_1x1x1(_vec1, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({10.0, 10.0, 10.0}), 10.0),
        _cells_2x2x2(_vec2, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({10.0, 10.0, 10.0}), 5.0),
        _cells_3x3x3(_vec3, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({10.0, 10.0, 10.0}), 3.0),
        _cells_11x4x4_nonZeroBoxMin(_vec4, std::array<double, 3>({2. / 3., 0.0, 0.0}),
                                    std::array<double, 3>({1., .125, .125}), 3. / 100.) {}
  ~CellBlock3DTest() override = default;

 protected:
  std::vector<std::array<double, 3>> getMesh(std::array<double, 3> start, std::array<double, 3> dx,
                                             std::array<int, 3> numParts) const;

  std::vector<autopas::FullParticleCell<autopas::MoleculeLJ>> _vec1, _vec2, _vec3, _vec4;
  autopas::CellBlock3D<autopas::FullParticleCell<autopas::MoleculeLJ>> _cells_1x1x1, _cells_2x2x2, _cells_3x3x3,
      _cells_11x4x4_nonZeroBoxMin;
};
