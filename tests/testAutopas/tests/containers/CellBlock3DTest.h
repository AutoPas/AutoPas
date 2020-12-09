/**
 * @file CellBlock3DTest.h
 * @author tchipev
 * @date 19.01.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/containers/CellBlock3D.h"
#include "testingHelpers/commonTypedefs.h"

namespace CellBlock3DTest {

class CellBlock3DTest : public AutoPasTestBase {
 public:
  CellBlock3DTest()
      : _cells_1x1x1(_vec1, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({10.0, 10.0, 10.0}), 10.0),
        _cells_1x1x1_cs2(_vec1_cs, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({10.0, 10.0, 10.0}),
                         5.0, 2.0),
        _cells_2x2x2(_vec2, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({10.0, 10.0, 10.0}), 5.0),
        _cells_2x2x2_cs05(_vec2_cs, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({10.0, 10.0, 10.0}),
                          10.0, 0.5),
        _cells_3x3x3(_vec3, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({10.0, 10.0, 10.0}), 3.0),
        _cells_11x4x4_nonZeroBoxMin(_vec4, std::array<double, 3>({2. / 3., 0.0, 0.0}),
                                    std::array<double, 3>({1., .125, .125}), 3. / 100.),
        _cells_19x19x19(_vec19, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({58.5, 58.5, 58.5}),
                        3.0) {}
  ~CellBlock3DTest() override = default;

  static std::vector<std::array<double, 3>> getMesh(std::array<double, 3> start, std::array<double, 3> dx,
                                                    std::array<int, 3> numParts);

 protected:
  std::vector<FMCell> _vec1, _vec1_cs, _vec2, _vec2_cs, _vec3, _vec4, _vec19;
  autopas::internal::CellBlock3D<FMCell> _cells_1x1x1, _cells_1x1x1_cs2, _cells_2x2x2, _cells_2x2x2_cs05, _cells_3x3x3,
      _cells_11x4x4_nonZeroBoxMin, _cells_19x19x19;
};

} // end namespace CellBlock3DTest
