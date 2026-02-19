/**
 * @file CellBlock3DTest.h
 * @author tchipev
 * @date 19.01.18
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/containers/CellBlock3D.h"
#include "testingHelpers/commonTypedefs.h"

class CellBlock3DTest : public AutoPasTestBase {
 public:
  void SetUp() override;

  static std::vector<std::array<double, 3>> getMesh(std::array<double, 3> start, std::array<double, 3> dx,
                                                    std::array<int, 3> numParts);

 protected:
  std::vector<FMCell> _vec1, _vec1_cs, _vec2, _vec2_cs, _vec3, _vec4, _vec19;

  // Use unique_ptr to delay construction of the CellBlock until logger exists
  std::unique_ptr<autopas::internal::CellBlock3D<FMCell>> _cells_1x1x1;
  std::unique_ptr<autopas::internal::CellBlock3D<FMCell>> _cells_1x1x1_cs2;
  std::unique_ptr<autopas::internal::CellBlock3D<FMCell>> _cells_2x2x2;
  std::unique_ptr<autopas::internal::CellBlock3D<FMCell>> _cells_2x2x2_cs05;
  std::unique_ptr<autopas::internal::CellBlock3D<FMCell>> _cells_3x3x3;
  std::unique_ptr<autopas::internal::CellBlock3D<FMCell>> _cells_11x4x4_nonZeroBoxMin;
  std::unique_ptr<autopas::internal::CellBlock3D<FMCell>> _cells_19x19x19;
};
