#pragma once
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "AutoPasTestBase.h"
#include "Objects/Objects.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayUtils.h"
#include "testingHelpers/commonTypedefs.h"

// testet die boxMin und boxMax funktionen für alle Objects
class ObjectsTest : public AutoPasTestBase {
 public:
  ObjectsTest()
      : AutoPasTestBase(),
        _CGrid{CubeGrid(particlesPerDim, 1., velocity, center, 0, 1.0, 1.0, 1.0)},
        _CGauss{CubeGauss(numParticles, boxlength, 5., 2., velocity, center, 0, 1.0, 1.0, 1.0)},
        _CUniform{CubeUniform(numParticles, boxlength, velocity, center, 0, 1.0, 1.0, 1.0)},
        _Sphere{(Sphere(center, 5, 1., velocity, 0, 1.0, 1.0, 1.0))} {}

 protected:
  std::array<double, 3> velocity = {0., 0., 0.};
  size_t numParticles = 100;
  std::array<double, 3> center = {0., 0., 0.};
  std::array<double, 3> boxlength = {10., 10., 10.};
  std::array<size_t, 3> particlesPerDim = {10, 10, 10};
  CubeGrid _CGrid;
  CubeGauss _CGauss;
  CubeUniform _CUniform;
  Sphere _Sphere;
};