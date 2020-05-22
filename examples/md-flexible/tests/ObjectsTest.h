/**
 * @file ObjectsTest.h
 * @author N. Fottner
 * @date 02/08/19
 */
#pragma once
#include <gtest/gtest.h>

#include <array>

#include "AutoPasTestBase.h"
#include "src/Objects/CubeGauss.h"
#include "src/Objects/CubeGrid.h"
#include "src/Objects/CubeUniform.h"
#include "src/Objects/Sphere.h"

// testet die boxMin und boxMax funktionen für alle Objects
class ObjectsTest : public AutoPasTestBase {
 public:
  ObjectsTest()
      : AutoPasTestBase(),
        _CGrid{velocity, 0, 1., 1., 1., particlesPerDim, 1, center},
        _CGauss{velocity, 0, 1., 1., 1., numParticles, boxlength, {5., 5., 5.}, {2., 2., 2.}, center},
        _CUniform{velocity, 0, 1., 1., 1., numParticles, boxlength, center},
        _Sphere{velocity, 0, 1., 1., 1., center, 5, 1.} {}

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