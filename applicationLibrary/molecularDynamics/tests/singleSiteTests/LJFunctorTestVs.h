/**
 * @file LJFunctorTestVs.h
 * @author F. Gratl
 * @date 20.03.20
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "LJFunctorTest.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"

template <class FuncType>
class LJFunctorTestVs : public LJFunctorTest {
 public:
  LJFunctorTestVs() : LJFunctorTest() {}

  constexpr static double cutoff{1.};
  constexpr static double epsilon{1.};
  constexpr static double sigma{1.};
};
