/**
 * @file TraversalComparison.h
 * @author humig
 * @date 12.07.19
 */

#pragma once

#include <gtest/gtest.h>
#include <cstdlib>
#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"
#include "testingHelpers/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * The tests in this class compare the calculated forces from all aos and soa traversals with a reference result.
 */
class TraversalComparison
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<std::tuple<autopas::ContainerOption, autopas::TraversalOption,
                                                      autopas::DataLayoutOption, autopas::Newton3Option>> {
 public:
  static void SetUpTestSuite();

 protected:
  static std::vector<std::array<double, 3>> calculateForces(autopas::ContainerOption containerOption,
                                                            autopas::TraversalOption traversalOption,
                                                            autopas::DataLayoutOption dataLayoutOption,
                                                            autopas::Newton3Option newton3Option,
                                                            unsigned long numMolecules, std::array<double, 3> boxMax);

  static std::array<double, 3> _boxMin;
  static std::vector<std::array<double, 3>> _boxMaxVector;
  static double _cutoff;

  static double _eps;
  static double _sig;
  static double _shift;
  static autopas::LJFunctor<Molecule, FMCell> _functor;

  static std::map<std::pair<int, std::array<double, 3>>, std::vector<std::array<double, 3>>> _forcesReference;

  static std::vector<int> _numParticlesVector;
};

std::array<double, 3> TraversalComparison::_boxMin{0, 0, 0};
std::vector<std::array<double, 3>> TraversalComparison::_boxMaxVector{{3, 3, 3}, {10, 10, 10}};
double TraversalComparison::_cutoff = 1.0;

double TraversalComparison::_eps = 1.0;
double TraversalComparison::_sig = 1.0;
double TraversalComparison::_shift = 0.0;
autopas::LJFunctor<Molecule, FMCell> TraversalComparison::_functor{_cutoff, _eps, _sig, _shift};

std::vector<int> TraversalComparison::_numParticlesVector{100, 1000, 2000};

std::map<std::pair<int, std::array<double, 3>>, std::vector<std::array<double, 3>>>
    TraversalComparison::_forcesReference{};