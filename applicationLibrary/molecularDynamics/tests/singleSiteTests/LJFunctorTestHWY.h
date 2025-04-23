/**
 * @file LJFunctorTestHWY.h
 * @author Luis Gall
 * @date 04/23/24
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/utils/SoA.h"
#include "molecularDynamicsLibrary/LJFunctorHWY.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "testingHelpers/commonTypedefs.h"

using VectorizationPattern = autopas::VectorizationPatternOption::Value;

using LJFunctorHWYTestingTuple =
    std::tuple<bool /*mixing*/, bool /*newton3*/, bool /*doDeleteSomeParticles*/, VectorizationPattern>;

class LJFunctorTestHWY : public AutoPasTestBase, public ::testing::WithParamInterface<LJFunctorHWYTestingTuple> {
 public:
  LJFunctorTestHWY() : AutoPasTestBase() {}

  constexpr static double _maxError = 1e-12;

  template <bool mixing>
  void testLJFunctorAVXvsLJFunctorHWYTwoCells(bool newton3, bool doDeleteSomeParticles, bool useUnalignedViews,
                                              VectorizationPattern pattern);

  template <bool mixing>
  void testLJFunctorAVXvsLJFunctorHWYOneCell(bool newton3, bool doDeleteSomeParticles, bool useUnalignedViews,
                                             VectorizationPattern pattern);

  template <bool mixing>
  void testLJFunctorAVXvsLJFunctorHWYVerlet(bool newton3, bool doDeleteSomeParticles);

  template <bool mixing>
  void testLJFunctorAVXvsLJFunctorHWYAoS(bool newton3, bool doDeleteSomeParticles);

  template <class SoAType>
  bool SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2);

  bool AoSParticlesEqual(FMCell &cell1, FMCell &cell2);

  bool particleEqual(Molecule &p1, Molecule &p2);

  constexpr static double _cutoff{6.};
  constexpr static double _skinPerTimestep{0.1};
  constexpr static unsigned int _rebuildFrequency{20};
  constexpr static double _interactionLengthSquare{(_cutoff + _skinPerTimestep * _rebuildFrequency) *
                                                   (_cutoff + _skinPerTimestep * _rebuildFrequency)};
  constexpr static double _epsilon{1.};
  constexpr static double _sigma{1.};
  const std::array<double, 3> _lowCorner{0., 0., 0.};
  const std::array<double, 3> _highCorner{6., 6., 6.};
};
