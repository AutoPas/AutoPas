/**
 * @file LJFunctorTestHWY.h
 * @author Luis Gall
 * @date 04/23/24
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/utils/SoA.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "testingHelpers/commonTypedefs.h"
#include "molecularDynamicsLibrary/LJFunctorHWY.h"

using VectorizationPattern = autopas::VectorizationPatternOption::Value;

using LJFunctorHWYTestingTuple = std::tuple<bool /*newton3*/, bool /*doDeleteSomeParticles*/, VectorizationPattern>;

class LJFunctorTestHWY : public AutoPasTestBase, public ::testing::WithParamInterface<LJFunctorHWYTestingTuple> {
 public:
  LJFunctorTestHWY() : AutoPasTestBase() {}

  constexpr static double _maxError = 1e-12;

  template <VectorizationPattern vecPattern = VectorizationPattern::p1xVec>
  void testLJFunctorAVXvsLJFunctorHWYTwoCells(bool newton3, bool doDeleteSomeParticles, bool useUnalignedViews);

  template <VectorizationPattern vecPattern = VectorizationPattern::p1xVec>
  void testLJFunctorAVXvsLJFunctorHWYOneCell(bool newton3, bool doDeleteSomeParticles, bool useUnalignedViews);

  template <VectorizationPattern vecPattern = VectorizationPattern::p1xVec>
  void testLJFunctorAVXvsLJFunctorHWYVerlet(bool newton3, bool doDeleteSomeParticles);

  void testLJFunctorAVXvsLJFunctorHWYAoS(bool newton3, bool doDeleteSomeParticles);

  template <class SoAType>
  bool SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2);

  bool AoSParticlesEqual(FMCell &cell1, FMCell &cell2);

  bool particleEqual(Particle &p1, Particle &p2);

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
