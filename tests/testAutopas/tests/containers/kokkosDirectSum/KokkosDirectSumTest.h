/**
 * @file LinkedCellsTest.h
 * @author lgaertner
 * @date 01.12.21
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ReferenceParticleCell.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/LinkedCellsReferences.h"
#include "autopas/particles/Particle.h"
#include "testingHelpers/commonTypedefs.h"

template <class TestingType>
class KokkosDirectSumTest : public AutoPasTestBase {
 public:
  using KokkosDirectSumType = typename TestingType::first_t;
  using keepListsValid = typename TestingType::second_t;

  KokkosDirectSumTest() : _kokkosDirectSum({0., 0., 0.}, {10., 10., 10.}, 1., 0.) {}

 private:
  KokkosDirectSumType _kokkosDirectSum;
  bool _keepListsValid{keepListsValid()};
};