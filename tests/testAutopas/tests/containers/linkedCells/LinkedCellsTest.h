/**
 * @file LinkedCellsTest.h
 * @author seckler
 * @date 27.04.18
 */

#pragma once

#include <autopas/containers/linkedCells/ReferenceLinkedCells.h>

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ReferenceParticleCell.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/particles/Particle.h"
#include "testingHelpers/commonTypedefs.h"

template <class LinkedCellsType>
class LinkedCellsTest : public AutoPasTestBase {
 public:
  LinkedCellsTest() : _linkedCells({0., 0., 0.}, {10., 10., 10.}, 1., 0., 1.) {}

 protected:
  LinkedCellsType _linkedCells;
};
