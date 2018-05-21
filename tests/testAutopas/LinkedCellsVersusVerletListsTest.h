/*
 * LinkedCellsVersusDirectSumTest.h
 *
 *  Created on: 21 May 2018
 *      Author: seckler
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopasIncludes.h"

class LinkedCellsVersusVerletListsTest : public AutoPasTestBase {
 public:
  LinkedCellsVersusVerletListsTest();

  ~LinkedCellsVersusVerletListsTest() override = default;

  std::array<double, 3> getBoxMin() const { return {0.0, 0.0, 0.0}; }

  std::array<double, 3> getBoxMax() const { return {3.0, 3.0, 3.0}; }

  double getCutoff() const { return 1.0; }

 protected:
  double fRand(double fMin, double fMax) const;

  std::array<double, 3> randomPosition(
      const std::array<double, 3> &boxMin,
      const std::array<double, 3> &boxMax) const;

  void fillContainerWithMolecules(
      unsigned long numMolecules,
      autopas::ParticleContainer<autopas::MoleculeLJ,
                                 autopas::FullParticleCell<autopas::MoleculeLJ>>
          &cont) const;

  void test(unsigned long numMolecules, double rel_err_tolerance);

  autopas::VerletLists<autopas::MoleculeLJ,
                     autopas::FullParticleCell<autopas::MoleculeLJ>>
      _verletLists;
  autopas::LinkedCells<autopas::MoleculeLJ,
                       autopas::FullParticleCell<autopas::MoleculeLJ>>
      _linkedCells;
};
