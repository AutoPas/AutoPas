/**
 * @file SingleCellIteratorTest.h
 * @author tchipev
 * @date 19.01.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "testingHelpers/commonTypedefs.h"

class SingleCellForEachTest : public AutoPasTestBase {
 public:
  SingleCellForEachTest() = default;

  void SetUp() override;

  void TearDown() override;

  ~SingleCellForEachTest() override = default;

  template <class Cell>
  void fillWithParticles(Cell *pc) {
    for (auto &m : _vecOfMolecules) {
      pc->addParticle(m);
    }
  }

  template <class Cell>
  void fillWithParticleReferences(Cell *pc) {
    for (auto &m : _vecOfMolecules) {
      pc->addParticleReference(&m);
    }
  }

  Molecule getMoleculeWithId(size_t id) {
    for (auto m : _vecOfMolecules) {
             if (m.getID() == id) return m;
         }
    autopas::utils::ExceptionHandler::exception("particle with ID " + std::to_string(id) + " not in _vecOfMolecules");
    return Molecule();
  }

 protected:
  // needs to be protected, because the test fixtures generate a derived class
  // for each unit test.
  std::vector<Molecule> _vecOfMolecules;

  std::array<double, 3> dummy{};

  template<typename Cell>
  void testCell(Cell cell, std::vector<size_t> &numMolecules, autopas::IteratorBehavior iteratorBehavior, std::array<double, 3> const lowerCorner, std::array<double, 3> const higherCorner);
};
