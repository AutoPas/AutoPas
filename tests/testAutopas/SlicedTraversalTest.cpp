/*
 * SlicedTraversalTest.cpp
 *
 *  Created on: 22 Jan 2018
 *      Author: gratl
 */

#include "SlicedTraversalTest.h"

typedef MockFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> MFunctor;
typedef autopas::CellFunctor<autopas::Particle,
                             autopas::FullParticleCell<autopas::Particle>, MFunctor, true, true> MCellFunctor;
typedef autopas::FullParticleCell<autopas::Particle> FPCell;


TEST(SlicedTraversalTest, testIsApplicableTooSmall) {
//  AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> autoPas;

//  autoPas.init({0,0,0},{1,1,1}, 1, autopas::linkedCells);

//  MFunctor functor;
//  MCellFunctor cellFunctor(&functor);

  std::vector<FPCell> cells;

  autopas::SlicedTraversal<FPCell, MCellFunctor> slicedTraversal(cells, {1,1,1}, nullptr);

  EXPECT_FALSE(slicedTraversal.isApplicable());
}

TEST(SlicedTraversalTest, testIsApplicableOk) {
  std::vector<FPCell> cells;

  autopas::SlicedTraversal<FPCell, MCellFunctor> slicedTraversal(cells, {5,5,5}, nullptr);

  EXPECT_TRUE(slicedTraversal.isApplicable());
}

TEST(SlicedTraversalTest, testIsApplicableOkOnlyOneDim) {

  std::vector<FPCell> cells;

  autopas::SlicedTraversal<FPCell, MCellFunctor> slicedTraversal(cells, {1,1,5}, nullptr);

  EXPECT_TRUE(slicedTraversal.isApplicable());
}