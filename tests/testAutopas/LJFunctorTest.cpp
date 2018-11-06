/**
 * @file LJFunctorTest.cpp
 * @author seckler
 * @date 06.11.18
 */

#include "LJFunctorTest.h"
#include <autopas/particles/MoleculeLJ.h>

typedef autopas::MoleculeLJ MolType;
typedef autopas::FullParticleCell<MolType> CellType;

TEST_F(LJFunctorTest, testAoSFunctorNoGlobalsNoN3) {
  double cutoff = 1.;
  double epsilon = 1.;
  double sigma = 1.;
  double shift = 0.1;
  autopas::LJFunctor<MolType, CellType> functor(cutoff, epsilon, sigma, shift);

  MolType p1({0., 0., 0.}, {0., 0., 0.}, 0);
  MolType p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1);

  functor.AoSFunctor(p1, p2, false);

  auto f1one = p1.getF();
  auto f2one = p2.getF();
  EXPECT_NEAR(f1one[0],-4547248.8989645941,1e-9);
  EXPECT_NEAR(f1one[1],-9094497.7979291882,1e-9);
  EXPECT_NEAR(f1one[2],-13641746.696893783,1e-9);

  EXPECT_DOUBLE_EQ(f2one[0],0);
  EXPECT_DOUBLE_EQ(f2one[1],0);
  EXPECT_DOUBLE_EQ(f2one[2],0);

  functor.AoSFunctor(p2, p1, false);
  auto f1two = p1.getF();
  auto f2two = p2.getF();
  EXPECT_NEAR(f1two[0],-4547248.8989645941,1e-9);
  EXPECT_NEAR(f1two[1],-9094497.7979291882,1e-9);
  EXPECT_NEAR(f1two[2],-13641746.696893783,1e-9);

  EXPECT_NEAR(f2two[0],+4547248.8989645941,1e-9);
  EXPECT_NEAR(f2two[1],+9094497.7979291882,1e-9);
  EXPECT_NEAR(f2two[2],+13641746.696893783,1e-9);
}
