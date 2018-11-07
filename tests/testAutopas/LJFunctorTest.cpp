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
  EXPECT_NEAR(f1one[0], -4547248.8989645941, 1e-9);
  EXPECT_NEAR(f1one[1], -9094497.7979291882, 1e-9);
  EXPECT_NEAR(f1one[2], -13641746.696893783, 1e-9);

  EXPECT_DOUBLE_EQ(f2one[0], 0);
  EXPECT_DOUBLE_EQ(f2one[1], 0);
  EXPECT_DOUBLE_EQ(f2one[2], 0);

  functor.AoSFunctor(p2, p1, false);
  auto f1two = p1.getF();
  auto f2two = p2.getF();
  EXPECT_NEAR(f1two[0], -4547248.8989645941, 1e-9);
  EXPECT_NEAR(f1two[1], -9094497.7979291882, 1e-9);
  EXPECT_NEAR(f1two[2], -13641746.696893783, 1e-9);

  EXPECT_NEAR(f2two[0], +4547248.8989645941, 1e-9);
  EXPECT_NEAR(f2two[1], +9094497.7979291882, 1e-9);
  EXPECT_NEAR(f2two[2], +13641746.696893783, 1e-9);
}

TEST_F(LJFunctorTest, testAoSFunctorNoGlobalsN3) {
  double cutoff = 1.;
  double epsilon = 1.;
  double sigma = 1.;
  double shift = 0.1;
  autopas::LJFunctor<MolType, CellType> functor(cutoff, epsilon, sigma, shift);

  MolType p1({0., 0., 0.}, {0., 0., 0.}, 0);
  MolType p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1);

  functor.AoSFunctor(p1, p2, true);

  auto f1one = p1.getF();
  auto f2one = p2.getF();
  EXPECT_NEAR(f1one[0], -4547248.8989645941, 1e-9);
  EXPECT_NEAR(f1one[1], -9094497.7979291882, 1e-9);
  EXPECT_NEAR(f1one[2], -13641746.696893783, 1e-9);

  EXPECT_NEAR(f2one[0], +4547248.8989645941, 1e-9);
  EXPECT_NEAR(f2one[1], +9094497.7979291882, 1e-9);
  EXPECT_NEAR(f2one[2], +13641746.696893783, 1e-9);

  functor.AoSFunctor(p2, p1, true);
  auto f1two = p1.getF();
  auto f2two = p2.getF();
  EXPECT_NEAR(f1two[0], 2 * -4547248.8989645941, 1e-9);
  EXPECT_NEAR(f1two[1], 2 * -9094497.7979291882, 1e-9);
  EXPECT_NEAR(f1two[2], 2 * -13641746.696893783, 1e-9);

  EXPECT_NEAR(f2two[0], 2 * +4547248.8989645941, 1e-9);
  EXPECT_NEAR(f2two[1], 2 * +9094497.7979291882, 1e-9);
  EXPECT_NEAR(f2two[2], 2 * +13641746.696893783, 1e-9);
}

TEST_F(LJFunctorTest, testAoSFunctorGlobalsThrowBad) {
  double cutoff = 1.;
  double epsilon = 1.;
  double sigma = 1.;
  double shift = 0.1;
  std::array<double, 3> lowCorner = {0., 0., 0.};
  std::array<double, 3> highCorner = {5., 5., 5.};
  bool duplicatedCalculation = true;
  typedef autopas::utils::ExceptionHandler::AutoPasException exception_type;
  {
    // throw if lowcorner == highcorner, but calculateglobals and duplicatedCalculation are true
    typedef autopas::LJFunctor<MolType, CellType, true> functortype;
    EXPECT_THROW(functortype functor(cutoff, epsilon, sigma, shift, lowCorner, {0., 0., 0.}, duplicatedCalculation),
                 exception_type);
  }

  autopas::LJFunctor<MolType, CellType, true> functor(cutoff, epsilon, sigma, shift, lowCorner, highCorner,
                                                      duplicatedCalculation);

  // getupot without postprocessing is not allowed
  EXPECT_THROW(functor.getUpot(), exception_type);
  EXPECT_THROW(functor.getVirial(), exception_type);

  EXPECT_NO_THROW(functor.resetGlobalValues());

  EXPECT_NO_THROW(functor.postProcessGlobalValues(true));
  EXPECT_NO_THROW(functor.resetGlobalValues());
  EXPECT_NO_THROW(functor.postProcessGlobalValues(true));
  // repeated postprocessing is not allowed
  EXPECT_THROW(functor.postProcessGlobalValues(true), exception_type);

  EXPECT_NO_THROW(functor.resetGlobalValues());
  EXPECT_NO_THROW(functor.postProcessGlobalValues(true));
}

TEST_F(LJFunctorTest, testAoSFunctorGlobalsInsideDuplicated) {
  double cutoff = 1.;
  double epsilon = 1.;
  double sigma = 1.;
  double shift = 0.1;
  std::array<double, 3> lowCorner = {0., 0., 0.};
  std::array<double, 3> highCorner = {5., 5., 5.};
  bool duplicatedCalculation = true;

  autopas::LJFunctor<MolType, CellType, true> functor(cutoff, epsilon, sigma, shift, lowCorner, highCorner,
                                                      duplicatedCalculation);

  MolType p1({0., 0., 0.}, {0., 0., 0.}, 0);
  MolType p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1);
  {
    // no newton 3
    functor.resetGlobalValues();

    functor.AoSFunctor(p1, p2, /*newton3*/ false);
    functor.AoSFunctor(p2, p1, /*newton3*/ false);

    functor.postProcessGlobalValues(/*newton3*/ false);

    double upot = functor.getUpot();
    double virial = functor.getVirial();

    EXPECT_NEAR(upot, 3178701.6514326506, 1.e-10);
    EXPECT_NEAR(virial, 6366148.4585504318, 1.e-10);
  }
  {
    // with newton 3
    functor.resetGlobalValues();

    functor.AoSFunctor(p1, p2, /*newton3*/ true);

    functor.postProcessGlobalValues(/*newton3*/ true);

    double upot = functor.getUpot();
    double virial = functor.getVirial();

    EXPECT_NEAR(upot, 3178701.6514326506, 1.e-10);
    EXPECT_NEAR(virial, 6366148.4585504318, 1.e-10);
  }
}

TEST_F(LJFunctorTest, testAoSFunctorGlobalsBoundaryDuplicated) {
  double cutoff = 1.;
  double epsilon = 1.;
  double sigma = 1.;
  double shift = 0.1;
  std::array<double, 3> lowCorner = {0., 0., 0.};
  std::array<double, 3> highCorner = {5., 5., 5.};
  bool duplicatedCalculation = true;

  autopas::LJFunctor<MolType, CellType, true> functor(cutoff, epsilon, sigma, shift, lowCorner, highCorner,
                                                      duplicatedCalculation);

  MolType p1({4.9, 0., 0.}, {0., 0., 0.}, 0);
  MolType p2({5.0, 0.2, 0.3}, {0., 0., 0.}, 1);
  {
    // no newton 3
    functor.resetGlobalValues();

    functor.AoSFunctor(p1, p2, /*newton3*/ false);
    functor.AoSFunctor(p2, p1, /*newton3*/ false);

    functor.postProcessGlobalValues(/*newton3*/ false);

    double upot = functor.getUpot();
    double virial = functor.getVirial();

    EXPECT_NEAR(upot, 3178701.65143265 * 0.5, 1.e-8);
    EXPECT_NEAR(virial, 6366148.4585504 * 0.5, 1.e-7);
  }
  {
    // with newton 3
    functor.resetGlobalValues();

    functor.AoSFunctor(p1, p2, /*newton3*/ true);

    functor.postProcessGlobalValues(/*newton3*/ true);

    double upot = functor.getUpot();
    double virial = functor.getVirial();

    EXPECT_NEAR(upot, 3178701.65143265 * 0.5, 1.e-8);
    EXPECT_NEAR(virial, 6366148.4585504 * 0.5, 1.e-7);
  }
}

TEST_F(LJFunctorTest, testAoSFunctorGlobalsInsideNonDuplicated) {
  double cutoff = 1.;
  double epsilon = 1.;
  double sigma = 1.;
  double shift = 0.1;
  std::array<double, 3> lowCorner = {0., 0., 0.};
  std::array<double, 3> highCorner = {5., 5., 5.};
  bool duplicatedCalculation = false;

  autopas::LJFunctor<MolType, CellType, true> functor(cutoff, epsilon, sigma, shift, lowCorner, highCorner,
                                                      duplicatedCalculation);

  MolType p1({0., 0., 0.}, {0., 0., 0.}, 0);
  MolType p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1);
  {
    // no newton 3
    functor.resetGlobalValues();

    functor.AoSFunctor(p1, p2, /*newton3*/ false);
    functor.AoSFunctor(p2, p1, /*newton3*/ false);

    functor.postProcessGlobalValues(/*newton3*/ false);

    double upot = functor.getUpot();
    double virial = functor.getVirial();

    EXPECT_NEAR(upot, 3178701.6514326506, 1.e-10);
    EXPECT_NEAR(virial, 6366148.4585504318, 1.e-10);
  }
  {
    // with newton 3
    functor.resetGlobalValues();

    functor.AoSFunctor(p1, p2, /*newton3*/ true);

    functor.postProcessGlobalValues(/*newton3*/ true);

    double upot = functor.getUpot();
    double virial = functor.getVirial();

    EXPECT_NEAR(upot, 3178701.6514326506, 1.e-10);
    EXPECT_NEAR(virial, 6366148.4585504318, 1.e-10);
  }
}

TEST_F(LJFunctorTest, testAoSFunctorGlobalsBoundaryNonDuplicated) {
  double cutoff = 1.;
  double epsilon = 1.;
  double sigma = 1.;
  double shift = 0.1;
  std::array<double, 3> lowCorner = {0., 0., 0.};
  std::array<double, 3> highCorner = {5., 5., 5.};
  bool duplicatedCalculation = false;

  autopas::LJFunctor<MolType, CellType, true> functor(cutoff, epsilon, sigma, shift, lowCorner, highCorner,
                                                      duplicatedCalculation);

  MolType p1({4.9, 0., 0.}, {0., 0., 0.}, 0);
  MolType p2({5.0, 0.2, 0.3}, {0., 0., 0.}, 1);
  {
    // no newton 3
    functor.resetGlobalValues();

    functor.AoSFunctor(p1, p2, /*newton3*/ false);
    functor.AoSFunctor(p2, p1, /*newton3*/ false);

    functor.postProcessGlobalValues(/*newton3*/ false);

    double upot = functor.getUpot();
    double virial = functor.getVirial();

    EXPECT_NEAR(upot, 3178701.65143265, 1.e-7);
    EXPECT_NEAR(virial, 6366148.4585504, 1.e-7);
  }
  {
    // with newton 3
    functor.resetGlobalValues();

    functor.AoSFunctor(p1, p2, /*newton3*/ true);

    functor.postProcessGlobalValues(/*newton3*/ true);

    double upot = functor.getUpot();
    double virial = functor.getVirial();

    EXPECT_NEAR(upot, 3178701.65143265, 1.e-7);
    EXPECT_NEAR(virial, 6366148.4585504, 1.e-7);
  }
}