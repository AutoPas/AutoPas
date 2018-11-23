/**
 * @file LJFunctorTest.cpp
 * @author seckler
 * @date 06.11.18
 */

#include "LJFunctorTest.h"
#include <autopas/particles/MoleculeLJ.h>
#include "testingHelpers/commonTypedefs.h"

void LJFunctorTest::testAoSNoGlobals(bool newton3) {
  autopas::LJFunctor<Molecule, FMCell> functor(cutoff, epsilon, sigma, shift);

  Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0);
  Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1);

  functor.AoSFunctor(p1, p2, newton3);

  auto f1one = p1.getF();
  auto f2one = p2.getF();
  EXPECT_NEAR(f1one[0], expectedForce[0], absDelta);
  EXPECT_NEAR(f1one[1], expectedForce[1], absDelta);
  EXPECT_NEAR(f1one[2], expectedForce[2], absDelta);

  if (newton3) {
    EXPECT_NEAR(f2one[0], -expectedForce[0], absDelta);
    EXPECT_NEAR(f2one[1], -expectedForce[1], absDelta);
    EXPECT_NEAR(f2one[2], -expectedForce[2], absDelta);
  } else {
    EXPECT_DOUBLE_EQ(f2one[0], 0);
    EXPECT_DOUBLE_EQ(f2one[1], 0);
    EXPECT_DOUBLE_EQ(f2one[2], 0);
  }

  functor.AoSFunctor(p2, p1, newton3);

  auto f1two = p1.getF();
  auto f2two = p2.getF();

  double factor = newton3 ? 2. : 1.;

  EXPECT_NEAR(f1two[0], factor * expectedForce[0], absDelta);
  EXPECT_NEAR(f1two[1], factor * expectedForce[1], absDelta);
  EXPECT_NEAR(f1two[2], factor * expectedForce[2], absDelta);

  EXPECT_NEAR(f2two[0], -factor * expectedForce[0], absDelta);
  EXPECT_NEAR(f2two[1], -factor * expectedForce[1], absDelta);
  EXPECT_NEAR(f2two[2], -factor * expectedForce[2], absDelta);
}

TEST_F(LJFunctorTest, testAoSFunctorNoGlobalsNoN3) {
  bool newton3 = false;
  testAoSNoGlobals(newton3);
}

TEST_F(LJFunctorTest, testAoSFunctorNoGlobalsN3) {
  bool newton3 = true;
  testAoSNoGlobals(newton3);
}

void LJFunctorTest::testSoANoGlobals(bool newton3, CellInteractionType interactionType) {
  autopas::LJFunctor<Molecule, FMCell> functor(cutoff, epsilon, sigma, shift);

  FMCell cell1, cell2;
  {
    Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0);
    cell1.addParticle(p1);
    Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1);
    switch (interactionType) {
      case CellInteractionType::own:
        cell1.addParticle(p2);
        break;
      case CellInteractionType::pair:
        cell2.addParticle(p2);
        break;
      default:
        FAIL();
    }
  }

  functor.SoALoader(cell1, cell1._particleSoABuffer);
  functor.SoALoader(cell2, cell2._particleSoABuffer);

  switch (interactionType) {
    case CellInteractionType::own:
      functor.SoAFunctor(cell1._particleSoABuffer, newton3);
      break;
    case CellInteractionType::pair:
      functor.SoAFunctor(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
      break;
  }
  functor.SoAExtractor(cell1, cell1._particleSoABuffer);
  functor.SoAExtractor(cell2, cell2._particleSoABuffer);

  auto f1one = cell1.begin()->getF();

  EXPECT_NEAR(f1one[0], expectedForce[0], absDelta);
  EXPECT_NEAR(f1one[1], expectedForce[1], absDelta);
  EXPECT_NEAR(f1one[2], expectedForce[2], absDelta);

  std::array<double, 3> f2one = {0., 0., 0.};
  switch (interactionType) {
    case CellInteractionType::own:
      f2one = (++cell1.begin())->getF();
      break;
    case CellInteractionType::pair:
      f2one = cell2.begin()->getF();
      break;
  }
  // if the interactiontype is own, then the forces of the second particle should always be calculated!
  if (newton3 or interactionType == CellInteractionType::own) {
    EXPECT_NEAR(f2one[0], -expectedForce[0], absDelta);
    EXPECT_NEAR(f2one[1], -expectedForce[1], absDelta);
    EXPECT_NEAR(f2one[2], -expectedForce[2], absDelta);
  } else {
    EXPECT_DOUBLE_EQ(f2one[0], 0);
    EXPECT_DOUBLE_EQ(f2one[1], 0);
    EXPECT_DOUBLE_EQ(f2one[2], 0);
  }

  if (interactionType == CellInteractionType::pair) {
    functor.SoALoader(cell1, cell1._particleSoABuffer);
    functor.SoALoader(cell2, cell2._particleSoABuffer);
    functor.SoAFunctor(cell2._particleSoABuffer, cell1._particleSoABuffer, newton3);
    functor.SoAExtractor(cell1, cell1._particleSoABuffer);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer);

    auto f1two = cell1.begin()->getF();
    auto f2two = cell2.begin()->getF();

    double factor = newton3 ? 2. : 1.;

    EXPECT_NEAR(f1two[0], factor * expectedForce[0], absDelta);
    EXPECT_NEAR(f1two[1], factor * expectedForce[1], absDelta);
    EXPECT_NEAR(f1two[2], factor * expectedForce[2], absDelta);

    EXPECT_NEAR(f2two[0], -factor * expectedForce[0], absDelta);
    EXPECT_NEAR(f2two[1], -factor * expectedForce[1], absDelta);
    EXPECT_NEAR(f2two[2], -factor * expectedForce[2], absDelta);
  }
}

TEST_F(LJFunctorTest, testSoAFunctorNoGlobalsNoN3Own) {
  bool newton3 = false;
  testSoANoGlobals(newton3, own);
}

TEST_F(LJFunctorTest, testSoAFunctorNoGlobalsN3Own) {
  bool newton3 = true;
  testSoANoGlobals(newton3, own);
}

TEST_F(LJFunctorTest, testSoAFunctorNoGlobalsNoN3Pair) {
  bool newton3 = false;
  testSoANoGlobals(newton3, pair);
}

TEST_F(LJFunctorTest, testSoAFunctorNoGlobalsN3Pair) {
  bool newton3 = true;
  testSoANoGlobals(newton3, pair);
}

TEST_F(LJFunctorTest, testFunctorGlobalsThrowBad) {
  bool duplicatedCalculation = true;
  typedef autopas::utils::ExceptionHandler::AutoPasException exception_type;
  {
    // throw if lowcorner == highcorner, but calculateglobals and duplicatedCalculation are true
    typedef autopas::LJFunctor<Molecule, FMCell, true> functortype;
    EXPECT_THROW(functortype functor(cutoff, epsilon, sigma, shift, lowCorner, {0., 0., 0.}, duplicatedCalculation),
                 exception_type);
  }

  autopas::LJFunctor<Molecule, FMCell, true> functor(cutoff, epsilon, sigma, shift, lowCorner, highCorner,
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

void LJFunctorTest::testAoSGlobals(LJFunctorTest::where_type where, bool newton3, bool duplicatedCalculation) {
  autopas::LJFunctor<Molecule, FMCell, true> functor(cutoff, epsilon, sigma, shift, lowCorner, highCorner,
                                                     duplicatedCalculation);
  double xOffset;
  double whereFactor;
  std::string where_str;
  switch (where) {
    case inside:
      xOffset = 0.;
      whereFactor = 1.;
      where_str = "inside";
      break;
    case boundary:
      xOffset = 4.9;
      // if there are no duplicated calculations all calculations count, therefore factor = 1
      // if there are duplicated calculations there shouldn't be only a partial (factor 0.5) contribution to the energy
      // if one particle is inside and one outside
      whereFactor = duplicatedCalculation ? 0.5 : 1;
      where_str = "boundary";
      break;
    case outside:
      xOffset = 5.0;
      // if there are no duplicated calculations all calculations count, therefore factor = 1
      // if there are duplicated calculations there shouldn't be any contribution to the energy if both particles are
      // outside
      whereFactor = duplicatedCalculation ? 0. : 1;
      where_str = "outside";
      break;
    default:
      throw "not in enum where_type";
  }
  Molecule p1({0. + xOffset, 0., 0.}, {0., 0., 0.}, 0);
  Molecule p2({0.1 + xOffset, 0.2, 0.3}, {0., 0., 0.}, 1);

  functor.resetGlobalValues();

  functor.AoSFunctor(p1, p2, newton3);
  if (not newton3) {
    functor.AoSFunctor(p2, p1, newton3);
  }
  functor.postProcessGlobalValues(newton3);

  double upot = functor.getUpot();
  double virial = functor.getVirial();

  EXPECT_NEAR(upot, whereFactor * expectedEnergy, absDelta)
      << "where: " << where_str << ", newton3: " << newton3 << ", duplicatedCalculation:" << duplicatedCalculation;
  EXPECT_NEAR(virial, whereFactor * expectedVirial, absDelta)
      << "where: " << where_str << ", newton3: " << newton3 << ", duplicatedCalculation:" << duplicatedCalculation;
}

TEST_F(LJFunctorTest, testAoSFunctorGlobals) {
  for (bool duplicatedCalculation : {false, true}) {
    for (where_type where : {inside, boundary, outside}) {
      for (bool newton3 : {false, true}) {
        testAoSGlobals(where, newton3, duplicatedCalculation);
      }
    }
  }
}

void LJFunctorTest::testSoAGlobals(LJFunctorTest::where_type where, bool newton3, bool duplicatedCalculation,
                                   CellInteractionType interactionType) {
  autopas::LJFunctor<Molecule, FMCell, true> functor(cutoff, epsilon, sigma, shift, lowCorner, highCorner,
                                                     duplicatedCalculation);
  double xOffset;
  double whereFactor;
  std::string where_str;
  switch (where) {
    case inside:
      xOffset = 0.;
      whereFactor = 1.;
      where_str = "inside";
      break;
    case boundary:
      xOffset = 4.9;
      // if there are no duplicated calculations all calculations count, therefore factor = 1
      // if there are duplicated calculations there shouldn't be only a partial (factor 0.5) contribution to the energy
      // if one particle is inside and one outside
      whereFactor = duplicatedCalculation ? 0.5 : 1;
      where_str = "boundary";
      break;
    case outside:
      xOffset = 5.0;
      // if there are no duplicated calculations all calculations count, therefore factor = 1
      // if there are duplicated calculations there shouldn't be any contribution to the energy if both particles are
      // outside
      whereFactor = duplicatedCalculation ? 0. : 1;
      where_str = "outside";
      break;
    default:
      throw "not in enum where_type";
  }
  FMCell cell1, cell2;
  {
    Molecule p1({0. + xOffset, 0., 0.}, {0., 0., 0.}, 0);
    cell1.addParticle(p1);
    Molecule p2({0.1 + xOffset, 0.2, 0.3}, {0., 0., 0.}, 1);
    switch (interactionType) {
      case CellInteractionType::own:
        cell1.addParticle(p2);
        break;
      case CellInteractionType::pair:
        cell2.addParticle(p2);
        break;
      default:
        FAIL();
    }
  }

  functor.resetGlobalValues();

  functor.SoALoader(cell1, cell1._particleSoABuffer);
  functor.SoALoader(cell2, cell2._particleSoABuffer);

  switch (interactionType) {
    case CellInteractionType::own:
      functor.SoAFunctor(cell1._particleSoABuffer, newton3);
      break;
    case CellInteractionType::pair:
      functor.SoAFunctor(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
      if (not newton3) {
        functor.SoAFunctor(cell2._particleSoABuffer, cell1._particleSoABuffer, newton3);
      }
      break;
  }
  functor.SoAExtractor(cell1, cell1._particleSoABuffer);
  functor.SoAExtractor(cell2, cell2._particleSoABuffer);

  functor.postProcessGlobalValues(newton3);

  double upot = functor.getUpot();
  double virial = functor.getVirial();

  EXPECT_NEAR(upot, whereFactor * expectedEnergy, absDelta)
      << "where: " << where_str << ", newton3: " << newton3 << ", duplicatedCalculation:" << duplicatedCalculation
      << ", interactionType: " << (interactionType == pair ? "pair" : "own");
  EXPECT_NEAR(virial, whereFactor * expectedVirial, absDelta)
      << "where: " << where_str << ", newton3: " << newton3 << ", duplicatedCalculation:" << duplicatedCalculation
      << ", interactionType: " << (interactionType == pair ? "pair" : "own");
}

TEST_F(LJFunctorTest, testSoAFunctorGlobalsOwn) {
  for (bool duplicatedCalculation : {false, true}) {
    // the own functor can only be called for inner or outside pairs! (if two particles lie in one cell they can be
    // either both inside the process or neither of them is)
    for (where_type where : {inside, outside}) {
      for (bool newton3 : {false, true}) {
        if (where == outside && not duplicatedCalculation) {
          // this case does not happen and it is expected to fail.
          continue;
        }
        testSoAGlobals(where, newton3, duplicatedCalculation, own);
      }
    }
  }
}

TEST_F(LJFunctorTest, testSoAFunctorGlobalsPair) {
  for (bool duplicatedCalculation : {false, true}) {
    for (where_type where : {inside, boundary, outside}) {
      for (bool newton3 : {false, true}) {
        testSoAGlobals(where, newton3, duplicatedCalculation, pair);
      }
    }
  }
}

TEST_F(LJFunctorTest, testAoSFunctorGlobalsOpenMPParallel) {
  bool duplicatedCalculation = false;
  bool newton3 = true;
  double multiparticlefactor = 2.;  // two particles, so factor 2
  double whereFactor = 1.;          // all inside, so factor 1
  std::string where_str = "inside";
  Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0);
  Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1);

  Molecule p3({0., 2., 0.}, {0., 0., 0.}, 0);
  Molecule p4({0.1, 2.2, 0.3}, {0., 0., 0.}, 1);

  autopas::LJFunctor<Molecule, FMCell, true> functor(cutoff, epsilon, sigma, shift, lowCorner, highCorner,
                                                     duplicatedCalculation);

  functor.resetGlobalValues();
  // This is a basic check for the global calculations, by checking the handling of two particle interactions in
  // parallel. If interactions are dangerous, archer will complain.
#if defined(AUTOPAS_OPENMP)
#pragma omp parallel
#endif
  {
#if defined(AUTOPAS_OPENMP)
#pragma omp sections
#endif
    {
#if defined(AUTOPAS_OPENMP)
#pragma omp section
#endif
      functor.AoSFunctor(p1, p2, newton3);
#if defined(AUTOPAS_OPENMP)
#pragma omp section
#endif
      functor.AoSFunctor(p3, p4, newton3);
    }
  }
  functor.postProcessGlobalValues(newton3);

  double upot = functor.getUpot();
  double virial = functor.getVirial();

  EXPECT_NEAR(upot, whereFactor * multiparticlefactor * expectedEnergy, absDelta)
      << "where: " << where_str << ", newton3: " << newton3 << ", duplicatedCalculation:" << duplicatedCalculation;
  EXPECT_NEAR(virial, whereFactor * multiparticlefactor * expectedVirial, absDelta)
      << "where: " << where_str << ", newton3: " << newton3 << ", duplicatedCalculation:" << duplicatedCalculation;
}