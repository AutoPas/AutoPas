/**
 * @file LJFunctorTest.cpp
 * @author seckler
 * @date 06.11.18
 */

#include "LJFunctorMixingTest.h"

std::array<double, 3> LJFunctorMixingTest::lennardForceCalculation(std::array<double, 3> x1, std::array<double, 3> x2) {
  double epsilon24 = PCL.mixing24E(0, 1);
  double ssigma = std::pow(PCL.mixingSS(0, 1), 0.5);
  std::array<double, 3> difference = autopas::ArrayMath::sub(x1, x2);
  double distance = L2Norm(difference);
  return autopas::ArrayMath::mulScalar(
      difference,
      (epsilon24) / (distance * distance) * (std::pow(ssigma / distance, 6) - 2 * std::pow(ssigma / distance, 12)));
}

void LJFunctorMixingTest::testAoSNoGlobals(bool newton3) {
  std::map<unsigned long, double> universalMap;
  for (unsigned long i = 0; i < 2; i++) {
    universalMap.emplace(i, 1.0);
  }
  ParticleClassLibrary PCL = ParticleClassLibrary(universalMap, universalMap, universalMap);
  autopas::LJFunctor<Molecule, FMCell> functor(cutoff, PCL, shift);

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

TEST_F(LJFunctorMixingTest, LennardApplication) {
  std::array<double, 3> p1({0., 0., 0.});
  std::array<double, 3> p2({0.1, 0.2, 0.3});
  std::array<double, 3> lennardO = lennardForceCalculation(p1, p2);
  std::cout.precision(17);
  std::cout << "output of calculation " << std::fixed << lennardO[0] << std::endl;
  std::cout << "output of calculation " << std::fixed << lennardO[1] << std::endl;
  std::cout << "output of calculation " << std::fixed << lennardO[2] << std::endl;
  for (auto i = 0; i < 3; i++) {
    EXPECT_NEAR(lennardO[i], expectedForce[i], absDelta);
  }
}
TEST_F(LJFunctorMixingTest, AoSNoGlobalNoN3) { testAoSNoGlobals(false); }
// TEST_F(LJFunctorMixingTest,AoSNoGlobalN3){
//    testAoSNoGlobals(true);
//}