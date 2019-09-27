/**
 * @file main.cpp
 * @date 14.09.19
 * @author Joachim Marin
 */

#include <cassert>
#include <complex>
#include <iostream>
#include <random>
#include "FmmParticle.h"
#include "Math3D.h"
#include "Octree.h"
#include "Operators.h"
#include "autopas/AutoPas.h"

void upwardPassRec(OctreeNode *node) {
  if (node->isLeaf()) {
    Operators::P2M(node);
  } else {
    for (int i = 0; i < 8; ++i) {
      auto child = node->getChild(i);
      upwardPassRec(child);
    }
    Operators::M2M(node);
  }
}

void upwardPass(Octree *tree) { upwardPassRec(tree->getRoot()); }

void downwardPassRec1(OctreeNode *node) {
  Operators::M2L(node);
  if (!node->isLeaf()) {
    for (int i = 0; i < 8; ++i) {
      auto child = node->getChild(i);
      downwardPassRec1(child);
    }
  }
}

void downwardPassRec2(OctreeNode *node) {
  Operators::L2L(node);
  if (!node->isLeaf()) {
    for (int i = 0; i < 8; ++i) {
      auto child = node->getChild(i);
      downwardPassRec2(child);
    }
  } else {
    Operators::L2P(node);
  }
}

void downwardPass(Octree *tree) {
  downwardPassRec1(tree->getRoot());
  downwardPassRec2(tree->getRoot());
}

void randomTest() {
  Octree tree(4, 2.0);

  auto vec = tree.getCell(tree.getHeight(), 2, 1, 3)->getCenter();
  assert(vec[0] == 5.0 && vec[1] == 3.0 && vec[2] == 7.0);

  double charge1 = 10;
  double charge2 = 20;

  auto corner1 = tree.getCell(tree.getHeight(), 0, 0, 0)->getContainer();
  auto corner2 = tree.getCell(tree.getHeight(), 3, 3, 3)->getContainer();

  std::random_device rd;
  std::mt19937 randomEngine(rd());
  std::uniform_real_distribution<double> random(0.0, 2.0);

  double x1, y1, z1;
  double x2, y2, z2;
  // x1 = y1 = z1 = 0.5;

  x1 = random(randomEngine);
  y1 = random(randomEngine);
  z1 = random(randomEngine);

  x2 = 6.0 + random(randomEngine);
  y2 = 6.0 + random(randomEngine);
  z2 = 6.0 + random(randomEngine);
  // x2 = y2 = z2 = 7.5;
  // x2 = 7;
  // y2 = 6.5;
  // z2 = 7.46;

  FmmParticle particle1({x1, y1, z1}, {0, 0, 0}, 0, charge1);
  FmmParticle particle2({x2, y2, z2}, {0, 0, 0}, 0, charge2);
  std::cout << "particle1.getR() = (" << particle1.getR()[0] << "," << particle1.getR()[1] << "," << particle1.getR()[2]
            << ")" << std::endl;
  std::cout << "particle2.getR() = (" << particle2.getR()[0] << "," << particle2.getR()[1] << "," << particle2.getR()[2]
            << ")" << std::endl;

  corner1->addParticle(particle1);
  corner2->addParticle(particle2);

  double distance = toSpherical(subtract(particle1.getR(), particle2.getR()))[0];
  std::cout << "distance = " << distance << std::endl;
  std::cout << "exact = " << charge2 / distance << std::endl;
  std::cout << "exact = " << charge1 / distance << std::endl;

  upwardPass(&tree);
  downwardPass(&tree);
}

bool almostEqual(double a, double b) { return std::abs(a - b) < 0.00001; }
bool almostEqual(std::complex<double> a, std::complex<double> b) {
  return std::abs(a.real() - b.real()) < 0.00001 && std::abs(a.imag() - b.imag()) < 0.00001;
}

double sphericalNormalization(int n) { return std::sqrt((2 * n + 1) / (4 * M_PI)); }

void test3DMath() {
  // Factorial
  assert(factorial(5) == 120);
  assert(factorial(4) == 24);

  // Double Factorial
  assert(doubleFactorial(9) == 945);
  assert(doubleFactorial(4) == 8);

  // Associated Legendre Polynomials
  assert(associatedLegendrePolynomial(0, 0, 0.7) == 1);
  assert(associatedLegendrePolynomial(0, 1, 0.4) == 0.4);
  assert(associatedLegendrePolynomial(2, 2, 0.3) == 2.73);

  // P3,4(0.2)
  std::cout << "P3,4(0.2)" << std::endl;
  double result = associatedLegendrePolynomial(3, 4, 0.2);
  double exact = -105.0 * 0.2 * std::pow(1.0 - 0.2 * 0.2, 3.0 / 2.0);
  std::cout << "result=" << result << std::endl;
  std::cout << "exact=" << exact << std::endl;
  assert(almostEqual(result, exact));

  // P1,3(0.6)
  std::cout << "P1,3(0.6)" << std::endl;
  result = associatedLegendrePolynomial(1, 3, 0.6);
  exact = -3.0 / 2.0 * (5.0 * 0.6 * 0.6 - 1.0) * std::pow(1 - 0.6 * 0.6, 0.5);
  std::cout << "result=" << result << std::endl;
  std::cout << "exact=" << exact << std::endl;
  assert(almostEqual(result, exact));

  std::complex<double> cResult;
  std::complex<double> cExact;

  std::random_device rd;
  std::mt19937 randomEngine(rd());
  std::uniform_real_distribution<double> random(0.0, M_2_PI);

  double theta = random(randomEngine);
  double phi = random(randomEngine);

  // Y0,0(theta, phi)
  std::cout << "Y0,0(" << theta << ", " << phi << ")" << std::endl;
  cResult = sphericalHarmonics(0, 0, theta, phi);
  cExact = 1.0 / 2.0 * std::sqrt(1.0 / M_PI) / sphericalNormalization(0);
  std::cout << "cResult=" << cResult << std::endl;
  std::cout << "cExact=" << cExact << std::endl;
  assert(almostEqual(cResult, cExact));

  // Y0,2(theta, phi)
  std::cout << "Y0,2(" << theta << ", " << phi << ")" << std::endl;
  cResult = sphericalHarmonics(0, 2, theta, phi);
  cExact = 1.0 / 4.0 * std::sqrt(5.0 / M_PI) * (3 * std::pow(std::cos(theta), 2) - 1) / sphericalNormalization(2);
  std::cout << "cResult=" << cResult << std::endl;
  std::cout << "cExact=" << cExact << std::endl;
  assert(almostEqual(cResult, cExact));

  // Y2,2(theta, phi)
  std::cout << "Y2,2(" << theta << ", " << phi << ")" << std::endl;
  cResult = sphericalHarmonics(2, 2, theta, phi);
  using namespace std::complex_literals;
  cExact = 1.0 / 4.0 * std::sqrt(7.5 / M_PI) * std::pow(std::sin(theta), 2) * std::exp(2.0 * 1i * phi) /
           sphericalNormalization(2);
  std::cout << "cResult=" << cResult << std::endl;
  std::cout << "cExact=" << cExact << std::endl;
  assert(almostEqual(cResult, cExact));

  // A2,5
  std::cout << "A2,5" << std::endl;
  result = getA(2, 5);
  exact = -1.0 / std::sqrt(30240);
  std::cout << "result=" << result << std::endl;
  std::cout << "exact=" << exact << std::endl;
  assert(almostEqual(result, exact));

  // A7,9
  std::cout << "A7,9" << std::endl;
  result = getA(7, 9);
  exact = -1.0 / std::sqrt(2 * factorial(16));
  std::cout << "result=" << result << std::endl;
  std::cout << "exact=" << exact << std::endl;
  assert(almostEqual(result, exact));

  // subtract
  std::cout << "subtract" << std::endl;
  auto vec = subtract({4, 3, 2}, {1, 4, 3});
  assert(vec[0] == 3 && vec[1] == -1 && vec[2] == -1);
}

int main(int argc, char **argv) {
  std::cout << "test" << std::endl;

  initMath();

  test3DMath();
  std::cout << "test 3DMath done" << std::endl << std::endl;

  // std::cout << factorial(4) << std::endl;

  // randomTest();

  Octree tree(4, 2.0);

  auto cell000 = tree.getCell(tree.getHeight(), 0, 0, 0)->getContainer();
  auto cell111 = tree.getCell(tree.getHeight(), 1, 1, 1)->getContainer();

  auto cell333 = tree.getCell(tree.getHeight(), 3, 3, 3)->getContainer();

  FmmParticle particle;

  // particle1
  particle = FmmParticle({0.5, 0.5, 0.5}, {0, 0, 0}, 1, 25.0);
  cell000->addParticle(particle);
  // particle2
  particle = FmmParticle({1.5, 0.5, 1.5}, {0, 0, 0}, 1, 50.0);
  cell000->addParticle(particle);
  // particle3
  particle = FmmParticle({2.5, 3.5, 3}, {0, 0, 0}, 1, 100.0);
  cell111->addParticle(particle);
  // particle4
  particle = FmmParticle({7.5, 6.5, 7}, {0, 0, 0}, 1, 500.0);
  cell333->addParticle(particle);

  // particle1:
  // interacts with 4
  // dist = sqrt(127.25) = 11.28
  // exact = 44.32

  // particle2:
  // interacts with 4
  // dist = sqrt(102.25) = 10.11
  // exact = 49.45

  // particle3:
  // interacts with 4
  // dist = sqrt(50) = 7.07
  // exact = 70.71

  // particle4:
  // interacts with 1,2,3
  // dist = sqrt(127.25) = 11.28
  // exact = 2.22
  // dist = sqrt(102.25) = 10.11
  // exact = 4.94
  // dist = sqrt(50) = 7.07
  // exact = 14.14
  // sum = 21.3

  upwardPass(&tree);
  downwardPass(&tree);

  // Printed Result:
  // result = 4.04145
  // result = 4.15036
  // result = 5.06055
  // result = 6.06218
  // result = 150.369

  // Expected Result:
  // (1/1/1):
  // p = 4.04       good
  //
  // (0.5/0.5/0.5)
  // p = 3.73       off by 0.42
  //
  // (1.5/1.5/1.5)
  // p = 4.41       off by 0.65
  //
  // (3,3,3)
  // p = 6.06       good
  //
  // (7,7,7)
  // p += 1.64
  // p += 2.22
  // p += 2.62
  // p += 144.30
  // p = 150.78     off by 0.41

  return EXIT_SUCCESS;
}