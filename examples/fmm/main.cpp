/**
 * @file main.cpp
 * @date 14.09.19
 * @author Joachim Marin
 */

#include <cassert>
#include <chrono>
#include <complex>
#include <iostream>
#include <random>
#include "FmmParticle.h"
#include "Math3D.h"
#include "Octree.h"
#include "Operators.h"
#include "autopas/AutoPas.h"

std::chrono::steady_clock::time_point lastTimePoint = std::chrono::steady_clock::now();

double measureTime() {
  std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
  double diff = std::chrono::duration_cast<std::chrono::milliseconds>(now - lastTimePoint).count();
  lastTimePoint = now;
  return diff;
}

void upwardPassRec(OctreeNode &node, Operators &op) {
  if (node.isLeaf()) {
    op.P2M(node);
  } else {
    for (int i = 0; i < 8; ++i) {
      auto child = node.getChild(i);
      upwardPassRec(*child, op);
    }
    op.M2M(node);
  }
}

void upwardPass(Octree &tree, Operators &op) {
  upwardPassRec(*tree.getRoot(), op);
  std::cout << "P2M and M2M took " << measureTime() << "ms" << std::endl;
}

void downwardPassRec1(OctreeNode &node, Operators &op) {
  op.M2L(node);
  if (!node.isLeaf()) {
    for (int i = 0; i < 8; ++i) {
      auto child = node.getChild(i);
      downwardPassRec1(*child, op);
    }
  }
}

void downwardPassRec2(OctreeNode &node, Operators &op) {
  op.L2L(node);
  if (!node.isLeaf()) {
    for (int i = 0; i < 8; ++i) {
      auto child = node.getChild(i);
      downwardPassRec2(*child, op);
    }
  } else {
    op.L2P(node);
  }
}

void downwardPass(Octree &tree, Operators &op) {
  downwardPassRec1(*tree.getRoot(), op);
  std::cout << "M2L took " << measureTime() << "ms" << std::endl;
  downwardPassRec2(*tree.getRoot(), op);
  std::cout << "L2L and L2P took " << measureTime() << "ms" << std::endl;
}

void randomTest(int orderOfExpansion) {
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

  double distance = Math3D::toSpherical(Math3D::subtract(particle1.getR(), particle2.getR()))[0];
  std::cout << "distance = " << distance << std::endl;
  std::cout << "exact = " << charge2 / distance << std::endl;
  std::cout << "exact = " << charge1 / distance << std::endl;

  Operators op(orderOfExpansion);
  upwardPass(tree, op);
  downwardPass(tree, op);
}

bool almostEqual(double a, double b) { return std::abs(a - b) < 0.01; }
bool almostEqual(std::complex<double> a, std::complex<double> b) {
  return std::abs(a.real() - b.real()) < 0.00001 && std::abs(a.imag() - b.imag()) < 0.00001;
}

double sphericalNormalization(int n) { return std::sqrt((2 * n + 1) / (4 * M_PI)); }

void test3DMath() {
  // Factorial
  assert(Math3D::factorial(5) == 120);
  assert(Math3D::factorial(4) == 24);

  // Double Factorial
  assert(Math3D::doubleFactorial(9) == 945);
  assert(Math3D::doubleFactorial(4) == 8);

  // Associated Legendre Polynomials
  assert(Math3D::associatedLegendrePolynomial(0, 0, 0.7) == 1);
  assert(Math3D::associatedLegendrePolynomial(0, 1, 0.4) == 0.4);
  assert(Math3D::associatedLegendrePolynomial(2, 2, 0.3) == 2.73);

  // P3,4(0.2)
  std::cout << "P3,4(0.2)" << std::endl;
  double result = Math3D::associatedLegendrePolynomial(3, 4, 0.2);
  double exact = -105.0 * 0.2 * std::pow(1.0 - 0.2 * 0.2, 3.0 / 2.0);
  std::cout << "result=" << result << std::endl;
  std::cout << "exact=" << exact << std::endl;
  assert(almostEqual(result, exact));

  // P1,3(0.6)
  std::cout << "P1,3(0.6)" << std::endl;
  result = Math3D::associatedLegendrePolynomial(1, 3, 0.6);
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
  cResult = Math3D::sphericalHarmonics(0, 0, theta, phi);
  cExact = 1.0 / 2.0 * std::sqrt(1.0 / M_PI) / sphericalNormalization(0);
  std::cout << "cResult=" << cResult << std::endl;
  std::cout << "cExact=" << cExact << std::endl;
  assert(almostEqual(cResult, cExact));

  // Y0,2(theta, phi)
  std::cout << "Y0,2(" << theta << ", " << phi << ")" << std::endl;
  cResult = Math3D::sphericalHarmonics(0, 2, theta, phi);
  cExact = 1.0 / 4.0 * std::sqrt(5.0 / M_PI) * (3 * std::pow(std::cos(theta), 2) - 1) / sphericalNormalization(2);
  std::cout << "cResult=" << cResult << std::endl;
  std::cout << "cExact=" << cExact << std::endl;
  assert(almostEqual(cResult, cExact));

  // Y2,2(theta, phi)
  std::cout << "Y2,2(" << theta << ", " << phi << ")" << std::endl;
  cResult = Math3D::sphericalHarmonics(2, 2, theta, phi);
  using namespace std::complex_literals;
  cExact = 1.0 / 4.0 * std::sqrt(7.5 / M_PI) * std::pow(std::sin(theta), 2) * std::exp(2.0 * 1i * phi) /
           sphericalNormalization(2);
  std::cout << "cResult=" << cResult << std::endl;
  std::cout << "cExact=" << cExact << std::endl;
  assert(almostEqual(cResult, cExact));

  // A2,5
  std::cout << "A2,5" << std::endl;
  result = Math3D::getA(2, 5);
  exact = -1.0 / std::sqrt(30240);
  std::cout << "result=" << result << std::endl;
  std::cout << "exact=" << exact << std::endl;
  assert(almostEqual(result, exact));

  // A7,9
  std::cout << "A7,9" << std::endl;
  result = Math3D::getA(7, 9);
  exact = -1.0 / std::sqrt(2 * Math3D::factorial(16));
  std::cout << "result=" << result << std::endl;
  std::cout << "exact=" << exact << std::endl;
  assert(almostEqual(result, exact));

  // subtract
  std::cout << "subtract" << std::endl;
  auto vec = Math3D::subtract({4, 3, 2}, {1, 4, 3});
  assert(vec[0] == 3 && vec[1] == -1 && vec[2] == -1);
}

int main(int argc, char **argv) {
  int orderOfExpansion = 8;

  if (argc == 2) {
    orderOfExpansion = static_cast<int>(std::strtol(argv[1], nullptr, 10));
  }
  std::cout << "orderOfExpansion = " << orderOfExpansion << std::endl;

  Math3D::initMath();

  // test3DMath();
  // std::cout << "test 3DMath done" << std::endl << std::endl;

  // std::cout << factorial(4) << std::endl;

  // randomTest(orderOfExpansion);

  measureTime();

  Octree tree(4, 2.0);

  auto cell000 = tree.getCell(tree.getHeight(), 0, 0, 0)->getContainer();
  auto cell111 = tree.getCell(tree.getHeight(), 1, 1, 1)->getContainer();

  auto cell333 = tree.getCell(tree.getHeight(), 3, 3, 3)->getContainer();

  // particle1
  auto particle1 = FmmParticle({0.5, 0.5, 0.5}, {0, 0, 0}, 1, 25.0);
  cell000->addParticle(particle1);
  // particle2
  auto particle2 = FmmParticle({1.5, 0.5, 1.5}, {0, 0, 0}, 1, 50.0);
  cell000->addParticle(particle2);
  // particle3
  auto particle3 = FmmParticle({2.5, 3.5, 3}, {0, 0, 0}, 1, 100.0);
  cell111->addParticle(particle3);
  // particle4
  auto particle4 = FmmParticle({7.5, 6.5, 7}, {0, 0, 0}, 1, 500.0);
  cell333->addParticle(particle4);

  std::cout << "Init took " << measureTime() << "ms" << std::endl;

  Operators op(orderOfExpansion);
  upwardPass(tree, op);
  downwardPass(tree, op);

  auto iter = cell000->begin();
  particle1 = *iter;
  ++iter;
  particle2 = *iter;
  iter = cell111->begin();
  particle3 = *iter;
  iter = cell333->begin();
  particle4 = *iter;

  std::cout << particle1.resultFMM << std::endl;
  std::cout << particle2.resultFMM << std::endl;
  std::cout << particle3.resultFMM << std::endl;
  std::cout << particle4.resultFMM << std::endl;

  assert(almostEqual(particle1.resultFMM, 44.32));
  assert(almostEqual(particle2.resultFMM, 49.45));
  assert(almostEqual(particle3.resultFMM, 70.71));
  assert(almostEqual(particle4.resultFMM, 21.3));

  /*std::cout << "fmm took " << measureTime() << "ms" << std::endl;

  int testSize = 1000000;
  Complex checksum = 0;
  for(int j = -testSize; j <= testSize; j++) {
    checksum += Math3D::powI(j) * static_cast<double>(j);
  }
  std::cout << checksum << std::endl;
  std::cout << "powI took " << measureTime() << "ms" << std::endl;
  checksum = 0;
  for(int j = -testSize; j <= testSize; j++) {
    using namespace std::complex_literals;
    checksum += std::pow(1i, j) * static_cast<double>(j);
  }
  std::cout << checksum << std::endl;
  std::cout << "std::pow took " << measureTime() << "ms" << std::endl;*/

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

  return EXIT_SUCCESS;
}