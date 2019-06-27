//
// Created by nicola on 22.05.19.
//

#include <gtest/gtest.h>
#include <math.h>
#include "../../../../examples/md-flexible/MDFlexParser.h"
#include "../../../../examples/md-flexible/PrintableMolecule.h"
#include "../../../../examples/md-flexible/TimeDiscretization.h"
#include "../../../../src/autopas/utils/ArrayMath.h"
#include "../../testingHelpers/GridGenerator.h"
#include "../../testingHelpers/RandomGenerator.h"
#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"

using namespace std;
using namespace autopas;

string arrayString(array<double, 3> input) {
  std::ostringstream os;
  for (double i : input) {
    os << i;
    os << " _ ";
  }
  std::string str(os.str());
  return str;
}

double L2Norm(std::array<double, 3> array) {
  double square_sum = 0;
  for (unsigned int i = 0; i < array.size(); i++) {
    square_sum += (array[i] * array[i]);
  }
  return sqrt(square_sum);
}

std::array<double, 3> lennardForceCalculation(std::array<double, 3> x1, std::array<double, 3> x2) {
  std::array<double, 3> difference = ArrayMath::sub(x1, x2);
  double distance = L2Norm(difference);
  double epsilon = 5;
  double sigma = 1;
  return ArrayMath::mulScalar(difference, (24 * epsilon) / (distance * distance) *
                                              (std::pow(sigma / distance, 6) - 2 * std::pow(sigma / distance, 12)));
}

TEST(TimeDTest, LennardJonesCalculation) {
  // dient zum vergleich zu den Forces die im Functor ausgerechnet werden
  std::array<double, 3> x1 = {1, 1, 1};
  std::array<double, 3> x2 = {1.5, 1.5, 1.5};
  cout << "lennardJonesFunction Calculation:  " << arrayString(lennardForceCalculation(x1, x2)) << endl;
  ASSERT_TRUE(true);
}

// Testet und visualisiert die Kräfte berechnungen und TimeDiscreatization Classe
TEST(TimeDTest, GeneralForceTest) {
  // auto container = autopas::LinkedCells<PrintableMolecule,ParticleCell<PrintableMolecule>>({0., 0., 0.},{10.,
  // 10., 10.},1);
  PrintableMolecule::setEpsilon(5.0);
  PrintableMolecule::setSigma(1.0);
  PrintableMolecule::setMass(1.0);
  auto *autoPas = new autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>(std::cout);
  double epsilon = 5.0;
  double sigma = 1.0;
  double cutoff = 2;
  array<double, 3> boxmin = {0., 0., 0.};
  array<double, 3> boxmax = {5., 5., 5.};
  PrintableMolecule::setEpsilon(epsilon);
  PrintableMolecule::setSigma(sigma);
  PrintableMolecule::setMass(1.0);
  autoPas->setBoxMax(boxmax);
  autoPas->setCutoff(cutoff);
  // erstmal auf linked cells testen
  autoPas->setAllowedContainers({autopas::ContainerOption::linkedCells});
  autoPas->init();
  PrintableMolecule p1({1., 1., 1.}, {0.5, 0.5, 0.5}, 0);
  autoPas->addParticle(p1);
  PrintableMolecule p2({1.5, 1.5, 1.5}, {0., 0.5, 0.}, 1);
  autoPas->addParticle(p2);
  //  for (auto iter = autoPas->getContainer()->begin(); iter.isValid(); ++iter) {
  //    cout << iter->toString() << endl;
  //    cout << "ParticleOldF= " << arrayString(iter->getOldf()) << endl;
  //  }
  double particleD = 0.01;
  int iterations = 0;
  // iterationen beginnend
  TimeDiscretization<decltype(autoPas)> td(particleD);
  auto *functor =
      new autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>, autopas::FunctorN3Modes::Both,
                             true>(cutoff, epsilon, sigma, 0.0, boxmin, boxmax, true);
  // domain vorbeireiten: -Force initialisieren
  autoPas->iteratePairwise(functor);
  // Dokumentation prints
  //  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  //  cout << "-----AFTER INITIALIZATION----" << endl;
  //  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  //  for (auto iter = autoPas->getContainer()->begin(); iter.isValid(); ++iter) {
  //    cout << iter->toString() << "  __END" << endl;
  //    cout << "ParticleOldF= " << arrayString(iter->getOldf()) << endl;
  //  }
  //  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  //  cout << "-------ITERATIONS START------" << endl;
  //  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;

  while (iterations < 10) {
    td.VSCalculateX(autoPas);
    autoPas->iteratePairwise(functor);
    td.VSCalculateV(autoPas);
    iterations++;
    //    for (auto iter = autoPas->getContainer()->begin(); iter.isValid(); ++iter) {
    //      cout << iter->toString() << endl;
    //      cout << "ParticleOldF= " << arrayString(iter->getOldf()) << endl;
    //    }
    //    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    //    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    //    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  }
  ASSERT_TRUE(true);
  delete functor;
  delete autoPas;
}

TEST(TimeDTest, CalcX) {
  PrintableMolecule::setEpsilon(5.0);
  PrintableMolecule::setSigma(1.0);
  PrintableMolecule::setMass(1.0);
  auto *autoPas = new autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>(std::cout);
  double epsilon = 5.0;
  double sigma = 1.0;
  double cutoff = 2;
  array<double, 3> boxmin = {0., 0., 0.};
  array<double, 3> boxmax = {5., 5., 5.};
  PrintableMolecule::setEpsilon(epsilon);
  PrintableMolecule::setSigma(sigma);
  PrintableMolecule::setMass(1.0);
  autoPas->setBoxMax(boxmax);
  autoPas->setCutoff(cutoff);
  // erstmal auf linked cells testen
  autoPas->setAllowedContainers({autopas::ContainerOption::linkedCells});
  autoPas->init();
  PrintableMolecule p1({1., 1., 1.}, {0.5, 0.5, 0.5}, 0);
  autoPas->addParticle(p1);
  PrintableMolecule p2({1.5, 1.5, 1.5}, {0., 0.5, 0.}, 1);
  autoPas->addParticle(p2);
  //  for (auto iter = autoPas->getContainer()->begin(); iter.isValid(); ++iter) {
  //    cout << iter->toString() << endl;
  //    cout << "ParticleOldF= " << arrayString(iter->getOldf()) << endl;
  //  }
  double particleD = 0.01;
  int iterations = 0;
  // iterationen beginnend
  TimeDiscretization<decltype(autoPas)> td(particleD);
  auto *functor =
      new autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>, autopas::FunctorN3Modes::Both,
                             true>(cutoff, epsilon, sigma, 0.0, boxmin, boxmax, true);
  // domain vorbeireiten: -Force initialisieren
  autoPas->iteratePairwise(functor);
  cout << "delta_t value =  " << particleD << endl;
  //  while (iterations < 10) {
  //    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  //    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  //    for (auto iter = autoPas->getContainer()->begin(); iter.isValid(); ++iter) {
  //      auto v = iter->getV();
  //      auto m = iter->getMass();
  //      auto f = iter->getF();
  //      iter->setOldf(f);
  //      cout << "Particle ID: " << iter->getID() << endl;
  //      cout << "initial Velocity: " << arrayString(v) << endl;
  //      v = autopas::ArrayMath::mulScalar(v, particleD);
  //      cout << "Velocity * delta_T= " << arrayString(v) << endl;
  //      cout << "initial F = " << arrayString(f) << endl;
  //      f = autopas::ArrayMath::mulScalar(f, (particleD * particleD / (2 * m)));
  //      cout << "F * delta² / 2*m = " << arrayString(f) << endl;
  //      cout << "Print old Positions:" << arrayString(iter->getR()) << endl;
  //      auto newR = autopas::ArrayMath::add(v, f);
  //      iter->addR(newR);
  //      cout << "Print new Positions: " << arrayString(iter->getR()) << endl;
  //      cout << endl;
  //    }
  //    iterations++;
  //  }
  ASSERT_TRUE(true);
}

TEST(TimeDTEst, CalcV) {
  PrintableMolecule::setEpsilon(5.0);
  PrintableMolecule::setSigma(1.0);
  PrintableMolecule::setMass(1.0);
  auto *autoPas = new autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>(std::cout);
  double epsilon = 5.0;
  double sigma = 1.0;
  double cutoff = 2;
  array<double, 3> boxmin = {0., 0., 0.};
  array<double, 3> boxmax = {5., 5., 5.};
  PrintableMolecule::setEpsilon(epsilon);
  PrintableMolecule::setSigma(sigma);
  PrintableMolecule::setMass(1.0);
  autoPas->setBoxMax(boxmax);
  autoPas->setCutoff(cutoff);
  // erstmal auf linked cells testen
  autoPas->setAllowedContainers({autopas::ContainerOption::linkedCells});
  autoPas->init();
  PrintableMolecule p1({1., 1., 1.}, {0.5, 0.5, 0.5}, 0);
  autoPas->addParticle(p1);
  PrintableMolecule p2({1.5, 1.5, 1.5}, {0., 0.5, 0.}, 1);
  autoPas->addParticle(p2);
  //  for (auto iter = autoPas->getContainer()->begin(); iter.isValid(); ++iter) {
  //    cout << iter->toString() << endl;
  //    cout << "ParticleOldF= " << arrayString(iter->getOldf()) << endl;
  //  }
  double particleD = 0.01;
  int iterations = 0;
  // iterationen beginnend
  TimeDiscretization<decltype(autoPas)> td(particleD);
  auto *functor =
      new autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>, autopas::FunctorN3Modes::Both,
                             true>(cutoff, epsilon, sigma, 0.0, boxmin, boxmax, true);
  // domain vorbeireiten: -Force initialisieren
  autoPas->iteratePairwise(functor);
  //  cout << "delta_t value =  " << particleD << endl;
  //  while (iterations < 10) {
  //    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  //    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  //
  //    for (auto iter = autoPas->getContainer()->begin(); iter.isValid(); ++iter) {
  //      auto m = iter->getMass();
  //      auto force = iter->getF();
  //      auto old_force = iter->getOldf();
  //      cout << "Particle ID: " << iter->getID() << endl;
  //      cout << "Old forces: " << arrayString(old_force) << endl;
  //      cout << "Current forces: " << arrayString(force) << endl;
  //      auto addedF = autopas::ArrayMath::add(force, old_force);
  //      cout << "OldF + Force =  " << arrayString(addedF) << endl;
  //      auto newV = autopas::ArrayMath::mulScalar(addedF, particleD / (2 * 1));
  //      cout << "Multiplied by delta_t and 2*m:" << endl << arrayString(newV) << endl;
  //      cout << "old Velocity= " << arrayString(iter->getV()) << endl;
  //      iter->addV(newV);
  //      cout << "new Velocity " << arrayString(iter->getV()) << endl;
  //      cout << endl;
  //    }
  //    iterations++;
  //  }
  ASSERT_TRUE(true);
}
