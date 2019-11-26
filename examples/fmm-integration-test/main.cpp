/**
 * @file main.cpp
 * @author Joachim Marin
 * @date 5.11.2019
 */

#include <iostream>
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "autopas/AutoPas.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/fastMultipoleMethod/FmmOperatorInterface.h"
#include "autopas/fastMultipoleMethod/PotentialOperators.h"
#include "autopas/particles/FmmParticle.h"
#include "autopas/utils/FmmMath.h"

using AutoPasCont = autopas::AutoPas<autopas::fmm::FmmParticle, autopas::FullParticleCell<autopas::fmm::FmmParticle>>;

int main(int argc, char **argv) {
  std::cout << std::endl;

  std::cout << "Init Test" << std::endl;

  AutoPasCont cont;
  cont.setAllowedContainers({autopas::ContainerOption::linkedCells});
  cont.setBoxMin({0, 0, 0});
  cont.setBoxMax({5, 4, 7});
  cont.init();

  RandomGenerator::fillWithParticles(cont, autopas::fmm::FmmParticle(), 30);
  autopas::fmm::FmmParticle part;
  /*part = autopas::fmm::FmmParticle({0.1, 0.1, 0.1}, {0, 0, 0},0, 1);
  cont.addParticle(part);
  part = autopas::fmm::FmmParticle({4.1, 0.1, 0.1}, {0, 0, 0},1, 1);
  cont.addParticle(part);

  part = autopas::fmm::FmmParticle({2, 0.1, 0.1}, {0, 0, 0},2, 1);
  cont.addParticle(part);*/

  std::cout << "Number of Particles: " << cont.getNumberOfParticles() << std::endl;

  auto fmmTree = cont.getFastMultipoleTree();

  autopas::utils::FmmMath<double, long>::initialize();
  autopas::utils::FmmMath<double, long> fmmMath;
  std::cout << "Factorial(5) = " << autopas::utils::FmmMath<double, long>::factorial(5) << std::endl;

  autopas::fmm::PotentialOperators<autopas::fmm::FmmParticle, autopas::FullParticleCell<autopas::fmm::FmmParticle>> op(
      5);
  op.RunFmm(*fmmTree, 5, cont);

  for (auto particle = cont.begin(); particle.isValid(); ++particle) {
    for (auto otherParticle = cont.begin(); otherParticle.isValid(); ++otherParticle) {
      if (particle->getID() != otherParticle->getID()) {
        double x = particle->getR()[0] - otherParticle->getR()[0];
        double y = particle->getR()[1] - otherParticle->getR()[1];
        double z = particle->getR()[2] - otherParticle->getR()[2];
        auto dist = 1.0 / std::sqrt(x * x + y * y + z * z);
        particle->resultExact += otherParticle->charge * dist;
      }
    }
  }

  bool isCorrect = true;
  for (auto particle = cont.begin(); particle.isValid(); ++particle) {
    std::cout << "[ID=" << particle->getID() << "] " << particle->getR()[0] << ", " << particle->getR()[1] << ", "
              << particle->getR()[2] << ", charge = " << particle->charge << std::endl;
    std::cout << "long range " << particle->longRange << std::endl;
    std::cout << "short range " << particle->shortRange << std::endl;
    std::cout << "resultFMM " << particle->resultFMM << std::endl;
    std::cout << "resultExact " << particle->resultExact << std::endl;
    double error = std::abs(particle->resultFMM - particle->resultExact);
    if (error > 0.01) {
      isCorrect = false;
    }
  }
  if (not isCorrect) {
    std::cerr << "wrong result" << std::endl;
  }

  std::cout << std::flush;
}
