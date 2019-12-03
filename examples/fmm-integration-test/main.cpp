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
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/FmmMath.h"
#include "autopas/utils/Timer.h"

using AutoPasCont = autopas::AutoPas<autopas::fmm::FmmParticle, autopas::FullParticleCell<autopas::fmm::FmmParticle>>;

int main(int argc, char **argv) {
  long orderOfExpansion = 5;
  double domainSize = 5;
  double particleDensity = 1;

  if (argc == 4) {
    orderOfExpansion = std::stol(argv[1]);
    domainSize = std::stod(argv[2]);
    particleDensity = std::stod(argv[3]);
  }

  long numberOfParticles = std::lround(domainSize * domainSize * domainSize * particleDensity);

  std::array<double, 3> boxMin = {0, 0, 0};
  std::array<double, 3> boxMax = {domainSize, domainSize, domainSize};

  // Print parameters.
  std::cout << "orderOfExpansion = " << orderOfExpansion << std::endl;
  std::cout << "domainSize = " << domainSize << std::endl;
  std::cout << "particleDensity = " << particleDensity << std::endl;
  std::cout << "numberOfParticles = " << numberOfParticles << std::endl;
  std::cout << "boxMin = " << autopas::utils::ArrayUtils::to_string(boxMin) << std::endl;
  std::cout << "boxMax = " << autopas::utils::ArrayUtils::to_string(boxMax) << std::endl;

  autopas::utils::Timer timer;
  timer.start();

  AutoPasCont cont;
  cont.setAllowedContainers({autopas::ContainerOption::linkedCells});
  cont.setBoxMin(boxMin);
  cont.setBoxMax(boxMax);
  cont.init();

  RandomGenerator::fillWithParticles(cont, autopas::fmm::FmmParticle(), cont.getBoxMin(), cont.getBoxMax(),
                                     numberOfParticles);

  autopas::fmm::PotentialOperators<autopas::fmm::FmmParticle, autopas::FullParticleCell<autopas::fmm::FmmParticle>> op;
  long initTime = timer.stop();
  std::cout << "Init: " << initTime << "us" << std::endl;
  timer.start();

  auto fmmTree = cont.getFastMultipoleTree();

  long fmmTreeTime = timer.stop();
  std::cout << "getFmmTree: " << fmmTreeTime << "us" << std::endl;

  std::cout << "Run Fmm..." << std::endl;
  timer.start();

  op.runFmm(*fmmTree, orderOfExpansion, cont);

  /*for (auto particle = cont.begin(); particle.isValid(); ++particle) {
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
    double error = std::abs(particle->resultFMM - particle->resultExact);
    if (error > 0.01) {
      isCorrect = false;

      std::cout << "[ID=" << particle->getID() << "] " << particle->getR()[0] << ", " << particle->getR()[1] << ", "
                << particle->getR()[2] << ", charge = " << particle->charge << std::endl;
      std::cout << "long range " << particle->longRange << std::endl;
      std::cout << "short range " << particle->shortRange << std::endl;
      std::cout << "resultFMM " << particle->resultFMM << std::endl;
      std::cout << "resultExact " << particle->resultExact << std::endl;
    }
  }
  std::cout << std::flush;
  if (not isCorrect) {
    std::cerr << "wrong result" << std::endl;
  }*/

  long runFmmTime = timer.stop();
  std::cout << "runFmm: " << runFmmTime << "us" << std::endl;

  long totalTime = initTime + fmmTreeTime + runFmmTime;

  double initPercentage = 100 * static_cast<double>(initTime) / static_cast<double>(totalTime);
  double treePercentage = 100 * static_cast<double>(fmmTreeTime) / static_cast<double>(totalTime);
  double fmmPercentage = 100 * static_cast<double>(runFmmTime) / static_cast<double>(totalTime);

  std::cout << "Init: " << initPercentage << "%" << std::endl;
  std::cout << "Tree: " << treePercentage << "%" << std::endl;
  std::cout << " Fmm: " << fmmPercentage << "%" << std::endl;

  std::cout << std::endl;
  std::cout << "totalTime (us) (ms)" << std::endl;

  std::cout << totalTime << std::endl;
  std::cout << totalTime / 1000 << std::endl;
}
