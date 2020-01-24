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
  long numberOfParticles = 1000;
  bool checkResults = false;

  if (argc == 5) {
    orderOfExpansion = std::stol(argv[1]);
    domainSize = (std::stod(argv[2]) + 0.5);
    numberOfParticles = std::stol(argv[3]);
    checkResults = std::stol(argv[4]) > 0;
  }

  std::array<double, 3> boxMin = {0, 0, 0};
  std::array<double, 3> boxMax = {domainSize, domainSize, domainSize};

  // Print parameters.
  std::cout << "orderOfExpansion = " << orderOfExpansion << std::endl;
  std::cout << "domainSize = " << domainSize << std::endl;
  std::cout << "numberOfParticles = " << numberOfParticles << std::endl;
  std::cout << "boxMin = " << autopas::utils::ArrayUtils::to_string(boxMin) << std::endl;
  std::cout << "boxMax = " << autopas::utils::ArrayUtils::to_string(boxMax) << std::endl;

  autopas::utils::Timer timer;
  timer.start();

  AutoPasCont cont;
  cont.setAllowedContainers({autopas::ContainerOption::linkedCells});
  cont.setBoxMin(boxMin);
  cont.setBoxMax(boxMax);
  cont.setCutoff(1.0);
  cont.setVerletSkin(0.0);
  cont.init();

  RandomGenerator::fillWithParticles(cont, autopas::fmm::FmmParticle(), cont.getBoxMin(), cont.getBoxMax(),
                                     numberOfParticles);

  std::random_device rd;
  std::default_random_engine randomEngine(42);
  std::uniform_real_distribution<double> random(0.0, 1.0);

  for (auto particle = cont.begin(); particle.isValid(); ++particle) {
    particle->charge = random(randomEngine) * 10.0;
  }

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

  if (checkResults) {
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

    double maxRelError = 0;
    double avgRelError = 0;
    double maxAbsError = 0;
    double avgAbsError = 0;

    // bool isCorrect = true;
    for (auto particle = cont.begin(); particle.isValid(); ++particle) {
      double absoluteError = std::abs(particle->resultFMM - particle->resultExact);
      double relativeError = absoluteError / std::abs(particle->resultExact);

      avgAbsError += absoluteError;
      avgRelError += relativeError;

      maxAbsError = std::max(maxAbsError, absoluteError);
      maxRelError = std::max(maxRelError, relativeError);
    }

    avgAbsError /= static_cast<double>(numberOfParticles);
    avgRelError /= static_cast<double>(numberOfParticles);

    std::cout << maxAbsError << "\t" << avgAbsError << "\t" << maxRelError << "\t" << avgRelError << "\t" << std::endl;
  }
}
