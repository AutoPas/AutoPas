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
  cont.setBoxMax({2, 3, 4});
  cont.init();

  RandomGenerator::fillWithParticles(cont, autopas::fmm::FmmParticle(), 10);

  std::cout << "Number of Particles: " << cont.getNumberOfParticles() << std::endl;

  auto fmmTree = cont.getFastMultipoleTree();

  autopas::utils::FmmMath<double, long>::initialize();
  autopas::utils::FmmMath<double, long> fmmMath;
  std::cout << "Factorial(5) = " << autopas::utils::FmmMath<double, long>::factorial(5) << std::endl;

  autopas::fmm::PotentialOperators<autopas::fmm::FmmParticle, autopas::FullParticleCell<autopas::fmm::FmmParticle>> op(
      8);
  op.RunFmm(*fmmTree, cont);

  std::cout << std::flush;
}
