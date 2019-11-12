/**
 * @file main.cpp
 * @author Joachim Marin
 * @date 5.11.2019
 */

#include <gtest/gtest.h>
#include <array>
#include <iostream>
#include "autopas/AutoPas.h"
#include "autopas/particles/Particle.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"

using AutoPasCont = autopas::AutoPas<autopas::Particle, autopas::FullParticleCell<autopas::Particle>>;

int main(int argc, char **argv) {
  std::cout << std::endl;

  std::cout << "Init Test" << std::endl;

  AutoPasCont cont;
  cont.setAllowedContainers({autopas::ContainerOption::linkedCells});
  cont.setBoxMin({0, 0, 0});
  cont.setBoxMax({2, 3, 4});
  cont.init();

  RandomGenerator::fillWithParticles(cont, autopas::Particle(), 10);

  std::cout << "Number of Particles: " << cont.getNumberOfParticles() << std::endl;

  auto fmmTree = cont.getFastMultipoleTree();

  std::cout << std::flush;
}
