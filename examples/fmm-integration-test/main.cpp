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

  auto cont = AutoPasCont();
  cont.setAllowedContainers(std::set<autopas::ContainerOption>{autopas::ContainerOption::linkedCells});
  cont.setBoxMin(std::array<double, 3>{0, 0, 0});
  cont.setBoxMax(std::array<double, 3>{2, 3, 4});
  cont.init();

  for (int i = 0; i < 10; ++i) {
    auto r = RandomGenerator::randomPosition(cont.getBoxMin(), cont.getBoxMax());
    auto v = std::array<double, 3>({0, 0, 0});
    auto part = autopas::Particle(r, v, i);
    cont.addParticle(part);
  }

  std::cout << "Number of Particles: " << cont.getNumberOfParticles() << std::endl;

  auto fmmTree = cont.getFastMultipoleTree();

  std::cout << std::flush;
}
