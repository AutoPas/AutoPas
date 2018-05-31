//
// Created by seckler on 19.01.18.
//

#include <array>
#include <iostream>
#include "autopasIncludes.h"
#include "sph/autopassph.h"
#include "utils/Timer.h"

template <class Container>
void measureContainer(Container *cont, int numParticles, int numIterations);

double fRand(double fMin, double fMax) {
  double f = static_cast<double>(rand()) / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

std::array<double, 3> randomPosition(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax) {
  std::array<double, 3> r{0, 0, 0};
  for (int d = 0; d < 3; ++d) {
    r[d] = fRand(boxMin[d], boxMax[d]);
  }
  return r;
}

void addParticles(
    autopas::LinkedCells<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> &sph_system,
    int numParticles) {
  // Place SPH particles

  srand(42);  // fixed seedpoint

  std::array<double, 3> boxMin(sph_system.getBoxMin()), boxMax(sph_system.getBoxMax());

  for (int i = 0; i < numParticles; ++i) {
    auto id = static_cast<unsigned long>(i);
    autopas::sph::SPHParticle particle(randomPosition(boxMin, boxMax), {0., 0., 0.}, id, 0.75, 0.012, 0.);
    // autopas::sph::SPHParticle ith(randomPosition(boxMin, boxMax), {0, 0, 0},
    // i++, 0.75, 0.012, 0. );
    sph_system.addParticle(particle);
  }

  //	for (auto it = cont->begin(); it.isValid(); ++it) {
  //		it->print();
  //	}

  for (auto part = sph_system.begin(); part.isValid(); ++part) {
    part->setMass(part->getMass() * boxMax[0] * boxMax[1] * boxMax[2] / (double)(numParticles));
  }
  // std::cout << "# of ptcls is... " << numParticles << std::endl;

  // Set the end time
}

int main(int argc, char *argv[]) {
  std::array<double, 3> boxMin({0., 0., 0.}), boxMax{};
  boxMax[0] = 0.15;
  boxMax[1] = boxMax[2] = boxMax[0] / 1.0;
  double cutoff = .03;

  autopas::LinkedCells<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> lcCont(
      boxMin, boxMax, cutoff);

  autopas::DirectSum<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> dirCont(
      boxMin, boxMax, cutoff);

  int numParticles = 16;
  int numIterations = 100000;
  int whichContainer = 0;
  int whichFunctor = 0;
  if (argc == 5) {
    numParticles = atoi(argv[1]);
    numIterations = atoi(argv[2]);
    whichContainer = atoi(argv[3]);
    whichFunctor = atoi(argv[4]);
  } else {
    exit(1);
  }

  addParticles(lcCont, numParticles);

  for (auto it = lcCont.begin(); it.isValid(); ++it) {
    dirCont.addParticle(*it);
  }

  if (whichContainer == 0) {
    measureContainer(&lcCont, numParticles, numIterations);
  } else if (whichContainer == 1) {
    measureContainer(&dirCont, numParticles, numIterations);
  }
}

template <class Container>
void measureContainer(Container *cont, int numParticles, int numIterations) {
  autopas::sph::SPHCalcDensityFunctor func;
  autopas::FlopCounterFunctor<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>>
      flopFunctor(cont->getCutoff());

  autopas::utils::Timer t;

  cont->iteratePairwiseAoS2(&flopFunctor);
  double flopsPerIteration = flopFunctor.getFlops(func.getNumFlopsPerKernelCall());

  t.start();
  for (int i = 0; i < numIterations; ++i) {
    cont->iteratePairwiseAoS2(&func);
  }
  double elapsedTime = t.stop();

  double flops = flopsPerIteration * numIterations;

  double MFUPS = numParticles * numIterations / elapsedTime * 1e-6;

  std::cout << numParticles << "\t" << numIterations << "\t" << elapsedTime / numIterations << "\t" << MFUPS;
  std::cout << "\t" << flops;
  std::cout << "\t" << flopFunctor.getHitRate();
  std::cout << "\t" << flops / elapsedTime * 1e-9 << std::endl;

  // std::cout << "measuring done" << std::endl;
}