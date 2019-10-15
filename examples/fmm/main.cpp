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
#include "NearFieldFunctor.h"
//#include "Octree.h"
#include "AdaptiveOctree.h"
#include "Operators.h"
#include "autopas/AutoPas.h"

std::chrono::steady_clock::time_point lastTimePoint = std::chrono::steady_clock::now();
int particleId = 42;

unsigned long long nearFieldCalculations = 0;
unsigned long long exactCalculations = 0;

double measureTime() {
  std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
  double diff = std::chrono::duration_cast<std::chrono::milliseconds>(now - lastTimePoint).count();
  lastTimePoint = now;
  return diff;
}

void upwardPassRec(AdaptiveOctreeNode &node, Operators &op) {
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

void upwardPass(AdaptiveOctree &tree, Operators &op) {
  upwardPassRec(*tree.getRoot(), op);
  std::cout << "P2M and M2M took " << measureTime() << "ms" << std::endl;
}

void downwardPassRec1(AdaptiveOctreeNode &node, Operators &op) {
  op.M2L(node);
  if (!node.isLeaf()) {
    for (int i = 0; i < 8; ++i) {
      auto child = node.getChild(i);
      downwardPassRec1(*child, op);
    }
  }
}

void downwardPassRec2(AdaptiveOctreeNode &node, Operators &op) {
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

void downwardPass(AdaptiveOctree &tree, Operators &op) {
  downwardPassRec1(*tree.getRoot(), op);
  std::cout << "M2L took " << measureTime() << "ms" << std::endl;
  downwardPassRec2(*tree.getRoot(), op);
  std::cout << "L2L and L2P took " << measureTime() << "ms" << std::endl;
}

void calculateNearField(AdaptiveOctreeNode &node) {
  if (node.isLeaf()) {
    for (auto particle =
             node.getTree()->getDomain()->getRegionIterator(node.getNodeMinCorner(), node.getNodeMaxCorner());
         particle.isValid(); ++particle) {
      for (auto neighbourNode : node.getNeighbourList()) {
        for (auto otherParticle = node.getTree()->getDomain()->getRegionIterator(neighbourNode->getNodeMinCorner(),
                                                                                 neighbourNode->getNodeMaxCorner());
             otherParticle.isValid(); ++otherParticle) {
          if (particle->getID() != otherParticle->getID()) {
            double x = particle->getR()[0] - otherParticle->getR()[0];
            double y = particle->getR()[1] - otherParticle->getR()[1];
            double z = particle->getR()[2] - otherParticle->getR()[2];
            auto dist = std::sqrt(x * x + y * y + z * z);
            otherParticle->resultFMM += particle->charge / dist;
            otherParticle->shortRange += particle->charge / dist;
            exactCalculations++;
          }
        }
      }
    }
  } else {
    for (int i = 0; i < 8; ++i) {
      calculateNearField(*node.getChild(i));
    }
  }
}

void addParticleToCont(AutoPasCont &cont, double x, double y, double z, double charge) {
  FmmParticle particle({x, y, z}, {0, 0, 0}, particleId, charge);
  particleId++;
  cont.addParticle(particle);
}

int main(int argc, char **argv) {
  int orderOfExpansion = 10;
  int maxParticlesPerNode = 32;
  int numberOfParticles = 1000;
  double errorTolerance = 0.001;

  if (argc == 4) {
    orderOfExpansion = static_cast<int>(std::strtol(argv[1], nullptr, 10));
    maxParticlesPerNode = static_cast<int>(std::strtol(argv[2], nullptr, 10));
    numberOfParticles = static_cast<int>(std::strtol(argv[3], nullptr, 10));
  }
  std::cout << "orderOfExpansion = " << orderOfExpansion << std::endl;
  std::cout << "numberOfParticles = " << numberOfParticles << std::endl;

  Math3D::initMath();

  measureTime();

  std::array<double, 3> lowCorner({0, 0, 0});
  std::array<double, 3> highCorner({8, 6, 4});
  AutoPasCont cont = AutoPasCont();

  cont.setAllowedContainers(std::set<autopas::ContainerOption>{autopas::ContainerOption::directSum});

  cont.setBoxMin(lowCorner);
  cont.setBoxMax(highCorner);

  cont.init();

  /*addParticleToCont(cont, 0.424, 0.5, 0.5, 10);
  addParticleToCont(cont, 1.114, 0.5, 0.5, 10);
  addParticleToCont(cont, 3.583, 0.5, 0.5, 10);
  addParticleToCont(cont, 2.143, 0.5, 0.5, 10);
  addParticleToCont(cont, 1.434, 0.5, 0.5, 10);

  addParticleToCont(cont, 7.7, 5.6, 3.8, 10);*/

  std::random_device rd;
  std::mt19937 randomEngine(rd());
  std::uniform_real_distribution<double> random(0.0, 1.0);

  // Add particles at random positions.
  for (int i = 0; i < numberOfParticles; ++i) {
    double x = random(randomEngine) * (highCorner[0] - lowCorner[0]) + lowCorner[0];
    double y = random(randomEngine) * (highCorner[1] - lowCorner[1]) + lowCorner[1];
    double z = random(randomEngine) * (highCorner[2] - lowCorner[2]) + lowCorner[2];
    double charge = random(randomEngine) * 1.0;
    addParticleToCont(cont, x, y, z, charge);
  }

  /*addParticleToCont(cont, 4.78628, 3.17841, 2.2938, 12.9846);
  addParticleToCont(cont, 3.11351, 5.81931, 0.747073, 44.769);
  addParticleToCont(cont, 1.16475, 4.12244, 3.60193, 62.506);
  addParticleToCont(cont, 7.59052, 1.39932, 3.65955, 75.8548);
  addParticleToCont(cont, 6.69388, 4.43934, 3.43072, 98.2199);
  addParticleToCont(cont, 6.66218, 0.471121, 1.16736, 48.0065);
  addParticleToCont(cont, 3.9152, 1.90016, 1.98089, 3.07885);
  addParticleToCont(cont, 5.41736, 1.65419, 3.5447, 26.1224);
  addParticleToCont(cont, 1.96786, 2.38841, 3.96098, 72.2801);
  addParticleToCont(cont, 6.36797, 0.528699, 3.6491, 47.3985);
  addParticleToCont(cont, 4.13174, 3.53585, 1.5189, 46.9957);
  addParticleToCont(cont, 3.74137, 4.02255, 3.33478, 4.53346);
  addParticleToCont(cont, 3.53866, 3.62115, 2.47588, 21.898);
  addParticleToCont(cont, 5.71515, 4.14601, 3.90504, 99.9418);
  addParticleToCont(cont, 3.8807, 0.573223, 3.25886, 9.29967);
  addParticleToCont(cont, 4.08636, 0.0616424, 0.117308, 25.6611);*/

  AdaptiveOctree tree = AdaptiveOctree(cont, maxParticlesPerNode, orderOfExpansion);

  std::cout << "Init took " << measureTime() << "ms" << std::endl;

  // FMM
  Operators op(orderOfExpansion);
  std::cout << "Operators done" << std::endl;
  upwardPass(tree, op);
  std::cout << "upwardPass done" << std::endl;
  downwardPass(tree, op);
  std::cout << "downwardPass done" << std::endl;

  // Near field:
  calculateNearField(*tree.getRoot());

  std::cout << "Near field calculation took " << measureTime() << "ms" << std::endl;

  // Calculate exact result by calculating all interactions directly.

  for (auto particle = cont.begin(); particle.isValid(); ++particle) {
    for (auto otherParticle = cont.begin(); otherParticle.isValid(); ++otherParticle) {
      if (particle->getID() != otherParticle->getID()) {
        double x = particle->getR()[0] - otherParticle->getR()[0];
        double y = particle->getR()[1] - otherParticle->getR()[1];
        double z = particle->getR()[2] - otherParticle->getR()[2];
        auto dist = std::sqrt(x * x + y * y + z * z);
        particle->resultExact += otherParticle->charge / dist;
        exactCalculations++;
      }
    }
  }

  std::cout << "exact took " << measureTime() << "ms" << std::endl;

  bool printParticles = false;

  // Check results.
  for (auto particle = cont.begin(); particle.isValid(); ++particle) {
    if (particle->resultExact != 0) {
      double error = std::abs(particle->resultFMM / particle->resultExact);

      if (std::abs(error - 1.0) > errorTolerance) {
        std::cout << particle->getR()[0] << ", " << particle->getR()[1] << ", " << particle->getR()[2]
                  << ", charge = " << particle->charge << std::endl;
        std::cout << "long range " << particle->longRange << std::endl;
        std::cout << "short range " << particle->shortRange << std::endl;
        std::cout << "resultFMM " << particle->resultFMM << std::endl;
        std::cout << "resultExact " << particle->resultExact << std::endl;

        printParticles = true;
      }
    }
  }
  printParticles = false;

  if (printParticles) {
    for (auto particle = cont.begin(); particle.isValid(); ++particle) {
      std::cout << particle->getR()[0] << ", " << particle->getR()[1] << ", " << particle->getR()[2]
                << ", charge = " << particle->charge << std::endl;
    }
  }
  return EXIT_SUCCESS;
}