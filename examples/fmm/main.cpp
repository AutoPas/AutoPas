/**
 * @file main.cpp
 * @date 14.09.19
 * @author Joachim Marin
 */

#include <chrono>
#include <iostream>
#include <random>
#include "AdaptiveOctree.h"
#include "FmmParticle.h"
#include "Math3D.h"
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

// P2M + M2M
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


// M2L
void downwardPassRec1(AdaptiveOctreeNode &node, Operators &op) {
  op.M2L(node);
  if (!node.isLeaf()) {
    for (int i = 0; i < 8; ++i) {
      auto child = node.getChild(i);
      downwardPassRec1(*child, op);
    }
  }
}


// L2L + L2P
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
    int tmp = 0;

    for (auto nearFieldNode : node.getNearFieldList()) {
      for (auto particle =
               node.getTree()->getDomain()->getRegionIterator(node.getNodeMinCorner(), node.getNodeMaxCorner());
           particle.isValid(); ++particle) {
        for (auto otherParticle = node.getTree()->getDomain()->getRegionIterator(nearFieldNode->getNodeMinCorner(),
                                                                                 nearFieldNode->getNodeMaxCorner());
             otherParticle.isValid(); ++otherParticle) {
          if (particle->getID() != otherParticle->getID()) {
            double x = particle->getR()[0] - otherParticle->getR()[0];
            double y = particle->getR()[1] - otherParticle->getR()[1];
            double z = particle->getR()[2] - otherParticle->getR()[2];
            auto dist = 1.0 / std::sqrt(x * x + y * y + z * z);
            particle->resultFMM += otherParticle->charge * dist;
            particle->shortRange += otherParticle->charge * dist;

            tmp++;
          }
        }
      }
    }
    nearFieldCalculations += tmp;

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
  int orderOfExpansion = 8;
  int maxParticlesPerNode = 16;
  int numberOfParticles = 200;
  double errorTolerance = 0.002;
  int minDepth = 0;
  int maxDepth = -1;

  if (argc == 6) {
    orderOfExpansion = static_cast<int>(std::strtol(argv[1], nullptr, 10));
    maxParticlesPerNode = static_cast<int>(std::strtol(argv[2], nullptr, 10));
    numberOfParticles = static_cast<int>(std::strtol(argv[3], nullptr, 10));
    minDepth = static_cast<int>(std::strtol(argv[4], nullptr, 10));
    maxDepth = static_cast<int>(std::strtol(argv[5], nullptr, 10));
  }
  std::cout << "orderOfExpansion = " << orderOfExpansion << std::endl;
  std::cout << "maxParticlesPerNode = " << maxParticlesPerNode << std::endl;
  std::cout << "numberOfParticles = " << numberOfParticles << std::endl;

  Math3D::initMath();

  measureTime();

  // Setup AutoPas object.
  std::array<double, 3> lowCorner({0, 0, 0});
  std::array<double, 3> highCorner({8, 6, 4});
  AutoPasCont cont = AutoPasCont();
  cont.setAllowedContainers(std::set<autopas::ContainerOption>{autopas::ContainerOption::directSum});
  cont.setBoxMin(lowCorner);
  cont.setBoxMax(highCorner);
  cont.init();

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

  AdaptiveOctree tree = AdaptiveOctree(cont, maxParticlesPerNode, orderOfExpansion, minDepth, maxDepth);

  std::cout << "currentMaxDepth = " << tree.currentMaxDepth << std::endl;
  std::cout << "numberOfNodes = " << tree.numberOfNodes << std::endl;
  std::cout << "numberOfLeaves = " << tree.numberOfLeaves << std::endl;
  std::cout << "near field nodes = " << tree.totalNearFieldNodes << " (total), "
            << (1.0 * tree.totalNearFieldNodes / tree.numberOfNodes) << " (average)" << std::endl;
  std::cout << "interaction nodes = " << tree.totalInteractionNodes << " (total), "
            << (1.0 * tree.totalInteractionNodes / tree.numberOfNodes) << " (average)" << std::endl;

  std::cout << "Init took " << measureTime() << "ms" << std::endl;

  // FMM
  Operators op(orderOfExpansion);
  std::cout << "Init Operators done" << std::endl;
  upwardPass(tree, op);
  std::cout << "UpwardPass done" << std::endl;
  downwardPass(tree, op);
  std::cout << "DownwardPass done" << std::endl;

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

  std::cout << "nearFieldCalculations = " << nearFieldCalculations << std::endl;
  std::cout << "exactCalculations = " << exactCalculations << std::endl;

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
      }
    }
  }

  return EXIT_SUCCESS;
}