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

void upwardPassRec(OctreeNode &node, Operators &op) {
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

void upwardPass(Octree &tree, Operators &op) {
  upwardPassRec(*tree.getRoot(), op);
  std::cout << "P2M and M2M took " << measureTime() << "ms" << std::endl;
}

void downwardPassRec1(OctreeNode &node, Operators &op) {
  op.M2L(node);
  if (!node.isLeaf()) {
    for (int i = 0; i < 8; ++i) {
      auto child = node.getChild(i);
      downwardPassRec1(*child, op);
    }
  }
}

void downwardPassRec2(OctreeNode &node, Operators &op) {
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

void downwardPass(Octree &tree, Operators &op) {
  downwardPassRec1(*tree.getRoot(), op);
  std::cout << "M2L took " << measureTime() << "ms" << std::endl;
  downwardPassRec2(*tree.getRoot(), op);
  std::cout << "L2L and L2P took " << measureTime() << "ms" << std::endl;
}

void addParticle(Octree &tree, const std::array<double, 3> &pos, double charge) {
  double x = pos[0];
  double y = pos[1];
  double z = pos[2];
  int cellX = static_cast<int>(x / tree.getCellSize());
  int cellY = static_cast<int>(y / tree.getCellSize());
  int cellZ = static_cast<int>(z / tree.getCellSize());
  FmmParticle particle({x, y, z}, {0, 0, 0}, particleId, charge);
  particleId++;

  tree.getCell(tree.getHeight(), cellX, cellY, cellZ)->getContainer()->addParticle(particle);
}

void generateParticleList(OctreeNode &node, std::vector<FmmParticle *> &particles) {
  if (node.isLeaf()) {
    for (auto iter = node.getContainer()->getRegionIterator(node.getLowCorner(), node.getHighCorner(),
                                                            autopas::IteratorBehavior::ownedOnly);
         iter.isValid(); ++iter) {
      particles.push_back(&*iter);
    }
  } else {
    for (int i = 0; i < 8; ++i) {
      generateParticleList(*node.getChild(i), particles);
    }
  }
}

void addParticlesFromNeighbours(OctreeNode &node) {
  if (node.isLeaf()) {
    // Check if cell contains particles.
    if (node.getContainer()
            ->getRegionIterator(node.getLowCorner(), node.getHighCorner(), autopas::IteratorBehavior::ownedOnly)
            .isValid()) {
      for (auto neighbour : *node.getNeighbourList()) {
        if (&node != neighbour && node.getContainer() != neighbour->getContainer()) {
          for (auto iter = neighbour->getContainer()->getRegionIterator(
                   neighbour->getLowCorner(), neighbour->getHighCorner(), autopas::IteratorBehavior::ownedOnly);
               iter.isValid(); ++iter) {
            node.getContainer()->addOrUpdateHaloParticle(*iter);
          }
        }
      }
    }
  } else {
    for (int i = 0; i < 8; ++i) {
      addParticlesFromNeighbours(*node.getChild(i));
    }
  }
}

void calculateNearField(OctreeNode &node /*, NearFieldFunctor &functor*/) {
  if (node.isLeaf()) {
    // node.getContainer()->iteratePairwise(&functor);

    for (auto part = node.getContainer()->getRegionIterator(node.getLowCorner(), node.getHighCorner(),
                                                            autopas::IteratorBehavior::ownedOnly);
         part.isValid(); ++part) {
      for (auto inter = node.getContainer()->getRegionIterator(node.getHaloLowCorner(), node.getHaloHighCorner());
           inter.isValid(); ++inter) {
        if (part->getID() != inter->getID()) {
          auto partPos = part->getR();
          auto interPos = inter->getR();
          double x = partPos[0] - interPos[0];
          double y = partPos[1] - interPos[1];
          double z = partPos[2] - interPos[2];
          auto dist = std::sqrt(x * x + y * y + z * z);
          double add = inter->charge / dist;
          part->resultFMM += add;
          part->shortRange += add;
          nearFieldCalculations++;
        }
      }
    }

  } else {
    for (int i = 0; i < 8; ++i) {
      calculateNearField(*node.getChild(i) /*, functor*/);
    }
  }
}

void addParticleToCont(AutoPasCont &cont, double x, double y, double z, double charge) {
  FmmParticle particle({x, y, z}, {0, 0, 0}, particleId, charge);
  particleId++;
  cont.addParticle(particle);
}

int main(int argc, char **argv) {
  int orderOfExpansion = 7;
  int treeSize = 8;
  int numberOfParticles = 10;
  double errorTolerance = 0.0001;

  if (argc == 4) {
    orderOfExpansion = static_cast<int>(std::strtol(argv[1], nullptr, 10));
    treeSize = static_cast<int>(std::strtol(argv[2], nullptr, 10));
    numberOfParticles = static_cast<int>(std::strtol(argv[3], nullptr, 10));
  }
  std::cout << "orderOfExpansion = " << orderOfExpansion << std::endl;
  std::cout << "treeSize = " << treeSize << std::endl;
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

  addParticleToCont(cont, 0.424, 0.5, 0.5, 10);
  addParticleToCont(cont, 1.114, 0.5, 0.5, 10);
  addParticleToCont(cont, 3.583, 0.5, 0.5, 10);
  addParticleToCont(cont, 2.143, 0.5, 0.5, 10);
  addParticleToCont(cont, 1.434, 0.5, 0.5, 10);

  AdaptiveOctree tree = AdaptiveOctree(cont, 4, 8);

  /*std::random_device rd;
  std::mt19937 randomEngine(rd());
  std::uniform_real_distribution<double> random(0.0, 1.0);

  // Add particles at random positions.
  for (int i = 0; i < numberOfParticles; ++i) {
    double x = random(randomEngine) * (highCorner[0] - lowCorner[0]) + lowCorner[0];
    double y = random(randomEngine) * (highCorner[1] - lowCorner[1]) + lowCorner[1];
    double z = random(randomEngine) * (highCorner[2] - lowCorner[2]) + lowCorner[2];
    double charge = random(randomEngine) * 100.0;
    FmmParticle particle({x, y, z}, {0, 0, 0}, particleId, charge);
    particleId++;
    cont.addParticle(particle);
  }*/



  /*Octree tree(treeSize, 2.0, 4);

  std::random_device rd;
  std::mt19937 randomEngine(rd());
  std::uniform_real_distribution<double> random(0.0, 1.0);

  // Add particles at random positions.
  for (int i = 0; i < numberOfParticles; ++i) {
    double x = random(randomEngine) * tree.getSize() * tree.getCellSize();
    double y = random(randomEngine) * tree.getSize() * tree.getCellSize();
    double z = random(randomEngine) * tree.getSize() * tree.getCellSize();
    double charge = random(randomEngine) * 100.0;
    addParticle(tree, {x, y, z}, charge);
  }

  std::cout << "Init took " << measureTime() << "ms" << std::endl;

  // FMM
  Operators op(orderOfExpansion);
  upwardPass(tree, op);
  downwardPass(tree, op);

  // Calculate near field:
  // NearFieldFunctor functor(tree.getSize() * tree.getCellSize() * 2.0);

  addParticlesFromNeighbours(*tree.getRoot());

  std::cout << "Near field setup took " << measureTime() << "ms" << std::endl;*/

  // calculateNearField(*tree.getRoot() /*, functor*/);

  /*
  std::cout << "Near field calculation took " << measureTime() << "ms" << std::endl;

  // Calculate exact result by calculating all interactions directly.
  // Create list of all particles.
  std::vector<FmmParticle *> particles(0);
  generateParticleList(*tree.getRoot(), particles);

  for (auto part : particles) {
    for (auto inter : particles) {
      if (part->getID() != inter->getID()) {
        auto distVec = Math3D::subtract(part->getR(), inter->getR());
        auto spherical = Math3D::toSpherical(distVec);
        auto dist = spherical[0];
        part->resultExact += inter->charge / dist;
        exactCalculations++;
      }
    }
  }

  std::cout << "exact took " << measureTime() << "ms" << std::endl;

  std::cout << "nearFieldCalculations = " << nearFieldCalculations << std::endl;
  std::cout << "exactCalculations = " << exactCalculations << std::endl;

  // Check results.
  for (auto part : particles) {
    if (part->resultExact != 0) {
      double error = std::abs(part->resultFMM / part->resultExact);

      if (std::abs(error - 1.0) > errorTolerance) {
        std::cout << part->getR()[0] << ", " << part->getR()[1] << ", " << part->getR()[2]
                  << ", charge = " << part->charge << std::endl;
        std::cout << "long range " << part->longRange << std::endl;
        std::cout << "short range " << part->shortRange << std::endl;
        std::cout << "resultFMM " << part->resultFMM << std::endl;
        std::cout << "resultExact " << part->resultExact << std::endl;

        assert(false);
      }
    }
  }*/

  return EXIT_SUCCESS;
}