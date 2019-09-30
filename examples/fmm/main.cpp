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
#include "Octree.h"
#include "Operators.h"
#include "autopas/AutoPas.h"

std::chrono::steady_clock::time_point lastTimePoint = std::chrono::steady_clock::now();
int particleId = 42;

int nearFieldCalculations = 0;
int exactCalculations = 0;

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
    for (auto iter = node.getContainer()->begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
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
    if (node.getContainer()->begin(autopas::IteratorBehavior::ownedOnly).isValid()) {
      for (auto neighbour : *node.getNeighbourList()) {
        if (&node != neighbour) {
          for (auto iter = neighbour->getContainer()->begin(autopas::IteratorBehavior::ownedOnly); iter.isValid();
               ++iter) {
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

    for (auto part = node.getContainer()->begin(autopas::IteratorBehavior::ownedOnly); part.isValid(); ++part) {
      for (auto inter = node.getContainer()->begin(); inter.isValid(); ++inter) {
        if (part->getID() != inter->getID()) {
          auto distVec = Math3D::subtract(part->getR(), inter->getR());
          auto spherical = Math3D::toSpherical(distVec);
          auto dist = spherical[0];
          part->resultFMM += inter->charge / dist;
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

int main(int argc, char **argv) {
  int orderOfExpansion = 7;
  int treeSize = 8;
  int numberOfParticles = 10;
  double errorTolerance = 0.01;

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

  Octree tree(treeSize, 2.0);

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
  calculateNearField(*tree.getRoot() /*, functor*/);

  std::cout << "Near field took " << measureTime() << "ms" << std::endl;

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
    double error = std::abs(part->resultFMM - part->resultExact);

    if (error > errorTolerance) {
      std::cout << part->getR()[0] << ", " << part->getR()[1] << ", " << part->getR()[2]
                << ", charge = " << part->charge << std::endl;
      std::cout << part->resultFMM << std::endl;
      std::cout << part->resultExact << std::endl;

      assert(false);
    }
  }

  return EXIT_SUCCESS;
}