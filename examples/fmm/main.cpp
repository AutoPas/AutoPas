/**
 * @file main.cpp
 * @date 14.09.19
 * @author Joachim Marin
 */

#include <getopt.h>
#include <exception>
#include <iomanip>
#include <iostream>
#include <random>
#include "AdaptiveOctree.h"
#include "FmmParticle.h"
#include "Math3D.h"
#include "Operators.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/Timer.h"

autopas::utils::Timer *timer;
int particleId = 0;

unsigned long long nearFieldCalculations = 0;
unsigned long long exactCalculations = 0;

double measureTime() {
  double ret = timer->stop();
  timer->start();
  return ret;
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
  std::cout << "P2M and M2M took " << measureTime() << "s" << std::endl;
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
  std::cout << "M2L took " << measureTime() << "s" << std::endl;
  downwardPassRec2(*tree.getRoot(), op);
  std::cout << "L2L and L2P took " << measureTime() << "s" << std::endl;
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
  // Default parameters
  int orderOfExpansion = 8;
  int maxParticlesPerNode = 64;
  int numberOfParticles = 500;
  double errorTolerance = 0.01;
  int minDepth = 0;
  int maxDepth = -1;

  bool displayHelp = false;
  int option, option_index;
  // clang-format off
  static struct option long_options[] = {{"help", no_argument, nullptr, 'h'},
                                         {"order", required_argument, nullptr, 'o'},
                                         {"particles-per-cell", required_argument, nullptr, 'p'},
                                         {"particles-total", required_argument, nullptr, 'n'},
                                         {"error-tolerance", required_argument, nullptr, 'e'},
                                         {"depth-min", required_argument, nullptr, 'm'},
                                         {"depth-max", required_argument, nullptr, 'M'},
                                         {nullptr, 0, nullptr, 0}};
  // clang-format on

  // Parse command line parameters.
  std::string strArg;
  while ((option = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    if (optarg != nullptr) strArg = optarg;
    transform(strArg.begin(), strArg.end(), strArg.begin(), ::tolower);
    switch (option) {
      case 'h': {
        displayHelp = true;
        break;
      }
      case 'o': {
        try {
          orderOfExpansion = std::stoi(strArg);
        } catch (const std::exception &) {
          std::cerr << "Error parsing order of expansion: " << optarg << std::endl;
          displayHelp = true;
        }
        break;
      }
      case 'p': {
        try {
          maxParticlesPerNode = std::stoi(strArg);
        } catch (const std::exception &) {
          std::cerr << "Error parsing maximum particles per cell: " << optarg << std::endl;
          displayHelp = true;
        }
        break;
      }
      case 'n': {
        try {
          numberOfParticles = std::stoi(strArg);
        } catch (const std::exception &) {
          std::cerr << "Error parsing total number of particles: " << optarg << std::endl;
          displayHelp = true;
        }
        break;
      }
      case 'e': {
        try {
          errorTolerance = std::stod(strArg);
        } catch (const std::exception &) {
          std::cerr << "Error parsing error tolerance:  " << optarg << std::endl;
          displayHelp = true;
        }
        break;
      }
      case 'm': {
        try {
          minDepth = std::stoi(strArg);
        } catch (const std::exception &) {
          std::cerr << "Error parsing minimum octree depth: " << optarg << std::endl;
          displayHelp = true;
        }
        break;
      }
      case 'M': {
        try {
          maxDepth = std::stoi(strArg);
        } catch (const std::exception &) {
          std::cerr << "Error parsing maximum octree depth: " << optarg << std::endl;
          displayHelp = true;
        }
        break;
      }
      default: {
        // error message handled by getopt
        displayHelp = true;
      }
    }
  }

  if (displayHelp) {
    std::cout << "Usage: " << argv[0] << std::endl;
    for (auto o : long_options) {
      if (o.name == nullptr) continue;
      std::cout << "    --" << std::setw(32) << std::left << o.name;
      if (o.has_arg) {
        std::cout << "option";
      }
      std::cout << std::endl;
    }

    return -1;
  }

  // Print parameters.
  std::cout << "orderOfExpansion = " << orderOfExpansion << std::endl;
  std::cout << "maxParticlesPerNode = " << maxParticlesPerNode << std::endl;
  std::cout << "numberOfParticles = " << numberOfParticles << std::endl;
  std::cout << "errorTolerance = " << errorTolerance << std::endl;
  std::cout << "minDepth = " << minDepth << std::endl;
  std::cout << "maxDepth = " << maxDepth << std::endl;

  // Start timer.
  auto tmp = autopas::utils::Timer();
  timer = &tmp;
  timer->start();

  // Initialize static fields of Math3D.
  Math3D::initialize();

  // Setup AutoPas object.
  std::array<double, 3> lowCorner({0, 0, 0});
  std::array<double, 3> highCorner({8, 8, 8});
  AutoPasCont cont = AutoPasCont();
  cont.setAllowedContainers(std::set<autopas::ContainerOption>{autopas::ContainerOption::linkedCells});
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

  // Create particle vector containing all particles in the domain. Will be used to calculate exact results.

  std::cout << "currentMaxDepth = " << tree.currentMaxDepth << std::endl;
  std::cout << "numberOfNodes = " << tree.numberOfNodes << std::endl;
  std::cout << "numberOfLeaves = " << tree.numberOfLeaves << std::endl;
  std::cout << "near field nodes = " << tree.totalNearFieldNodes << " (total), "
            << (1.0 * tree.totalNearFieldNodes / tree.numberOfNodes) << " (average)" << std::endl;
  std::cout << "interaction nodes = " << tree.totalInteractionNodes << " (total), "
            << (1.0 * tree.totalInteractionNodes / tree.numberOfNodes) << " (average)" << std::endl;

  std::cout << "Init took " << measureTime() << "s" << std::endl;

  // FMM
  Operators op(orderOfExpansion);
  std::cout << "Init Operators done" << std::endl;
  upwardPass(tree, op);
  std::cout << "UpwardPass done" << std::endl;
  downwardPass(tree, op);
  std::cout << "DownwardPass done" << std::endl;

  // Near field:
  calculateNearField(*tree.getRoot());

  std::cout << "Near field calculation took " << measureTime() << "s" << std::endl;

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

  std::cout << "exact took " << measureTime() << "s" << std::endl;

  std::cout << "nearFieldCalculations = " << nearFieldCalculations << std::endl;
  std::cout << "exactCalculations = " << exactCalculations << std::endl;

  bool correctResult = true;

  // Check results.
  for (auto particle = cont.begin(); particle.isValid(); ++particle) {
    if (particle->resultExact != 0) {
      int id = particle->getID();
      double error = std::abs(particle->resultFMM / particle->resultExact);

      if (std::abs(error - 1.0) > errorTolerance) {
        std::cout << "[ID=" << id << "] " << particle->getR()[0] << ", " << particle->getR()[1] << ", "
                  << particle->getR()[2] << ", charge = " << particle->charge << std::endl;
        std::cout << "long range " << particle->longRange << std::endl;
        std::cout << "short range " << particle->shortRange << std::endl;
        std::cout << "resultFMM " << particle->resultFMM << std::endl;
        std::cout << "resultExact " << particle->resultExact << std::endl;
        correctResult = false;
      }
    }
  }
  if (correctResult) {
    return EXIT_SUCCESS;
  } else {
    std::cout << "At least 1 result of the fast multipole method is wrong." << std::endl;
    return 1;
  }
}