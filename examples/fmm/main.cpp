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
#include "../../tools/autopasTools/generators/GaussianGenerator.h"
#include "../../tools/autopasTools/generators/RandomGenerator.h"
#include "AdaptiveOctree.h"
#include "FmmParticle.h"
#include "Math3D.h"
#include "Operators.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/Timer.h"

int particleId = 0;

unsigned long long nearFieldCalculations = 0;
unsigned long long exactCalculations = 0;

double getRatio(long a, long b) { return static_cast<double>(a) / static_cast<double>(b); }

// P2M
void upwardPassP2M(AdaptiveOctreeNode &node, Operators &op) {
  if (node.isLeaf()) {
    op.P2M(node);
  } else {
    for (int i = 0; i < 8; ++i) {
      auto child = node.getChild(i);
      upwardPassP2M(*child, op);
    }
  }
}

// M2M
void upwardPassM2M(AdaptiveOctreeNode &node, Operators &op) {
  if (not node.isLeaf()) {
    for (int i = 0; i < 8; ++i) {
      auto child = node.getChild(i);
      upwardPassM2M(*child, op);
    }
    op.M2M(node);
  }
}

// M2L
void downwardPassM2L(AdaptiveOctreeNode &node, Operators &op) {
  op.M2L(node);
  if (!node.isLeaf()) {
    for (int i = 0; i < 8; ++i) {
      auto child = node.getChild(i);
      downwardPassM2L(*child, op);
    }
  }
}

// L2L
void downwardPassL2L(AdaptiveOctreeNode &node, Operators &op) {
  op.L2L(node);
  if (!node.isLeaf()) {
    for (int i = 0; i < 8; ++i) {
      auto child = node.getChild(i);
      downwardPassL2L(*child, op);
    }
  }
}

// L2P
void downwardPassL2P(AdaptiveOctreeNode &node, Operators &op) {
  if (!node.isLeaf()) {
    for (int i = 0; i < 8; ++i) {
      auto child = node.getChild(i);
      downwardPassL2P(*child, op);
    }
  } else {
    op.L2P(node);
  }
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
  bool checkResults = false;
  bool uniform = true;

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
                                         {"check", no_argument, nullptr, 'c'},
                                         {"non-uniform", no_argument, nullptr, 'u'},
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
      case 'c': {
        checkResults = true;
        break;
      }
      case 'u': {
        uniform = false;
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
  auto timer = autopas::utils::Timer();
  timer.start();

  // Initialize static fields of Math3D.
  Math3D::initialize();

  // Setup AutoPas object.
  std::array<double, 3> lowCorner({0, 0, 0});
  std::array<double, 3> highCorner({8, 8, 8});
  AutoPasCont cont = AutoPasCont();
  cont.setAllowedContainers(std::set<autopas::ContainerOption>{autopas::ContainerOption::linkedCells});
  cont.setBoxMin(lowCorner);
  cont.setBoxMax(highCorner);
  cont.setCutoff(1);
  cont.setVerletSkin(0);
  cont.init();

  std::random_device rd;
  std::default_random_engine randomEngine(42);
  std::uniform_real_distribution<double> random(0.0, 1.0);

  // Add particles at random positions.
  if (uniform) {
    FmmParticle defParticle({0, 0, 0}, {0, 0, 0}, 0, 1);
    autopasTools::generators::RandomGenerator::fillWithParticles(cont, defParticle, cont.getBoxMin(), cont.getBoxMax(),
                                                                 numberOfParticles);
  } else {
    FmmParticle defParticle1({0, 0, 0}, {0, 0, 0}, 0, 1);
    FmmParticle defParticle2({0, 0, 0}, {0, 0, 0}, numberOfParticles / 2, 1);
    autopasTools::generators::GaussianGenerator::fillWithParticles(cont, cont.getBoxMin(), cont.getBoxMax(),
                                                                   numberOfParticles / 2, defParticle1, {2, 2, 2});
    autopasTools::generators::GaussianGenerator::fillWithParticles(cont, cont.getBoxMin(), cont.getBoxMax(),
                                                                   numberOfParticles / 2, defParticle2, {6, 6, 6});
  }
  for (auto particle = cont.begin(); particle.isValid(); ++particle) {
    double charge = random(randomEngine) * 10.0;
    particle->charge = charge;
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

  long timeInit = timer.stop();
  std::cout << "Init took " << (timeInit / 1000) << "ms" << std::endl;
  timer.start();

  autopas::utils::Timer timerFmm;
  timerFmm.start();

  // FMM
  Operators op(orderOfExpansion);
  std::cout << "Init Operators done" << std::endl;

  upwardPassP2M(*tree.getRoot(), op);
  long timeP2M = timer.stop();

  std::cout << "P2M took " << (timeP2M / 1000) << "ms" << std::endl;

  timer.start();
  upwardPassM2M(*tree.getRoot(), op);
  long timeM2M = timer.stop();

  std::cout << "M2M took " << (timeM2M / 1000) << "ms" << std::endl;

  std::cout << "UpwardPass done" << std::endl;

  timer.start();
  downwardPassM2L(*tree.getRoot(), op);
  long timeM2L = timer.stop();
  std::cout << "M2L took " << (timeM2L / 1000) << "ms" << std::endl;

  timer.start();
  downwardPassL2L(*tree.getRoot(), op);
  long timeL2L = timer.stop();
  std::cout << "L2L took " << (timeL2L / 1000) << "ms" << std::endl;

  timer.start();
  downwardPassL2P(*tree.getRoot(), op);
  long timeL2P = timer.stop();

  std::cout << "L2P took " << (timeL2P / 1000) << "ms" << std::endl;

  std::cout << "DownwardPass done" << std::endl;

  // Near field:
  timer.start();
  calculateNearField(*tree.getRoot());
  long timeNear = timer.stop();

  std::cout << "Near field calculation took " << (timeNear / 1000) << "ms" << std::endl;

  long fmmTime = timeP2M + timeM2M + timeM2L + timeL2L + timeL2P + timeNear;
  long timeFar = fmmTime - timeNear;

  std::cout << std::endl << std::endl;
  std::cout << (fmmTime / 1000) << std::endl << std::endl;

  std::cout << (timeFar / 1000) << "\t" << (timeNear / 1000) << "\t" << getRatio(timeFar, fmmTime) << "\t"
            << getRatio(timeNear, fmmTime) << std::endl;

  std::cout << getRatio(timeP2M, timeFar) << "\t" << getRatio(timeM2M, timeFar) << "\t" << getRatio(timeM2L, timeFar)
            << "\t" << getRatio(timeL2L, timeFar) << "\t" << getRatio(timeL2P, timeFar) << "\t" << std::endl;

  if (checkResults) {
    // Calculate exact result by calculating all interactions directly.
    double maxAbsError = 0;
    double totalAbsError = 0;
    double maxRelError = 0;
    double totalRelError = 0;

    timer.start();
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
    long timeExact = timer.stop();

    std::cout << "exact took " << (timeExact / 1000) << "ms" << std::endl;

    std::cout << "nearFieldCalculations = " << nearFieldCalculations << std::endl;
    std::cout << "exactCalculations = " << exactCalculations << std::endl;

    bool correctResult = true;

    // Check results.
    for (auto particle = cont.begin(); particle.isValid(); ++particle) {
      if (particle->resultExact != 0) {
        int id = particle->getID();

        double absoluteError = std::abs(particle->resultFMM - particle->resultExact);

        double relativeError = absoluteError / std::abs(particle->resultExact);

        totalAbsError += absoluteError;
        maxAbsError = std::max(maxAbsError, absoluteError);

        totalRelError += relativeError;
        maxRelError = std::max(maxRelError, relativeError);

        if (relativeError > errorTolerance) {
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
    double avgAbsError = totalAbsError / numberOfParticles;
    double avgRelError = totalRelError / numberOfParticles;
    std::cout << "errors:" << std::endl
              << maxAbsError << "\t" << avgAbsError << "\t" << maxRelError << "\t" << avgRelError << std::endl;
    if (correctResult) {
      return EXIT_SUCCESS;
    } else {
      std::cout << "At least 1 result of the fast multipole method is wrong." << std::endl;
      return 1;
    }
  }
  return EXIT_SUCCESS;
}