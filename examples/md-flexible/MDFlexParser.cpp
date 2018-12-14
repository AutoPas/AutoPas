/**
 * @file MDFlexParser.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include "MDFlexParser.h"

bool MDFlexParser::parseInput(int argc, char **argv) {
  bool displayHelp = false;
  int option, option_index;
  static struct option long_options[] = {{"box-length", required_argument, nullptr, 'b'},
                                         {"container", required_argument, nullptr, 'c'},
                                         {"container-selector-strategy", required_argument, nullptr, 'k'},
                                         {"cutoff", required_argument, nullptr, 'C'},
                                         {"distribution-mean", required_argument, nullptr, 'm'},
                                         {"distribution-stddeviation", required_argument, nullptr, 'z'},
                                         {"data-layout", required_argument, nullptr, 'd'},
                                         {"functor", required_argument, nullptr, 'f'},
                                         {"help", no_argument, nullptr, 'h'},
                                         {"iterations", required_argument, nullptr, 'i'},
                                         {"no-flops", no_argument, nullptr, 'F'},
                                         {"particles-generator", required_argument, nullptr, 'g'},
                                         {"particles-per-dimension", required_argument, nullptr, 'n'},
                                         {"particles-total", required_argument, nullptr, 'N'},
                                         {"particle-spacing", required_argument, nullptr, 's'},
                                         {"traversal", required_argument, nullptr, 't'},
                                         {"traversal-selector-strategy", required_argument, nullptr, 'T'},
                                         {"tuning-interval", required_argument, nullptr, 'I'},
                                         {"tuning-samples", required_argument, nullptr, 'S'},
                                         {"log-level", required_argument, nullptr, 'l'},
                                         {"verlet-rebuild-frequency", required_argument, nullptr, 'v'},
                                         {"verlet-skin-radius", required_argument, nullptr, 'r'},
                                         {"vtk", required_argument, nullptr, 'w'},
                                         {nullptr, 0, nullptr, 0}};  // needed to signal the end of the array
  string strArg;
  int bla = 0;
  while ((option = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    ++bla;
    if (optarg != nullptr) strArg = optarg;
    transform(strArg.begin(), strArg.end(), strArg.begin(), ::tolower);
    switch (option) {
      case 'b': {
        try {
          boxLength = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing box length: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'c': {
        // delete default argument
        containerOptions.clear();
        if (strArg.find("direct") != string::npos or strArg.find("ds") != string::npos) {
          containerOptions.push_back(autopas::directSum);
        }
        if (strArg.find("linked") != string::npos or strArg.find("lc") != string::npos) {
          containerOptions.push_back(autopas::linkedCells);
        }
        if (strArg.find("verlet") != string::npos or strArg.find("vl") != string::npos) {
          containerOptions.push_back(autopas::verletLists);
        }
        if (strArg.find("vcells") != string::npos) {
          containerOptions.push_back(autopas::verletListsCells);
        }
        if (strArg.find("vcluster") != string::npos) {
          containerOptions.push_back(autopas::verletClusterLists);
        }
        if (containerOptions.empty()) {
          cerr << "Unknown container option: " << strArg << endl;
          cerr << "Please use 'DirectSum', 'LinkedCells', 'VerletLists', 'vcells' or 'vcluster'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'C': {
        try {
          cutoff = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing cutoff Radius: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'd': {
        if (strArg.find("aos") != string::npos) {
          dataLayoutOption = autopas::aos;
        } else if (strArg.find("soa") != string::npos) {
          dataLayoutOption = autopas::soa;
        } else {
          cerr << "Unknown data layout : " << strArg << endl;
          cerr << "Please use 'AoS' or 'SoA'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'f': {
        if (strArg.find("lj") != string::npos || strArg.find("lennard-jones") != string::npos) {
          functorOption = lj12_6;
        } else {
          cerr << "Unknown functor : " << strArg << endl;
          cerr << "Please use 'Lennard-Jones', you have no options here :P" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'F': {
        measureFlops = false;
        break;
      }
      case 'g': {
        if (strArg.find("grid") != string::npos) {
          generatorOption = GeneratorOption::grid;
        } else if (strArg.find("uni") != string::npos) {
          generatorOption = GeneratorOption::uniform;
        } else if (strArg.find("gaus") != string::npos) {
          generatorOption = GeneratorOption::gaussian;
        } else {
          cerr << "Unknown generator : " << strArg << endl;
          cerr << "Please use 'Grid' or 'Gaussian'" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'h': {
        displayHelp = true;
        break;
      }
      case 'i': {
        try {
          iterations = stoul(strArg);
          if (iterations < 1) {
            cerr << "IterationNumber of iterations has to be a positive integer!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing number of iterations: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'I': {
        try {
          tuningInterval = (unsigned int)stoul(strArg);
          if (tuningInterval < 1) {
            cerr << "Tuning interval has to be a positive integer!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing tuning interval: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'k': {
        if (strArg.find("abs") != string::npos) {
          containerSelectorStrategy = autopas::SelectorStrategy::fastestAbs;
        } else if (strArg.find("mea") != string::npos) {
          containerSelectorStrategy = autopas::SelectorStrategy::fastestMean;
        } else if (strArg.find("med") != string::npos) {
          containerSelectorStrategy = autopas::SelectorStrategy::fastestMedian;
        } else {
          cerr << "Unknown Container Selector Strategy: " << strArg << endl;
          cerr << "Please use 'fastestAbs', 'fastestMean' or 'fastestMedian'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'l': {
        switch (strArg[0]) {
          case 't': {
            logLevel = spdlog::level::trace;
            break;
          }
          case 'd': {
            logLevel = spdlog::level::debug;
            break;
          }
          case 'i': {
            logLevel = spdlog::level::info;
            break;
          }
          case 'w': {
            logLevel = spdlog::level::warn;
            break;
          }
          case 'e': {
            logLevel = spdlog::level::err;
            break;
          }
          case 'c': {
            logLevel = spdlog::level::critical;
            break;
          }
          case 'o': {
            logLevel = spdlog::level::off;
            break;
          }
          default: {
            cerr << "Unknown Log Level : " << strArg << endl;
            cerr << "Please use 'trace', 'debug', 'info', 'warning', 'error', 'critical' or 'off'." << endl;
            displayHelp = true;
          }
        }
        break;
      }
      case 'm': {
        try {
          distributionMean = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing distribution mean: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'n': {
        try {
          particlesPerDim = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing number of particles per dimension: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'N': {
        try {
          particlesTotal = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing total number of particles: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 's': {
        try {
          particleSpacing = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing separation of particles: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'S': {
        try {
          tuningSamples = (unsigned int)stoul(strArg);
          if (tuningSamples < 1) {
            cerr << "Tuning samples has to be a positive integer!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing number of tuning samples: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 't': {
        traversalOptions.clear();
        if (strArg.find("c08") != string::npos) {
          traversalOptions.push_back(autopas::TraversalOptions::c08);
        }
        if (strArg.find("c01") != string::npos) {
          traversalOptions.push_back(autopas::TraversalOptions::c01);
        }
        if (strArg.find("c18") != string::npos) {
          traversalOptions.push_back(autopas::TraversalOptions::c18);
        }
        if (strArg.find("sli") != string::npos) {
          traversalOptions.push_back(autopas::TraversalOptions::sliced);
        }
        if (strArg.find("dir") != string::npos) {
          traversalOptions.push_back(autopas::TraversalOptions::directSumTraversal);
        }
        if (traversalOptions.empty()) {
          cerr << "Unknown Traversal : " << strArg << endl;
          cerr << "Please use 'c08', 'c01', 'c18' or 'sliced'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'T': {
        if (strArg.find("abs") != string::npos) {
          traversalSelectorStrategy = autopas::SelectorStrategy::fastestAbs;
        } else if (strArg.find("mea") != string::npos) {
          traversalSelectorStrategy = autopas::SelectorStrategy::fastestMean;
        } else if (strArg.find("med") != string::npos) {
          traversalSelectorStrategy = autopas::SelectorStrategy::fastestMedian;
        } else {
          cerr << "Unknown Traversal Selector Strategy: " << strArg << endl;
          cerr << "Please use 'fastestAbs', 'fastestMean' or 'fastestMedian'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'v': {
        try {
          verletRebuildFrequency = (unsigned int)stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet-rebuild-frequency: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'w': {
        writeVTK = strArg;
        break;
      }
      case 'r': {
        try {
          verletSkinRadius = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet-skin-radius: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'z': {
        try {
          distributionStdDev = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing distribution standard deviation: " << optarg << endl;
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
    cout << "Usage: " << argv[0] << endl;
    for (auto o : long_options) {
      if (o.name == nullptr) continue;
      cout << "    --" << setw(valueOffset + 2) << left << o.name;
      if (o.has_arg) cout << "option";
      cout << endl;
    }
    return false;
  }
  return true;
}

std::string selectorStrategyToString(autopas::SelectorStrategy selectorStrategy) {
  std::string retString;
  switch (selectorStrategy) {
    case autopas::SelectorStrategy::fastestAbs: {
      retString = "Fastest Absolute Value";
      break;
    }
    case autopas::SelectorStrategy::fastestMean: {
      retString = "Fastest Mean Value";
      break;
    }
    case autopas::SelectorStrategy::fastestMedian: {
      retString = "Fastest Median Value";
      break;
    }
  }
  return retString;
}

void MDFlexParser::printConfig() {
  constexpr size_t valueOffset = 32;
  cout << setw(valueOffset) << left << "Container"
       << ":  ";
  for (auto &op : containerOptions) {
    switch (op) {
      case autopas::ContainerOptions::directSum: {
        cout << "DirectSum, ";
        break;
      }
      case autopas::ContainerOptions::linkedCells: {
        cout << "LinkedCells, ";
        break;
      }
      case autopas::ContainerOptions::verletLists: {
        cout << "VerletLists, ";
        break;
      }
      case autopas::ContainerOptions::verletListsCells: {
        cout << "VerletListsCells, ";
        break;
      }
      case autopas::ContainerOptions::verletClusterLists: {
        cout << "VerletClusterLists, ";
        break;
      }
    }
  }
  // deletes last comma
  cout << "\b\b  " << endl;

  // if verlet lists are in the container options print verlet config data
  if (find(containerOptions.begin(), containerOptions.end(), autopas::ContainerOptions::verletLists) !=
      containerOptions.end()) {
    cout << setw(valueOffset) << left << "Verlet rebuild frequency"
         << ":  " << verletRebuildFrequency << endl;

    cout << setw(valueOffset) << left << "Verlet skin radius"
         << ":  " << verletSkinRadius << endl;
  }

  if (containerOptions.size() > 1) {
    cout << setw(valueOffset) << left << "Container Selector Strategy"
         << ":  " << selectorStrategyToString(containerSelectorStrategy) << endl;
  }

  cout << setw(valueOffset) << left << "Data Layout"
       << ":  ";
  switch (dataLayoutOption) {
    case autopas::DataLayoutOption::aos: {
      cout << "Array-of-Structures" << endl;
      break;
    }
    case autopas::DataLayoutOption::soa: {
      cout << "Structure-of-Arrays" << endl;
      break;
    }
  }

  cout << setw(valueOffset) << left << "Functor"
       << ":  ";
  switch (functorOption) {
    case FunctorOption::lj12_6: {
      cout << "Lennard-Jones (12-6)" << endl;
      break;
    }
  }

  cout << setw(valueOffset) << left << "Cutoff radius"
       << ":  " << cutoff << endl;

  cout << setw(valueOffset) << left << "Particle Generator"
       << ":  ";
  switch (generatorOption) {
    case GeneratorOption::grid: {
      cout << "Grid generator" << endl;
      cout << setw(valueOffset) << left << "Particle spacing"
           << ":  " << particleSpacing << endl;

      cout << "Particles" << endl;
      cout << setw(valueOffset) << left << "  per dimension"
           << ":  " << particlesPerDim << endl;
      cout << setw(valueOffset) << left << "  total"
           << ":  " << (particlesPerDim * particlesPerDim * particlesPerDim) << endl;
      break;
    }
    case GeneratorOption::gaussian: {
      cout << "Gaussian generator" << endl;
      cout << setw(valueOffset) << left << "Box length"
           << ":  " << boxLength << endl;
      cout << setw(valueOffset) << left << "Distribution mean"
           << ":  " << distributionMean << endl;
      cout << setw(valueOffset) << left << "Distribution standard deviation"
           << ":  " << distributionStdDev << endl;

      cout << "Particles" << endl;
      cout << setw(valueOffset) << left << "  total"
           << ":  " << particlesTotal << endl;
      break;
    }
    case GeneratorOption::uniform: {
      cout << "Uniform generator" << endl;
      cout << setw(valueOffset) << left << "Box length"
           << ":  " << boxLength << endl;
      cout << "Particles" << endl;
      cout << setw(valueOffset) << left << "  total"
           << ":  " << particlesTotal << endl;
      break;
    }
  }

  cout << setw(valueOffset) << left << "Allowed traversals"
       << ":  ";
  for (auto &t : traversalOptions) {
    switch (t) {
      case autopas::TraversalOptions::c08: {
        cout << "c08, ";
        break;
      }
      case autopas::TraversalOptions::sliced: {
        cout << "sliced, ";
        break;
      }
      case autopas::TraversalOptions::c18: {
        cout << "c18, ";
        break;
      }
      case autopas::TraversalOptions::c01: {
        cout << "c01, ";
        break;
      }
      case autopas::TraversalOptions::directSumTraversal: {
        cout << "direct sum, ";
        break;
      }
      default:
        break;
    }
  }
  // deletes last comma
  cout << "\b\b  " << endl;

  if (traversalOptions.size() > 1) {
    cout << setw(valueOffset) << left << "Traversal Selector Strategy"
         << ":  " << selectorStrategyToString(traversalSelectorStrategy) << endl;
  }

  cout << setw(valueOffset) << left << "Iterations"
       << ":  " << iterations << endl;
  cout << setw(valueOffset) << left << "Tuning Interval"
       << ":  " << tuningInterval << endl;
  cout << setw(valueOffset) << left << "Tuning Samples"
       << ":  " << tuningSamples << endl;
}

std::vector<autopas::ContainerOptions> MDFlexParser::getContainerOptions() const { return containerOptions; }

double MDFlexParser::getCutoff() const { return cutoff; }

autopas::DataLayoutOption MDFlexParser::getDataLayoutOption() const { return dataLayoutOption; }

MDFlexParser::FunctorOption MDFlexParser::getFunctorOption() const { return functorOption; }

size_t MDFlexParser::getIterations() const { return iterations; }

size_t MDFlexParser::getParticlesPerDim() const { return particlesPerDim; }

double MDFlexParser::getParticleSpacing() const { return particleSpacing; }

size_t MDFlexParser::getParticlesTotal() const { return particlesTotal; }

const vector<autopas::TraversalOptions> &MDFlexParser::getTraversalOptions() const { return traversalOptions; }

unsigned int MDFlexParser::getVerletRebuildFrequency() const { return verletRebuildFrequency; }

double MDFlexParser::getVerletSkinRadius() const { return verletSkinRadius; }

MDFlexParser::GeneratorOption MDFlexParser::getGeneratorOption() const { return generatorOption; }

double MDFlexParser::getDistributionMean() const { return distributionMean; }

double MDFlexParser::getDistributionStdDev() const { return distributionStdDev; }

string MDFlexParser::getWriteVTK() const { return writeVTK; }

double MDFlexParser::getBoxLength() {
  if (boxLength == -1) boxLength = ceil(2 * distributionMean);
  return boxLength;
}

bool MDFlexParser::getMeasureFlops() const { return measureFlops; }

unsigned int MDFlexParser::getTuningInterval() const { return tuningInterval; }

unsigned int MDFlexParser::getTuningSamples() const { return tuningSamples; }

autopas::Logger::LogLevel MDFlexParser::getLogLevel() const { return logLevel; }

autopas::SelectorStrategy MDFlexParser::getTraversalSelectorStrategy() const { return traversalSelectorStrategy; }

autopas::SelectorStrategy MDFlexParser::getContainerSelectorStrategy() const { return containerSelectorStrategy; }
// spdlog::level::level_enum MDFlexParser::getLogLevel() const { return logLevel; }
