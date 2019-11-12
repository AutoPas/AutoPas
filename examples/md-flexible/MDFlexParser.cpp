/**
 * @file MDFlexParser.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include "MDFlexParser.h"
#include "autopas/utils/StringUtils.h"

template <class T>
std::string iterableToString(T arr) {
  std::ostringstream ss;
  for (auto a : arr) {
    ss << a.to_string() << ", ";
  }
  // remove last comma by removing last two characters.
  auto str = ss.str();
  return arr.empty() ? "" : str.substr(0, str.size() - 2);
}

bool MDFlexParser::parseInput(int argc, char **argv) {
  bool displayHelp = false;
  int option, option_index;
  static struct option long_options[] = {{"box-length", required_argument, nullptr, 'b'},
                                         {"container", required_argument, nullptr, 'c'},
                                         {"selector-strategy", required_argument, nullptr, 'y'},
                                         {"cutoff", required_argument, nullptr, 'C'},
                                         {"cell-size-factor", required_argument, nullptr, 'a'},
                                         {"distribution-mean", required_argument, nullptr, 'm'},
                                         {"distribution-stddeviation", required_argument, nullptr, 'z'},
                                         {"data-layout", required_argument, nullptr, 'd'},
                                         {"functor", required_argument, nullptr, 'f'},
                                         {"help", no_argument, nullptr, 'h'},
                                         {"iterations", required_argument, nullptr, 'i'},
                                         {"no-flops", no_argument, nullptr, 'F'},
                                         {"newton3", required_argument, nullptr, '3'},
                                         {"particles-generator", required_argument, nullptr, 'g'},
                                         {"particles-per-dimension", required_argument, nullptr, 'n'},
                                         {"particles-total", required_argument, nullptr, 'N'},
                                         {"particle-spacing", required_argument, nullptr, 's'},
                                         {"traversal", required_argument, nullptr, 't'},
                                         {"tuning-acquisition-function", required_argument, nullptr, 'A'},
                                         {"tuning-interval", required_argument, nullptr, 'I'},
                                         {"tuning-samples", required_argument, nullptr, 'S'},
                                         {"tuning-max-evidence", required_argument, nullptr, 'E'},
                                         {"tuning-strategy", required_argument, nullptr, 'T'},
                                         {"log-level", required_argument, nullptr, 'l'},
                                         {"log-file", required_argument, nullptr, 'L'},
                                         {"verlet-rebuild-frequency", required_argument, nullptr, 'v'},
                                         {"verlet-cluster-size", required_argument, nullptr, 'q'},
                                         {"precision", required_argument, nullptr, 'p'},
                                         {"verlet-skin-radius", required_argument, nullptr, 'r'},
                                         {"vtk", required_argument, nullptr, 'w'},
                                         {nullptr, 0, nullptr, 0}};  // needed to signal the end of the array
  std::string strArg;
  while ((option = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    if (optarg != nullptr) strArg = optarg;
    transform(strArg.begin(), strArg.end(), strArg.begin(), ::tolower);
    switch (option) {
      case '3': {
        newton3Options = autopas::Newton3Option::parseOptions(strArg);
        if (newton3Options.empty()) {
          cerr << "Unknown Newton3 option: " << strArg << endl;
          cerr << "Please use 'enabled' or 'disabled'!" << endl;
          displayHelp = true;
        }
        break;
      }
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
        // overwrite default argument
        containerOptions = autopas::ContainerOption::parseOptions(strArg);
        if (containerOptions.empty()) {
          cerr << "Unknown container option: " << strArg << endl;
          cerr << "Please use 'DirectSum', 'LinkedCells', 'VerletLists', 'VCells' or 'VCluster'!" << endl;
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
      case 'a': {
        cellSizeFactors = autopas::utils::StringUtils::parseNumberSet(strArg);
        if (cellSizeFactors->isEmpty()) {
          cerr << "Error parsing cell size factors: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'd': {
        dataLayoutOptions = autopas::DataLayoutOption::parseOptions(strArg);
        if (dataLayoutOptions.empty()) {
          cerr << "Unknown data layouts: " << strArg << endl;
          cerr << "Please use 'AoS' or 'SoA'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'f': {
        if (strArg.find("avx") != string::npos) {
          functorOption = lj12_6_AVX;
        } else if (strArg.find("lj") != string::npos || strArg.find("lennard-jones") != string::npos) {
          if (strArg.find("lobal") != string::npos) {
            functorOption = lj12_6_Globals;
          } else {
            functorOption = lj12_6;
          }
        } else {
          cerr << "Unknown functor: " << strArg << endl;
          cerr << "Please use 'Lennard-Jones' or 'Lennard-Jones-AVX'" << endl;
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
          cerr << "Unknown generator: " << strArg << endl;
          cerr << "Please use 'Grid' or 'Gaussian'" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'p': {
        if (strArg.find("float") != string::npos) {
          precisionOption = PrecisionOption::FP32;
        } else if (strArg.find("double") != string::npos) {
          precisionOption = PrecisionOption::FP64;
        } else {
          cerr << "Unknown precision: " << strArg << endl;
          cerr << "Please use 'float' or 'double'" << endl;
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
      case 'y': {
        auto parsedStrategies = autopas::SelectorStrategyOption::parseOptions(strArg);
        if (parsedStrategies.size() != 1) {
          cerr << "Please provide exactly one selector strategy. Parsed strategies are: "
               << iterableToString(parsedStrategies) << endl;
          displayHelp = true;
        }
        selectorStrategy = *parsedStrategies.begin();
        if (selectorStrategy == autopas::SelectorStrategyOption()) {
          cerr << "Unknown Selector Strategy: " << strArg << endl;
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
            cerr << "Unknown Log Level: " << strArg << endl;
            cerr << "Please use 'trace', 'debug', 'info', 'warning', 'error', 'critical' or 'off'." << endl;
            displayHelp = true;
          }
        }
        break;
      }
      case 'L': {
        logFileName = strArg;
        autopas::TuningStrategyOption::parseOptions(strArg);
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
      case 'E': {
        try {
          tuningMaxEvidence = (unsigned int)stoul(strArg);
          if (tuningMaxEvidence < 1) {
            cerr << "Tuning max evidence has to be a positive integer!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing number of tuning max evidence: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'A': {
        auto parsedAcquisitionFunctions = autopas::AcquisitionFunctionOption::parseOptions(strArg);
        if (parsedAcquisitionFunctions.size() != 1) {
          cerr << "Please provide exactly one acquisition function. Parsed strategies are: "
               << iterableToString(parsedAcquisitionFunctions) << endl;
          displayHelp = true;
        }
        tuningAcqFun = *parsedAcquisitionFunctions.begin();
        if (tuningAcqFun == autopas::AcquisitionFunctionOption()) {
          cerr << "Unknown acquisition function: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 't': {
        traversalOptions = autopas::TraversalOption::parseOptions(strArg);
        if (traversalOptions.empty()) {
          cerr << "Unknown Traversal: " << strArg << endl;
          cerr << "Please use 'c08', 'c01', 'c18', 'sliced' or 'direct'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'T': {
        auto parsedTuningStrategies = autopas::TuningStrategyOption::parseOptions(strArg);
        if (parsedTuningStrategies.size() != 1) {
          cerr << "Please provide exactly one tuning strategy. Parsed strategies are: "
               << iterableToString(parsedTuningStrategies) << endl;
          displayHelp = true;
        }
        tuningStrategyOption = *parsedTuningStrategies.begin();
        if (tuningStrategyOption == autopas::TuningStrategyOption()) {
          cerr << "Unknown Tuning Strategy: " << strArg << endl;
          cerr << "Please use 'full-search', 'bayesian-search' or 'active-harmony'!" << endl;
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
      case 'q': {
        try {
          verletClusterSize = (unsigned int)stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet-cluster-size: " << optarg << endl;
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

void MDFlexParser::printConfig() {
  constexpr size_t valueOffset = 32;
  cout << setw(valueOffset) << left << "Container"
       << ":  " << iterableToString(containerOptions) << endl;

  // if verlet lists are in the container options print verlet config data
  if (iterableToString(containerOptions).find("erlet") != std::string::npos) {
    cout << setw(valueOffset) << left << "Verlet rebuild frequency"
         << ":  " << verletRebuildFrequency << endl;

    cout << setw(valueOffset) << left << "Verlet skin radius"
         << ":  " << verletSkinRadius << endl;

    cout << setw(valueOffset) << left << "Verlet cluster size"
         << ":  " << verletClusterSize << endl;
  }

  if (containerOptions.size() > 1 or traversalOptions.size() > 1 or dataLayoutOptions.size() > 1) {
    cout << setw(valueOffset) << left << "Selector Strategy"
         << ":  " << selectorStrategy.to_string() << endl;
  }

  cout << setw(valueOffset) << left << "Data Layout"
       << ":  " << iterableToString(dataLayoutOptions) << endl;

  cout << setw(valueOffset) << left << "Functor"
       << ":  ";
  switch (functorOption) {
    case FunctorOption::lj12_6: {
      cout << "Lennard-Jones (12-6)" << endl;
      break;
    }
    case FunctorOption::lj12_6_AVX: {
      cout << "Lennard-Jones (12-6) AVX intrinsics" << endl;
      break;
    }
    case FunctorOption::lj12_6_Globals: {
      cout << "Lennard-Jones-globals (12-6)" << endl;
      break;
    }
  }

  cout << setw(valueOffset) << left << "Newton3"
       << ":  " << iterableToString(newton3Options) << endl;

  cout << setw(valueOffset) << left << "Cutoff radius"
       << ":  " << cutoff << endl;

  cout << setw(valueOffset) << left << "Cell size factor"
       << ":  " << static_cast<std::string>(*cellSizeFactors) << endl;

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
  switch (precisionOption) {
    case MDFlexParser::PrecisionOption::FP32: {
      cout << setw(valueOffset) << left << "Particle precision"
           << ":  "
           << "float" << endl;
      break;
    }
    case MDFlexParser::PrecisionOption::FP64: {
      cout << setw(valueOffset) << left << "Particle precision"
           << ":  "
           << "double" << endl;
      break;
    }
  }

  cout << setw(valueOffset) << left << "Allowed traversals"
       << ":  " << iterableToString(traversalOptions) << endl;

  cout << setw(valueOffset) << left << "Iterations"
       << ":  " << iterations << endl;
  cout << setw(valueOffset) << left << "Tuning Strategy"
       << ":  " << tuningStrategyOption.to_string() << endl;
  cout << setw(valueOffset) << left << "Tuning Interval"
       << ":  " << tuningInterval << endl;
  cout << setw(valueOffset) << left << "Tuning Samples"
       << ":  " << tuningSamples << endl;
  cout << setw(valueOffset) << left << "Tuning Max evidence"
       << ":  " << tuningMaxEvidence << endl;
  cout << setw(valueOffset) << left << "Tuning Acquisition Function"
       << ":  " << tuningAcqFun.to_string() << endl;
}

std::set<autopas::ContainerOption> MDFlexParser::getContainerOptions() const { return containerOptions; }

double MDFlexParser::getCutoff() const { return cutoff; }

const autopas::NumberSet<double> &MDFlexParser::getCellSizeFactors() const { return *cellSizeFactors; }

std::set<autopas::DataLayoutOption> MDFlexParser::getDataLayoutOptions() const { return dataLayoutOptions; }

MDFlexParser::FunctorOption MDFlexParser::getFunctorOption() const { return functorOption; }

size_t MDFlexParser::getIterations() const { return iterations; }

size_t MDFlexParser::getParticlesPerDim() const { return particlesPerDim; }

double MDFlexParser::getParticleSpacing() const { return particleSpacing; }

size_t MDFlexParser::getParticlesTotal() const { return particlesTotal; }

const std::set<autopas::TraversalOption> &MDFlexParser::getTraversalOptions() const { return traversalOptions; }

unsigned int MDFlexParser::getVerletRebuildFrequency() const { return verletRebuildFrequency; }

unsigned int MDFlexParser::getVerletClusterSize() const { return verletClusterSize; };

double MDFlexParser::getVerletSkinRadius() const { return verletSkinRadius; }

MDFlexParser::GeneratorOption MDFlexParser::getGeneratorOption() const { return generatorOption; }

MDFlexParser::PrecisionOption MDFlexParser::getPrecisionOption() const { return precisionOption; }

double MDFlexParser::getDistributionMean() const { return distributionMean; }

double MDFlexParser::getDistributionStdDev() const { return distributionStdDev; }

std::string MDFlexParser::getWriteVTK() const { return writeVTK; }

double MDFlexParser::getBoxLength() {
  if (boxLength == -1) boxLength = ceil(2 * distributionMean);
  return boxLength;
}

bool MDFlexParser::getMeasureFlops() const { return measureFlops; }

unsigned int MDFlexParser::getTuningInterval() const { return tuningInterval; }

unsigned int MDFlexParser::getTuningSamples() const { return tuningSamples; }

unsigned int MDFlexParser::getTuningMaxEvidence() const { return tuningMaxEvidence; }

autopas::AcquisitionFunctionOption MDFlexParser::getAcquisitionFunctionOption() const { return tuningAcqFun; }

autopas::Logger::LogLevel MDFlexParser::getLogLevel() const { return logLevel; }

autopas::SelectorStrategyOption MDFlexParser::getSelectorStrategy() const { return selectorStrategy; }

std::set<autopas::Newton3Option> MDFlexParser::getNewton3Options() const { return newton3Options; }

const string &MDFlexParser::getLogFileName() const { return logFileName; }

autopas::TuningStrategyOption MDFlexParser::getTuningStrategyOption() const { return tuningStrategyOption; }
