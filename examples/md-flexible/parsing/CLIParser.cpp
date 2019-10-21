/**
 * @file CLIParser.cpp
 * @author F. Gratl
 * @date 10/18/19
 */

#include "CLIParser.h"

bool CLIParser::parseInput(int argc, char **argv, MDFlexConfig &config) {
  using namespace std;
  bool displayHelp = false;
  int option, option_index;
  static struct option long_options[] = {{"yaml-filename", required_argument, nullptr, 'Y'},
                                         {"box-min", required_argument, nullptr, 'k'},
                                         {"box-max", required_argument, nullptr, 'K'},
                                         {"container", required_argument, nullptr, 'c'},
                                         {"cutoff", required_argument, nullptr, 'C'},
                                         {"cell-size-factor", required_argument, nullptr, 'a'},
                                         {"data-layout", required_argument, nullptr, 'd'},
                                         {"distribution-mean", required_argument, nullptr, 'm'},
                                         {"distribution-stddeviation", required_argument, nullptr, 'z'},
                                         {"deltaT", required_argument, nullptr, 'D'},
                                         {"functor", required_argument, nullptr, 'f'},
                                         {"help", no_argument, nullptr, 'h'},
                                         {"iterations", required_argument, nullptr, 'i'},
                                         {"no-flops", no_argument, nullptr, 'F'},
                                         {"newton3", required_argument, nullptr, '3'},
                                         {"particles-generator", required_argument, nullptr, 'g'},
                                         {"particles-per-dimension", required_argument, nullptr, 'n'},
                                         {"particles-total", required_argument, nullptr, 'N'},
                                         {"particle-spacing", required_argument, nullptr, 's'},
                                         {"periodic", required_argument, nullptr, 'p'},
                                         {"selector-strategy", required_argument, nullptr, 'y'},
                                         {"thermostat", required_argument, nullptr, 'u'},
                                         {"traversal", required_argument, nullptr, 't'},
                                         {"tuning-interval", required_argument, nullptr, 'I'},
                                         {"tuning-samples", required_argument, nullptr, 'S'},
                                         {"tuning-max-evidence", required_argument, nullptr, 'E'},
                                         {"tuning-strategy", required_argument, nullptr, 'T'},
                                         {"log-level", required_argument, nullptr, 'l'},
                                         {"log-file", required_argument, nullptr, 'L'},
                                         {"verlet-rebuild-frequency", required_argument, nullptr, 'v'},
                                         {"verlet-skin-radius", required_argument, nullptr, 'r'},
                                         {"vtk-filename", required_argument, nullptr, 'w'},
                                         {"vtk-write-frequency", required_argument, nullptr, 'Z'},
                                         {nullptr, no_argument, nullptr, 0}};  // needed to signal the end of the array
  // reset getopt to scan from the start of argv
  optind = 1;
  string strArg;
  while ((option = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    if (optarg != nullptr) strArg = optarg;
    transform(strArg.begin(), strArg.end(), strArg.begin(), ::tolower);
    switch (option) {
      case '3': {
        config.newton3Options = autopas::utils::StringUtils::parseNewton3Options(strArg, false);
        if (config.newton3Options.empty()) {
          cerr << "Unknown Newton3 option: " << strArg << endl;
          cerr << "Please use 'enabled' or 'disabled'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'a': {
        config.cellSizeFactors = autopas::utils::StringUtils::parseNumberSet(strArg);
        if (config.cellSizeFactors->isEmpty()) {
          cerr << "Error parsing cell size factors: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'c': {
        // overwrite default argument
        config.containerOptions = autopas::utils::StringUtils::parseContainerOptions(strArg, false);
        if (config.containerOptions.empty()) {
          cerr << "Unknown container option: " << strArg << endl;
          cerr << "Please use 'DirectSum', 'LinkedCells', 'VerletLists', 'VCells' or 'VCluster'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'C': {
        try {
          config.cutoff = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing cutoff Radius: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'D': {
        try {
          config.deltaT = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing epsilon value: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'd': {
        config.dataLayoutOptions = autopas::utils::StringUtils::parseDataLayout(strArg);
        if (config.dataLayoutOptions.empty()) {
          cerr << "Unknown data layouts: " << strArg << endl;
          cerr << "Please use 'AoS' or 'SoA'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'E': {
        try {
          config.tuningMaxEvidence = (unsigned int)stoul(strArg);
          if (config.tuningMaxEvidence < 1) {
            cerr << "Tuning max evidence has to be a positive integer!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing number of tuning max evidence: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'f': {
        if (strArg.find("avx") != string::npos) {
          config.functorOption = MDFlexConfig::FunctorOption::lj12_6_AVX;
        } else if (strArg.find("lj") != string::npos || strArg.find("lennard-jones") != string::npos) {
          config.functorOption = MDFlexConfig::FunctorOption::lj12_6;
        } else {
          cerr << "Unknown functor: " << strArg << endl;
          cerr << "Please use 'Lennard-Jones' or 'Lennard-Jones-AVX'" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'F': {
        config.measureFlops = false;
        break;
      }
      case 'g': {
        if (strArg.find("grid") != string::npos) {
          config.generatorOption = MDFlexConfig::GeneratorOption::grid;
        } else if (strArg.find("uni") != string::npos) {
          config.generatorOption = MDFlexConfig::GeneratorOption::uniform;
        } else if (strArg.find("gaus") != string::npos) {
          config.generatorOption = MDFlexConfig::GeneratorOption::gaussian;
        } else {
          cerr << "Unknown generator: " << strArg << endl;
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
          config.iterations = stoul(strArg);
          if (config.iterations < 1) {
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
          config.tuningInterval = (unsigned int)stoul(strArg);
          if (config.tuningInterval < 1) {
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
        try {
          config.boxMin = autopas::utils::StringUtils::parseBoxOption(strArg);
        } catch (const exception &) {
          cerr << "Error parsing BoxMinOption: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'K': {
        try {
          config.boxMax = autopas::utils::StringUtils::parseBoxOption(strArg);
        } catch (const exception &) {
          cerr << "Error parsing BoxMaxOption: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'l': {
        switch (strArg[0]) {
          case 't': {
            config.logLevel = autopas::Logger::LogLevel::trace;
            break;
          }
          case 'd': {
            config.logLevel = autopas::Logger::LogLevel::debug;
            break;
          }
          case 'i': {
            config.logLevel = autopas::Logger::LogLevel::info;
            break;
          }
          case 'w': {
            config.logLevel = autopas::Logger::LogLevel::warn;
            break;
          }
          case 'e': {
            config.logLevel = autopas::Logger::LogLevel::err;
            break;
          }
          case 'c': {
            config.logLevel = autopas::Logger::LogLevel::critical;
            break;
          }
          case 'o': {
            config.logLevel = autopas::Logger::LogLevel::off;
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
        config.logFileName = strArg;
        break;
      }
      case 'm': {
        try {
          config.distributionMean = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing distribution mean: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'n': {
        try {
          config.particlesPerDim = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing number of particles per dimension: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'P': {
        try {
          config.defaultParticlesTotal = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing total number of particles: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'p': {
        config.periodic = true;
        break;
      }
      case 'r': {
        try {
          config.verletSkinRadius = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet-skin-radius: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'S': {
        try {
          config.tuningSamples = (unsigned int)stoul(strArg);
          if (config.tuningSamples < 1) {
            cerr << "Tuning samples has to be a positive integer!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing number of tuning samples: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 's': {
        try {
          config.particleSpacing = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing separation of particles: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 't': {
        config.traversalOptions = autopas::utils::StringUtils::parseTraversalOptions(strArg);
        if (config.traversalOptions.empty()) {
          cerr << "Unknown Traversal: " << strArg << endl;
          cerr << "Please use 'c08', 'c01', 'c18', 'sliced' or 'direct'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'T': {
        config.tuningStrategyOption = autopas::utils::StringUtils::parseTuningStrategyOption(strArg);
        if (config.tuningStrategyOption == autopas::TuningStrategyOption(-1)) {
          cerr << "Unknown Tuning Strategy: " << strArg << endl;
          cerr << "Please use 'full-search' or 'bayesian-search'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'u': {
        config.thermostat = autopas::utils::StringUtils::parseBoolOption(strArg);
        break;
      }
      case 'v': {
        try {
          config.verletRebuildFrequency = (unsigned int)stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet-rebuild-frequency: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'w': {
        config.VTKFileName = strArg;
        break;
      }
      case 'y': {
        config.selectorStrategy = autopas::utils::StringUtils::parseSelectorStrategy(strArg);
        if (config.selectorStrategy == autopas::SelectorStrategyOption(-1)) {
          cerr << "Unknown Selector Strategy: " << strArg << endl;
          cerr << "Please use 'fastestAbs', 'fastestMean' or 'fastestMedian'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'Y': {
        // already parsed in CLIParser::yamlFilePresent
        break;
      }
      case 'z': {
        config.vtkWriteFrequency = stoul(strArg);
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
      if (o.name == nullptr) {
        continue;
      }
      cout << "    --" << setw(config.valueOffset + 2) << left << o.name;
      if (o.has_arg) {
        cout << "option";
      }
      cout << endl;
    }
    return false;
  }
  return true;
}

bool CLIParser::yamlFilePresent(int argc, char **argv, MDFlexConfig &config) {
  int option, optionIndex;
  // suppress error messages since we only want to look if the yaml option is there
  auto opterrBefore = opterr;
  std::vector<std::string> argvBefore(argv, argv + argc);
  opterr = 0;
  static struct option longOptions[] = {{"yaml-filename", required_argument, nullptr, 'Y'},
                                        {nullptr, 0, nullptr, 0}};  // needed to signal the end of the array
  std::string strArg;
  // Yaml Parsing file parameter must be set before all other Options
  option = getopt_long(argc, argv, "", longOptions, &optionIndex);
  while ((option = getopt_long(argc, argv, "", longOptions, &optionIndex)) != -1) {
    if (option == 'Y') {
      config.yamlFilename = optarg;
      return true;
    }
  }

  opterr = opterrBefore;
  return false;
}

