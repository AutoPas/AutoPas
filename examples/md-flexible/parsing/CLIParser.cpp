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
  static struct option long_options[] = {{"help", no_argument, nullptr, 'h'},
                                         {MDFlexConfig::boxLengthStr, required_argument, nullptr, 'b'},
                                         {MDFlexConfig::cellSizeFactorsStr, required_argument, nullptr, 'a'},
                                         {MDFlexConfig::containerOptionsStr, required_argument, nullptr, 'c'},
                                         {MDFlexConfig::cutoffStr, required_argument, nullptr, 'C'},
                                         {MDFlexConfig::dataLayoutOptionsStr, required_argument, nullptr, 'd'},
                                         {MDFlexConfig::deltaTStr, required_argument, nullptr, 'D'},
                                         {MDFlexConfig::distributionMeanStr, required_argument, nullptr, 'm'},
                                         {MDFlexConfig::distributionStdDevStr, required_argument, nullptr, 'z'},
                                         {MDFlexConfig::functorOptionStr, required_argument, nullptr, 'f'},
                                         {MDFlexConfig::generatorOptionStr, required_argument, nullptr, 'g'},
                                         {MDFlexConfig::iterationsStr, required_argument, nullptr, 'i'},
                                         {MDFlexConfig::logFileNameStr, required_argument, nullptr, 'L'},
                                         {MDFlexConfig::logLevelStr, required_argument, nullptr, 'l'},
                                         {MDFlexConfig::measureFlopsStr, no_argument, nullptr, 'F'},
                                         {MDFlexConfig::newton3OptionsStr, required_argument, nullptr, '3'},
                                         {MDFlexConfig::particlesPerDimStr, required_argument, nullptr, 'n'},
                                         {MDFlexConfig::particlesSpacingStr, required_argument, nullptr, 's'},
                                         {MDFlexConfig::particlesTotalStr, required_argument, nullptr, 'N'},
                                         {MDFlexConfig::periodicStr, required_argument, nullptr, 'p'},
                                         {MDFlexConfig::selectorStrategyStr, required_argument, nullptr, 'y'},
                                         {MDFlexConfig::thermostatStr, required_argument, nullptr, 'u'},
                                         {MDFlexConfig::traversalOptionsStr, required_argument, nullptr, 't'},
                                         {MDFlexConfig::tuningIntervalStr, required_argument, nullptr, 'I'},
                                         {MDFlexConfig::tuningMaxEvidenceStr, required_argument, nullptr, 'E'},
                                         {MDFlexConfig::tuningSamplesStr, required_argument, nullptr, 'S'},
                                         {MDFlexConfig::tuningStrategyOptionsStr, required_argument, nullptr, 'T'},
                                         {MDFlexConfig::verletRebuildFrequencyStr, required_argument, nullptr, 'v'},
                                         {MDFlexConfig::verletSkinRadiusStr, required_argument, nullptr, 'r'},
                                         {MDFlexConfig::vtkFileNameStr, required_argument, nullptr, 'w'},
                                         {MDFlexConfig::vtkWriteFrequencyStr, required_argument, nullptr, 'W'},
                                         {MDFlexConfig::yamlFilenameStr, required_argument, nullptr, 'Y'},
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
      case 'b': {
        try {
          config.boxLength = stod(strArg);
          if (config.boxLength < 0) {
            cerr << "Box length has to be a positive (floating point) number!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing number of the box length: " << optarg << endl;
          displayHelp = true;
        }
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
        } else if (strArg.find("sp") != string::npos) {
          config.generatorOption = MDFlexConfig::GeneratorOption::sphere;
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
      case 'N': {
        try {
          config.particlesTotal = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing total number of particles: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'P': {
        try {
          config.particlesTotal = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing total number of particles: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'p': {
        try {
          config.periodic = autopas::utils::StringUtils::parseBoolOption(strArg);
        } catch (const exception &) {
          cerr << "Error parsing whether there should be periodic boundary conditions: " << strArg << endl;
          displayHelp = true;
        }
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
        config.useThermostat = autopas::utils::StringUtils::parseBoolOption(strArg);
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
        config.vtkFileName = strArg;
        break;
      }
      case 'W': {
        try {
          config.vtkWriteFrequency = (size_t)stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing vtk write frequency: " << optarg << endl;
          displayHelp = true;
        }
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
        try {
          config.vtkWriteFrequency = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet-rebuild-frequency: " << optarg << endl;
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

  // only create objects if nothing was set by a yaml file
  if (config.cubeGaussObjects.empty() and config.cubeGridObjects.empty() and config.cubeUniformObjects.empty() and
      config.sphereObjects.empty()) {
    switch (config.generatorOption) {
      case MDFlexConfig::GeneratorOption::grid: {
        CubeGrid grid({config.particlesPerDim, config.particlesPerDim, config.particlesPerDim}, config.particleSpacing,
                      {0, 0, 0}, {0, 0, 0}, 0, 1, 1, 1);
        config.cubeGridObjects.push_back(grid);
        break;
      }
      case MDFlexConfig::GeneratorOption::gaussian: {
        CubeGauss cubeGauss(config.particlesTotal, {config.boxLength, config.boxLength, config.boxLength},
                            config.distributionMean, config.distributionStdDev, {0, 0, 0}, {0, 0, 0}, 0, 1, 1, 1);
        config.cubeGaussObjects.push_back(cubeGauss);
        break;
      }
      case MDFlexConfig::GeneratorOption::uniform: {
        CubeUniform cubeUniform(config.particlesTotal, {config.boxLength, config.boxLength, config.boxLength},
                                {0, 0, 0}, {0, 0, 0}, 0, 1, 1, 1);
        config.cubeUniformObjects.push_back(cubeUniform);
        break;
      }
      case MDFlexConfig::GeneratorOption::sphere: {
        auto centerOfBox = config.particlesPerDim / 2.;
        Sphere sphere({centerOfBox, centerOfBox, centerOfBox}, centerOfBox, config.particleSpacing, {0, 0, 0}, 0, 1, 1,
                      1);
        config.sphereObjects.push_back(sphere);
        break;
      }
    }
  }

  if (displayHelp) {
    cout << "Usage: " << argv[0] << endl;
    for (auto o : long_options) {
      if (o.name == nullptr) {
        continue;
      }
      cout << "    --" << setw(MDFlexConfig::valueOffset + 2) << left << o.name;
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
  static struct option longOptions[] = {{MDFlexConfig::yamlFilenameStr, required_argument, nullptr, 'Y'},
                                        {nullptr, 0, nullptr, 0}};  // needed to signal the end of the array
  std::string strArg;
    optind =1;
    // Yaml Parsing file parameter must be set before all other Options
//  option = getopt_long(argc, argv, "", longOptions, &optionIndex);
  while ((option = getopt_long(argc, argv, "", longOptions, &optionIndex)) != -1) {
    if (option == 'Y') {
      config.yamlFilename = optarg;
      return true;
    }
  }

  opterr = opterrBefore;
  return false;
}
