/**
 * @file CLIParser.cpp
 * @author F. Gratl
 * @date 18.10.2019
 */

#include "CLIParser.h"

#include <sys/stat.h>

#include "autopas/utils/StringUtils.h"

bool CLIParser::parseInput(int argc, char **argv, MDFlexConfig &config) {
  using namespace std;
  bool displayHelp = false;
  int option, option_index;
  static struct option long_options[] = {{"help", no_argument, nullptr, 'h'},
                                         {MDFlexConfig::newton3OptionsStr, required_argument, nullptr, '3'},
                                         {MDFlexConfig::checkpointfileStr, required_argument, nullptr, '4'},
                                         {MDFlexConfig::acquisitionFunctionOptionStr, required_argument, nullptr, 'A'},
                                         {MDFlexConfig::cellSizeFactorsStr, required_argument, nullptr, 'a'},
                                         {MDFlexConfig::boxLengthStr, required_argument, nullptr, 'b'},
                                         {MDFlexConfig::containerOptionsStr, required_argument, nullptr, 'c'},
                                         {MDFlexConfig::cutoffStr, required_argument, nullptr, 'C'},
                                         {MDFlexConfig::dataLayoutOptionsStr, required_argument, nullptr, 'd'},
                                         {MDFlexConfig::deltaTStr, required_argument, nullptr, 'D'},
                                         {MDFlexConfig::dontCreateEndConfigStr, no_argument, nullptr, 'e'},
                                         {MDFlexConfig::tuningMaxEvidenceStr, required_argument, nullptr, 'E'},
                                         {MDFlexConfig::functorOptionStr, required_argument, nullptr, 'f'},
                                         {MDFlexConfig::dontMeasureFlopsStr, no_argument, nullptr, 'F'},
                                         {MDFlexConfig::generatorOptionStr, required_argument, nullptr, 'g'},
                                         {MDFlexConfig::iterationsStr, required_argument, nullptr, 'i'},
                                         {MDFlexConfig::tuningIntervalStr, required_argument, nullptr, 'I'},
                                         {MDFlexConfig::logLevelStr, required_argument, nullptr, 'l'},
                                         {MDFlexConfig::logFileNameStr, required_argument, nullptr, 'L'},
                                         {MDFlexConfig::distributionMeanStr, required_argument, nullptr, 'm'},
                                         {MDFlexConfig::maxTuningPhasesWithoutTestStr, required_argument, nullptr, 'M'},
                                         {MDFlexConfig::particlesPerDimStr, required_argument, nullptr, 'n'},
                                         {MDFlexConfig::particlesTotalStr, required_argument, nullptr, 'N'},
                                         {MDFlexConfig::relativeOptimumRangeStr, required_argument, nullptr, 'o'},
                                         {MDFlexConfig::periodicStr, required_argument, nullptr, 'p'},
                                         {MDFlexConfig::tuningPhasesStr, required_argument, nullptr, 'P'},
                                         {MDFlexConfig::verletClusterSizeStr, required_argument, nullptr, 'q'},
                                         {MDFlexConfig::verletSkinRadiusStr, required_argument, nullptr, 'r'},
                                         {MDFlexConfig::particlesSpacingStr, required_argument, nullptr, 's'},
                                         {MDFlexConfig::tuningSamplesStr, required_argument, nullptr, 'S'},
                                         {MDFlexConfig::traversalOptionsStr, required_argument, nullptr, 't'},
                                         {MDFlexConfig::tuningStrategyOptionsStr, required_argument, nullptr, 'T'},
                                         {MDFlexConfig::thermostatStr, required_argument, nullptr, 'u'},
                                         {MDFlexConfig::verletRebuildFrequencyStr, required_argument, nullptr, 'v'},
                                         {MDFlexConfig::vtkFileNameStr, required_argument, nullptr, 'w'},
                                         {MDFlexConfig::vtkWriteFrequencyStr, required_argument, nullptr, 'W'},
                                         {MDFlexConfig::selectorStrategyStr, required_argument, nullptr, 'y'},
                                         {MDFlexConfig::yamlFilenameStr, required_argument, nullptr, 'Y'},
                                         {MDFlexConfig::distributionStdDevStr, required_argument, nullptr, 'z'},
                                         {nullptr, no_argument, nullptr, 0}};  // needed to signal the end of the array
  // reset getopt to scan from the start of argv
  optind = 1;
  string strArg;
  while ((option = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    if (optarg != nullptr) strArg = optarg;
    transform(strArg.begin(), strArg.end(), strArg.begin(), ::tolower);
    switch (option) {
      case '3': {
        config.newton3Options = autopas::Newton3Option::parseOptions(strArg);
        if (config.newton3Options.empty()) {
          cerr << "Unknown Newton3 option: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case '4': {
        // already parsed in CLIParser::inputFilesPresent
        break;
      }
      case 'A': {
        auto parsedOptions = autopas::AcquisitionFunctionOption::parseOptions(strArg);
        if (parsedOptions.size() != 1) {
          cerr << "Pass exactly one tuning acquisition function." << endl
               << "Passed: " << strArg << endl
               << "Parsed: " << autopas::utils::ArrayUtils::to_string(parsedOptions) << endl;
          displayHelp = true;
        }
        config.acquisitionFunctionOption = *parsedOptions.begin();
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
        break;
      }
      case 'c': {
        config.containerOptions = autopas::ContainerOption::parseOptions(strArg);
        if (config.containerOptions.empty()) {
          cerr << "Unknown container option: " << strArg << endl;
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
        config.dataLayoutOptions = autopas::DataLayoutOption::parseOptions(strArg);
        if (config.dataLayoutOptions.empty()) {
          cerr << "Unknown data layouts: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'e': {
        config.dontCreateEndConfig = false;
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
        } else if (strArg.find("glob") != string::npos) {
          config.functorOption = MDFlexConfig::FunctorOption::lj12_6_Globals;
        } else if (strArg.find("lj") != string::npos || strArg.find("lennard-jones") != string::npos) {
          config.functorOption = MDFlexConfig::FunctorOption::lj12_6;
        } else {
          cerr << "Unknown functor: " << strArg << endl;
          cerr << "Please use 'Lennard-Jones', 'Lennard-Jones-With-Globals' or 'Lennard-Jones-AVX'" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'F': {
        config.dontMeasureFlops = false;
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
          auto mean = stod(strArg);
          config.distributionMean = {mean, mean, mean};
        } catch (const exception &) {
          cerr << "Error parsing distribution mean: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'M': {
        try {
          config.maxTuningPhasesWithoutTest = (unsigned int)stoul(strArg);
          if (config.maxTuningPhasesWithoutTest < 1) {
            cerr << "Max tuning phases without test has to be positive!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing max tuning phases without test: " << optarg << endl;
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
      case 'o': {
        try {
          config.relativeOptimumRange = (double)stoul(strArg);
          if (config.relativeOptimumRange < 1) {
            cerr << "Relative optimum range has to be greater or equal one!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing relative optimum range: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'P': {
        try {
          config.tuningPhases = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing number of tuning phases: " << strArg << endl;
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
      case 'q': {
        try {
          config.verletClusterSize = (unsigned int)stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet cluster size: " << optarg << endl;
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
        config.traversalOptions = autopas::TraversalOption::parseOptions(strArg);
        if (config.traversalOptions.empty()) {
          cerr << "Unknown Traversal: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'T': {
        auto parsedOptions = autopas::TuningStrategyOption::parseOptions(strArg);
        if (parsedOptions.size() != 1) {
          cerr << "Pass exactly one tuning strategy option." << endl
               << "Passed: " << strArg << endl
               << "Parsed: " << autopas::utils::ArrayUtils::to_string(parsedOptions) << endl;
          displayHelp = true;
        }
        config.tuningStrategyOption = *parsedOptions.begin();
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
        auto parsedOptions = autopas::SelectorStrategyOption::parseOptions(strArg);
        if (parsedOptions.size() != 1) {
          cerr << "Pass exactly one selector strategy option." << endl
               << "Passed: " << strArg << endl
               << "Parsed: " << autopas::utils::ArrayUtils::to_string(parsedOptions) << endl;
          displayHelp = true;
        }
        config.selectorStrategy = *parsedOptions.begin();
        break;
      }
      case 'Y': {
        // already parsed in CLIParser::inputFilesPresent
        break;
      }
      case 'z': {
        try {
          auto stdDev = stod(strArg);
          config.distributionStdDev = {stdDev, stdDev, stdDev};
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

  // only create objects if nothing was set by a yaml file and there was no checkpoint
  if (config.checkpointfile.empty() and config.cubeGaussObjects.empty() and config.cubeGridObjects.empty() and
      config.cubeUniformObjects.empty() and config.sphereObjects.empty()) {
    // common settings for any object type:
    unsigned int typeID = 0;
    double epsilon = 1.;
    double sigma = 1.;
    double mass = 1.;
    std::array<double, 3> bottomLeftCorner = {0, 0, 0};
    std::array<double, 3> velocity = {0, 0, 0};

    switch (config.generatorOption) {
      case MDFlexConfig::GeneratorOption::grid: {
        CubeGrid grid(velocity, typeID, epsilon, sigma, mass,
                      {config.particlesPerDim, config.particlesPerDim, config.particlesPerDim}, config.particleSpacing,
                      bottomLeftCorner);
        config.cubeGridObjects.push_back(grid);
        break;
      }
      case MDFlexConfig::GeneratorOption::gaussian: {
        CubeGauss cubeGauss(velocity, typeID, epsilon, sigma, mass, config.particlesTotal,
                            {config.boxLength, config.boxLength, config.boxLength}, config.distributionMean,
                            config.distributionStdDev, bottomLeftCorner);
        config.cubeGaussObjects.push_back(cubeGauss);
        break;
      }
      case MDFlexConfig::GeneratorOption::uniform: {
        CubeUniform cubeUniform(velocity, typeID, epsilon, sigma, mass, config.particlesTotal,
                                {config.boxLength, config.boxLength, config.boxLength}, bottomLeftCorner);
        config.cubeUniformObjects.push_back(cubeUniform);
        break;
      }
      case MDFlexConfig::GeneratorOption::sphere: {
        auto centerOfBox = config.particlesPerDim / 2.;
        Sphere sphere(velocity, typeID, epsilon, sigma, mass, {centerOfBox, centerOfBox, centerOfBox}, centerOfBox,
                      config.particleSpacing);
        config.sphereObjects.push_back(sphere);
        break;
      }
    }
  }

  if (displayHelp) {
    // filter out null values and copy rest in more sane data structure
    std::vector<std::pair<std::string, bool>> options;
    for (auto &o : long_options) {
      if (o.name != nullptr) {
        options.emplace_back(std::make_pair(o.name, o.has_arg));
      }
    }

    // by default sort sorts by first member of pair
    std::sort(std::begin(options), std::end(options));

    // print everything
    cout << "Usage: " << argv[0] << endl;
    for (auto &o : options) {
      cout << "    --" << setw(MDFlexConfig::valueOffset + 2) << left << o.first;
      if (o.second) {
        cout << "option";
      }
      cout << endl;
    }
    return false;
  }
  return true;
}

// anonymous namespace to hide helper function
namespace {

/**
 * Checks if a file with the given path exists.
 * @param filename
 * @return True iff the file exists.
 */
bool checkFileExists(const std::string &filename) {
  struct stat buffer;
  return (stat(filename.c_str(), &buffer) == 0);
}

}  // namespace

void CLIParser::inputFilesPresent(int argc, char **argv, MDFlexConfig &config) {
  int option = 0, optionIndex = 0;
  // suppress error messages since we only want to look if the yaml option is there
  auto opterrBefore = opterr;
  opterr = 0;
  static struct option longOptions[] = {{MDFlexConfig::checkpointfileStr, required_argument, nullptr, 'C'},
                                        {MDFlexConfig::yamlFilenameStr, required_argument, nullptr, 'Y'},
                                        {nullptr, 0, nullptr, 0}};  // needed to signal the end of the array
  std::string strArg;
  optind = 1;

  // search all cli parameters for input file options
  while ((option = getopt_long(argc, argv, "", longOptions, &optionIndex)) != -1) {
    switch (option) {
      case 'C':
        config.checkpointfile = optarg;
        if (not checkFileExists(optarg)) {
          throw std::runtime_error("CLIParser::inputFilesPresent: Checkpoint-File " + config.checkpointfile +
                                   " not found!");
        }
        break;
      case 'Y':
        config.yamlFilename = optarg;
        if (not checkFileExists(optarg)) {
          throw std::runtime_error("CLIParser::inputFilesPresent: Yaml-File " + config.yamlFilename + " not found!");
        }
        break;
      default: {
        // do nothing
      }
    }
  }

  opterr = opterrBefore;
}
