/**
 * @file CLIParser.cpp
 * @author F. Gratl
 * @date 18.10.2019
 */

#include "CLIParser.h"

#include <sys/stat.h>

#include <any>
#include <fstream>

MDFlexParser::exitCodes MDFlexParser::CLIParser::parseInput(int argc, char **argv, MDFlexConfig &config) {
  using namespace std;

  // the following, shorter version does not work with icpc 2019.4.243. Error:
  // error: class template name must be a placeholder for the complete type being initialized
  // (not for a component of that type)
  //
  // static const std::tuple relevantOptions{
  //
  // therefore workaround with make_tuple and auto
  static const auto relevantOptions{std::make_tuple(
      config.newton3Options, config.checkpointfile, config.acquisitionFunctionOption, config.cellSizeFactors,
      config.boxLength, config.containerOptions, config.cutoff, config.dataLayoutOptions, config.deltaT,
      config.dontCreateEndConfig, config.tuningMaxEvidence, config.functorOption, config.dontMeasureFlops,
      config.generatorOption, config.iterations, config.tuningInterval, config.logLevel, config.logFileName,
      config.distributionMean, config.maxTuningPhasesWithoutTest, config.particlesPerDim, config.particlesTotal,
      config.relativeOptimumRange, config.periodic, config.tuningPhases, config.verletClusterSize,
      config.verletSkinRadius, config.particleSpacing, config.tuningSamples, config.traversalOptions,
      config.tuningStrategyOption, config.useThermostat, config.verletRebuildFrequency, config.vtkFileName,
      config.vtkWriteFrequency, config.selectorStrategy, config.yamlFilename, config.distributionStdDev,
      MDFlexConfig::MDFlexOption<std::string, 'Z'>("", "zsh-completions", false, "Generate completions file for zsh."),
      MDFlexConfig::MDFlexOption<std::string, 'h'>("", "help", false, "Display this message."))};

  constexpr auto relevantOptionsSize = std::tuple_size_v<decltype(relevantOptions)>;

  // sanity check that all getopt chars are unique. Brackets for scoping.
  {
    // map tracking mappings of getopt chars to strings
    std::map<char, std::string> getoptCharsToName;
    // look for clashes by checking if getopt chars are in the map and otherwise add them
    autopas::utils::TupleUtils::for_each(relevantOptions, [&](auto &opt) {
      if (auto iterAtClash = getoptCharsToName.find(opt.getoptChar); iterAtClash != getoptCharsToName.end()) {
        throw std::runtime_error("CLIParser::parseInput: the following options share the same getopt char!\n" +
                                 opt.name + " : " + opt.getoptChar + "\n" + iterAtClash->second + " : " +
                                 iterAtClash->first);
      } else {
        getoptCharsToName.insert({opt.getoptChar, opt.name});
      }
    });
  }

  // create data structure for options that getopt can use
  std::vector<struct option> long_options;
  // reserve space for all relevant options and terminal field
  long_options.reserve(relevantOptionsSize + 1);

  autopas::utils::TupleUtils::for_each(relevantOptions,
                                       [&](auto &elem) { long_options.push_back(elem.toGetoptOption()); });

  // needed to signal the end of the array
  long_options.push_back({nullptr, no_argument, nullptr, 0});

  // reset getopt to scan from the start of argv
  optind = 1;
  bool displayHelp = false;
  for (int cliOption = 0, cliOptionIndex = 0;
       (cliOption = getopt_long(argc, argv, "", long_options.data(), &cliOptionIndex)) != -1;) {
    string strArg;
    if (optarg != nullptr) strArg = optarg;
    transform(strArg.begin(), strArg.end(), strArg.begin(), ::tolower);
    switch (cliOption) {
      case '3': {
        config.newton3Options.value = autopas::Newton3Option::parseOptions(strArg);
        if (config.newton3Options.value.empty()) {
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
        config.acquisitionFunctionOption.value = *parsedOptions.begin();
        break;
      }
      case 'a': {
        config.cellSizeFactors.value = autopas::utils::StringUtils::parseNumberSet(strArg);
        if (config.cellSizeFactors.value->isEmpty()) {
          cerr << "Error parsing cell size factors: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'b': {
        try {
          config.boxLength.value = stod(strArg);
          if (config.boxLength.value < 0) {
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
        config.containerOptions.value = autopas::ContainerOption::parseOptions(strArg);
        if (config.containerOptions.value.empty()) {
          cerr << "Unknown container option: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'C': {
        try {
          config.cutoff.value = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing cutoff Radius: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'D': {
        try {
          config.deltaT.value = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing epsilon value: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'd': {
        config.dataLayoutOptions.value = autopas::DataLayoutOption::parseOptions(strArg);
        if (config.dataLayoutOptions.value.empty()) {
          cerr << "Unknown data layouts: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'e': {
        config.dontCreateEndConfig.value = false;
        break;
      }
      case 'E': {
        try {
          config.tuningMaxEvidence.value = (unsigned int)stoul(strArg);
          if (config.tuningMaxEvidence.value < 1) {
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
          config.functorOption.value = MDFlexConfig::FunctorOption::lj12_6_AVX;
        } else if (strArg.find("glob") != string::npos) {
          config.functorOption.value = MDFlexConfig::FunctorOption::lj12_6_Globals;
        } else if (strArg.find("lj") != string::npos || strArg.find("lennard-jones") != string::npos) {
          config.functorOption.value = MDFlexConfig::FunctorOption::lj12_6;
        } else {
          cerr << "Unknown functor: " << strArg << endl;
          cerr << "Please use 'Lennard-Jones', 'Lennard-Jones-With-Globals' or 'Lennard-Jones-AVX'" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'F': {
        config.dontMeasureFlops.value = false;
        break;
      }
      case 'g': {
        if (strArg.find("grid") != string::npos) {
          config.generatorOption.value = MDFlexConfig::GeneratorOption::grid;
        } else if (strArg.find("uni") != string::npos) {
          config.generatorOption.value = MDFlexConfig::GeneratorOption::uniform;
        } else if (strArg.find("gaus") != string::npos) {
          config.generatorOption.value = MDFlexConfig::GeneratorOption::gaussian;
        } else if (strArg.find("sp") != string::npos) {
          config.generatorOption.value = MDFlexConfig::GeneratorOption::sphere;
        } else {
          cerr << "Unknown generator: " << strArg << endl;
          cerr << "Please use 'Grid' or 'Gaussian'" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'h': {
        printHelpMessage(std::cout, argv[0], relevantOptions);
        return MDFlexParser::exitCodes::helpFlagFound;
      }
      case 'Z': {
        // generate the completions file and do nothing else
        createZSHCompletionFile(relevantOptions);
        return MDFlexParser::exitCodes::completionsFlagFound;
      }
      case 'i': {
        try {
          config.iterations.value = stoul(strArg);
          if (config.iterations.value < 1) {
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
          config.tuningInterval.value = (unsigned int)stoul(strArg);
          if (config.tuningInterval.value < 1) {
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
            config.logLevel.value = autopas::Logger::LogLevel::trace;
            break;
          }
          case 'd': {
            config.logLevel.value = autopas::Logger::LogLevel::debug;
            break;
          }
          case 'i': {
            config.logLevel.value = autopas::Logger::LogLevel::info;
            break;
          }
          case 'w': {
            config.logLevel.value = autopas::Logger::LogLevel::warn;
            break;
          }
          case 'e': {
            config.logLevel.value = autopas::Logger::LogLevel::err;
            break;
          }
          case 'c': {
            config.logLevel.value = autopas::Logger::LogLevel::critical;
            break;
          }
          case 'o': {
            config.logLevel.value = autopas::Logger::LogLevel::off;
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
        config.logFileName.value = strArg;
        break;
      }
      case 'm': {
        try {
          auto mean = stod(strArg);
          config.distributionMean.value = {mean, mean, mean};
        } catch (const exception &) {
          cerr << "Error parsing distribution mean: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'M': {
        try {
          config.maxTuningPhasesWithoutTest.value = (unsigned int)stoul(strArg);
          if (config.maxTuningPhasesWithoutTest.value < 1) {
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
          config.particlesPerDim.value = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing number of particles per dimension: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'N': {
        try {
          config.particlesTotal.value = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing total number of particles: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'o': {
        try {
          config.relativeOptimumRange.value = (double)stoul(strArg);
          if (config.relativeOptimumRange.value < 1) {
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
          config.tuningPhases.value = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing number of tuning phases: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'p': {
        try {
          config.periodic.value = autopas::utils::StringUtils::parseBoolOption(strArg);
        } catch (const exception &) {
          cerr << "Error parsing whether there should be periodic boundary conditions: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'q': {
        try {
          config.verletClusterSize.value = (unsigned int)stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet cluster size: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'r': {
        try {
          config.verletSkinRadius.value = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet-skin-radius: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'S': {
        try {
          config.tuningSamples.value = (unsigned int)stoul(strArg);
          if (config.tuningSamples.value < 1) {
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
          config.particleSpacing.value = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing separation of particles: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case 't': {
        config.traversalOptions.value = autopas::TraversalOption::parseOptions(strArg);
        if (config.traversalOptions.value.empty()) {
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
        config.tuningStrategyOption.value = *parsedOptions.begin();
        break;
      }
      case 'u': {
        config.useThermostat.value = autopas::utils::StringUtils::parseBoolOption(strArg);
        break;
      }
      case 'v': {
        try {
          config.verletRebuildFrequency.value = (unsigned int)stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet-rebuild-frequency: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'w': {
        config.vtkFileName.value = strArg;
        break;
      }
      case 'W': {
        try {
          config.vtkWriteFrequency.value = (size_t)stoul(strArg);
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
        config.selectorStrategy.value = *parsedOptions.begin();
        break;
      }
      case 'Y': {
        // already parsed in CLIParser::inputFilesPresent
        break;
      }
      case 'z': {
        try {
          auto stdDev = stod(strArg);
          config.distributionStdDev.value = {stdDev, stdDev, stdDev};
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
  if (config.checkpointfile.value.empty() and config.cubeGaussObjects.empty() and config.cubeGridObjects.empty() and
      config.cubeUniformObjects.empty() and config.sphereObjects.empty()) {
    // common settings for any object type:
    unsigned int typeID = 0;
    double epsilon = 1.;
    double sigma = 1.;
    double mass = 1.;
    std::array<double, 3> bottomLeftCorner = {0, 0, 0};
    std::array<double, 3> velocity = {0, 0, 0};

    switch (config.generatorOption.value) {
      case MDFlexConfig::GeneratorOption::grid: {
        CubeGrid grid(velocity, typeID, epsilon, sigma, mass,
                      {config.particlesPerDim.value, config.particlesPerDim.value, config.particlesPerDim.value},
                      config.particleSpacing.value, bottomLeftCorner);
        config.cubeGridObjects.push_back(grid);
        break;
      }
      case MDFlexConfig::GeneratorOption::gaussian: {
        CubeGauss cubeGauss(velocity, typeID, epsilon, sigma, mass, config.particlesTotal.value,
                            {config.boxLength.value, config.boxLength.value, config.boxLength.value},
                            config.distributionMean.value, config.distributionStdDev.value, bottomLeftCorner);
        config.cubeGaussObjects.push_back(cubeGauss);
        break;
      }
      case MDFlexConfig::GeneratorOption::uniform: {
        CubeUniform cubeUniform(velocity, typeID, epsilon, sigma, mass, config.particlesTotal.value,
                                {config.boxLength.value, config.boxLength.value, config.boxLength.value},
                                bottomLeftCorner);
        config.cubeUniformObjects.push_back(cubeUniform);
        break;
      }
      case MDFlexConfig::GeneratorOption::sphere: {
        auto centerOfBox = config.particlesPerDim.value / 2.;
        Sphere sphere(velocity, typeID, epsilon, sigma, mass, {centerOfBox, centerOfBox, centerOfBox}, centerOfBox,
                      config.particleSpacing.value);
        config.sphereObjects.push_back(sphere);
        break;
      }
    }
  }

  if (displayHelp) {
    printHelpMessage(std::cout, argv[0], relevantOptions);
    return MDFlexParser::exitCodes::parsingError;
  }
  return MDFlexParser::exitCodes::success;
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

void MDFlexParser::CLIParser::inputFilesPresent(int argc, char **argv, MDFlexConfig &config) {
  // suppress error messages since we only want to look if the yaml option is there
  auto opterrBefore = opterr;
  opterr = 0;
  static struct option longOptions[] = {config.checkpointfile.toGetoptOption(),
                                        config.yamlFilename.toGetoptOption(),
                                        {nullptr, 0, nullptr, 0}};  // needed to signal the end of the array
  optind = 1;

  // search all cli parameters for input file options
  for (int cliOption = 0, cliOptionIndex = 0;
       (cliOption = getopt_long(argc, argv, "", longOptions, &cliOptionIndex)) != -1;) {
    std::string strArg;
    switch (cliOption) {
      case 'K':
        config.checkpointfile.value = optarg;
        if (not checkFileExists(optarg)) {
          throw std::runtime_error("CLIParser::inputFilesPresent: Checkpoint-File " + config.checkpointfile.value +
                                   " not found!");
        }
        break;
      case 'Y':
        config.yamlFilename.value = optarg;
        if (not checkFileExists(optarg)) {
          throw std::runtime_error("CLIParser::inputFilesPresent: Yaml-File " + config.yamlFilename.value +
                                   " not found!");
        }
        break;
      default: {
        // do nothing
      }
    }
  }

  opterr = opterrBefore;
}
