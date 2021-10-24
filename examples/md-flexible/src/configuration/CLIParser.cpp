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

  // utility options
  // getoptChars for all other options are line numbers so use negative numbers here to avoid clashes
  // also do not use -1 because it is used by getopt to signal that there is no cli option
  const static auto helpOption{MDFlexConfig::MDFlexOption<std::string, -2>("", "help", false, "Display this message.")};
  const static auto zshCompletionsOption{
      MDFlexConfig::MDFlexOption<std::string, -3>("", "zsh-completions", false, "Generate completions file for zsh.")};

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
      config.dontCreateEndConfig, config.dontShowProgressBar, config.tuningMaxEvidence,
      config.extrapolationMethodOption, config.evidenceFirstPrediction, config.functorOption, config.dontMeasureFlops,
      config.generatorOption, config.iterations, config.tuningInterval, config.logLevel, config.logFileName,
      config.distributionMean, config.maxTuningPhasesWithoutTest, config.particlesPerDim, config.particlesTotal,
      config.relativeOptimumRange, config.relativeBlacklistRange, config.periodic, config.tuningPhases,
      config.verletClusterSize, config.verletSkinRadius, config.particleSpacing, config.tuningSamples,
      config.traversalOptions, config.tuningStrategyOption, config.mpiStrategyOption, config.useThermostat,
      config.verletRebuildFrequency, config.vtkFileName, config.vtkWriteFrequency, config.selectorStrategy,
      config.yamlFilename, config.distributionStdDev, config.globalForce, config.loadBalancer, zshCompletionsOption,
      helpOption)};

  constexpr auto relevantOptionsSize = std::tuple_size_v<decltype(relevantOptions)>;

  // sanity check that all getopt chars are unique. Brackets for scoping.
  {
    // map tracking mappings of getopt chars to strings
    std::map<int, std::string> getoptCharsToName;
    // look for clashes by checking if getopt chars are in the map and otherwise add them
    autopas::utils::TupleUtils::for_each(relevantOptions, [&](auto &opt) {
      if (auto iterAtClash = getoptCharsToName.find(opt.getoptChar); iterAtClash != getoptCharsToName.end()) {
        throw std::runtime_error("CLIParser::parseInput: the following options share the same getopt char!\n" +
                                 opt.name + " : " + std::to_string(opt.getoptChar) + "\n" + iterAtClash->second +
                                 " : " + std::to_string(iterAtClash->first));
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
      case decltype(config.newton3Options)::getoptChar: {
        config.newton3Options.value = autopas::Newton3Option::parseOptions(strArg);
        if (config.newton3Options.value.empty()) {
          cerr << "Unknown Newton3 option: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.checkpointfile)::getoptChar: {
        // already parsed in CLIParser::inputFilesPresent
        break;
      }
      case decltype(config.acquisitionFunctionOption)::getoptChar: {
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
      case decltype(config.cellSizeFactors)::getoptChar: {
        config.cellSizeFactors.value = autopas::utils::StringUtils::parseNumberSet(strArg);
        if (config.cellSizeFactors.value->isEmpty()) {
          cerr << "Error parsing cell size factors: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.boxLength)::getoptChar: {
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
      case decltype(config.containerOptions)::getoptChar: {
        config.containerOptions.value = autopas::ContainerOption::parseOptions(strArg);
        if (config.containerOptions.value.empty()) {
          cerr << "Unknown container option: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.cutoff)::getoptChar: {
        try {
          config.cutoff.value = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing cutoff Radius: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.deltaT)::getoptChar: {
        try {
          config.deltaT.value = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing epsilon value: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.dataLayoutOptions)::getoptChar: {
        config.dataLayoutOptions.value = autopas::DataLayoutOption::parseOptions(strArg);
        if (config.dataLayoutOptions.value.empty()) {
          cerr << "Unknown data layouts: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.dontCreateEndConfig)::getoptChar: {
        config.dontCreateEndConfig.value = false;
        break;
      }
      case decltype(config.dontShowProgressBar)::getoptChar: {
        config.dontShowProgressBar.value = true;
        break;
      }
      case decltype(config.tuningMaxEvidence)::getoptChar: {
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
      case decltype(config.extrapolationMethodOption)::getoptChar: {
        auto parsedOptions = autopas::ExtrapolationMethodOption::parseOptions(strArg);
        if (parsedOptions.size() != 1) {
          cerr << "Pass exactly one extrapolation method option." << endl
               << "Passed: " << strArg << endl
               << "Parsed: " << autopas::utils::ArrayUtils::to_string(parsedOptions) << endl;
          displayHelp = true;
        }
        config.extrapolationMethodOption.value = *parsedOptions.begin();
        break;
      }
      case decltype(config.evidenceFirstPrediction)::getoptChar: {
        try {
          config.evidenceFirstPrediction.value = (unsigned int)stoul(strArg);
          if (config.evidenceFirstPrediction.value < 2) {
            cerr << "The number of evidence for the first prediction has to be at least two!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing max tuning phases without test: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.functorOption)::getoptChar: {
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
      case decltype(config.dontMeasureFlops)::getoptChar: {
        config.dontMeasureFlops.value = false;
        break;
      }
      case decltype(config.generatorOption)::getoptChar: {
        if (strArg.find("grid") != string::npos) {
          config.generatorOption.value = MDFlexConfig::GeneratorOption::grid;
        } else if (strArg.find("uni") != string::npos) {
          config.generatorOption.value = MDFlexConfig::GeneratorOption::uniform;
        } else if (strArg.find("gaus") != string::npos) {
          config.generatorOption.value = MDFlexConfig::GeneratorOption::gaussian;
        } else if (strArg.find("sp") != string::npos) {
          config.generatorOption.value = MDFlexConfig::GeneratorOption::sphere;
        } else if (strArg.find("cl") != string::npos) {
          config.generatorOption.value = MDFlexConfig::GeneratorOption::closestPacked;
        } else {
          cerr << "Unknown generator: " << strArg << endl;
          cerr << "Please use 'Grid' or 'Gaussian'" << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(helpOption)::getoptChar: {
        printHelpMessage(std::cout, argv[0], relevantOptions);
        return MDFlexParser::exitCodes::helpFlagFound;
      }
      case decltype(zshCompletionsOption)::getoptChar: {
        // generate the completions file and do nothing else
        createZSHCompletionFile(relevantOptions);
        return MDFlexParser::exitCodes::completionsFlagFound;
      }
      case decltype(config.iterations)::getoptChar: {
        try {
          config.iterations.value = stoul(strArg);
          if (config.iterations.value < 1) {
            cerr << "Number of iterations has to be a positive integer!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing number of iterations: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.tuningInterval)::getoptChar: {
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
      case decltype(config.logLevel)::getoptChar: {
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
      case decltype(config.logFileName)::getoptChar: {
        config.logFileName.value = strArg;
        break;
      }
      case decltype(config.distributionMean)::getoptChar: {
        try {
          auto mean = stod(strArg);
          config.distributionMean.value = {mean, mean, mean};
        } catch (const exception &) {
          cerr << "Error parsing distribution mean: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.maxTuningPhasesWithoutTest)::getoptChar: {
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
      case decltype(config.particlesPerDim)::getoptChar: {
        try {
          config.particlesPerDim.value = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing number of particles per dimension: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.particlesTotal)::getoptChar: {
        try {
          config.particlesTotal.value = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing total number of particles: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.relativeOptimumRange)::getoptChar: {
        try {
          config.relativeOptimumRange.value = stod(strArg);
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
      case decltype(config.relativeBlacklistRange)::getoptChar: {
        try {
          config.relativeBlacklistRange.value = stod(strArg);
          if (config.relativeBlacklistRange.value < 1 && config.relativeBlacklistRange.value != 0) {
            cerr << "Relative range for blacklist range has to be greater or equal one or has to be zero!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing relative range for blacklist: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.tuningPhases)::getoptChar: {
        try {
          config.tuningPhases.value = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing number of tuning phases: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.periodic)::getoptChar: {
        try {
          config.periodic.value = autopas::utils::StringUtils::parseBoolOption(strArg);
        } catch (const exception &) {
          cerr << "Error parsing whether there should be periodic boundary conditions: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.verletClusterSize)::getoptChar: {
        try {
          config.verletClusterSize.value = (unsigned int)stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet cluster size: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.verletSkinRadius)::getoptChar: {
        try {
          config.verletSkinRadius.value = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet-skin-radius: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.tuningSamples)::getoptChar: {
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
      case decltype(config.particleSpacing)::getoptChar: {
        try {
          config.particleSpacing.value = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing separation of particles: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.traversalOptions)::getoptChar: {
        config.traversalOptions.value = autopas::TraversalOption::parseOptions(strArg);
        if (config.traversalOptions.value.empty()) {
          cerr << "Unknown Traversal: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.loadEstimatorOptions)::getoptChar: {
        config.loadEstimatorOptions.value = autopas::LoadEstimatorOption::parseOptions(strArg);
        if (config.loadEstimatorOptions.value.empty()) {
          cerr << "Unknown load estimator: " << strArg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.tuningStrategyOption)::getoptChar: {
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
      case decltype(config.mpiStrategyOption)::getoptChar: {
        auto parsedOptions = autopas::MPIStrategyOption::parseOptions(strArg);
        if (parsedOptions.size() != 1) {
          cerr << "Pass exactly one mpi strategy option. AutoPas cannot switch between several." << endl
               << "Passed: " << strArg << endl
               << "Parsed: " << autopas::utils::ArrayUtils::to_string(parsedOptions) << endl;
          displayHelp = true;
        }
        config.mpiStrategyOption.value = *parsedOptions.begin();
        break;
      }
      case decltype(config.useThermostat)::getoptChar: {
        config.useThermostat.value = autopas::utils::StringUtils::parseBoolOption(strArg);
        break;
      }
      case decltype(config.verletRebuildFrequency)::getoptChar: {
        try {
          config.verletRebuildFrequency.value = (unsigned int)stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet-rebuild-frequency: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.vtkFileName)::getoptChar: {
        config.vtkFileName.value = strArg;
        break;
      }
      case decltype(config.vtkWriteFrequency)::getoptChar: {
        try {
          config.vtkWriteFrequency.value = (size_t)stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing vtk write frequency: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.selectorStrategy)::getoptChar: {
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
      case decltype(config.yamlFilename)::getoptChar: {
        // already parsed in CLIParser::inputFilesPresent
        break;
      }
      case decltype(config.distributionStdDev)::getoptChar: {
        try {
          auto stdDev = stod(strArg);
          config.distributionStdDev.value = {stdDev, stdDev, stdDev};
        } catch (const exception &) {
          cerr << "Error parsing distribution standard deviation: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.globalForce)::getoptChar: {
        // when passing via cli global force can only be in z-direction. For fancier forces use yaml input.
        try {
          auto force = stod(strArg);
          config.globalForce.value = {0, 0, force};
        } catch (const exception &) {
          cerr << "Error parsing global force: " << optarg << endl;
          cerr << "Expecting one double as force along the z-axis." << endl;
          displayHelp = true;
        }
        break;
      }
      case decltype(config.loadBalancer)::getoptChar: {
        auto parsedOptions = LoadBalancerOption::parseOptions(strArg);
        if (parsedOptions.size() != 1) {
          cerr << "Pass exactly one load balancer option." << endl
               << "Passed: " << strArg << endl
               << "Parsed: " << autopas::utils::ArrayUtils::to_string(parsedOptions) << endl;
          displayHelp = true;
        }
        config.loadBalancer.value = *parsedOptions.begin();
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
      config.cubeUniformObjects.empty() and config.sphereObjects.empty() and config.cubeClosestPackedObjects.empty()) {
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
      case MDFlexConfig::GeneratorOption::closestPacked: {
        CubeClosestPacked cubeClosestPacked(velocity, typeID, epsilon, sigma, mass, config.particleSpacing.value,
                                            {config.boxLength.value, config.boxLength.value, config.boxLength.value},
                                            bottomLeftCorner);
        config.cubeClosestPackedObjects.push_back(cubeClosestPacked);
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
      case decltype(config.checkpointfile)::getoptChar:
        config.checkpointfile.value = optarg;
        if (not checkFileExists(optarg)) {
          throw std::runtime_error("CLIParser::inputFilesPresent: Checkpoint-File " + config.checkpointfile.value +
                                   " not found!");
        }
        break;
      case decltype(config.yamlFilename)::getoptChar:
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
