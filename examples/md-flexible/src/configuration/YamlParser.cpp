/**
 * @file YamlParser.cpp
 * @author N. Fottner, D. Martin
 * @date 15.07.2019, 11.04.2023
 */
#include "YamlParser.h"

#include <string>

void MDFlexParser::YamlParser::throwObjectParseException(const char *key, const char *exp) {
  std::string msg;
  msg.append("Error parsing ");
  msg.append(key);
  msg.append(". Make sure that key \"");
  msg.append(key);
  msg.append("\" exists and has the expected value: ");
  msg.append(exp);
  throw std::runtime_error(msg);
}

template <typename T>
T MDFlexParser::YamlParser::parseObjectValue(YAML::iterator &it, const char *key, const char *exp) {
  T value;
  try {
    value = it->second[key].as<T>();
  } catch (const std::exception &e) {
    throwObjectParseException(key, exp);
  }
  return value;
}

template <typename T, size_t S>
std::array<T, S> MDFlexParser::YamlParser::parseObjectValue(YAML::iterator &it, const char *key, const char *exp) {
  std::array<T, S> value;
  try {
    YAML::Node n = it->second[key];
    for (int i = 0; i < S; i++) {
      value[i] = n[i].as<T>();
    }

  } catch (const std::exception &e) {
    throwObjectParseException(key, exp);
  }
  return value;
}

std::array<double, 3> MDFlexParser::YamlParser::parseVelocity(MDFlexConfig &config, YAML::iterator &it) {
  return parseObjectValue<double, 3>(it, MDFlexConfig::velocityStr, "Three doubles as YAML-sequence");
}

unsigned long MDFlexParser::YamlParser::parseParticleType(MDFlexConfig &config, YAML::iterator &it) {
  return parseObjectValue<unsigned long>(it, MDFlexConfig::particleTypeStr, "Unsigned Integer");
}

double MDFlexParser::YamlParser::parseEpsilon(MDFlexConfig &config, YAML::iterator &it) {
  return parseObjectValue<double>(it, config.epsilonMap.name.c_str(), "Double");
}

double MDFlexParser::YamlParser::parseSigma(MDFlexConfig &config, YAML::iterator &it) {
  return parseObjectValue<double>(it, config.sigmaMap.name.c_str(), "Double");
}

double MDFlexParser::YamlParser::parseMass(MDFlexConfig &config, YAML::iterator &it) {
  return parseObjectValue<double>(it, config.massMap.name.c_str(), "Double");
}

double MDFlexParser::YamlParser::parseParticleSpacing(MDFlexConfig &config, YAML::iterator &it) {
  return parseObjectValue<double>(it, config.particleSpacing.name.c_str(), "Double");
}

std::array<unsigned long, 3> MDFlexParser::YamlParser::parseParticlesPerDim(MDFlexConfig &config, YAML::iterator &it) {
  return parseObjectValue<unsigned long, 3>(it, config.particlesPerDim.name.c_str(),
                                            "Three unsigned integers as YAML-sequence");
}

std::array<double, 3> MDFlexParser::YamlParser::parseBottomLeftCorner(MDFlexConfig &config, YAML::iterator &it) {
  return parseObjectValue<double, 3>(it, MDFlexConfig::bottomLeftBackCornerStr, "Three doubles as YAML-sequence");
}

std::array<double, 3> MDFlexParser::YamlParser::parseDistrMean(MDFlexConfig &config, YAML::iterator &it) {
  return parseObjectValue<double, 3>(it, config.distributionMean.name.c_str(), "Three doubles as YAML-sequence");
}

std::array<double, 3> MDFlexParser::YamlParser::parseDistrStdDev(MDFlexConfig &config, YAML::iterator &it) {
  return parseObjectValue<double, 3>(it, config.distributionStdDev.name.c_str(), "Three doubles as YAML-sequence");
}

size_t MDFlexParser::YamlParser::parseNumParticles(MDFlexConfig &config, YAML::iterator &it) {
  return parseObjectValue<size_t>(it, MDFlexConfig::particlesPerObjectStr, "Unsigned Integer");
}

std::array<double, 3> MDFlexParser::YamlParser::parseBoxLength(MDFlexConfig &config, YAML::iterator &it) {
  return parseObjectValue<double, 3>(it, config.boxLength.name.c_str(), "Three doubles as YAML-sequence");
}

std::array<double, 3> MDFlexParser::YamlParser::parseCenter(MDFlexConfig &config, YAML::iterator &it) {
  return parseObjectValue<double, 3>(it, MDFlexConfig::sphereCenterStr, "Three doubles as YAML-sequence");
}

double MDFlexParser::YamlParser::parseRadius(MDFlexConfig &config, YAML::iterator &it) {
  return parseObjectValue<double>(it, MDFlexConfig::sphereRadiusStr, "Double");
}

CubeGrid MDFlexParser::YamlParser::parseCubeGridObject(MDFlexConfig &config, YAML::iterator &it) {
  std::array<double, 3> velocity = parseVelocity(config, it);
  unsigned long particleType = parseParticleType(config, it);
  double epsilon = parseEpsilon(config, it);
  double sigma = parseSigma(config, it);
  double mass = parseMass(config, it);
  std::array<unsigned long, 3> particlesPerDim = parseParticlesPerDim(config, it);
  double particleSpacing = parseParticleSpacing(config, it);
  std::array<double, 3> bottomLeftCorner = parseBottomLeftCorner(config, it);

  CubeGrid cubeGrid(velocity, particleType, epsilon, sigma, mass, particlesPerDim, particleSpacing, bottomLeftCorner);
  return cubeGrid;
}

CubeUniform MDFlexParser::YamlParser::parseCubeUniformObject(MDFlexConfig &config, YAML::iterator &it) {
  std::array<double, 3> velocity = parseVelocity(config, it);
  unsigned long particleType = parseParticleType(config, it);
  double epsilon = parseEpsilon(config, it);
  double sigma = parseSigma(config, it);
  double mass = parseMass(config, it);
  size_t numParticles = parseNumParticles(config, it);
  std::array<double, 3> boxLength = parseBoxLength(config, it);
  std::array<double, 3> bottomLeftCorner = parseBottomLeftCorner(config, it);

  CubeUniform cubeUniform(velocity, particleType, epsilon, sigma, mass, numParticles, boxLength, bottomLeftCorner);
  return cubeUniform;
}

CubeGauss MDFlexParser::YamlParser::parseCubeGaussObject(MDFlexConfig &config, YAML::iterator &it) {
  std::array<double, 3> velocity = parseVelocity(config, it);
  unsigned long particleType = parseParticleType(config, it);
  double epsilon = parseEpsilon(config, it);
  double sigma = parseSigma(config, it);
  double mass = parseMass(config, it);
  size_t numParticles = parseNumParticles(config, it);
  std::array<double, 3> boxLength = parseBoxLength(config, it);
  std::array<double, 3> distributionMean = parseDistrMean(config, it);
  std::array<double, 3> distributionStdDev = parseDistrStdDev(config, it);
  std::array<double, 3> bottomLeftCorner = parseBottomLeftCorner(config, it);

  CubeGauss cubeGauss(velocity, particleType, epsilon, sigma, mass, numParticles, boxLength, distributionMean,
                      distributionStdDev, bottomLeftCorner);
  return cubeGauss;
}

Sphere MDFlexParser::YamlParser::parseSphereObject(MDFlexConfig &config, YAML::iterator &it) {
  std::array<double, 3> velocity = parseVelocity(config, it);
  unsigned long particleType = parseParticleType(config, it);
  double epsilon = parseEpsilon(config, it);
  double sigma = parseSigma(config, it);
  double mass = parseMass(config, it);
  std::array<double, 3> sphereCenter = parseCenter(config, it);
  double sphereRadius = parseRadius(config, it);
  double particleSpacing = parseParticleSpacing(config, it);

  Sphere sphere(velocity, particleType, epsilon, sigma, mass, sphereCenter, sphereRadius, particleSpacing);
  return sphere;
}

CubeClosestPacked MDFlexParser::YamlParser::parseCubeClosestPacked(MDFlexConfig &config, YAML::iterator &it) {
  std::array<double, 3> velocity = parseVelocity(config, it);
  unsigned long particleType = parseParticleType(config, it);
  double epsilon = parseEpsilon(config, it);
  double sigma = parseSigma(config, it);
  double mass = parseMass(config, it);
  double particleSpacing = parseParticleSpacing(config, it);
  std::array<double, 3> boxLength = parseBoxLength(config, it);
  std::array<double, 3> bottomLeftCorner = parseBottomLeftCorner(config, it);

  CubeClosestPacked cubeClosestPacked(velocity, particleType, epsilon, sigma, mass, particleSpacing, boxLength,
                                      bottomLeftCorner);

  return cubeClosestPacked;
}

bool MDFlexParser::YamlParser::parseYamlFile(MDFlexConfig &config) {
  /*
  Global variables used to print the expected input and a description of the parameter if an error occurs while
  parsing. Yaml mark is used to identify the current line of the error.
  */
  std::string expected;
  std::string description;
  YAML::Mark m;

  YAML::Node node = YAML::LoadFile(config.yamlFilename.value);

  // We iterate over all keys to identify known/unknown parameters.
  for (auto itemIterator = node.begin(); itemIterator != node.end(); ++itemIterator) {
    std::string key;
    try {
      key = itemIterator->first.as<std::string>();
      m = node[key].Mark();

      if (key == config.containerOptions.name) {
        expected = "YAML-sequence of possible values.";
        description = config.containerOptions.description;

        config.containerOptions.value =
            autopas::ContainerOption::parseOptions(autopas::utils::ArrayUtils::to_string(node[key], ", ", {"", ""}));

        if (config.containerOptions.value.empty()) {
          throw YamlParserException("Parsed container list is empty. You used possibly an unknown container option.");
        }

      } else if (key == config.boxMin.name) {
        expected = "YAML-sequence of three floats. Example: [0, 0, 0].";
        description = config.boxMin.description;

        auto tmpNode = node[key];
        config.boxMin.value = {tmpNode[0].as<double>(), tmpNode[1].as<double>(), tmpNode[2].as<double>()};
      } else if (key == config.boxMax.name) {
        expected = "YAML-sequence of three floats. Example: [42, 42, 42].";
        description = config.boxMax.description;

        auto tmpNode = node[key];
        config.boxMax.value = {tmpNode[0].as<double>(), tmpNode[1].as<double>(), tmpNode[2].as<double>()};
      } else if (key == config.subdivideDimension.name) {
        expected = "YAML-sequence of three ints in [0, 1].";
        description = config.subdivideDimension.description;

        auto tmpNode = node[key];
        config.subdivideDimension.value = {tmpNode[0].as<bool>(), tmpNode[1].as<bool>(), tmpNode[2].as<bool>()};
      } else if (key == config.loadBalancingInterval.name) {
        expected = "Unsigned Integer";
        description = config.loadBalancingInterval.description;

        int tmp = node[key].as<int>();
        if (tmp < 0) {
          throw YamlParserException("Load balancing interval must be a positive integer.");
        }

        config.loadBalancingInterval.value = tmp;
      } else if (key == config.selectorStrategy.name) {
        expected = "Exactly one selector strategy out of the possible values.";
        description = config.selectorStrategy.description;

        std::set<autopas::options::SelectorStrategyOption> parsedOptions;
        if (node[key].IsSequence()) {
          if (node[key].size() != 1) {
            throw YamlParserException("Pass Exactly one selector strategy.");
          }
          parsedOptions = autopas::SelectorStrategyOption::parseOptions(
              autopas::utils::ArrayUtils::to_string(node[key], "", {"", ""}));
        } else {
          parsedOptions = autopas::SelectorStrategyOption::parseOptions(node[key].as<std::string>());
        }
        config.selectorStrategy.value = *parsedOptions.begin();

      } else if (key == config.boundaryOption.name) {
        expected = "YAML-sequence of three possible values.";
        description = config.boundaryOption.description;

        auto tmpNode = node[key];
        config.boundaryOption.value = {options::BoundaryTypeOption::parseOptionExact(tmpNode[0].as<std::string>()),
                                       options::BoundaryTypeOption::parseOptionExact(tmpNode[1].as<std::string>()),
                                       options::BoundaryTypeOption::parseOptionExact(tmpNode[2].as<std::string>())};
      } else if (key == config.cutoff.name) {
        expected = "Positive floating point value > 0.";
        description = config.cutoff.description;

        double tmp = node[key].as<double>();
        if (tmp <= 0) {
          throw YamlParserException("Cutoff has to be > 0!");
        }

        config.cutoff.value = tmp;
      } else if (key == config.cellSizeFactors.name) {
        expected = "YAML-sequence of floats.";
        description = config.cellSizeFactors.description;

        config.cellSizeFactors.value = autopas::utils::StringUtils::parseNumberSet(
            autopas::utils::ArrayUtils::to_string(node[key], ", ", {"", ""}));

        if (config.cellSizeFactors.value->isEmpty()) {
          throw YamlParserException("Parsed cell-size-factor-list is empty.");
        }
      } else if (key == config.dataLayoutOptions.name) {
        expected = "YAML-sequence of possible values.";
        description = config.dataLayoutOptions.description;

        config.dataLayoutOptions.value =
            autopas::DataLayoutOption::parseOptions(autopas::utils::ArrayUtils::to_string(node[key], ", ", {"", ""}));
        if (config.dataLayoutOptions.value.empty()) {
          throw YamlParserException("Parsed data-layouts-list is empty.");
        }
      } else if (key == config.functorOption.name) {
        expected = "One of the possible values.";
        description = config.functorOption.description;

        auto strArg = node[key].as<std::string>();
        transform(strArg.begin(), strArg.end(), strArg.begin(), ::tolower);
        if (strArg.find("avx") != std::string::npos) {
          config.functorOption.value = MDFlexConfig::FunctorOption::lj12_6_AVX;
        } else if (strArg.find("sve") != std::string::npos) {
          config.functorOption.value = MDFlexConfig::FunctorOption::lj12_6_SVE;
        } else if (strArg.find("glob") != std::string::npos) {
          config.functorOption.value = MDFlexConfig::FunctorOption::lj12_6_Globals;
        } else if (strArg.find("lj") != std::string::npos or strArg.find("lennard-jones") != std::string::npos) {
          config.functorOption.value = MDFlexConfig::FunctorOption::lj12_6;
        } else {
          throw YamlParserException("Unrecognized functor!");
        }
      } else if (key == config.iterations.name) {
        expected = "Unsigned Integer > 0";
        description = config.iterations.description;

        long tmp = node[key].as<long>();
        if (tmp < 1) {
          throw YamlParserException("The number of iterations has to be a positive integer > 0.");
        }
        config.iterations.value = tmp;

      } else if (key == config.tuningPhases.name) {
        expected = "Unsigned Integer";
        description = config.tuningPhases.description;

        long tmp = node[key].as<long>();
        if (tmp < 0) {
          throw YamlParserException("The number of tuning phases has to be a positive integer.");
        }

        config.tuningPhases.value = tmp;
      } else if (key == config.dontMeasureFlops.name) {
        expected = "Boolean Value";
        description = config.dontMeasureFlops.description;

        // "not" needed because of semantics
        config.dontMeasureFlops.value = not node[key].as<bool>();
      } else if (key == config.dontCreateEndConfig.name) {
        expected = "Boolean Value";
        description = config.dontCreateEndConfig.description;

        // "not" needed because of semantics
        config.dontCreateEndConfig.value = not node[key].as<bool>();
      } else if (key == config.dontShowProgressBar.name) {
        expected = "Boolean Value";
        description = config.dontShowProgressBar.description;

        config.dontShowProgressBar.value = node[key].as<bool>();
      } else if (key == config.newton3Options.name) {
        expected = "YAML-sequence of possible values.";
        description = config.newton3Options.description;

        config.newton3Options.value =
            autopas::Newton3Option::parseOptions(autopas::utils::ArrayUtils::to_string(node[key], ", ", {"", ""}));
        if (config.newton3Options.value.empty()) {
          throw YamlParserException("Unknown Newton3 option!");
        }
      } else if (key == config.deltaT.name) {
        expected = "Positive floating point value.";
        description = config.deltaT.description;

        config.deltaT.value = node[key].as<double>();
      } else if (key == config.traversalOptions.name) {
        expected = "YAML-sequence of possible values.";
        description = config.traversalOptions.description;

        config.traversalOptions.value =
            autopas::TraversalOption::parseOptions(autopas::utils::ArrayUtils::to_string(node[key], ", ", {"", ""}));

        if (config.traversalOptions.value.empty()) {
          throw YamlParserException("Parsed traversal-list is empty. Maybe you used an unknown option.");
        }

      } else if (key == config.loadEstimatorOptions.name) {
        expected = "YAML-sequence of possible values.";
        description = config.loadEstimatorOptions.description;

        config.loadEstimatorOptions.value = autopas::LoadEstimatorOption::parseOptions(
            autopas::utils::ArrayUtils::to_string(node[key], ", ", {"", ""}));

        if (config.loadEstimatorOptions.value.empty()) {
          throw YamlParserException("Parsed load-estimator-list is empty. Maybe you used an unknown option.");
        }

      } else if (key == config.tuningInterval.name) {
        expected = "Unsigned Integer";
        description = config.tuningInterval.description;

        int tmp = node[key].as<int>();
        if (tmp < 1) {
          throw YamlParserException("Tuning interval has to be a positive integer!");
        }

        config.tuningInterval.value = tmp;
      } else if (key == config.tuningSamples.name) {
        expected = "Unsigned Integer >= 1";
        description = config.tuningSamples.description;

        int tmp = node[key].as<int>();
        if (tmp < 1) {
          throw YamlParserException("Tuning samples has to be a positive integer!");
        }

        config.tuningSamples.value = tmp;
      } else if (key == config.tuningMaxEvidence.name) {
        expected = "Unsigned Integer >= 1";
        description = config.tuningMaxEvidence.description;

        int tmp = node[key].as<int>();
        if (tmp < 1) {
          throw YamlParserException("Tuning max evidence has to be a positive integer >= 1!");
        }

        config.tuningMaxEvidence.value = tmp;
      } else if (key == config.relativeOptimumRange.name) {
        expected = "Floating point value >= 1";
        description = config.relativeOptimumRange.description;

        double tmp = node[key].as<double>();
        if (tmp < 1.0) {
          throw YamlParserException("Relative optimum range has to be greater or equal one!");
        }

        config.relativeOptimumRange.value = tmp;
      } else if (key == config.maxTuningPhasesWithoutTest.name) {
        expected = "Unsigned Integer";
        description = config.maxTuningPhasesWithoutTest.description;

        int tmp = node[key].as<int>();
        if (tmp < 1) {
          throw YamlParserException("Max tuning phases without test has to be positive!");
        }

        config.maxTuningPhasesWithoutTest.value = tmp;
      } else if (key == config.relativeBlacklistRange.name) {
        expected = "Floating point value >= 1 or 0";
        description = config.relativeBlacklistRange.description;

        double tmp = node[key].as<double>();
        if (tmp < 1.0 and tmp != 0.0) {
          throw YamlParserException(
              "Relative range for blacklist range has to be greater or equal one or has to be zero!");
        }

        config.relativeBlacklistRange.value = tmp;
      } else if (key == config.evidenceFirstPrediction.name) {
        expected = "Unsigned Integer >= 2";
        description = config.evidenceFirstPrediction.description;

        int tmp = node[key].as<int>();
        if (tmp < 2) {
          throw YamlParserException("The number of evidence for the first prediction has to be at least two!");
        }

        config.evidenceFirstPrediction.value = tmp;
      } else if (key == config.extrapolationMethodOption.name) {
        expected = "Exactly one extrapolation method out of the possible values.";
        description = config.extrapolationMethodOption.description;

        std::set<autopas::options::ExtrapolationMethodOption> parsedOptions;
        if (node[key].IsSequence()) {
          if (node[key].size() != 1) {
            throw YamlParserException("Pass exactly one extrapolation method!");
          }
          parsedOptions = autopas::ExtrapolationMethodOption::parseOptions(
              autopas::utils::ArrayUtils::to_string(node[key], "", {"", ""}));
        } else {
          parsedOptions = autopas::ExtrapolationMethodOption::parseOptions(node[key].as<std::string>());
        }
        config.extrapolationMethodOption.value = *parsedOptions.begin();

      } else if (key == config.tuningStrategyOption.name) {
        expected = "Exactly one tuning strategy option out of the possible values.";
        description = config.tuningStrategyOption.description;

        std::set<autopas::options::TuningStrategyOption> parsedOptions;
        if (node[key].IsSequence()) {
          if (node[key].size() != 1) {
            throw YamlParserException("Pass Exactly one tuning strategy!");
          }
          parsedOptions = autopas::TuningStrategyOption::parseOptions(
              autopas::utils::ArrayUtils::to_string(node[key], "", {"", ""}));
        } else {
          parsedOptions = autopas::TuningStrategyOption::parseOptions(node[key].as<std::string>());
        }
        config.tuningStrategyOption.value = *parsedOptions.begin();
      } else if (key == config.mpiStrategyOption.name) {
        expected = "Exactly one MPI strategy option out of the possible values.";
        description = config.mpiStrategyOption.description;

        std::set<autopas::options::MPIStrategyOption> parsedOptions;
        if (node[key].IsSequence()) {
          if (node[key].size() != 1) {
            throw YamlParserException("Pass exactly one MPI strategy!");
          }
          parsedOptions =
              autopas::MPIStrategyOption::parseOptions(autopas::utils::ArrayUtils::to_string(node[key], "", {"", ""}));
        } else {
          parsedOptions = autopas::MPIStrategyOption::parseOptions(node[key].as<std::string>());
        }
        config.mpiStrategyOption.value = *parsedOptions.begin();
      } else if (key == config.MPITuningMaxDifferenceForBucket.name) {
        expected = "Floating-point Value";
        description = config.MPITuningMaxDifferenceForBucket.description;

        config.MPITuningMaxDifferenceForBucket.value = node[key].as<double>();
      } else if (key == config.MPITuningWeightForMaxDensity.name) {
        expected = "Floating-point Value";
        description = config.MPITuningWeightForMaxDensity.description;

        config.MPITuningWeightForMaxDensity.value = node[key].as<double>();
      } else if (key == config.acquisitionFunctionOption.name) {
        expected = "Exactly one acquisition function option out of the possible values.";
        description = config.acquisitionFunctionOption.description;

        std::set<autopas::options::AcquisitionFunctionOption> parsedOptions;
        if (node[key].IsSequence()) {
          if (node[key].size() != 1) {
            throw YamlParserException("Pass Exactly one acquisition function option!");
          }
          parsedOptions = autopas::AcquisitionFunctionOption::parseOptions(
              autopas::utils::ArrayUtils::to_string(node[key], "", {"", ""}));
        } else {
          parsedOptions = autopas::AcquisitionFunctionOption::parseOptions(node[key].as<std::string>());
        }
        config.acquisitionFunctionOption.value = *parsedOptions.begin();
      } else if (key == config.logLevel.name) {
        expected = "Log level out of the possible values.";
        description = config.logLevel.description;

        auto strArg = node[key].as<std::string>();
        switch (std::tolower(strArg[0])) {
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
            throw YamlParserException("Unknown Log Level parsed!");
          }
        }
      } else if (key == config.checkpointfile.name) {
        expected = "String";
        description = config.checkpointfile.description;

        config.checkpointfile.value = node[key].as<std::string>();
        if (config.checkpointfile.value.empty()) {
          throw YamlParserException("Parsed checkpoint filename is empty!");
        }
      } else if (key == config.logFileName.name) {
        expected = "String";
        description = config.logFileName.description;

        config.logFileName.value = node[key].as<std::string>();
        if (config.logFileName.value.empty()) {
          throw YamlParserException("Parsed log filename is empty!");
        }
      } else if (key == config.verletRebuildFrequency.name) {
        expected = "Unsigned Integer >= 1";
        description = config.verletRebuildFrequency.description;

        int tmp = node[key].as<int>();
        if (tmp < 1) {
          throw YamlParserException("Verlet rebuild frequency has to be a positive integer >= 1!");
        }

        config.verletRebuildFrequency.value = tmp;
      } else if (key == config.verletSkinRadiusPerTimestep.name) {
        expected = "Positive floating-point value.";
        description = config.verletSkinRadiusPerTimestep.description;

        config.verletSkinRadiusPerTimestep.value = node[key].as<double>();
      } else if (key == config.fastParticlesThrow.name) {
        expected = "Boolean Value";
        description = config.fastParticlesThrow.description;

        config.fastParticlesThrow.value = node[key].as<bool>();
      } else if (key == config.verletClusterSize.name) {
        expected = "Unsigned Integer";
        description = config.verletClusterSize.description;

        int tmp = node[key].as<int>();
        if (tmp < 0) {
          throw YamlParserException("Verlet cluster size has to be a positive integer!");
        }

        config.verletClusterSize.value = tmp;
      } else if (key == config.vtkFileName.name) {
        expected = "String";
        description = config.vtkFileName.description;

        config.vtkFileName.value = node[key].as<std::string>();
        if (config.vtkFileName.value.empty()) {
          throw YamlParserException("Parsed VTK filename is empty!");
        }
      } else if (key == config.vtkWriteFrequency.name) {
        expected = "Unsigned Integer >= 1";
        description = config.vtkWriteFrequency.description;

        int tmp = node[key].as<int>();
        if (tmp < 1) {
          throw YamlParserException("VTK write frequency has to be a positive integer >= 1!");
        }

        config.vtkWriteFrequency.value = (size_t)tmp;
      } else if (key == config.globalForce.name) {
        expected = "YAML-sequence of three floats. Example: [0, 0, -9.81].";
        description = config.globalForce.description;

        config.globalForce.value = {node[key][0].as<double>(), node[key][1].as<double>(), node[key][2].as<double>()};
      } else if (key == MDFlexConfig::objectsStr) {
        expected = "See AllOptions.yaml for examples.";
        description = "";

        // remove default objects
        config.cubeGridObjects.clear();
        config.cubeGaussObjects.clear();
        config.cubeUniformObjects.clear();
        config.sphereObjects.clear();
        config.cubeClosestPackedObjects.clear();
        config.epsilonMap.value.clear();
        config.sigmaMap.value.clear();
        config.massMap.value.clear();

        int objID = 0;
        const char *generatorName;

        for (auto objectIterator = node[MDFlexConfig::objectsStr].begin();
             objectIterator != node[MDFlexConfig::objectsStr].end(); ++objectIterator) {
          try {
            if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeGridObjectsStr) {
              generatorName = MDFlexConfig::cubeGridObjectsStr;
              for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
                objID = std::distance(objectIterator->second.begin(), it);
                CubeGrid cubeGrid = parseCubeGridObject(config, it);

                config.cubeGridObjects.emplace_back(cubeGrid);
                config.addParticleType(it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                       it->second[config.epsilonMap.name].as<double>(),
                                       it->second[config.sigmaMap.name].as<double>(),
                                       it->second[config.massMap.name].as<double>());
              }
            } else if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeGaussObjectsStr) {
              generatorName = MDFlexConfig::cubeGaussObjectsStr;
              for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
                objID = std::distance(objectIterator->second.begin(), it);
                CubeGauss cubeGauss = parseCubeGaussObject(config, it);

                config.cubeGaussObjects.emplace_back(cubeGauss);
                config.addParticleType(it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                       it->second[config.epsilonMap.name].as<double>(),
                                       it->second[config.sigmaMap.name].as<double>(),
                                       it->second[config.massMap.name].as<double>());
              }
            } else if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeUniformObjectsStr) {
              generatorName = MDFlexConfig::cubeUniformObjectsStr;
              for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
                objID = std::distance(objectIterator->second.begin(), it);
                CubeUniform cubeUniform = parseCubeUniformObject(config, it);

                config.cubeUniformObjects.emplace_back(cubeUniform);
                config.addParticleType(it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                       it->second[config.epsilonMap.name].as<double>(),
                                       it->second[config.sigmaMap.name].as<double>(),
                                       it->second[config.massMap.name].as<double>());
              }
            } else if (objectIterator->first.as<std::string>() == MDFlexConfig::sphereObjectsStr) {
              generatorName = MDFlexConfig::sphereObjectsStr;
              for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
                objID = std::distance(objectIterator->second.begin(), it);
                Sphere sphere = parseSphereObject(config, it);

                config.sphereObjects.emplace_back(sphere);
                config.addParticleType(it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                       it->second[config.epsilonMap.name].as<double>(),
                                       it->second[config.sigmaMap.name].as<double>(),
                                       it->second[config.massMap.name].as<double>());
              }
            } else if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeClosestPackedObjectsStr) {
              generatorName = MDFlexConfig::cubeClosestPackedObjectsStr;
              for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
                objID = std::distance(objectIterator->second.begin(), it);
                CubeClosestPacked cubeClosestPacked = parseCubeClosestPacked(config, it);

                config.cubeClosestPackedObjects.emplace_back(cubeClosestPacked);
                config.addParticleType(it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                       it->second[config.epsilonMap.name].as<double>(),
                                       it->second[config.sigmaMap.name].as<double>(),
                                       it->second[config.massMap.name].as<double>());
              }
            } else {
              std::cerr << "YamlParser: Unrecognized generator \"" << objectIterator->first.as<std::string>()
                        << "\" used." << std::endl;
              return false;
            }
          } catch (const std::exception &e) {
            std::cerr << "YamlParser: Error parsing " << generatorName << " object with ID " << objID << "."
                      << std::endl
                      << "Message: " << e.what() << std::endl
                      << "See AllOptions.yaml for examples." << std::endl;
            return false;
          }
        }
      } else if (key == config.useThermostat.name) {
        expected = "See AllOptions.yaml for examples.";
        description = config.useThermostat.description;

        config.useThermostat.value = true;

        m = node[key][config.initTemperature.name].Mark();
        expected = "Floating-Point Value";
        description = config.initTemperature.description;
        config.initTemperature.value = node[key][config.initTemperature.name].as<double>();

        m = node[key][config.thermostatInterval.name].Mark();
        expected = "Unsigned Integer > 0";
        description = config.thermostatInterval.description;

        int tmp = node[key][config.thermostatInterval.name].as<size_t>();
        if (tmp <= 1) {
          throw YamlParserException("thermostatInterval has to be > 0!");
        }
        config.thermostatInterval.value = node[key][config.thermostatInterval.name].as<size_t>();

        m = node[key][config.targetTemperature.name].Mark();
        expected = "Floating-Point Value";
        description = config.targetTemperature.description;
        config.targetTemperature.value = node[key][config.targetTemperature.name].as<double>();

        m = node[key][config.deltaTemp.name].Mark();
        expected = "Floating-Point Value";
        description = config.deltaTemp.description;
        config.deltaTemp.value = node[key][config.deltaTemp.name].as<double>();

        m = node[key][config.addBrownianMotion.name].Mark();
        expected = "Boolean Value";
        description = config.addBrownianMotion.description;
        config.addBrownianMotion.value = node[key][config.addBrownianMotion.name].as<bool>();

      } else if (key == config.loadBalancer.name) {
        expected = "YAML-sequence of possible values.";
        description = config.loadBalancer.description;

        std::set<LoadBalancerOption> parsedOptions;
        if (node[key].IsSequence()) {
          if (node[key].size() != 1) {
            throw YamlParserException("Pass Exactly one load balancer option!");
          }
          parsedOptions =
              LoadBalancerOption::parseOptions(autopas::utils::ArrayUtils::to_string(node[key], "", {"", ""}));
        } else {
          parsedOptions = LoadBalancerOption::parseOptions(node[key].as<std::string>());
        }
        config.loadBalancer.value = *parsedOptions.begin();
      } else {
        std::cerr << "Unrecognized option in input YAML: " + key << std::endl;
        // return false;
      }
    } catch (const YAML::Exception &e) {
      // We do not use e.mark, as this don't provides the correct line number in some cases. Use the mark from above;
      std::cerr << "Error while parsing the YAML-file in line " << (m.line + 1) << " at column " << m.column
                << ", key: " << key << std::endl
                << "Expected: " << expected << std::endl
                << "Parameter description: " << description << std::endl;
      return false;
    } catch (const YamlParserException &e) {
      std::cerr << "Error while parsing the YAML-file in line " << (m.line + 1) << " at column " << m.column
                << std::endl
                << "Incorrect input-parameter for key " << key << ": " << e.what() << std::endl
                << "Expected: " << expected << std::endl
                << "Parameter description: " << description << std::endl;
      return false;
    } catch (const std::exception &e) {
      std::cerr << "Error while parsing the YAML-file in line " << (m.line + 1) << " at column " << m.column
                << ", key: " << key << std::endl
                << "Message: " << e.what() << std::endl
                << "Expected: " << expected << std::endl
                << "Parameter description: " << description << std::endl;
      return false;
    }
  }

  return true;
}
