/**
 * @file YamlParser.cpp
 * @author N. Fottner, D. Martin
 * @date 15.07.2019, 11.04.2023
 */
#include "YamlParser.h"

const std::string MDFlexParser::YamlParser::parseSequenceOneElementExpected(const YAML::Node node,
                                                                            const std::string &errMsg) {
  if (node.IsSequence()) {
    if (node.size() != 1) {
      throw std::runtime_error(errMsg);
    }
    return autopas::utils::ArrayUtils::to_string(node, "", {"", ""});
  } else {
    return node.as<std::string>();
  }
}

const std::string MDFlexParser::YamlParser::makeErrorMsg(const YAML::Mark &mark, const std::string &key,
                                                         const std::string &errorMsg, const std::string &expected,
                                                         const std::string &description) {
  std::stringstream ss;
  ss << "YamlParser: Parsing error in line " << (mark.line + 1) << " at column " << mark.column << ", key: " << key
     << std::endl
     << "Message: " << errorMsg << std::endl
     << "Expected: " << expected << std::endl
     << "Parameter description: " << description << std::endl;
  return ss.str();
}

const CubeGrid MDFlexParser::YamlParser::parseCubeGridObject(const MDFlexConfig &config, const YAML::Node node,
                                                             std::vector<std::string> &objectErrors) {
  const auto velocity = parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::velocityStr, objectErrors);
  const auto particleType = parseComplexTypeValueSingle<unsigned long>(node, MDFlexConfig::particleTypeStr, objectErrors);
  const auto particlesPerDim =
      parseComplexTypeValueSequence<unsigned long, 3>(node, config.particlesPerDim.name, objectErrors);
  const auto particleSpacing = parseComplexTypeValueSingle<double>(node, config.particleSpacing.name.c_str(), objectErrors);
  const auto bottomLeftCorner =
      parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::bottomLeftBackCornerStr, objectErrors);

  const CubeGrid cubeGrid(velocity, particleType, particlesPerDim, particleSpacing, bottomLeftCorner);
  return cubeGrid;
}

const CubeUniform MDFlexParser::YamlParser::parseCubeUniformObject(const MDFlexConfig &config, const YAML::Node node,
                                                                   std::vector<std::string> &objectErrors) {
  const auto velocity = parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::velocityStr, objectErrors);
  const auto particleType = parseComplexTypeValueSingle<unsigned long>(node, MDFlexConfig::particleTypeStr, objectErrors);
  const auto numParticles = parseComplexTypeValueSingle<size_t>(node, MDFlexConfig::particlesPerObjectStr, objectErrors);
  const auto boxLength = parseComplexTypeValueSequence<double, 3>(node, config.boxLength.name, objectErrors);
  const auto bottomLeftCorner =
      parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::bottomLeftBackCornerStr, objectErrors);

  const CubeUniform cubeUniform(velocity, particleType, numParticles, boxLength, bottomLeftCorner);
  return cubeUniform;
}

const CubeGauss MDFlexParser::YamlParser::parseCubeGaussObject(const MDFlexConfig &config, const YAML::Node node,
                                                               std::vector<std::string> &objectErrors) {
  const auto velocity = parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::velocityStr, objectErrors);
  const auto particleType = parseComplexTypeValueSingle<unsigned long>(node, MDFlexConfig::particleTypeStr, objectErrors);
  const auto numParticles = parseComplexTypeValueSingle<size_t>(node, MDFlexConfig::particlesPerObjectStr, objectErrors);
  const auto boxLength = parseComplexTypeValueSequence<double, 3>(node, config.boxLength.name, objectErrors);
  const auto distributionMean = parseComplexTypeValueSequence<double, 3>(node, config.distributionMean.name, objectErrors);
  const auto distributionStdDev =
      parseComplexTypeValueSequence<double, 3>(node, config.distributionStdDev.name, objectErrors);
  const auto bottomLeftCorner =
      parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::bottomLeftBackCornerStr, objectErrors);

  const CubeGauss cubeGauss(velocity, particleType, numParticles, boxLength, distributionMean, distributionStdDev, bottomLeftCorner);
  return cubeGauss;
}

const Sphere MDFlexParser::YamlParser::parseSphereObject(const MDFlexConfig &config, const YAML::Node node,
                                                         std::vector<std::string> &objectErrors) {
  const auto velocity = parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::velocityStr, objectErrors);
  const auto particleType = parseComplexTypeValueSingle<unsigned long>(node, MDFlexConfig::particleTypeStr, objectErrors);
  const auto sphereCenter = parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::sphereCenterStr, objectErrors);
  const auto sphereRadius = parseComplexTypeValueSingle<double>(node, MDFlexConfig::sphereRadiusStr, objectErrors);
  const auto particleSpacing = parseComplexTypeValueSingle<double>(node, config.particleSpacing.name.c_str(), objectErrors);

  const Sphere sphere(velocity, particleType, sphereCenter, sphereRadius, particleSpacing);
  return sphere;
}

const CubeClosestPacked MDFlexParser::YamlParser::parseCubeClosestPacked(const MDFlexConfig &config,
                                                                         const YAML::Node node,
                                                                         std::vector<std::string> &objectErrors) {
  const auto velocity = parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::velocityStr, objectErrors);
  const auto particleType = parseComplexTypeValueSingle<unsigned long>(node, MDFlexConfig::particleTypeStr, objectErrors);
  const auto particleSpacing = parseComplexTypeValueSingle<double>(node, config.particleSpacing.name.c_str(), objectErrors);
  const auto boxLength = parseComplexTypeValueSequence<double, 3>(node, config.boxLength.name, objectErrors);
  const auto bottomLeftCorner =
      parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::bottomLeftBackCornerStr, objectErrors);

  const CubeClosestPacked cubeClosestPacked(velocity, particleType, particleSpacing, boxLength, bottomLeftCorner);

  return cubeClosestPacked;
}

bool MDFlexParser::YamlParser::parseYamlFile(MDFlexConfig &config) {
  /*
  Global variables used to print the expected input and a description of the parameter if an error occurs while
  parsing. Yaml mark is used to identify the current line of the error.
  */
  std::string expected;
  std::string description;
  YAML::Mark mark;
  std::vector<std::string> errors;

  const auto node = YAML::LoadFile(config.yamlFilename.value);

  // We iterate over all keys to identify known/unknown parameters.
  for (auto itemIterator = node.begin(); itemIterator != node.end(); ++itemIterator) {
    std::string key;
    try {
      key = itemIterator->first.as<std::string>();
      mark = node[key].Mark();

      if (key == config.containerOptions.name) {
        expected = "YAML-sequence of possible values.";
        description = config.containerOptions.description;

        config.containerOptions.value =
            autopas::ContainerOption::parseOptions(autopas::utils::ArrayUtils::to_string(node[key], ", ", {"", ""}));

        if (config.containerOptions.value.empty()) {
          throw std::runtime_error("Parsed container list is empty. You used possibly an unknown container option.");
        }

      } else if (key == config.boxMin.name) {
        expected = "YAML-sequence of three floats. Example: [0, 0, 0].";
        description = config.boxMin.description;

        config.boxMin.value = {node[key][0].as<double>(), node[key][1].as<double>(), node[key][2].as<double>()};
      } else if (key == config.boxMax.name) {
        expected = "YAML-sequence of three floats. Example: [42, 42, 42].";
        description = config.boxMax.description;

        config.boxMax.value = {node[key][0].as<double>(), node[key][1].as<double>(), node[key][2].as<double>()};
      } else if (key == config.subdivideDimension.name) {
        expected = "YAML-sequence of three booleans.";
        description = config.subdivideDimension.description;

        config.subdivideDimension.value = {node[key][0].as<bool>(), node[key][1].as<bool>(), node[key][2].as<bool>()};
      } else if (key == config.loadBalancingInterval.name) {
        expected = "Unsigned Integer";
        description = config.loadBalancingInterval.description;

        config.loadBalancingInterval.value = node[key].as<int>();
        if (config.loadBalancingInterval.value < 0) {
          throw std::runtime_error("Load balancing interval must be a positive integer.");
        }
      } else if (key == config.selectorStrategy.name) {
        expected = "Exactly one selector strategy out of the possible values.";
        description = config.selectorStrategy.description;

        const auto parsedOptions = autopas::SelectorStrategyOption::parseOptions(
            parseSequenceOneElementExpected(node[key], "Pass Exactly one selector strategy!"));

        config.selectorStrategy.value = *parsedOptions.begin();

      } else if (key == config.boundaryOption.name) {
        expected = "YAML-sequence of three possible values.";
        description = config.boundaryOption.description;

        config.boundaryOption.value = {options::BoundaryTypeOption::parseOptionExact(node[key][0].as<std::string>()),
                                       options::BoundaryTypeOption::parseOptionExact(node[key][1].as<std::string>()),
                                       options::BoundaryTypeOption::parseOptionExact(node[key][2].as<std::string>())};
      } else if (key == config.cutoff.name) {
        expected = "Positive floating point value > 0.";
        description = config.cutoff.description;

        config.cutoff.value = node[key].as<double>();
        if (config.cutoff.value <= 0) {
          throw std::runtime_error("Cutoff has to be > 0!");
        }
      } else if (key == config.cellSizeFactors.name) {
        expected = "YAML-sequence of floats.";
        description = config.cellSizeFactors.description;

        config.cellSizeFactors.value = autopas::utils::StringUtils::parseNumberSet(
            autopas::utils::ArrayUtils::to_string(node[key], ", ", {"", ""}));

        if (config.cellSizeFactors.value->isEmpty()) {
          throw std::runtime_error("Parsed cell-size-factor-list is empty.");
        }
      } else if (key == config.dataLayoutOptions.name) {
        expected = "YAML-sequence of possible values.";
        description = config.dataLayoutOptions.description;

        config.dataLayoutOptions.value =
            autopas::DataLayoutOption::parseOptions(autopas::utils::ArrayUtils::to_string(node[key], ", ", {"", ""}));
        if (config.dataLayoutOptions.value.empty()) {
          throw std::runtime_error("Parsed data-layouts-list is empty.");
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
          throw std::runtime_error("Unrecognized functor!");
        }
      } else if (key == config.iterations.name) {
        expected = "Unsigned Integer > 0";
        description = config.iterations.description;

        config.iterations.value = node[key].as<long>();
        if (config.iterations.value < 1) {
          throw std::runtime_error("The number of iterations has to be a positive integer > 0.");
        }

      } else if (key == config.tuningPhases.name) {
        expected = "Unsigned Integer";
        description = config.tuningPhases.description;

        config.tuningPhases.value = node[key].as<long>();
        if (config.tuningPhases.value < 0) {
          throw std::runtime_error("The number of tuning phases has to be a positive integer.");
        }
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
          throw std::runtime_error("Unknown Newton3 option!");
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
          throw std::runtime_error("Parsed traversal-list is empty. Maybe you used an unknown option.");
        }

      } else if (key == config.loadEstimatorOptions.name) {
        expected = "YAML-sequence of possible values.";
        description = config.loadEstimatorOptions.description;

        config.loadEstimatorOptions.value = autopas::LoadEstimatorOption::parseOptions(
            autopas::utils::ArrayUtils::to_string(node[key], ", ", {"", ""}));

        if (config.loadEstimatorOptions.value.empty()) {
          throw std::runtime_error("Parsed load-estimator-list is empty. Maybe you used an unknown option.");
        }

      } else if (key == config.tuningInterval.name) {
        expected = "Unsigned Integer";
        description = config.tuningInterval.description;

        config.tuningInterval.value = node[key].as<int>();
        if (config.tuningInterval.value < 1) {
          throw std::runtime_error("Tuning interval has to be a positive integer!");
        }
      } else if (key == config.tuningSamples.name) {
        expected = "Unsigned Integer >= 1";
        description = config.tuningSamples.description;

        config.tuningSamples.value = node[key].as<int>();
        if (config.tuningSamples.value < 1) {
          throw std::runtime_error("Tuning samples has to be a positive integer!");
        }
      } else if (key == config.tuningMaxEvidence.name) {
        expected = "Unsigned Integer >= 1";
        description = config.tuningMaxEvidence.description;

        config.tuningMaxEvidence.value = node[key].as<int>();
        if (config.tuningMaxEvidence.value < 1) {
          throw std::runtime_error("Tuning max evidence has to be a positive integer >= 1!");
        }
      } else if (key == config.relativeOptimumRange.name) {
        expected = "Floating point value >= 1";
        description = config.relativeOptimumRange.description;

        config.relativeOptimumRange.value = node[key].as<double>();
        if (config.relativeOptimumRange.value < 1.0) {
          throw std::runtime_error("Relative optimum range has to be greater or equal one!");
        }
      } else if (key == config.maxTuningPhasesWithoutTest.name) {
        expected = "Unsigned Integer";
        description = config.maxTuningPhasesWithoutTest.description;

        config.maxTuningPhasesWithoutTest.value = node[key].as<int>();
        if (config.maxTuningPhasesWithoutTest.value < 1) {
          throw std::runtime_error("Max tuning phases without test has to be positive!");
        }
      } else if (key == config.relativeBlacklistRange.name) {
        expected = "Floating point value >= 1 or 0";
        description = config.relativeBlacklistRange.description;

        config.relativeBlacklistRange.value = node[key].as<double>();
        if (config.relativeBlacklistRange.value < 1.0 and config.relativeBlacklistRange.value != 0.0) {
          throw std::runtime_error(
              "Relative range for blacklist range has to be greater or equal one or has to be zero!");
        }
      } else if (key == config.evidenceFirstPrediction.name) {
        expected = "Unsigned Integer >= 2";
        description = config.evidenceFirstPrediction.description;

        config.evidenceFirstPrediction.value = node[key].as<int>();
        if (config.evidenceFirstPrediction.value < 2) {
          throw std::runtime_error("The number of evidence for the first prediction has to be at least two!");
        }
      } else if (key == config.extrapolationMethodOption.name) {
        expected = "Exactly one extrapolation method out of the possible values.";
        description = config.extrapolationMethodOption.description;

        const auto parsedOptions = autopas::ExtrapolationMethodOption::parseOptions(
            parseSequenceOneElementExpected(node[key], "Pass exactly one extrapolation method!"));

        config.extrapolationMethodOption.value = *parsedOptions.begin();

      } else if (key == config.tuningStrategyOption.name) {
        expected = "Exactly one tuning strategy option out of the possible values.";
        description = config.tuningStrategyOption.description;

        const auto parsedOptions = autopas::TuningStrategyOption::parseOptions(
            parseSequenceOneElementExpected(node[key], "Pass Exactly one tuning strategy!"));

        config.tuningStrategyOption.value = *parsedOptions.begin();
      } else if (key == config.mpiStrategyOption.name) {
        expected = "Exactly one MPI strategy option out of the possible values.";
        description = config.mpiStrategyOption.description;

        const auto parsedOptions = autopas::MPIStrategyOption::parseOptions(
            parseSequenceOneElementExpected(node[key], "Pass exactly one MPI strategy!"));

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

        const auto parsedOptions = autopas::AcquisitionFunctionOption::parseOptions(
            parseSequenceOneElementExpected(node[key], "Pass Exactly one acquisition function option!"));

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
            throw std::runtime_error("Unknown Log Level parsed!");
          }
        }
      } else if (key == config.checkpointfile.name) {
        expected = "String";
        description = config.checkpointfile.description;

        config.checkpointfile.value = node[key].as<std::string>();
        if (config.checkpointfile.value.empty()) {
          throw std::runtime_error("Parsed checkpoint filename is empty!");
        }
      } else if (key == config.logFileName.name) {
        expected = "String";
        description = config.logFileName.description;

        config.logFileName.value = node[key].as<std::string>();
        if (config.logFileName.value.empty()) {
          throw std::runtime_error("Parsed log filename is empty!");
        }
      } else if (key == config.verletRebuildFrequency.name) {
        expected = "Unsigned Integer >= 1";
        description = config.verletRebuildFrequency.description;

        config.verletRebuildFrequency.value = node[key].as<int>();

        if (config.verletRebuildFrequency.value < 1) {
          throw std::runtime_error("Verlet rebuild frequency has to be a positive integer >= 1!");
        }
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

        config.verletClusterSize.value = node[key].as<int>();
        if (config.verletClusterSize.value < 0) {
          throw std::runtime_error("Verlet cluster size has to be a positive integer!");
        }
      } else if (key == config.vtkFileName.name) {
        expected = "String";
        description = config.vtkFileName.description;

        config.vtkFileName.value = node[key].as<std::string>();
        if (config.vtkFileName.value.empty()) {
          throw std::runtime_error("Parsed VTK filename is empty!");
        }
      } else if (key == config.vtkWriteFrequency.name) {
        expected = "Unsigned Integer >= 1";
        description = config.vtkWriteFrequency.description;

        config.vtkWriteFrequency.value = node[key].as<size_t>();
        if (config.vtkWriteFrequency.value < 1) {
          throw std::runtime_error("VTK write frequency has to be a positive integer >= 1!");
        }
      } else if (key == config.globalForce.name) {
        expected = "YAML-sequence of three floats. Example: [0, 0, -9.81].";
        description = config.globalForce.description;

        config.globalForce.value = {node[key][0].as<double>(), node[key][1].as<double>(), node[key][2].as<double>()};
      } else if (key == MDFlexConfig::siteStr) {
        expected = "See AllOptions.yaml for examples.";
        description = "";

        // remove default objects
        config.epsilonMap.value.clear();
        config.sigmaMap.value.clear();
        config.massMap.value.clear();

        int siteID = 0;
        std::vector<std::string> siteErrors;

        auto pushSiteError = [&](const std::string &error) {
          std::stringstream ss;
          ss << "YamlParser: Error parsing site with ID " << siteID << "." << std::endl
             << "Message: " << error << std::endl
             << "See AllOptions.yaml for examples." << std::endl;
          errors.push_back(ss.str());
        };

        for (auto siteIterator = node[MDFlexConfig::siteStr].begin(); siteIterator != node[MDFlexConfig::siteStr].end(); ++siteIterator) {
          siteErrors.clear();
          siteID = std::distance(node[MDFlexConfig::siteStr].begin(), siteIterator);

          const auto epsilon = parseComplexTypeValueSingle<double>(siteIterator->second, config.epsilonMap.name.c_str(), siteErrors);
          const auto sigma = parseComplexTypeValueSingle<double>(siteIterator->second, config.sigmaMap.name.c_str(), siteErrors);
          const auto mass = parseComplexTypeValueSingle<double>(siteIterator->second, config.massMap.name.c_str(), siteErrors);

          config.addSiteType(siteID, epsilon, sigma, mass);
        }
      } else if (key == MDFlexConfig::moleculesStr) {
        // todo throw error if momentOfInertia with zero element is used (physically nonsense + breaks the quaternion update)
        expected = "See AllOptions.yaml for examples.";
        description = "";

        // remove default objects
        config.molToSiteIdMap.clear();
        config.molToSitePosMap.clear();
        config.momentOfInertiaMap.clear();

        int molID = 0;
        std::vector<std::string> molErrors;

        auto pushMolError = [&](const std::string &error) {
          std::stringstream ss;
          ss << "YamlParser: Error parsing multi-site molecule with ID " << molID << "." << std::endl
             << "Message: " << error << std::endl
             << "See AllOptions.yaml for examples." << std::endl;
          errors.push_back(ss.str());
        };

#if MD_FLEXIBLE_MODE==MULTISITE
        for (auto molIterator = node[MDFlexConfig::moleculesStr].begin(); molIterator != node[MDFlexConfig::moleculesStr].end(); ++molIterator) {
          molErrors.clear();
          molID = std::distance(node[MDFlexConfig::moleculesStr].begin(), molIterator);

          const auto molToSiteId = parseComplexTypeValueSequence<unsigned long>(molIterator->second, MDFlexConfig::moleculeToSiteIdStr, molErrors);
          const auto molToSitePos = parseComplexTypeValueSequence<std::array<double, 3>>(molIterator->second, MDFlexConfig::moleculeToSitePosStr, molErrors);
          const auto momentOfInertia = parseComplexTypeValueSequence<double, 3>(molIterator->second, MDFlexConfig::momentOfInertiaStr, molErrors);

          if (molToSiteId.size() != molToSitePos.size()) {
            std::stringstream ss;
            ss << "Lengths of " << MDFlexConfig::moleculeToSiteIdStr << " and " << MDFlexConfig::moleculeToSitePosStr << "do not match!";
            molErrors.push_back(ss.str());
          }

          config.addMolType(molID, molToSiteId, molToSitePos, momentOfInertia);
        }

        for_each(molErrors.begin(), molErrors.end(), pushMolError);
#else
        AutoPasLog(WARN, "Multi-Site Molecule information has been provided, however md-flexible has been compiled without Multi-Site support.");
#endif
      }else if (key == MDFlexConfig::objectsStr) {
        expected = "See AllOptions.yaml for examples.";
        description = "";

        // remove default objects
        config.cubeGridObjects.clear();
        config.cubeGaussObjects.clear();
        config.cubeUniformObjects.clear();
        config.sphereObjects.clear();
        config.cubeClosestPackedObjects.clear();

        int objID = 0;
        std::string generatorName;
        std::vector<std::string> objectErrors;

        auto pushObjectError = [&](const std::string &error) {
          std::stringstream ss;
          ss << "YamlParser: Error parsing " << generatorName << " object with ID " << objID << "." << std::endl
             << "Message: " << error << std::endl
             << "See AllOptions.yaml for examples." << std::endl;
          errors.push_back(ss.str());
        };

        for (auto objectIterator = node[MDFlexConfig::objectsStr].begin();
             objectIterator != node[MDFlexConfig::objectsStr].end(); ++objectIterator) {
          if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeGridObjectsStr) {
            generatorName = MDFlexConfig::cubeGridObjectsStr;
            for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
              objectErrors.clear();
              objID = std::distance(objectIterator->second.begin(), it);
              const auto cubeGrid = parseCubeGridObject(config, it->second, objectErrors);

              config.cubeGridObjects.emplace_back(cubeGrid);
              for_each(objectErrors.begin(), objectErrors.end(), pushObjectError);
            }
          } else if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeGaussObjectsStr) {
            generatorName = MDFlexConfig::cubeGaussObjectsStr;
            for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
              objectErrors.clear();
              objID = std::distance(objectIterator->second.begin(), it);
              const auto cubeGauss = parseCubeGaussObject(config, it->second, objectErrors);

              config.cubeGaussObjects.emplace_back(cubeGauss);
              for_each(objectErrors.begin(), objectErrors.end(), pushObjectError);
            }
          } else if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeUniformObjectsStr) {
            generatorName = MDFlexConfig::cubeUniformObjectsStr;
            for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
              objectErrors.clear();
              objID = std::distance(objectIterator->second.begin(), it);
              const auto cubeUniform = parseCubeUniformObject(config, it->second, objectErrors);

              config.cubeUniformObjects.emplace_back(cubeUniform);
              for_each(objectErrors.begin(), objectErrors.end(), pushObjectError);
            }
          } else if (objectIterator->first.as<std::string>() == MDFlexConfig::sphereObjectsStr) {
            generatorName = MDFlexConfig::sphereObjectsStr;
            for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
              objectErrors.clear();
              objID = std::distance(objectIterator->second.begin(), it);
              const auto sphere = parseSphereObject(config, it->second, objectErrors);

              config.sphereObjects.emplace_back(sphere);
              for_each(objectErrors.begin(), objectErrors.end(), pushObjectError);
            }
          } else if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeClosestPackedObjectsStr) {
            generatorName = MDFlexConfig::cubeClosestPackedObjectsStr;
            for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
              objectErrors.clear();
              objID = std::distance(objectIterator->second.begin(), it);
              const auto cubeClosestPacked = parseCubeClosestPacked(config, it->second, objectErrors);

              config.cubeClosestPackedObjects.emplace_back(cubeClosestPacked);
              for_each(objectErrors.begin(), objectErrors.end(), pushObjectError);
            }
          } else {
            std::stringstream ss;
            ss << "YamlParser: Unrecognized generator \"" << objectIterator->first.as<std::string>() << "\" used."
               << std::endl;
            errors.push_back(ss.str());
          }
        }
      } else if (key == config.useThermostat.name) {
        expected = "See AllOptions.yaml for examples.";
        description = config.useThermostat.description;

        config.useThermostat.value = true;

        // Parse initTemperature.
        mark = node[key][config.initTemperature.name].Mark();
        expected = "Floating-Point Value";
        description = config.initTemperature.description;
        try {
          config.initTemperature.value = node[key][config.initTemperature.name].as<double>();
        } catch (const std::exception &e) {
          errors.push_back(makeErrorMsg(mark, key, e.what(), expected, description));
        }

        // Parse thermostatInterval.
        mark = node[key][config.thermostatInterval.name].Mark();
        expected = "Unsigned Integer > 0";
        description = config.thermostatInterval.description;
        try {
          config.thermostatInterval.value = node[key][config.thermostatInterval.name].as<size_t>();
          if (config.thermostatInterval.value <= 1) {
            throw std::runtime_error("thermostatInterval has to be > 0!");
          }
        } catch (const std::exception &e) {
          errors.push_back(makeErrorMsg(mark, key, e.what(), expected, description));
        }

        // Parse targetTemperature.
        mark = node[key][config.targetTemperature.name].Mark();
        expected = "Floating-Point Value";
        description = config.targetTemperature.description;
        try {
          config.targetTemperature.value = node[key][config.targetTemperature.name].as<double>();
        } catch (const std::exception &e) {
          errors.push_back(makeErrorMsg(mark, key, e.what(), expected, description));
        }

        // Parse deltaTemp.
        mark = node[key][config.deltaTemp.name].Mark();
        expected = "Floating-Point Value";
        description = config.deltaTemp.description;
        try {
          config.deltaTemp.value = node[key][config.deltaTemp.name].as<double>();
        } catch (const std::exception &e) {
          errors.push_back(makeErrorMsg(mark, key, e.what(), expected, description));
        }

        // Parse addBrownianMotion.
        mark = node[key][config.addBrownianMotion.name].Mark();
        expected = "Boolean Value";
        description = config.addBrownianMotion.description;
        try {
          config.addBrownianMotion.value = node[key][config.addBrownianMotion.name].as<bool>();
        } catch (const std::exception &e) {
          errors.push_back(makeErrorMsg(mark, key, e.what(), expected, description));
        }

      } else if (key == config.loadBalancer.name) {
        expected = "YAML-sequence of possible values.";
        description = config.loadBalancer.description;

        const auto parsedOptions = LoadBalancerOption::parseOptions(
            parseSequenceOneElementExpected(node[key], "Pass Exactly one load balancer option!"));

        config.loadBalancer.value = *parsedOptions.begin();
      } else {
        std::stringstream ss;
        ss << "YamlParser: Unrecognized option in input YAML: " + key << std::endl;
        errors.push_back(ss.str());
      }
    } catch (const std::exception &e) {
      errors.push_back(makeErrorMsg(mark, key, e.what(), expected, description));
    }
  }

  if (!errors.empty()) {
    for (std::string &err : errors) {
      std::cerr << err << std::endl << std::endl;
    }
    return false;
  }

  return true;
}