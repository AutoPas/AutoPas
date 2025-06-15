/**
 * @file YamlParser.cpp
 * @author N. Fottner, D. Martin
 * @date 15.07.2019, 11.04.2023
 */
#include "YamlParser.h"

#include "autopas/options/TuningMetricOption.h"

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
  const auto particleType =
      parseComplexTypeValueSingle<unsigned long>(node, MDFlexConfig::particleTypeStr, objectErrors);
  const auto particlesPerDim =
      parseComplexTypeValueSequence<unsigned long, 3>(node, config.particlesPerDim.name, objectErrors);
  const auto particleSpacing =
      parseComplexTypeValueSingle<double>(node, config.particleSpacing.name.c_str(), objectErrors);
  const auto bottomLeftCorner =
      parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::bottomLeftBackCornerStr, objectErrors);

  const CubeGrid cubeGrid(velocity, particleType, particlesPerDim, particleSpacing, bottomLeftCorner);
  return cubeGrid;
}

const CubeUniform MDFlexParser::YamlParser::parseCubeUniformObject(const MDFlexConfig &config, const YAML::Node node,
                                                                   std::vector<std::string> &objectErrors) {
  const auto velocity = parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::velocityStr, objectErrors);
  const auto particleType =
      parseComplexTypeValueSingle<unsigned long>(node, MDFlexConfig::particleTypeStr, objectErrors);
  const auto numParticles =
      parseComplexTypeValueSingle<size_t>(node, MDFlexConfig::particlesPerObjectStr, objectErrors);
  const auto boxLength = parseComplexTypeValueSequence<double, 3>(node, config.boxLength.name, objectErrors);
  const auto bottomLeftCorner =
      parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::bottomLeftBackCornerStr, objectErrors);

  const CubeUniform cubeUniform(velocity, particleType, numParticles, boxLength, bottomLeftCorner);
  return cubeUniform;
}

const CubeGauss MDFlexParser::YamlParser::parseCubeGaussObject(const MDFlexConfig &config, const YAML::Node node,
                                                               std::vector<std::string> &objectErrors) {
  const auto velocity = parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::velocityStr, objectErrors);
  const auto particleType =
      parseComplexTypeValueSingle<unsigned long>(node, MDFlexConfig::particleTypeStr, objectErrors);
  const auto numParticles =
      parseComplexTypeValueSingle<size_t>(node, MDFlexConfig::particlesPerObjectStr, objectErrors);
  const auto boxLength = parseComplexTypeValueSequence<double, 3>(node, config.boxLength.name, objectErrors);
  const auto distributionMean =
      parseComplexTypeValueSequence<double, 3>(node, config.distributionMean.name, objectErrors);
  const auto distributionStdDev =
      parseComplexTypeValueSequence<double, 3>(node, config.distributionStdDev.name, objectErrors);
  const auto bottomLeftCorner =
      parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::bottomLeftBackCornerStr, objectErrors);

  const CubeGauss cubeGauss(velocity, particleType, numParticles, boxLength, distributionMean, distributionStdDev,
                            bottomLeftCorner);
  return cubeGauss;
}

const Sphere MDFlexParser::YamlParser::parseSphereObject(const MDFlexConfig &config, const YAML::Node node,
                                                         std::vector<std::string> &objectErrors) {
  const auto velocity = parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::velocityStr, objectErrors);
  const auto particleType =
      parseComplexTypeValueSingle<unsigned long>(node, MDFlexConfig::particleTypeStr, objectErrors);
  const auto sphereCenter = parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::sphereCenterStr, objectErrors);
  const auto sphereRadius = parseComplexTypeValueSingle<double>(node, MDFlexConfig::sphereRadiusStr, objectErrors);
  const auto particleSpacing =
      parseComplexTypeValueSingle<double>(node, config.particleSpacing.name.c_str(), objectErrors);

  const Sphere sphere(velocity, particleType, sphereCenter, sphereRadius, particleSpacing);
  return sphere;
}

const CubeClosestPacked MDFlexParser::YamlParser::parseCubeClosestPacked(const MDFlexConfig &config,
                                                                         const YAML::Node node,
                                                                         std::vector<std::string> &objectErrors) {
  const auto velocity = parseComplexTypeValueSequence<double, 3>(node, MDFlexConfig::velocityStr, objectErrors);
  const auto particleType =
      parseComplexTypeValueSingle<unsigned long>(node, MDFlexConfig::particleTypeStr, objectErrors);
  const auto particleSpacing =
      parseComplexTypeValueSingle<double>(node, config.particleSpacing.name.c_str(), objectErrors);
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
      } else if (key == config.cutoffFactorElectrostatics.name) {
        expected = "Positive floating point value > 0.";
        description = config.cutoffFactorElectrostatics.description;

        config.cutoffFactorElectrostatics.value = node[key].as<double>();
        if (config.cutoffFactorElectrostatics.value <= 0) {
          throw std::runtime_error("Cutoff factor has to be > 0!");
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
      } else if (key == config.dataLayoutOptions3B.name) {
        expected = "YAML-sequence of possible values.";
        description = config.dataLayoutOptions3B.description;

        config.dataLayoutOptions3B.value =
            autopas::DataLayoutOption::parseOptions(autopas::utils::ArrayUtils::to_string(node[key], ", ", {"", ""}));
        if (config.dataLayoutOptions3B.value.empty()) {
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
        } else if (strArg.find("lj") != std::string::npos or strArg.find("lennard-jones") != std::string::npos) {
          config.functorOption.value = MDFlexConfig::FunctorOption::lj12_6;
        } else {
          throw std::runtime_error("Unrecognized pairwise functor!");
        }
        config.addInteractionType(autopas::InteractionTypeOption::pairwise);
      } else if (key == config.functorOption3B.name) {
        expected = "One of the possible values.";
        description = config.functorOption3B.description;

        auto strArg = node[key].as<std::string>();
        transform(strArg.begin(), strArg.end(), strArg.begin(), ::tolower);
        if (strArg.find("at") != std::string::npos or strArg.find("axilrod-teller") != std::string::npos) {
          config.functorOption3B.value = MDFlexConfig::FunctorOption3B::at;
        } else {
          throw std::runtime_error("Unrecognized triwise functor!");
        }
        config.addInteractionType(autopas::InteractionTypeOption::triwise);
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
      } else if (key == config.newton3Options3B.name) {
        expected = "YAML-sequence of possible values.";
        description = config.newton3Options3B.description;

        config.newton3Options3B.value =
            autopas::Newton3Option::parseOptions(autopas::utils::ArrayUtils::to_string(node[key], ", ", {"", ""}));
        if (config.newton3Options3B.value.empty()) {
          throw std::runtime_error("Unknown Newton3 option!");
        }
      } else if (key == config.deltaT.name) {
        expected = "Positive floating point value.";
        description = config.deltaT.description;

        config.deltaT.value = node[key].as<double>();
      } else if (key == config.energySensorOption.name) {
        expected = "Exactly one energy sensor out of the possible options.";
        description = config.energySensorOption.description;
        const auto parsedOptions = autopas::EnergySensorOption::parseOptions(
            parseSequenceOneElementExpected(node[key], "Pass Exactly one energy sensor!"));
        config.energySensorOption.value = *parsedOptions.begin();
      } else if (key == config.pauseSimulationDuringTuning.name) {
        expected = "Boolean Value";
        description = config.pauseSimulationDuringTuning.description;

        config.pauseSimulationDuringTuning.value = node[key].as<bool>();
      } else if (key == config.sortingThreshold.name) {
        expected = "Unsigned Integer >= 0.";
        description = config.sortingThreshold.description;

        config.sortingThreshold.value = node[key].as<size_t>();
      } else if (key == config.respaStepSize.name) {
        expected = "Unsigned Integer >= 1.";
        description = config.respaStepSize.description;

        config.respaStepSize.value = node[key].as<size_t>();
      } else if (key == config.rotationalAnalysisLagSteps.name) {
        expected = "Unsigned Integer >= 0.";
        description = config.rotationalAnalysisLagSteps.description;

        config.rotationalAnalysisLagSteps.value = node[key].as<size_t>();
      } else if (key == config.rotationalAnalysisStepInterval.name) {
        expected = "Unsigned Integer >= 0.";
        description = config.rotationalAnalysisStepInterval.description;

        config.rotationalAnalysisStepInterval.value = node[key].as<size_t>();
      } else if (key == config.rotationalAnalysisFilename.name) {
        expected = "String.";
        description = config.rotationalAnalysisFilename.description;

        config.rotationalAnalysisFilename.value = node[key].as<std::string>();
      } else if (key == config.rotationalAnalysisOutputFolder.name) {
        expected = "String.";
        description = config.rotationalAnalysisOutputFolder.description;

        config.rotationalAnalysisOutputFolder.value = node[key].as<std::string>();
      } else if (key == config.rotationalAnalysisStartIteration.name) {
        expected = "Unsigned Integer >= 0.";
        description = config.rotationalAnalysisStartIteration.description;

        config.rotationalAnalysisStartIteration.value = node[key].as<size_t>();
      } else if (key == config.rotationalAnalysisEndIteration.name) {
        expected = "Unsigned Integer >= 0.";
        description = config.rotationalAnalysisEndIteration.description;

        config.rotationalAnalysisEndIteration.value = node[key].as<size_t>();
      } else if (key == config.useApproxForceRespa.name) {
        expected = "Boolean Value.";
        description = config.useApproxForceRespa.description;

        config.useApproxForceRespa.value = node[key].as<bool>();
      } else if (key == config.useSecondAutpasInstance.name) {
        expected = "Boolean Value.";
        description = config.useSecondAutpasInstance.description;

        config.useSecondAutpasInstance.value = node[key].as<bool>();
      } else if (key == config.multiMultisiteModelsRespa.name) {  //
        expected = "Boolean Value.";
        description = config.multiMultisiteModelsRespa.description;

        config.multiMultisiteModelsRespa.value = node[key].as<bool>();
      } else if (key == config.respaMoleculeTypes.name) {
        expected = "YAML-sequence of possible values.";
        description = config.respaMoleculeTypes.description;

        std::vector<std::string> errors;
        config.respaMoleculeTypes.value =
            parseComplexTypeValueSequence<unsigned long>(node, "respa-molecule-types", errors);
        if (not errors.empty()) {
          throw std::runtime_error(errors[0]);
        }
        if (config.respaMoleculeTypes.value.empty()) {
          throw std::runtime_error("Parsed respa-molecule-type-list is empty.");
        }
      } else if (key == config.respaDistanceClassMode.name) {
        expected = "YAML-sequence of possible values.";
        description = config.respaDistanceClassMode.description;

        auto strArg = node[key].as<std::string>();

        config.respaDistanceClassMode.value = autopas::DistanceClassOption::parseOptionExact(strArg);
      } else if (key == config.traversalOptions.name) {
        expected = "YAML-sequence of possible values.";
        description = config.traversalOptions.description;

        config.traversalOptions.value =
            autopas::TraversalOption::parseOptions(autopas::utils::ArrayUtils::to_string(node[key], ", ", {"", ""}));

        if (config.traversalOptions.value.empty()) {
          throw std::runtime_error("Parsed traversal-list is empty. Maybe you used an unknown option.");
        }

      } else if (key == config.traversalOptions3B.name) {
        expected = "YAML-sequence of possible values.";
        description = config.traversalOptions3B.description;

        config.traversalOptions3B.value =
            autopas::TraversalOption::parseOptions(autopas::utils::ArrayUtils::to_string(node[key], ", ", {"", ""}));

        if (config.traversalOptions3B.value.empty()) {
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
      } else if (key == config.earlyStoppingFactor.name) {
        expected = "Floating point value > 1";
        description = config.earlyStoppingFactor.description;

        config.earlyStoppingFactor.value = node[key].as<double>();
        if (config.earlyStoppingFactor.value <= 1) {
          throw std::runtime_error("EarlyStoppingFactor has to be greater than 1!");
        }
      } else if (key == config.useLOESSSmoothening.name) {
        expected = "Boolean Value";
        description = config.useLOESSSmoothening.description;

        config.useLOESSSmoothening.value = node[key].as<bool>();
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

      } else if (key == config.tuningStrategyOptions.name) {
        expected = "List of tuning strategies that will be applied in the given order.";
        description = config.tuningStrategyOptions.description;

        config.tuningStrategyOptions.value =
            autopas::TuningStrategyOption::parseOptions<std::vector<autopas::TuningStrategyOption>>(
                autopas::utils::ArrayUtils::to_string(node[key], ", ", {"", ""}));
      } else if (key == config.tuningMetricOption.name) {
        expected = "Exactly one tuning metric option out of the possible values.";
        description = config.tuningMetricOption.description;

        const auto parsedOptions = autopas::TuningMetricOption::parseOptions(
            parseSequenceOneElementExpected(node[key], "Pass Exactly one tuning metric!"));

        config.tuningMetricOption.value = *parsedOptions.begin();
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
      } else if (key == config.ruleFilename.name) {
        expected = "String";
        description = config.ruleFilename.description;

        config.ruleFilename.value = node[key].as<std::string>();
        if (config.ruleFilename.value.empty()) {
          throw std::runtime_error("Parsed rule filename is empty!");
        }
      } else if (key == config.fuzzyRuleFilename.name) {
        expected = "String";
        description = config.fuzzyRuleFilename.description;

        config.fuzzyRuleFilename.value = node[key].as<std::string>();
        if (config.fuzzyRuleFilename.value.empty()) {
          throw std::runtime_error("Parsed rule filename is empty!");
        }
      } else if (key == config.verletRebuildFrequency.name) {
        expected = "Unsigned Integer >= 1";
        description = config.verletRebuildFrequency.description;

        config.verletRebuildFrequency.value = node[key].as<int>();

        if (config.verletRebuildFrequency.value < 1) {
          throw std::runtime_error("Verlet rebuild frequency has to be a positive integer >= 1!");
        }
      } else if (key == config.verletSkinRadius.name) {
        expected = "Positive floating-point value.";
        description = config.verletSkinRadius.description;

        config.verletSkinRadius.value = node[key].as<double>();
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
      } else if (key == config.vtkOutputFolder.name) {
        expected = "String";
        description = config.vtkOutputFolder.description;

        config.vtkOutputFolder.value = node[key].as<std::string>();
        if (config.vtkOutputFolder.value.empty()) {
          throw std::runtime_error("Parsed VTK output folder name is empty");
        }
      } else if (key == config.vtkWriteFrequency.name) {
        expected = "Unsigned Integer >= 1";
        description = config.vtkWriteFrequency.description;

        config.vtkWriteFrequency.value = node[key].as<size_t>();
        if (config.vtkWriteFrequency.value < 1) {
          throw std::runtime_error("VTK write frequency has to be a positive integer >= 1!");
        }
      } else if (key == config.statisticsOutputFolder.name) {
        expected = "String";
        description = config.statisticsOutputFolder.description;

        config.statisticsOutputFolder.value = node[key].as<std::string>();
        if (config.statisticsOutputFolder.value.empty()) {
          throw std::runtime_error("Parsed statistics output folder name is empty");
        }
      } else if (key == config.statisticsOutputFilename.name) {
        expected = "String";
        description = config.statisticsOutputFilename.description;

        config.statisticsOutputFilename.value = node[key].as<std::string>();
        if (config.statisticsOutputFilename.value.empty()) {
          throw std::runtime_error("Parsed statistics filename is empty!");
        }
      } else if (key == config.rdfOutputFolder.name) {
        expected = "String";
        description = config.rdfOutputFolder.description;

        config.rdfOutputFolder.value = node[key].as<std::string>();
        if (config.rdfOutputFolder.value.empty()) {
          throw std::runtime_error("Parsed RDF output folder name is empty");
        }
      } else if (key == config.rdfFileName.name) {
        expected = "String";
        description = config.rdfFileName.description;

        config.rdfFileName.value = node[key].as<std::string>();
        if (config.rdfFileName.value.empty()) {
          throw std::runtime_error("Parsed RDF filename is empty!");
        }
      } else if (key == config.rdfNumBins.name) {
        expected = "Unsigned Integer >= 1";
        description = config.rdfNumBins.description;

        config.rdfNumBins.value = node[key].as<size_t>();
        if (config.rdfNumBins.value < 1) {
          throw std::runtime_error("RDF num bins has to be a positive integer >= 1!");
        }
      } else if (key == config.rdfCaptureFreuency.name) {
        expected = "Unsigned Integer >= 1";
        description = config.rdfCaptureFreuency.description;

        config.rdfCaptureFreuency.value = node[key].as<size_t>();
        if (config.rdfCaptureFreuency.value < 1) {
          throw std::runtime_error("RDF write frequency has to be a positive integer >= 1!");
        }
      } else if (key == config.rdfStartIteration.name) {
        expected = "Unsigned Integer >= 0";
        description = config.rdfStartIteration.description;

        config.rdfStartIteration.value = node[key].as<size_t>();
        if (config.rdfStartIteration.value < 0) {
          throw std::runtime_error("RDF start iteration has to be a positive integer >= 0!");
        }
      } else if (key == config.rdfEndIteration.name) {
        expected = "Unsigned Integer >= 0";
        description = config.rdfEndIteration.description;

        config.rdfEndIteration.value = node[key].as<size_t>();
        if (config.rdfEndIteration.value < 0) {
          throw std::runtime_error("RDF end iteration has to be a positive integer >= 0!");
        }
      } else if (key == config.rdfRadius.name) {
        expected = "Positive floating-point value.";
        description = config.rdfRadius.description;

        config.rdfRadius.value = node[key].as<double>();
        if (config.rdfRadius.value <= 0) {
          throw std::runtime_error("RDF radius has to be a positive value > 0!");
        }
      } else if (key == config.rdfGuardArea.name) {
        expected = "Positive floating-point value.";
        description = config.rdfGuardArea.description;

        config.odfGuardArea.value = node[key].as<double>();
        if (config.odfGuardArea.value < 0) {
          throw std::runtime_error("ODF guard area has to be a positive value >= 0!");
        }
      } else if (key == config.odfOutputFolder.name) {
        expected = "String";
        description = config.odfOutputFolder.description;

        config.odfOutputFolder.value = node[key].as<std::string>();
        if (config.odfOutputFolder.value.empty()) {
          throw std::runtime_error("Parsed ODF output folder name is empty");
        }
      } else if (key == config.odfFileName.name) {
        expected = "String";
        description = config.odfFileName.description;

        config.odfFileName.value = node[key].as<std::string>();
        if (config.odfFileName.value.empty()) {
          throw std::runtime_error("Parsed ODF filename is empty!");
        }
      } else if (key == config.odfNumBins.name) {
        expected = "Unsigned Integer >= 1";
        description = config.odfNumBins.description;

        config.odfNumBins.value = node[key].as<size_t>();
        if (config.odfNumBins.value < 1) {
          throw std::runtime_error("ODF num bins has to be a positive integer >= 1!");
        }
      } else if (key == config.odfCaptureFreuency.name) {
        expected = "Unsigned Integer >= 1";
        description = config.odfCaptureFreuency.description;

        config.odfCaptureFreuency.value = node[key].as<size_t>();
        if (config.odfCaptureFreuency.value < 1) {
          throw std::runtime_error("ODF write frequency has to be a positive integer >= 1!");
        }
      } else if (key == config.odfStartIteration.name) {
        expected = "Unsigned Integer >= 0";
        description = config.odfStartIteration.description;

        config.odfStartIteration.value = node[key].as<size_t>();
        if (config.odfStartIteration.value < 0) {
          throw std::runtime_error("ODF start iteration has to be a positive integer >= 0!");
        }
      } else if (key == config.odfEndIteration.name) {
        expected = "Unsigned Integer >= 0";
        description = config.odfEndIteration.description;

        config.odfEndIteration.value = node[key].as<size_t>();
        if (config.odfEndIteration.value < 0) {
          throw std::runtime_error("ODF end iteration has to be a positive integer >= 0!");
        }
      } else if (key == config.odfRadius.name) {
        expected = "Positive floating-point value.";
        description = config.odfRadius.description;

        config.odfRadius.value = node[key].as<double>();
        if (config.odfRadius.value <= 0) {
          throw std::runtime_error("ODF radius has to be a positive value > 0!");
        }
      } else if (key == config.odfGuardArea.name) {
        expected = "Positive floating-point value.";
        description = config.odfGuardArea.description;

        config.odfGuardArea.value = node[key].as<double>();
        if (config.odfGuardArea.value < 0) {
          throw std::runtime_error("ODF guard area has to be a positive value >= 0!");
        }
      } else if (key == config.ibiEquilibrateIterations.name) {
        expected = "Unsigned Integer >= 0";
        description = config.ibiEquilibrateIterations.description;

        config.ibiEquilibrateIterations.value = node[key].as<size_t>();
        if (config.ibiEquilibrateIterations.value < 0) {
          throw std::runtime_error("IBI equilibrate iterations has to be a positive integer >= 0!");
        }
      } else if (key == config.ibiConvergenceThreshold.name) {
        expected = "double > 0 and <= 1";
        description = config.ibiConvergenceThreshold.description;

        config.ibiConvergenceThreshold.value = node[key].as<double>();
        if (config.ibiConvergenceThreshold.value <= 0 or config.ibiConvergenceThreshold.value > 1.0) {
          throw std::runtime_error("IBI convergence threshold has to be a positive double > 0 and <= 1!");
        }
      } else if (key == config.ibiUpdateAlpha.name) {
        expected = "double > 0";
        description = config.ibiUpdateAlpha.description;

        config.ibiUpdateAlpha.value = node[key].as<double>();
        if (config.ibiUpdateAlpha.value < 0) {
          throw std::runtime_error("IBI update alpha has to be a positive double > 0!");
        }
      } else if (key == config.lutOutputFolder.name) {
        expected = "String";
        description = config.lutOutputFolder.description;

        config.lutOutputFolder.value = node[key].as<std::string>();
        if (config.lutOutputFolder.value.empty()) {
          throw std::runtime_error("Parsed Lookup table output folder is empty");
        }
      } else if (key == config.lutInputFile.name) {
        expected = "String";
        description = config.lutInputFile.description;

        config.lutInputFile.value = node[key].as<std::string>();
        std::cout << "parsing file " << config.lutInputFile.value << std::endl;
        if (config.lutInputFile.value.empty()) {
          throw std::runtime_error("Parsed Lookup table input file is empty");
        }
      } else if (key == config.lutFileName.name) {
        expected = "String";
        description = config.lutFileName.description;

        config.lutFileName.value = node[key].as<std::string>();
        if (config.lutFileName.value.empty()) {
          throw std::runtime_error("Parsed Lookup table file name is empty");
        }
      } else if (key == config.useTuningLogger.name) {
        expected = "Boolean Value";
        description = config.useTuningLogger.description;

        config.useTuningLogger.value = node[config.useTuningLogger.name].as<bool>();
      } else if (key == config.outputSuffix.name) {
        expected = "String";
        description = config.outputSuffix.description;

        config.outputSuffix.value = node[config.outputSuffix.name].as<std::string>();
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
        config.nuMap.value.clear();
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

        for (auto siteIterator = node[MDFlexConfig::siteStr].begin(); siteIterator != node[MDFlexConfig::siteStr].end();
             ++siteIterator) {
          siteErrors.clear();
          siteID = std::distance(node[MDFlexConfig::siteStr].begin(), siteIterator);

          // Check Coulomb parameters
          const auto charge = parseComplexTypeValueSingle<double>(siteIterator->second, config.chargeMap.name.c_str(),
                                                                  siteErrors, false);
          const auto coulombEpsilon = parseComplexTypeValueSingle<double>(
              siteIterator->second, config.coulombEpsilonMap.name.c_str(), siteErrors, false);

          config.addCoulombParametersToSite(siteID, coulombEpsilon, charge);

          const auto mass =
              parseComplexTypeValueSingle<double>(siteIterator->second, config.massMap.name.c_str(), siteErrors);

          config.addSiteType(siteID, mass);
          // Check LJ parameters
          const auto epsilon = parseComplexTypeValueSingle<double>(siteIterator->second, config.epsilonMap.name.c_str(),
                                                                   siteErrors, false);
          const auto sigma = parseComplexTypeValueSingle<double>(siteIterator->second, config.sigmaMap.name.c_str(),
                                                                 siteErrors, false);
          config.addLJParametersToSite(siteID, epsilon, sigma);

          // Check Axilrod-Teller parameter
          const auto nu =
              parseComplexTypeValueSingle<double>(siteIterator->second, config.nuMap.name.c_str(), siteErrors, false);
          config.addATParametersToSite(siteID, nu);
        }
      } else if (key == MDFlexConfig::moleculesStr) {
        // todo throw error if momentOfInertia with zero element is used (physically nonsense + breaks the quaternion
        // update)
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

#if MD_FLEXIBLE_MODE == MULTISITE
        for (auto molIterator = node[MDFlexConfig::moleculesStr].begin();
             molIterator != node[MDFlexConfig::moleculesStr].end(); ++molIterator) {
          molErrors.clear();
          molID = std::distance(node[MDFlexConfig::moleculesStr].begin(), molIterator);

          const auto molToSiteId = parseComplexTypeValueSequence<unsigned long>(
              molIterator->second, MDFlexConfig::moleculeToSiteIdStr, molErrors);
          const auto molToSitePos = parseComplexTypeValueSequence<std::array<double, 3>>(
              molIterator->second, MDFlexConfig::moleculeToSitePosStr, molErrors);
          const auto momentOfInertia = parseComplexTypeValueSequence<double, 3>(
              molIterator->second, MDFlexConfig::momentOfInertiaStr, molErrors);

          if (molToSiteId.size() != molToSitePos.size()) {
            std::stringstream ss;
            ss << "Lengths of " << MDFlexConfig::moleculeToSiteIdStr << " and " << MDFlexConfig::moleculeToSitePosStr
               << "do not match!";
            molErrors.push_back(ss.str());
          }

          config.addMolType(molID, molToSiteId, molToSitePos, momentOfInertia);
        }

        std::for_each(molErrors.begin(), molErrors.end(), pushMolError);
#else
        AutoPasLog(WARN,
                   "Multi-Site Molecule information has been provided, however md-flexible has been compiled without "
                   "Multi-Site support.");
#endif
      } else if (key == MDFlexConfig::objectsStr) {
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
              std::for_each(objectErrors.begin(), objectErrors.end(), pushObjectError);
            }
          } else if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeGaussObjectsStr) {
            generatorName = MDFlexConfig::cubeGaussObjectsStr;
            for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
              objectErrors.clear();
              objID = std::distance(objectIterator->second.begin(), it);
              const auto cubeGauss = parseCubeGaussObject(config, it->second, objectErrors);

              config.cubeGaussObjects.emplace_back(cubeGauss);
              std::for_each(objectErrors.begin(), objectErrors.end(), pushObjectError);
            }
          } else if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeUniformObjectsStr) {
            generatorName = MDFlexConfig::cubeUniformObjectsStr;
            for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
              objectErrors.clear();
              objID = std::distance(objectIterator->second.begin(), it);
              const auto cubeUniform = parseCubeUniformObject(config, it->second, objectErrors);

              config.cubeUniformObjects.emplace_back(cubeUniform);
              std::for_each(objectErrors.begin(), objectErrors.end(), pushObjectError);
            }
          } else if (objectIterator->first.as<std::string>() == MDFlexConfig::sphereObjectsStr) {
            generatorName = MDFlexConfig::sphereObjectsStr;
            for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
              objectErrors.clear();
              objID = std::distance(objectIterator->second.begin(), it);
              const auto sphere = parseSphereObject(config, it->second, objectErrors);

              config.sphereObjects.emplace_back(sphere);
              std::for_each(objectErrors.begin(), objectErrors.end(), pushObjectError);
            }
          } else if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeClosestPackedObjectsStr) {
            generatorName = MDFlexConfig::cubeClosestPackedObjectsStr;
            for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
              objectErrors.clear();
              objID = std::distance(objectIterator->second.begin(), it);
              const auto cubeClosestPacked = parseCubeClosestPacked(config, it->second, objectErrors);

              config.cubeClosestPackedObjects.emplace_back(cubeClosestPacked);
              std::for_each(objectErrors.begin(), objectErrors.end(), pushObjectError);
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
          if (config.thermostatInterval.value < 1) {
            throw std::runtime_error("thermostatInterval has to be > 0");
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

#ifndef MD_FLEXIBLE_ENABLE_ALLLBL
        if (config.loadBalancer.value == LoadBalancerOption::all) {
          errors.push_back(makeErrorMsg(mark, key,
                                        "The input file requests ALL but md-flexible was not compiled with ALL.",
                                        expected, description));
        }
#endif
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
