/**
 * @file YamlParser.cpp
 * @author N. Fottner
 * @date 15.07.2019
 */
#include "YamlParser.h"

#include <string>

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

        auto tmpNode = node[config.boxMax.name];
        config.boxMax.value = {tmpNode[0].as<double>(), tmpNode[1].as<double>(), tmpNode[2].as<double>()};
      } else if (key == config.subdivideDimension.name) {
        expected = "YAML-sequence of three ints in [0, 1].";
        description = config.subdivideDimension.description;

        auto tmpNode = node[config.subdivideDimension.name];
        config.subdivideDimension.value = {tmpNode[0].as<bool>(), tmpNode[1].as<bool>(), tmpNode[2].as<bool>()};
      } else if (key == config.loadBalancingInterval.name) {
        expected = "Unsigned Integer";
        description = config.loadBalancingInterval.description;

        int tmp = node[config.loadBalancingInterval.name].as<int>();
        if (tmp < 0) {
          throw YamlParserException("Load balancing interval must be a positive integer.");
        }

        config.loadBalancingInterval.value = tmp;
      } else if (key == config.selectorStrategy.name) {
        expected = "Exactly one selector strategy out of the possible values.";
        description = config.selectorStrategy.description;

        std::set<autopas::options::SelectorStrategyOption> parsedOptions;
        if (node[config.selectorStrategy.name].IsSequence()) {
          if (node[config.selectorStrategy.name].size() != 1) {
            throw YamlParserException("Pass Exactly one selector strategy.");
          }
          parsedOptions = autopas::SelectorStrategyOption::parseOptions(
              autopas::utils::ArrayUtils::to_string(node[config.selectorStrategy.name], "", {"", ""}));
        } else {
          parsedOptions =
              autopas::SelectorStrategyOption::parseOptions(node[config.selectorStrategy.name].as<std::string>());
        }
        config.selectorStrategy.value = *parsedOptions.begin();

      } else if (key == config.boundaryOption.name) {
        expected = "YAML-sequence of three possible values.";
        description = config.boundaryOption.description;

        auto tmpNode = node[config.boundaryOption.name];
        config.boundaryOption.value = {options::BoundaryTypeOption::parseOptionExact(tmpNode[0].as<std::string>()),
                                       options::BoundaryTypeOption::parseOptionExact(tmpNode[1].as<std::string>()),
                                       options::BoundaryTypeOption::parseOptionExact(tmpNode[2].as<std::string>())};
      } else if (key == config.cutoff.name) {
        expected = "Positive floating point value.";
        description = config.cutoff.description;

        config.cutoff.value = node[config.cutoff.name].as<double>();
      } else if (key == config.cellSizeFactors.name) {
        expected = "YAML-sequence of floats.";
        description = config.cellSizeFactors.description;

        config.cellSizeFactors.value = autopas::utils::StringUtils::parseNumberSet(
            autopas::utils::ArrayUtils::to_string(node[config.cellSizeFactors.name], ", ", {"", ""}));

        if (config.cellSizeFactors.value->isEmpty()) {
          throw YamlParserException("Parsed cell-size-factor-list is empty.");
        }
      } else if (key == config.dataLayoutOptions.name) {
        expected = "YAML-sequence of possible values.";
        description = config.dataLayoutOptions.description;

        config.dataLayoutOptions.value = autopas::DataLayoutOption::parseOptions(
            autopas::utils::ArrayUtils::to_string(node[config.dataLayoutOptions.name], ", ", {"", ""}));
        if (config.dataLayoutOptions.value.empty()) {
          throw YamlParserException("Parsed data-layouts-list is empty.");
        }
      } else if (key == config.functorOption.name) {
        expected = "One of the possible values.";
        description = config.functorOption.description;

        auto strArg = node[config.functorOption.name].as<std::string>();
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
        expected = "Unsigned Integer";
        description = config.iterations.description;

        long tmp = node[config.iterations.name].as<long>();
        if (tmp < 1) {
          throw YamlParserException("The number of iterations has to be a positive integer.");
        }
        config.iterations.value = tmp;

      } else if (key == config.tuningPhases.name) {
        expected = "Unsigned Integer";
        description = config.tuningPhases.description;

        long tmp = node[config.tuningPhases.name].as<long>();
        if (tmp < 0) {
          throw YamlParserException("The number of tuning phases has to be a positive integer.");
        }

        config.tuningPhases.value = tmp;
      } else if (key == config.dontMeasureFlops.name) {
        expected = "Boolean Value";
        description = config.dontMeasureFlops.description;

        // "not" needed because of semantics
        config.dontMeasureFlops.value = not node[config.dontMeasureFlops.name].as<bool>();
      } else if (key == config.dontCreateEndConfig.name) {
        expected = "Boolean Value";
        description = config.dontCreateEndConfig.description;

        // "not" needed because of semantics
        config.dontCreateEndConfig.value = not node[config.dontCreateEndConfig.name].as<bool>();
      } else if (key == config.dontShowProgressBar.name) {
        expected = "Boolean Value";
        description = config.dontShowProgressBar.description;

        config.dontShowProgressBar.value = node[config.dontShowProgressBar.name].as<bool>();
      } else if (key == config.newton3Options.name) {
        expected = "YAML-sequence of possible values.";
        description = config.newton3Options.description;

        config.newton3Options.value = autopas::Newton3Option::parseOptions(
            autopas::utils::ArrayUtils::to_string(node[config.newton3Options.name], ", ", {"", ""}));
        if (config.newton3Options.value.empty()) {
          throw YamlParserException("Unknown Newton3 option!");
        }
      } else if (key == config.deltaT.name) {
        expected = "Positive floating point value.";
        description = config.deltaT.description;

        config.deltaT.value = node[config.deltaT.name].as<double>();
      } else if (key == config.traversalOptions.name) {
        expected = "YAML-sequence of possible values.";
        description = config.traversalOptions.description;

        config.traversalOptions.value = autopas::TraversalOption::parseOptions(
            autopas::utils::ArrayUtils::to_string(node[config.traversalOptions.name], ", ", {"", ""}));

        if (config.traversalOptions.value.empty()) {
          throw YamlParserException("Parsed traversal-list is empty. Maybe you used an unknown option.");
        }

      } else if (key == config.loadEstimatorOptions.name) {
        expected = "YAML-sequence of possible values.";
        description = config.loadEstimatorOptions.description;

        config.loadEstimatorOptions.value = autopas::LoadEstimatorOption::parseOptions(
            autopas::utils::ArrayUtils::to_string(node[config.loadEstimatorOptions.name], ", ", {"", ""}));

        if (config.loadEstimatorOptions.value.empty()) {
          throw YamlParserException("Parsed load-estimator-list is empty. Maybe you used an unknown option.");
        }

      } else if (key == config.tuningInterval.name) {
        expected = "Unsigned Integer";
        description = config.tuningInterval.description;

        int tmp = node[config.tuningInterval.name].as<int>();
        if (tmp < 1) {
          throw YamlParserException("Tuning interval has to be a positive integer!");
        }

        config.tuningInterval.value = tmp;
      } else if (key == config.tuningSamples.name) {
        expected = "Unsigned Integer >= 1";
        description = config.tuningSamples.description;

        int tmp = node[config.tuningSamples.name].as<int>();
        if (tmp < 1) {
          throw YamlParserException("Tuning samples has to be a positive integer!");
        }

        config.tuningSamples.value = tmp;
      } else if (key == config.tuningMaxEvidence.name) {
        expected = "Unsigned Integer";
        description = config.tuningMaxEvidence.description;

        int tmp = node[config.tuningMaxEvidence.name].as<int>();
        if (tmp < 1) {
          throw YamlParserException("Tuning max evidence has to be a positive integer!");
        }

        config.tuningMaxEvidence.value = tmp;
      } else if (key == config.relativeOptimumRange.name) {
        expected = "Floating point value >= 1";
        description = config.relativeOptimumRange.description;

        double tmp = node[config.relativeOptimumRange.name].as<double>();
        if (tmp < 1.0) {
          throw YamlParserException("Relative optimum range has to be greater or equal one!");
        }

        config.relativeOptimumRange.value = tmp;
      } else if (key == config.maxTuningPhasesWithoutTest.name) {
        expected = "Unsigned Integer";
        description = config.maxTuningPhasesWithoutTest.description;

        int tmp = node[config.maxTuningPhasesWithoutTest.name].as<int>();
        if (tmp < 1) {
          throw YamlParserException("Max tuning phases without test has to be positive!");
        }

        config.maxTuningPhasesWithoutTest.value = tmp;
      } else if (key == config.relativeBlacklistRange.name) {
        expected = "Floating point value >= 1 or 0";
        description = config.relativeBlacklistRange.description;

        double tmp = node[config.relativeBlacklistRange.name].as<double>();
        if (tmp < 1.0 and tmp != 0.0) {
          throw YamlParserException(
              "Relative range for blacklist range has to be greater or equal one or has to be zero!");
        }

        config.relativeBlacklistRange.value = tmp;
      } else if (key == config.evidenceFirstPrediction.name) {
        expected = "Unsigned Integer >= 2";
        description = config.evidenceFirstPrediction.description;

        int tmp = node[config.evidenceFirstPrediction.name].as<int>();
        if (tmp < 2) {
          throw YamlParserException("The number of evidence for the first prediction has to be at least two!");
        }

        config.evidenceFirstPrediction.value = tmp;
      } else if (key == config.extrapolationMethodOption.name) {
        expected = "Exactly one extrapolation method out of the possible values.";
        description = config.extrapolationMethodOption.description;

        std::set<autopas::options::ExtrapolationMethodOption> parsedOptions;
        if (node[config.extrapolationMethodOption.name].IsSequence()) {
          if (node[config.extrapolationMethodOption.name].size() != 1) {
            throw YamlParserException("Pass exactly one extrapolation method!");
          }
          parsedOptions = autopas::ExtrapolationMethodOption::parseOptions(
              autopas::utils::ArrayUtils::to_string(node[config.extrapolationMethodOption.name], "", {"", ""}));
        } else {
          parsedOptions = autopas::ExtrapolationMethodOption::parseOptions(
              node[config.extrapolationMethodOption.name].as<std::string>());
        }
        config.extrapolationMethodOption.value = *parsedOptions.begin();

      } else if (key == config.tuningStrategyOption.name) {
        expected = "Exactly one tuning strategy option out of the possible values.";
        description = config.tuningStrategyOption.description;

        std::set<autopas::options::TuningStrategyOption> parsedOptions;
        if (node[config.tuningStrategyOption.name].IsSequence()) {
          if (node[config.tuningStrategyOption.name].size() != 1) {
            throw YamlParserException("Pass Exactly one tuning strategy!");
          }
          parsedOptions = autopas::TuningStrategyOption::parseOptions(
              autopas::utils::ArrayUtils::to_string(node[config.tuningStrategyOption.name], "", {"", ""}));
        } else {
          parsedOptions =
              autopas::TuningStrategyOption::parseOptions(node[config.tuningStrategyOption.name].as<std::string>());
        }
        config.tuningStrategyOption.value = *parsedOptions.begin();
      } else if (key == config.mpiStrategyOption.name) {
        expected = "Exactly one MPI strategy option out of the possible values.";
        description = config.mpiStrategyOption.description;

        std::set<autopas::options::MPIStrategyOption> parsedOptions;
        if (node[config.mpiStrategyOption.name].IsSequence()) {
          if (node[config.mpiStrategyOption.name].size() != 1) {
            throw YamlParserException("Pass exactly one MPI strategy!");
          }
          parsedOptions = autopas::MPIStrategyOption::parseOptions(
              autopas::utils::ArrayUtils::to_string(node[config.mpiStrategyOption.name], "", {"", ""}));
        } else {
          parsedOptions =
              autopas::MPIStrategyOption::parseOptions(node[config.mpiStrategyOption.name].as<std::string>());
        }
        config.mpiStrategyOption.value = *parsedOptions.begin();
      } else if (key == config.MPITuningMaxDifferenceForBucket.name) {
        expected = "Floating-point Value";
        description = config.MPITuningMaxDifferenceForBucket.description;

        config.MPITuningMaxDifferenceForBucket.value = node[config.MPITuningMaxDifferenceForBucket.name].as<double>();
      } else if (key == config.MPITuningWeightForMaxDensity.name) {
        expected = "Floating-point Value";
        description = config.MPITuningWeightForMaxDensity.description;

        config.MPITuningWeightForMaxDensity.value = node[config.MPITuningWeightForMaxDensity.name].as<double>();
      } else if (key == config.acquisitionFunctionOption.name) {
        expected = "Exactly one acquisition function option out of the possible values.";
        description = config.acquisitionFunctionOption.description;

        std::set<autopas::options::AcquisitionFunctionOption> parsedOptions;
        if (node[config.acquisitionFunctionOption.name].IsSequence()) {
          if (node[config.acquisitionFunctionOption.name].size() != 1) {
            throw YamlParserException("Pass Exactly one acquisition function option!");
          }
          parsedOptions = autopas::AcquisitionFunctionOption::parseOptions(
              autopas::utils::ArrayUtils::to_string(node[config.acquisitionFunctionOption.name], "", {"", ""}));
        } else {
          parsedOptions = autopas::AcquisitionFunctionOption::parseOptions(
              node[config.acquisitionFunctionOption.name].as<std::string>());
        }
        config.acquisitionFunctionOption.value = *parsedOptions.begin();
      } else if (key == config.logLevel.name) {
        expected = "Log level out of the possible values.";
        description = config.logLevel.description;

        auto strArg = node[config.logLevel.name].as<std::string>();
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

        config.checkpointfile.value = node[config.checkpointfile.name].as<std::string>();
        if (config.checkpointfile.value.empty()) {
          throw YamlParserException("Parsed checkpoint filename is empty!");
        }
      } else if (key == config.logFileName.name) {
        expected = "String";
        description = config.logFileName.description;

        config.logFileName.value = node[config.logFileName.name].as<std::string>();
        if (config.logFileName.value.empty()) {
          throw YamlParserException("Parsed log filename is empty!");
        }
      } else if (key == config.verletRebuildFrequency.name) {
        expected = "Unsigned Integer";
        description = config.verletRebuildFrequency.description;

        int tmp = node[config.verletRebuildFrequency.name].as<int>();
        if (tmp < 0) {
          throw YamlParserException("Verlet rebuild frequency has to be a positive integer!");
        }

        config.verletRebuildFrequency.value = tmp;
      } else if (key == config.verletSkinRadiusPerTimestep.name) {
        expected = "Positive floating-point value.";
        description = config.verletSkinRadiusPerTimestep.description;

        config.verletSkinRadiusPerTimestep.value = node[config.verletSkinRadiusPerTimestep.name].as<double>();
      } else if (key == config.fastParticlesThrow.name) {
        expected = "Boolean Value";
        description = config.fastParticlesThrow.description;

        config.fastParticlesThrow.value = node[config.fastParticlesThrow.name].as<bool>();
      } else if (key == config.verletClusterSize.name) {
        expected = "Unsigned Integer";
        description = config.verletClusterSize.description;

        int tmp = node[config.verletClusterSize.name].as<int>();
        if (tmp < 0) {
          throw YamlParserException("Verlet cluster size has to be a positive integer!");
        }

        config.verletClusterSize.value = tmp;
      } else if (key == config.vtkFileName.name) {
        expected = "String";
        description = config.vtkFileName.description;

        config.vtkFileName.value = node[config.vtkFileName.name].as<std::string>();
        if (config.vtkFileName.value.empty()) {
          throw YamlParserException("Parsed VTK filename is empty!");
        }
      } else if (key == config.vtkWriteFrequency.name) {
        expected = "Unsigned Integer >= 1";
        description = config.vtkWriteFrequency.description;

        int tmp = node[config.vtkWriteFrequency.name].as<int>();
        if (tmp < 1) {
          throw YamlParserException("VTK write frequency has to be a positive integer >= 1!");
        }

        config.vtkWriteFrequency.value = (size_t)tmp;
      } else if (key == config.globalForce.name) {
        expected = "YAML-sequence of three floats. Example: [0, 0, -9.81].";
        description = config.globalForce.description;

        config.globalForce.value = {node[config.globalForce.name][0].as<double>(),
                                    node[config.globalForce.name][1].as<double>(),
                                    node[config.globalForce.name][2].as<double>()};
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

        for (auto objectIterator = node[MDFlexConfig::objectsStr].begin();
             objectIterator != node[MDFlexConfig::objectsStr].end(); ++objectIterator) {
          if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeGridObjectsStr) {
            try {
              for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
                CubeGrid cubeGrid({it->second[MDFlexConfig::velocityStr][0].as<double>(),
                                   it->second[MDFlexConfig::velocityStr][1].as<double>(),
                                   it->second[MDFlexConfig::velocityStr][2].as<double>()},
                                  it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                  it->second[config.epsilonMap.name].as<double>(),
                                  it->second[config.sigmaMap.name].as<double>(),
                                  it->second[config.massMap.name].as<double>(),
                                  {it->second[config.particlesPerDim.name][0].as<unsigned long>(),
                                   it->second[config.particlesPerDim.name][1].as<unsigned long>(),
                                   it->second[config.particlesPerDim.name][2].as<unsigned long>()},
                                  it->second[config.particleSpacing.name].as<double>(),
                                  {it->second[MDFlexConfig::bottomLeftBackCornerStr][0].as<double>(),
                                   it->second[MDFlexConfig::bottomLeftBackCornerStr][1].as<double>(),
                                   it->second[MDFlexConfig::bottomLeftBackCornerStr][2].as<double>()});

                config.cubeGridObjects.emplace_back(cubeGrid);
                config.addParticleType(it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                       it->second[config.epsilonMap.name].as<double>(),
                                       it->second[config.sigmaMap.name].as<double>(),
                                       it->second[config.massMap.name].as<double>());
              }
            } catch (const std::exception &e) {
              // std::cerr << e.what() << std::endl;
              std::cerr << "YamlParser: Error parsing one of " << MDFlexConfig::cubeGridObjectsStr << " objects. "
                        << "See AllOptions.yaml for examples." << std::endl;
              return false;
            }
          } else if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeGaussObjectsStr) {
            try {
              for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
                CubeGauss cubeGauss({it->second[MDFlexConfig::velocityStr][0].as<double>(),
                                     it->second[MDFlexConfig::velocityStr][1].as<double>(),
                                     it->second[MDFlexConfig::velocityStr][2].as<double>()},
                                    it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                    it->second[config.epsilonMap.name].as<double>(),
                                    it->second[config.sigmaMap.name].as<double>(),
                                    it->second[config.massMap.name].as<double>(),
                                    it->second[MDFlexConfig::particlesPerObjectStr].as<size_t>(),
                                    {it->second[config.boxLength.name][0].as<double>(),
                                     it->second[config.boxLength.name][1].as<double>(),
                                     it->second[config.boxLength.name][2].as<double>()},
                                    {it->second[config.distributionMean.name][0].as<double>(),
                                     it->second[config.distributionMean.name][1].as<double>(),
                                     it->second[config.distributionMean.name][2].as<double>()},
                                    {it->second[config.distributionStdDev.name][0].as<double>(),
                                     it->second[config.distributionStdDev.name][1].as<double>(),
                                     it->second[config.distributionStdDev.name][2].as<double>()},
                                    {it->second[MDFlexConfig::bottomLeftBackCornerStr][0].as<double>(),
                                     it->second[MDFlexConfig::bottomLeftBackCornerStr][1].as<double>(),
                                     it->second[MDFlexConfig::bottomLeftBackCornerStr][2].as<double>()});
                config.cubeGaussObjects.emplace_back(cubeGauss);
                config.addParticleType(it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                       it->second[config.epsilonMap.name].as<double>(),
                                       it->second[config.sigmaMap.name].as<double>(),
                                       it->second[config.massMap.name].as<double>());
              }
            } catch (const std::exception &e) {
              std::cerr << "YamlParser: Error parsing one of " << MDFlexConfig::cubeGaussObjectsStr << " objects. "
                        << "See AllOptions.yaml for examples." << std::endl;
              return false;
            }
          } else if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeUniformObjectsStr) {
            try {
              for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
                CubeUniform cubeUniform({it->second[MDFlexConfig::velocityStr][0].as<double>(),
                                         it->second[MDFlexConfig::velocityStr][1].as<double>(),
                                         it->second[MDFlexConfig::velocityStr][2].as<double>()},
                                        it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                        it->second[config.epsilonMap.name].as<double>(),
                                        it->second[config.sigmaMap.name].as<double>(),
                                        it->second[config.massMap.name].as<double>(),
                                        it->second[MDFlexConfig::particlesPerObjectStr].as<size_t>(),
                                        {it->second[config.boxLength.name][0].as<double>(),
                                         it->second[config.boxLength.name][1].as<double>(),
                                         it->second[config.boxLength.name][2].as<double>()},
                                        {it->second[MDFlexConfig::bottomLeftBackCornerStr][0].as<double>(),
                                         it->second[MDFlexConfig::bottomLeftBackCornerStr][1].as<double>(),
                                         it->second[MDFlexConfig::bottomLeftBackCornerStr][2].as<double>()});
                config.cubeUniformObjects.emplace_back(cubeUniform);
                config.addParticleType(it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                       it->second[config.epsilonMap.name].as<double>(),
                                       it->second[config.sigmaMap.name].as<double>(),
                                       it->second[config.massMap.name].as<double>());
              }
            } catch (const std::exception &e) {
              std::cerr << "YamlParser: Error parsing one of " << MDFlexConfig::cubeUniformObjectsStr << " objects. "
                        << "See AllOptions.yaml for examples." << std::endl;
              return false;
            }
          } else if (objectIterator->first.as<std::string>() == MDFlexConfig::sphereObjectsStr) {
            try {
              for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
                Sphere sphere({it->second[MDFlexConfig::velocityStr][0].as<double>(),
                               it->second[MDFlexConfig::velocityStr][1].as<double>(),
                               it->second[MDFlexConfig::velocityStr][2].as<double>()},
                              it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                              it->second[config.epsilonMap.name].as<double>(),
                              it->second[config.sigmaMap.name].as<double>(),
                              it->second[config.massMap.name].as<double>(),
                              {it->second[MDFlexConfig::sphereCenterStr][0].as<double>(),
                               it->second[MDFlexConfig::sphereCenterStr][1].as<double>(),
                               it->second[MDFlexConfig::sphereCenterStr][2].as<double>()},
                              it->second[MDFlexConfig::sphereRadiusStr].as<int>(),
                              it->second[config.particleSpacing.name].as<double>());
                config.sphereObjects.emplace_back(sphere);
                config.addParticleType(it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                       it->second[config.epsilonMap.name].as<double>(),
                                       it->second[config.sigmaMap.name].as<double>(),
                                       it->second[config.massMap.name].as<double>());
              }
            } catch (const std::exception &e) {
              std::cerr << "YamlParser: Error parsing one of " << MDFlexConfig::sphereObjectsStr << " objects. "
                        << "See AllOptions.yaml for examples." << std::endl;
              return false;
            }
          } else if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeClosestPackedObjectsStr) {
            try {
              for (auto it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
                CubeClosestPacked cubeClosestPacked(
                    {it->second[MDFlexConfig::velocityStr][0].as<double>(),
                     it->second[MDFlexConfig::velocityStr][1].as<double>(),
                     it->second[MDFlexConfig::velocityStr][2].as<double>()},
                    it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                    it->second[config.epsilonMap.name].as<double>(), it->second[config.sigmaMap.name].as<double>(),
                    it->second[config.massMap.name].as<double>(), it->second[config.particleSpacing.name].as<double>(),
                    {it->second[config.boxLength.name][0].as<double>(),
                     it->second[config.boxLength.name][1].as<double>(),
                     it->second[config.boxLength.name][2].as<double>()},
                    {it->second[MDFlexConfig::bottomLeftBackCornerStr][0].as<double>(),
                     it->second[MDFlexConfig::bottomLeftBackCornerStr][1].as<double>(),
                     it->second[MDFlexConfig::bottomLeftBackCornerStr][2].as<double>()});
                config.cubeClosestPackedObjects.emplace_back(cubeClosestPacked);
                config.addParticleType(it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                       it->second[config.epsilonMap.name].as<double>(),
                                       it->second[config.sigmaMap.name].as<double>(),
                                       it->second[config.massMap.name].as<double>());
              }
            } catch (const std::exception &e) {
              std::cerr << "YamlParser: Error parsing one of " << MDFlexConfig::cubeClosestPackedObjectsStr
                        << " objects. "
                        << "See AllOptions.yaml for examples." << std::endl;
              return false;
            }
          } else {
            std::cerr << "YamlParser: Unrecognized generator \"" << objectIterator->first.as<std::string>()
                      << "\" used." << std::endl;
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
        config.initTemperature.value = node[config.useThermostat.name][config.initTemperature.name].as<double>();

        m = node[key][config.thermostatInterval.name].Mark();
        expected = "Unsigned Integer";
        description = config.thermostatInterval.description;
        config.thermostatInterval.value = node[config.useThermostat.name][config.thermostatInterval.name].as<size_t>();

        m = node[key][config.targetTemperature.name].Mark();
        expected = "Floating-Point Value";
        description = config.targetTemperature.description;
        config.targetTemperature.value = node[config.useThermostat.name][config.targetTemperature.name].as<double>();

        m = node[key][config.deltaTemp.name].Mark();
        expected = "Floating-Point Value";
        description = config.deltaTemp.description;
        config.deltaTemp.value = node[config.useThermostat.name][config.deltaTemp.name].as<double>();

        m = node[key][config.addBrownianMotion.name].Mark();
        expected = "Boolean Value";
        description = config.addBrownianMotion.description;
        config.addBrownianMotion.value = node[config.useThermostat.name][config.addBrownianMotion.name].as<bool>();

      } else if (key == config.loadBalancer.name) {
        expected = "YAML-sequence of possible values.";
        description = config.loadBalancer.description;

        std::set<LoadBalancerOption> parsedOptions;
        if (node[config.loadBalancer.name].IsSequence()) {
          if (node[config.loadBalancer.name].size() != 1) {
            throw YamlParserException("Pass Exactly one load balancer option!");
          }
          parsedOptions = LoadBalancerOption::parseOptions(
              autopas::utils::ArrayUtils::to_string(node[config.loadBalancer.name], "", {"", ""}));
        } else {
          parsedOptions = LoadBalancerOption::parseOptions(node[config.loadBalancer.name].as<std::string>());
        }
        config.loadBalancer.value = *parsedOptions.begin();
      } else {
        std::cerr << "Unrecognized option in input YAML: " + key << std::endl;
        // return false;
      }
    } catch (const YAML::Exception &e) {
      // We do not use e.mark, as this don't provides the correct line number in some cases. Use the mark from above;
      std::cerr << "Error while parsing the YAML-file in line " << (m.line + 1) << " at column " << m.column
                << std::endl
                << "Expected: " << expected << std::endl
                << "Parameter description: " << description << std::endl;
      return false;
    } catch (const YamlParserException &e) {
      std::cerr << "Incorrect input-parameter for key " << key << ": " << e.what() << std::endl
                << "Message: " << e.what() << std::endl
                << "Expected: " << expected << std::endl
                << "Parameter description: " << description << std::endl;
      return false;
    } catch (const std::exception &e) {
      std::cerr << "Error while parsing the YAML-file in line " << (m.line + 1) << " at column " << m.column
                << std::endl
                << "Message: " << e.what() << std::endl
                << "Expected: " << expected << std::endl
                << "Parameter description: " << description << std::endl;
      return false;
    }
  }

  return true;
}
