/**
 * @file YamlParser.cpp
 * @author N. Fottner
 * @date 15.07.2019
 */

#include "YamlParser.h"

#include <sys/stat.h>

#include <ostream>

namespace YamlParser {

bool parseYamlFile(MDFlexConfig &config) {
  YAML::Node node = YAML::LoadFile(config.yamlFilename);

  if (node[MDFlexConfig::containerOptionsStr]) {
    config.containerOptions = autopas::ContainerOption::parseOptions(
        autopas::utils::ArrayUtils::to_string(node[MDFlexConfig::containerOptionsStr], ", ", {"", ""}));
  }
  if (node[MDFlexConfig::boxMinStr]) {
    auto tmpNode = node[MDFlexConfig::boxMinStr];
    config.boxMin = {tmpNode[0].as<double>(), tmpNode[1].as<double>(), tmpNode[2].as<double>()};
  }
  if (node[MDFlexConfig::boxMaxStr]) {
    auto tmpNode = node[MDFlexConfig::boxMaxStr];
    config.boxMax = {tmpNode[0].as<double>(), tmpNode[1].as<double>(), tmpNode[2].as<double>()};
  }
  if (node[MDFlexConfig::selectorStrategyStr]) {
    auto parsedOptions =
        autopas::SelectorStrategyOption::parseOptions(node[MDFlexConfig::selectorStrategyStr].as<std::string>());
    if (parsedOptions.size() != 1) {
      throw std::runtime_error("YamlParser::parseYamlFile: Pass exactly one selector strategy option!");
    }
    config.selectorStrategy = *parsedOptions.begin();
  }
  if (node[MDFlexConfig::periodicStr]) {
    config.periodic = node[MDFlexConfig::periodicStr].as<bool>();
  }
  if (node[MDFlexConfig::cutoffStr]) {
    config.cutoff = node[MDFlexConfig::cutoffStr].as<double>();
  }
  if (node[MDFlexConfig::cellSizeFactorsStr]) {
    config.cellSizeFactors = autopas::utils::StringUtils::parseNumberSet(
        autopas::utils::ArrayUtils::to_string(node[MDFlexConfig::cellSizeFactorsStr], ", ", {"", ""}));
  }
  if (node[MDFlexConfig::dataLayoutOptionsStr]) {
    config.dataLayoutOptions = autopas::DataLayoutOption::parseOptions(
        autopas::utils::ArrayUtils::to_string(node[MDFlexConfig::dataLayoutOptionsStr], ", ", {"", ""}));
  }
  if (node[MDFlexConfig::functorOptionStr]) {
    auto strArg = node[MDFlexConfig::functorOptionStr].as<std::string>();
    if (strArg.find("avx") != std::string::npos) {
      config.functorOption = MDFlexConfig::FunctorOption::lj12_6_AVX;
    } else if (strArg.find("glob") != std::string::npos) {
      config.functorOption = MDFlexConfig::FunctorOption::lj12_6_Globals;
    } else if (strArg.find("lj") != std::string::npos || strArg.find("lennard-jones") != std::string::npos) {
      config.functorOption = MDFlexConfig::FunctorOption::lj12_6;
    }
  }
  if (node[MDFlexConfig::iterationsStr]) {
    config.iterations = node[MDFlexConfig::iterationsStr].as<unsigned long>();
  }
  if (node[MDFlexConfig::tuningPhasesStr]) {
    config.tuningPhases = node[MDFlexConfig::tuningPhasesStr].as<unsigned long>();
  }
  if (node[MDFlexConfig::dontMeasureFlopsStr]) {
    // "not" needed because of semantics
    config.dontMeasureFlops = not node[MDFlexConfig::dontMeasureFlopsStr].as<bool>();
  }
  if (node[MDFlexConfig::dontCreateEndConfigStr]) {
    // "not" needed because of semantics
    config.dontCreateEndConfig = not node[MDFlexConfig::dontCreateEndConfigStr].as<bool>();
  }
  if (node[MDFlexConfig::newton3OptionsStr]) {
    config.newton3Options = autopas::Newton3Option::parseOptions(
        autopas::utils::ArrayUtils::to_string(node[MDFlexConfig::newton3OptionsStr], ", ", {"", ""}));
  }
  if (node[MDFlexConfig::deltaTStr]) {
    config.deltaT = node[MDFlexConfig::deltaTStr].as<double>();
  }
  if (node[MDFlexConfig::traversalOptionsStr]) {
    config.traversalOptions = autopas::TraversalOption::parseOptions(
        autopas::utils::ArrayUtils::to_string(node[MDFlexConfig::traversalOptionsStr], ", ", {"", ""}));
  }
  if (node[MDFlexConfig::tuningIntervalStr]) {
    config.tuningInterval = node[MDFlexConfig::tuningIntervalStr].as<unsigned int>();
  }
  if (node[MDFlexConfig::tuningSamplesStr]) {
    config.tuningSamples = node[MDFlexConfig::tuningSamplesStr].as<unsigned int>();
  }
  if (node[MDFlexConfig::tuningMaxEvidenceStr]) {
    config.tuningMaxEvidence = node[MDFlexConfig::tuningMaxEvidenceStr].as<unsigned int>();
  }
  if (node[MDFlexConfig::relativeOptimumRangeStr]) {
    config.relativeOptimumRange = node[MDFlexConfig::relativeOptimumRangeStr].as<double>();
  }
  if (node[MDFlexConfig::maxTuningPhasesWithoutTestStr]) {
    config.maxTuningPhasesWithoutTest = node[MDFlexConfig::maxTuningPhasesWithoutTestStr].as<unsigned int>();
  }
  if (node[MDFlexConfig::tuningStrategyOptionsStr]) {
    auto parsedOptions =
        autopas::TuningStrategyOption::parseOptions(node[MDFlexConfig::tuningStrategyOptionsStr].as<std::string>());
    if (parsedOptions.size() != 1) {
      throw std::runtime_error("YamlParser::parseYamlFile: Pass exactly one tuning strategy option!");
    }
    config.tuningStrategyOption = *parsedOptions.begin();
  }
  if (node[MDFlexConfig::logLevelStr]) {
    auto strArg = node[MDFlexConfig::logLevelStr].as<std::string>();
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
    }
  }
  if (node[MDFlexConfig::checkpointfileStr]) {
    config.checkpointfile = node[MDFlexConfig::checkpointfileStr].as<std::string>();
  }
  if (node[MDFlexConfig::logFileNameStr]) {
    config.logFileName = node[MDFlexConfig::logFileNameStr].as<std::string>();
  }
  if (node[MDFlexConfig::verletRebuildFrequencyStr]) {
    config.verletRebuildFrequency = node[MDFlexConfig::verletRebuildFrequencyStr].as<unsigned int>();
  }
  if (node[MDFlexConfig::verletSkinRadiusStr]) {
    config.verletSkinRadius = node[MDFlexConfig::verletSkinRadiusStr].as<double>();
  }
  if (node[MDFlexConfig::verletClusterSizeStr]) {
    config.verletClusterSize = node[MDFlexConfig::verletClusterSizeStr].as<unsigned int>();
  }
  if (node[MDFlexConfig::vtkFileNameStr]) {
    config.vtkFileName = node[MDFlexConfig::vtkFileNameStr].as<std::string>();
  }
  if (node[MDFlexConfig::vtkWriteFrequencyStr]) {
    config.vtkWriteFrequency = node[MDFlexConfig::vtkWriteFrequencyStr].as<size_t>();
  }
  if (node[MDFlexConfig::objectsStr]) {
    // remove default objects
    config.cubeGridObjects.clear();
    config.cubeGaussObjects.clear();
    config.cubeUniformObjects.clear();
    config.sphereObjects.clear();
    config.epsilonMap.clear();
    config.sigmaMap.clear();
    config.massMap.clear();

    for (YAML::const_iterator objectIterator = node[MDFlexConfig::objectsStr].begin();
         objectIterator != node[MDFlexConfig::objectsStr].end(); ++objectIterator) {
      if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeGridObjectsStr) {
        for (YAML::const_iterator it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
          CubeGrid cubeGrid({it->second[MDFlexConfig::velocityStr][0].as<double>(),
                             it->second[MDFlexConfig::velocityStr][1].as<double>(),
                             it->second[MDFlexConfig::velocityStr][2].as<double>()},
                            it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                            it->second[MDFlexConfig::epsilonStr].as<double>(),
                            it->second[MDFlexConfig::sigmaStr].as<double>(),
                            it->second[MDFlexConfig::massStr].as<double>(),
                            {it->second[MDFlexConfig::particlesPerDimStr][0].as<unsigned long>(),
                             it->second[MDFlexConfig::particlesPerDimStr][1].as<unsigned long>(),
                             it->second[MDFlexConfig::particlesPerDimStr][2].as<unsigned long>()},
                            it->second[MDFlexConfig::particlesSpacingStr].as<double>(),
                            {it->second[MDFlexConfig::bottomLeftBackCornerStr][0].as<double>(),
                             it->second[MDFlexConfig::bottomLeftBackCornerStr][1].as<double>(),
                             it->second[MDFlexConfig::bottomLeftBackCornerStr][2].as<double>()});

          config.cubeGridObjects.emplace_back(cubeGrid);
          config.addParticleType(it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                 it->second[MDFlexConfig::epsilonStr].as<double>(),
                                 it->second[MDFlexConfig::sigmaStr].as<double>(),
                                 it->second[MDFlexConfig::massStr].as<double>());
        }
        continue;
      }
      if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeGaussObjectsStr) {
        for (YAML::const_iterator it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
          CubeGauss cubeGauss({it->second[MDFlexConfig::velocityStr][0].as<double>(),
                               it->second[MDFlexConfig::velocityStr][1].as<double>(),
                               it->second[MDFlexConfig::velocityStr][2].as<double>()},
                              it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                              it->second[MDFlexConfig::epsilonStr].as<double>(),
                              it->second[MDFlexConfig::sigmaStr].as<double>(),
                              it->second[MDFlexConfig::massStr].as<double>(),
                              it->second[MDFlexConfig::particlesPerObjectStr].as<size_t>(),
                              {it->second[MDFlexConfig::boxLengthStr][0].as<double>(),
                               it->second[MDFlexConfig::boxLengthStr][1].as<double>(),
                               it->second[MDFlexConfig::boxLengthStr][2].as<double>()},
                              {it->second[MDFlexConfig::distributionMeanStr][0].as<double>(),
                               it->second[MDFlexConfig::distributionMeanStr][1].as<double>(),
                               it->second[MDFlexConfig::distributionMeanStr][2].as<double>()},
                              {it->second[MDFlexConfig::distributionStdDevStr][0].as<double>(),
                               it->second[MDFlexConfig::distributionStdDevStr][1].as<double>(),
                               it->second[MDFlexConfig::distributionStdDevStr][2].as<double>()},
                              {it->second[MDFlexConfig::bottomLeftBackCornerStr][0].as<double>(),
                               it->second[MDFlexConfig::bottomLeftBackCornerStr][1].as<double>(),
                               it->second[MDFlexConfig::bottomLeftBackCornerStr][2].as<double>()});
          config.cubeGaussObjects.emplace_back(cubeGauss);
          config.addParticleType(it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                 it->second[MDFlexConfig::epsilonStr].as<double>(),
                                 it->second[MDFlexConfig::sigmaStr].as<double>(),
                                 it->second[MDFlexConfig::massStr].as<double>());
        }
        continue;
      }
      if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeUniformObjectsStr) {
        for (YAML::const_iterator it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
          CubeUniform cubeUniform({it->second[MDFlexConfig::velocityStr][0].as<double>(),
                                   it->second[MDFlexConfig::velocityStr][1].as<double>(),
                                   it->second[MDFlexConfig::velocityStr][2].as<double>()},
                                  it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                  it->second[MDFlexConfig::epsilonStr].as<double>(),
                                  it->second[MDFlexConfig::sigmaStr].as<double>(),
                                  it->second[MDFlexConfig::massStr].as<double>(),
                                  it->second[MDFlexConfig::particlesPerObjectStr].as<size_t>(),
                                  {it->second[MDFlexConfig::boxLengthStr][0].as<double>(),
                                   it->second[MDFlexConfig::boxLengthStr][1].as<double>(),
                                   it->second[MDFlexConfig::boxLengthStr][2].as<double>()},
                                  {it->second[MDFlexConfig::bottomLeftBackCornerStr][0].as<double>(),
                                   it->second[MDFlexConfig::bottomLeftBackCornerStr][1].as<double>(),
                                   it->second[MDFlexConfig::bottomLeftBackCornerStr][2].as<double>()});
          config.cubeUniformObjects.emplace_back(cubeUniform);
          config.addParticleType(it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                 it->second[MDFlexConfig::epsilonStr].as<double>(),
                                 it->second[MDFlexConfig::sigmaStr].as<double>(),
                                 it->second[MDFlexConfig::massStr].as<double>());
        }
        continue;
      }
      if (objectIterator->first.as<std::string>() == MDFlexConfig::sphereObjectsStr) {
        for (YAML::const_iterator it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
          Sphere sphere({it->second[MDFlexConfig::velocityStr][0].as<double>(),
                         it->second[MDFlexConfig::velocityStr][1].as<double>(),
                         it->second[MDFlexConfig::velocityStr][2].as<double>()},
                        it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                        it->second[MDFlexConfig::epsilonStr].as<double>(),
                        it->second[MDFlexConfig::sigmaStr].as<double>(), it->second[MDFlexConfig::massStr].as<double>(),
                        {it->second[MDFlexConfig::sphereCenterStr][0].as<double>(),
                         it->second[MDFlexConfig::sphereCenterStr][1].as<double>(),
                         it->second[MDFlexConfig::sphereCenterStr][2].as<double>()},
                        it->second[MDFlexConfig::sphereRadiusStr].as<int>(),
                        it->second[MDFlexConfig::particlesSpacingStr].as<double>());
          config.sphereObjects.emplace_back(sphere);
          config.addParticleType(it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                 it->second[MDFlexConfig::epsilonStr].as<double>(),
                                 it->second[MDFlexConfig::sigmaStr].as<double>(),
                                 it->second[MDFlexConfig::massStr].as<double>());
        }
        continue;
      }
    }
  }
  if (node[MDFlexConfig::thermostatStr]) {
    config.useThermostat = true;

    config.initTemperature = node[MDFlexConfig::thermostatStr][MDFlexConfig::initTemperatureStr].as<double>();
    config.thermostatInterval = node[MDFlexConfig::thermostatStr][MDFlexConfig::thermostatIntervalStr].as<size_t>();
    config.targetTemperature = node[MDFlexConfig::thermostatStr][MDFlexConfig::targetTemperatureStr].as<double>();
    config.deltaTemp = node[MDFlexConfig::thermostatStr][MDFlexConfig::deltaTempStr].as<double>();
    config.addBrownianMotion = node[MDFlexConfig::thermostatStr][MDFlexConfig::addBrownianMotionStr].as<bool>();
  }
  return true;
}

}  // namespace YamlParser