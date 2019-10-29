/**
 * @file YamlParser.cpp
 * @author N. Fottner
 * @date 15/7/19
 */

#include "YamlParser.h"
#include <sys/stat.h>

bool YamlParser::checkFileExists(const std::string &filename) {
  struct stat buffer;
  return (stat(filename.c_str(), &buffer) == 0);
}

bool YamlParser::parseYamlFile(MDFlexConfig &config) {
  if (not checkFileExists(config.yamlFilename)) {
    throw std::runtime_error("YamlParser::parseYamlFile: File " + config.yamlFilename + " not found!");
  }

  YAML::Node node = YAML::LoadFile(config.yamlFilename);

  if (node[MDFlexConfig::containerOptionsStr]) {
    config.containerOptions = autopas::utils::StringUtils::parseContainerOptions(
        node[MDFlexConfig::containerOptionsStr].as<std::string>(), false);
  }
  if (node[MDFlexConfig::selectorStrategyStr]) {
    config.selectorStrategy =
        autopas::utils::StringUtils::parseSelectorStrategy(node[MDFlexConfig::selectorStrategyStr].as<std::string>());
  }
  if (node[MDFlexConfig::periodicStr]) {
    config.periodic = node[MDFlexConfig::periodicStr].as<bool>();
  }
  if (node[MDFlexConfig::cutoffStr]) {
    config.cutoff = node[MDFlexConfig::cutoffStr].as<double>();
  }
  if (node[MDFlexConfig::cellSizeFactorsStr]) {
    config.cellSizeFactors =
        autopas::utils::StringUtils::parseNumberSet(node[MDFlexConfig::cellSizeFactorsStr].as<std::string>());
  }
  if (node[MDFlexConfig::dataLayoutOptionsStr]) {
    config.dataLayoutOptions =
        autopas::utils::StringUtils::parseDataLayout(node[MDFlexConfig::dataLayoutOptionsStr].as<std::string>());
  }
  if (node[MDFlexConfig::functorOptionStr]) {
    auto strArg = node[MDFlexConfig::functorOptionStr].as<std::string>();
    if (strArg.find("avx") != std::string::npos) {
      config.functorOption = MDFlexConfig::FunctorOption::lj12_6_AVX;
    } else if (strArg.find("lj") != std::string::npos || strArg.find("lennard-jones") != std::string::npos) {
      config.functorOption = MDFlexConfig::FunctorOption::lj12_6;
    }
  }
  if (node[MDFlexConfig::iterationsStr]) {
    config.iterations = node[MDFlexConfig::iterationsStr].as<unsigned long>();
  }
  if (node[MDFlexConfig::measureFlopsStr]) {
    config.measureFlops = node[MDFlexConfig::measureFlopsStr].as<bool>();
  }
  if (node[MDFlexConfig::newton3OptionsStr]) {
    config.newton3Options =
        autopas::utils::StringUtils::parseNewton3Options(node[MDFlexConfig::newton3OptionsStr].as<std::string>());
  }
  if (node[MDFlexConfig::deltaTStr]) {
    config.deltaT = node[MDFlexConfig::deltaTStr].as<double>();
  }
  if (node[MDFlexConfig::traversalOptionsStr]) {
    config.traversalOptions =
        autopas::utils::StringUtils::parseTraversalOptions(node[MDFlexConfig::traversalOptionsStr].as<std::string>());
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
  if (node[MDFlexConfig::tuningStrategyOptionsStr]) {
    config.tuningStrategyOption = autopas::utils::StringUtils::parseTuningStrategyOption(
        node[MDFlexConfig::tuningStrategyOptionsStr].as<std::string>());
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
    config.checkpointfile = node["checkpointFile"].as<std::string>();
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

    for (YAML::const_iterator objectIterator = node[MDFlexConfig::objectsStr].begin();
         objectIterator != node[MDFlexConfig::objectsStr].end(); ++objectIterator) {
      if (objectIterator->first.as<std::string>() == MDFlexConfig::cubeGridObjectsStr) {
        for (YAML::const_iterator it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
          CubeGrid cubeGrid({it->second[MDFlexConfig::particlesPerDimStr][0].as<unsigned long>(),
                             it->second[MDFlexConfig::particlesPerDimStr][1].as<unsigned long>(),
                             it->second[MDFlexConfig::particlesPerDimStr][2].as<unsigned long>()},
                            it->second[MDFlexConfig::particlesSpacingStr].as<double>(),
                            //                            {0.,0.,0.},
                            {it->second[MDFlexConfig::bottomLeftBackCornerStr][0].as<double>(),
                             it->second[MDFlexConfig::bottomLeftBackCornerStr][1].as<double>(),
                             it->second[MDFlexConfig::bottomLeftBackCornerStr][2].as<double>()},
                            {it->second[MDFlexConfig::velocityStr][0].as<double>(),
                             it->second[MDFlexConfig::velocityStr][1].as<double>(),
                             it->second[MDFlexConfig::velocityStr][2].as<double>()},
                            it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                            it->second[MDFlexConfig::epsilonStr].as<double>(),
                            it->second[MDFlexConfig::sigmaStr].as<double>(),
                            it->second[MDFlexConfig::massStr].as<double>());

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
          CubeGauss cubeGauss(it->second[MDFlexConfig::particlesPerObjectStr].as<size_t>(),
                              {it->second[MDFlexConfig::boxLengthStr][0].as<double>(),
                               it->second[MDFlexConfig::boxLengthStr][1].as<double>(),
                               it->second[MDFlexConfig::boxLengthStr][2].as<double>()},
                              it->second[MDFlexConfig::distributionMeanStr].as<double>(),
                              it->second[MDFlexConfig::distributionStdDevStr].as<double>(),
                              {it->second[MDFlexConfig::bottomLeftBackCornerStr][0].as<double>(),
                               it->second[MDFlexConfig::bottomLeftBackCornerStr][1].as<double>(),
                               it->second[MDFlexConfig::bottomLeftBackCornerStr][2].as<double>()},
                              {it->second[MDFlexConfig::velocityStr][0].as<double>(),
                               it->second[MDFlexConfig::velocityStr][1].as<double>(),
                               it->second[MDFlexConfig::velocityStr][2].as<double>()},
                              it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                              it->second[MDFlexConfig::epsilonStr].as<double>(),
                              it->second[MDFlexConfig::sigmaStr].as<double>(),
                              it->second[MDFlexConfig::massStr].as<double>());
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
          CubeUniform cubeUniform(it->second[MDFlexConfig::particlesPerObjectStr].as<size_t>(),
                                  {it->second[MDFlexConfig::boxLengthStr][0].as<double>(),
                                   it->second[MDFlexConfig::boxLengthStr][1].as<double>(),
                                   it->second[MDFlexConfig::boxLengthStr][2].as<double>()},
                                  {it->second[MDFlexConfig::bottomLeftBackCornerStr][0].as<double>(),
                                   it->second[MDFlexConfig::bottomLeftBackCornerStr][1].as<double>(),
                                   it->second[MDFlexConfig::bottomLeftBackCornerStr][2].as<double>()},
                                  {it->second[MDFlexConfig::velocityStr][0].as<double>(),
                                   it->second[MDFlexConfig::velocityStr][1].as<double>(),
                                   it->second[MDFlexConfig::velocityStr][2].as<double>()},
                                  it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                                  it->second[MDFlexConfig::epsilonStr].as<double>(),
                                  it->second[MDFlexConfig::sigmaStr].as<double>(),
                                  it->second[MDFlexConfig::massStr].as<double>());
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
          Sphere sphere({it->second[MDFlexConfig::sphereCenterStr][0].as<double>(),
                         it->second[MDFlexConfig::sphereCenterStr][1].as<double>(),
                         it->second[MDFlexConfig::sphereCenterStr][2].as<double>()},
                        it->second[MDFlexConfig::sphereRadiusStr].as<int>(),
                        it->second[MDFlexConfig::particlesSpacingStr].as<double>(),
                        {it->second[MDFlexConfig::velocityStr][0].as<double>(),
                         it->second[MDFlexConfig::velocityStr][1].as<double>(),
                         it->second[MDFlexConfig::velocityStr][2].as<double>()},
                        it->second[MDFlexConfig::particleTypeStr].as<unsigned long>(),
                        it->second[MDFlexConfig::epsilonStr].as<double>(),
                        it->second[MDFlexConfig::sigmaStr].as<double>(),
                        it->second[MDFlexConfig::massStr].as<double>());
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
    config.useCurrentTempForBrownianMotion =
        node[MDFlexConfig::thermostatStr][MDFlexConfig::useCurrentTempForBrownianMotionStr].as<bool>();
  }
  return true;
}
