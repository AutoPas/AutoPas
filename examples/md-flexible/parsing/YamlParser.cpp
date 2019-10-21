/**
 * @file YamlParser.cpp
 * @author N. Fottner
 * @date 15/7/19
 */

#include "YamlParser.h"

bool YamlParser::parseYamlFile(MDFlexConfig &config) {
  using namespace autopas;
  YAML::Node node = YAML::LoadFile(config.yamlFilename);

  if (node["container"]) {
    config.containerOptions =
        autopas::utils::StringUtils::parseContainerOptions(node["container"].as<std::string>(), false);
  }
  if (node["selector-strategy"]) {
    config.selectorStrategy =
        autopas::utils::StringUtils::parseSelectorStrategy(node["selector-strategy"].as<std::string>());
  }
  if (node["periodic-boundaries"]) {
    config.periodic = node["periodic-boundaries"].as<bool>();
  }
  if (node["cutoff"]) {
    config.cutoff = node["cutoff"].as<double>();
  }
  if (node["cell-Size-Factor"]) {
    config.cellSizeFactors = autopas::utils::StringUtils::parseNumberSet(node["cell-Size-Factor"].as<std::string>());
  }
  if (node["data-layout"]) {
    config.dataLayoutOptions = autopas::utils::StringUtils::parseDataLayout(node["data-layout"].as<std::string>());
  }
  if (node["functor"]) {
    auto strArg = node["functor"].as<std::string>();
    if (strArg.find("avx") != std::string::npos) {
      config.functorOption = MDFlexConfig::FunctorOption::lj12_6_AVX;
    } else if (strArg.find("lj") != std::string::npos || strArg.find("lennard-jones") != std::string::npos) {
      config.functorOption = MDFlexConfig::FunctorOption::lj12_6;
    }
  }
  if (node["iterations"]) {
    config.iterations = node["iterations"].as<unsigned long>();
  }
  if (node["no-flops"]) {
    config.measureFlops = node["iterations"].as<bool>();
  }
  if (node["newton3"]) {
    config.newton3Options = autopas::utils::StringUtils::parseNewton3Options(node["newton3"].as<std::string>());
  }
  if (node["deltaT"]) {
    config.deltaT = node["deltaT"].as<double>();
  }
  if (node["traversal"]) {
    config.traversalOptions = autopas::utils::StringUtils::parseTraversalOptions(node["traversal"].as<std::string>());
  }
  if (node["tuning-interval"]) {
    config.tuningInterval = node["tuning-interval"].as<unsigned int>();
  }
  if (node["tuning-samples"]) {
    config.tuningSamples = node["tuning-samples"].as<unsigned int>();
  }
  if (node["tuning-max-evidence"]) {
    config.tuningMaxEvidence = node["tuning-max-evidence"].as<unsigned int>();
  }
  if (node["tuning-strategy"]) {
    config.tuningStrategyOption =
        autopas::utils::StringUtils::parseTuningStrategyOption(node["tuning-strategy"].as<std::string>());
  }
  if (node["log-level"]) {
    auto strArg = node["log-level"].as<std::string>();
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
  if (node["log-file"]) {
    config.logFileName = node["log-file"].as<std::string>();
  }
  if (node["verlet-rebuild-frequency"]) {
    config.verletRebuildFrequency = node["verlet-rebuild-frequency"].as<unsigned int>();
  }
  if (node["verlet-skin-radius"]) {
    config.verletSkinRadius = node["verlet-skin-radius"].as<double>();
  }
  if (node["vtk-filename"]) {
    config.VTKFileName = node["vtk-filename"].as<std::string>();
  }
  if (node["vtk-write-frequency"]) {
    config.vtkWriteFrequency = node["vtk-write-frequency"].as<size_t>();
  }
  if (node["Objects"]) {
    // remove default objects
    config.cubeGridObjects.clear();
    config.cubeGaussObjects.clear();
    config.cubeUniformObjects.clear();
    config.sphereObjects.clear();

    for (YAML::const_iterator objectIterator = node["Objects"].begin(); objectIterator != node["Objects"].end();
         ++objectIterator) {
      if (objectIterator->first.as<std::string>() == "CubeGrid") {
        for (YAML::const_iterator it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
          CubeGrid cubeGrid(
              {it->second["particles-per-Dim"][0].as<unsigned long>(),
               it->second["particles-per-Dim"][1].as<unsigned long>(),
               it->second["particles-per-Dim"][2].as<unsigned long>()},
              it->second["particleSpacing"].as<double>(),
              {it->second["bottomLeftCorner"][0].as<double>(), it->second["bottomLeftCorner"][1].as<double>(),
               it->second["bottomLeftCorner"][2].as<double>()},
              {it->second["velocity"][0].as<double>(), it->second["velocity"][1].as<double>(),
               it->second["velocity"][2].as<double>()},
              it->second["particle-type"].as<unsigned long>(), it->second["particle-epsilon"].as<double>(),
              it->second["particle-sigma"].as<double>(), it->second["particle-mass"].as<double>());
          config.cubeGridObjects.emplace_back(cubeGrid);
          config.addParticleType(it->second["particle-type"].as<unsigned long>(),
                                 it->second["particle-epsilon"].as<double>(), it->second["particle-sigma"].as<double>(),
                                 it->second["particle-mass"].as<double>());
        }
        continue;
      }
      if (objectIterator->first.as<std::string>() == "CubeGauss") {
        for (YAML::const_iterator it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
          CubeGauss cubeGauss(
              it->second["numberOfParticles"].as<size_t>(),
              {it->second["box-length"][0].as<double>(), it->second["box-length"][1].as<double>(),
               it->second["box-length"][2].as<double>()},
              it->second["distribution-mean"].as<double>(), it->second["distribution-stddev"].as<double>(),
              {it->second["bottomLeftCorner"][0].as<double>(), it->second["bottomLeftCorner"][1].as<double>(),
               it->second["bottomLeftCorner"][2].as<double>()},
              {it->second["velocity"][0].as<double>(), it->second["velocity"][1].as<double>(),
               it->second["velocity"][2].as<double>()},
              it->second["particle-type"].as<unsigned long>(), it->second["particle-epsilon"].as<double>(),
              it->second["particle-sigma"].as<double>(), it->second["particle-mass"].as<double>());
          config.cubeGaussObjects.emplace_back(cubeGauss);
          config.addParticleType(it->second["particle-type"].as<unsigned long>(),
                                 it->second["particle-epsilon"].as<double>(), it->second["particle-sigma"].as<double>(),
                                 it->second["particle-mass"].as<double>());
        }
        continue;
      }
      if (objectIterator->first.as<std::string>() == "CubeUniform") {
        for (YAML::const_iterator it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
          CubeUniform cubeUniform(
              it->second["numberOfParticles"].as<size_t>(),
              {it->second["box-length"][0].as<double>(), it->second["box-length"][1].as<double>(),
               it->second["box-length"][2].as<double>()},

              {it->second["bottomLeftCorner"][0].as<double>(), it->second["bottomLeftCorner"][1].as<double>(),
               it->second["bottomLeftCorner"][2].as<double>()},
              {it->second["velocity"][0].as<double>(), it->second["velocity"][1].as<double>(),
               it->second["velocity"][2].as<double>()},
              it->second["particle-type"].as<unsigned long>(), it->second["particle-epsilon"].as<double>(),
              it->second["particle-sigma"].as<double>(), it->second["particle-mass"].as<double>());
          config.cubeUniformObjects.emplace_back(cubeUniform);
          config.addParticleType(it->second["particle-type"].as<unsigned long>(),
                                 it->second["particle-epsilon"].as<double>(), it->second["particle-sigma"].as<double>(),
                                 it->second["particle-mass"].as<double>());
        }
        continue;
      }
      if (objectIterator->first.as<std::string>() == "Sphere") {
        for (YAML::const_iterator it = objectIterator->second.begin(); it != objectIterator->second.end(); ++it) {
          Sphere sphere({it->second["center"][0].as<double>(), it->second["center"][1].as<double>(),
                         it->second["center"][2].as<double>()},
                        it->second["radius"].as<int>(), it->second["particleSpacing"].as<double>(),
                        {it->second["velocity"][0].as<double>(), it->second["velocity"][1].as<double>(),
                         it->second["velocity"][2].as<double>()},
                        it->second["particle-type"].as<unsigned long>(), it->second["particle-epsilon"].as<double>(),
                        it->second["particle-sigma"].as<double>(), it->second["particle-mass"].as<double>());
          config.sphereObjects.emplace_back(sphere);
          config.addParticleType(it->second["particle-type"].as<unsigned long>(),
                                 it->second["particle-epsilon"].as<double>(), it->second["particle-sigma"].as<double>(),
                                 it->second["particle-mass"].as<double>());
        }
        continue;
      }
    }
  }
  if (node["Thermostat"]) {
    config.thermostat = true;

    YAML::const_iterator iterNode = node["Thermostat"].begin();
    config.initializeThermostat = iterNode->second.as<bool>();
    ++iterNode;
    config.initTemperature = iterNode->second.as<double>();
    ++iterNode;
    config.numberOfTimesteps = iterNode->second.as<size_t>();
    if (iterNode != node["Thermostat"].end()) {  // if target value is specified
      config.thermoTarget = true;
      ++iterNode;
      config.targetTemperature = iterNode->second["targetTemperature"].as<double>();
      config.deltaTemp = iterNode->second["deltaTemp"].as<double>();
    }
  }
  return true;
}
