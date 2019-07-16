//
// Created by nicola on 15.07.19.
//

#include "YamlParser.h"

YamlParser::YamlParser(string filename) {
  YAML::Node config = YAML::LoadFile(filename);

  if (config["box-length"]) {
    this->boxLength = config["box-length"].as<int>();
  }
  if (config["container"]) {
    this->containerOptions =
        autopas::utils::StringUtils::parseContainerOptions(config["container"].as<std::string>(), false);
  }
  if (config["selector-strategy"]) {
    this->selectorStrategy =
        autopas::utils::StringUtils::parseSelectorStrategy(config["selector-strategy"].as<std::string>());
  }
  if (config["cutoff"]) {
    this->cutoff = config["cutoff"].as<double>();
  }
  if (config["cell-Size-Factor"]) {
    this->cellSizeFactors = autopas::utils::StringUtils::parseNumberSet(config["cell-Size-Factor"].as<std::string>());
  }
  if (config["distribution-mean"]) {
    this->distributionMean = config["distribution-mean"].as<double>();
  }
  if (config["distriubtion-stddeviation"]) {
    this->distributionStdDev = config["distribution-stddeviataion"].as<double>();
  }
  if (config["data-layout"]) {
    this->dataLayoutOptions = autopas::utils::StringUtils::parseDataLayout(config["data-layout"].as<std::string>());
  }
  if (config["functor"]) {
    auto strArg = config["functor"].as<std::string>();
    if (strArg.find("avx") != string::npos) {
      this->functorOption = lj12_6_AVX;
    } else if (strArg.find("lj") != string::npos || strArg.find("lennard-jones") != string::npos) {
      this->functorOption = lj12_6;
    }
  }
  if (config["iterations"]) {
    this->iterations = config["iterations"].as<unsigned long>();
  }
  if (config["no-flops"]) {
    this->measureFlops = config["iterations"].as<bool>();
  }
  if (config["newton3"]) {
    this->newton3Options = autopas::utils::StringUtils::parseNewton3Options(config["newton3"].as<std::string>());
  }
  if (config["delta_t"]) {
    this->delta_t = config["delta_t"].as<double>();
  }
  if (config["epsilon"]) {
    this->epsilon = config["epsilon"].as<double>();
  }
  if (config["sigma"]) {
    this->sigma = config["sigma"].as<double>();
  }
  if (config["particle-mass"]) {
    this->mass = config["mass"].as<double>();
  }
  if (config["particles-generator"]) {
    auto strArg = config["particles-generator"].as<std::string>();
    if (strArg.find("grid") != string::npos) {
      generatorOption = GeneratorOption::grid;
    } else if (strArg.find("uni") != string::npos) {
      generatorOption = GeneratorOption::uniform;
    } else
      (strArg.find("gaus") != string::npos);
    { generatorOption = GeneratorOption::gaussian; }
  }
  if (config["particles-per-dimension"]) {
    this->particlesPerDim = config["particles-per-dimension"].as<unsigned long>();
  }
  if (config["particles-total"]) {
    this->particlesTotal = config["particles-total"].as<unsigned long>();
  }
  if (config["particle-spacing"]) {
    this->particleSpacing = config["particle-spacing"].as<double>();
  }
  if (config["traversal"]) {
    this->traversalOptions = autopas::utils::StringUtils::parseTraversalOptions(config["traversal"].as<std::string>());
  }
  if (config["tuning-interval"]) {
    this->tuningInterval = config["tuning-interval"].as<unsigned int>();
  }
  if (config["tuning-samples"]) {
    this->tuningSamples = config["tuning-samples"].as<unsigned int>();
  }
  if (config["tuning-max-evidence"]) {
    this->tuningMaxEvidence = config["tuning-max-evidence"].as<unsigned int>();
  }
  if (config["tuning-strategy"]) {
    this->tuningStrategyOption =
        autopas::utils::StringUtils::parseTuningStrategyOption(config["tuning-strategy"].as<std::string>());
  }
  if (config["log-level"]) {
    auto strArg = config["log-level"].as<std::string>();
    switch (strArg[0]) {
      case 't': {
        logLevel = spdlog::level::trace;
        break;
      }
      case 'd': {
        logLevel = spdlog::level::debug;
        break;
      }
      case 'i': {
        logLevel = spdlog::level::info;
        break;
      }
      case 'w': {
        logLevel = spdlog::level::warn;
        break;
      }
      case 'e': {
        logLevel = spdlog::level::err;
        break;
      }
      case 'c': {
        logLevel = spdlog::level::critical;
        break;
      }
      case 'o': {
        logLevel = spdlog::level::off;
        break;
      }
    }
  }
  if (config["log-file"]) {
    this->logFileName = config["log-file"].as<std::string>();
  }
  if (config["verlet-rebuild-frequency"]) {
    this->verletRebuildFrequency = config["verlet-rebuild-frequency"].as<unsigned int>();
  }
  if (config["verlet-skin-radius"]) {
    this->verletSkinRadius = config["verlet-skin-raduis"].as<double>();
  }
  if (config["vtk"]) {
    this->writeVTK = config["vtk"].as<std::string>();
  }
}