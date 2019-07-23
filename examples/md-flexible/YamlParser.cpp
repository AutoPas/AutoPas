//
// Created by nicola on 15.07.19.
//

#include "YamlParser.h"

void YamlParser::parseInput(string &filename) {
  YAML::Node config = YAML::LoadFile(filename);

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

  if (config["Objects"]) {
    for (YAML::const_iterator it = config["Objects"].begin(); it != config["Objects"].end(); ++it) {
      if (it->first.as<std::string>() == "CubeGrid") {
        CubeGrid C({it->second["particles-per-Dim"][0].as<unsigned long>(),
                    it->second["particles-per-Dim"][1].as<unsigned long>(),
                    it->second["particles-per-Dim"][2].as<unsigned long>()},
                   it->second["particleSpacing"].as<double>(),
                   {it->second["velocity"][0].as<double>(), it->second["velocity"][1].as<double>(),
                    it->second["velocity"][2].as<double>()});
        CubeGridObjects.emplace_back(C);
        continue;
      }
      if (it->first.as<std::string>() == "CubeGauss") {
        CubeGauss C({it->second["box-length"][0].as<double>(), it->second["box-length"][1].as<double>(),
                     it->second["box-length"][2].as<double>()},
                    it->second["numberOfParticles"].as<size_t>(), it->second["distribution-mean"].as<double>(),
                    it->second["distribution-stddeviation"].as<double>(),
                    {it->second["velocity"][0].as<double>(), it->second["velocity"][1].as<double>(),
                     it->second["velocity"][2].as<double>()});
        CubeGaussObjects.emplace_back(C);
        continue;
      }
      if (it->first.as<std::string>() == "CubeUniform") {
        CubeUniform C({it->second["box-length"][0].as<double>(), it->second["box-length"][1].as<double>(),
                       it->second["box-length"][2].as<double>()},
                      it->second["numberOfParticles"].as<size_t>(),
                      {it->second["velocity"][0].as<double>(), it->second["velocity"][1].as<double>(),
                       it->second["velocity"][2].as<double>()});
        CubeUniformObjects.emplace_back(C);
        continue;
      }
      if (it->first.as<std::string>() == "Sphere") {
        Sphere S({it->second["center"][0].as<double>(), it->second["center"][1].as<double>(),
                  it->second["center"][2].as<double>()},
                 it->second["radius"].as<int>(), it->second["particleSpacing"].as<double>(),
                 it->second["firstId"].as<unsigned long>(),
                 {it->second["velocity"][0].as<double>(), it->second["velocity"][1].as<double>(),
                  it->second["velocity"][2].as<double>()});
        SphereObjects.emplace_back(S);
        continue;
      }
    }
  }
}

template <class T>
std::string iterableToString(T arr) {
  std::ostringstream ss;
  for (auto a : arr) {
    ss << autopas::utils::StringUtils::to_string(a) << ", ";
  }
  // deletes last comma
  ss << "\b\b";
  return ss.str();
}

void YamlParser::printConfig() {
  constexpr size_t valueOffset = 32;
  cout << setw(valueOffset) << left << "Container"
       << ":  " << iterableToString(containerOptions) << endl;

  // if verlet lists are in the container options print verlet config data
  if (iterableToString(containerOptions).find("erlet") != std::string::npos) {
    cout << setw(valueOffset) << left << "Verlet rebuild frequency"
         << ":  " << verletRebuildFrequency << endl;

    cout << setw(valueOffset) << left << "Verlet skin radius"
         << ":  " << verletSkinRadius << endl;
  }

  if (containerOptions.size() > 1 or traversalOptions.size() > 1 or dataLayoutOptions.size() > 1) {
    cout << setw(valueOffset) << left << "Selector Strategy"
         << ":  " << autopas::utils::StringUtils::to_string(selectorStrategy) << endl;
  }

  cout << setw(valueOffset) << left << "Data Layout"
       << ":  " << iterableToString(dataLayoutOptions) << endl;

  cout << setw(valueOffset) << left << "Functor"
       << ":  ";
  switch (functorOption) {
    case FunctorOption::lj12_6: {
      cout << "Lennard-Jones (12-6)" << endl;
      break;
    }
    case FunctorOption::lj12_6_AVX: {
      cout << "Lennard-Jones (12-6) AVX intrinsics" << endl;
      break;
    }
  }

  cout << setw(valueOffset) << left << "Newton3"
       << ":  " << iterableToString(newton3Options) << endl;

  cout << setw(valueOffset) << left << "Cutoff radius"
       << ":  " << cutoff << endl;

  cout << setw(valueOffset) << left << "Cell size factor"
       << ":  " << static_cast<std::string>(*cellSizeFactors) << endl;

  cout << setw(valueOffset) << left << "Object Generation:" << endl;
  int i=1;
  for (auto c : CubeGridObjects) {
      cout << "-Cube Grid Nr " << i <<  ":  " << endl;
    c.printConfig();
    i++;
  }
  i=1;
  for (auto c : CubeGaussObjects) {
      cout << "-Cube Gauss Nr" << i <<  ":  " << endl;
    c.printConfig();
    i++;
  }
  i=1;
  for (auto c : CubeUniformObjects) {
      cout << "-Cube Uniform Nr " << i <<  ":  " << endl;
    c.printConfig();
    i++;
  }
  i=1;
  for (auto c : SphereObjects) {
    cout << "-Sphere Nr " << i <<  ":  " << endl;
    c.printConfig();
    i++;
  }

  cout << setw(valueOffset) << left << "Allowed traversals"
       << ":  " << iterableToString(traversalOptions) << endl;
  cout << setw(valueOffset) << left << "Particles Mass"
       << ":  " << mass << endl;
  cout << setw(valueOffset) << left
       << "Particles Epsilon"  //@todo verändern wenn verschieden ParticleType in der Simulation sind
       << ":  " << epsilon << endl;
  cout << setw(valueOffset) << left
       << "Particles Sigma"  //@todo verändern wenn verschieden ParticleType in der Simulation sind
       << ":  " << sigma << endl;
  cout << setw(valueOffset) << left << "delta_t"
       << ":  " << delta_t << endl;
  cout << setw(valueOffset) << left << "Iterations"  // iterations * delta_t = time_end;
       << ":  " << iterations << endl;
  cout << setw(valueOffset) << left << "Tuning Strategy"
       << ":  " << autopas::utils::StringUtils::to_string(tuningStrategyOption) << endl;
  cout << setw(valueOffset) << left << "Tuning Interval"
       << ":  " << tuningInterval << endl;
  cout << setw(valueOffset) << left << "Tuning Samples"
       << ":  " << tuningSamples << endl;
  cout << setw(valueOffset) << left << "Tuning Max evidence"
       << ":  " << tuningMaxEvidence << endl;
}

size_t YamlParser::particlesTotal() {
  size_t particlesTotal = 0;
  for (auto e : CubeGridObjects) {
    particlesTotal += e.getParticlesTotal();
  }
  for (auto e : CubeGaussObjects) {
    particlesTotal += e.getNumParticles();
  }
  for (auto e : CubeUniformObjects) {
    particlesTotal += e.getNumParticles();
  }
  for (auto e : SphereObjects) {
    particlesTotal += e.particlesTotal();
  }
  return particlesTotal;
}

void YamlParser::calcAutopasBox() {
    std::vector<double> XCoordinates;
    std::vector<double> YCoordinates;
    std::vector<double> ZCoordinates;

    for (auto e : CubeGridObjects) {
//get array of (vector of x), (vector of y), (vector of z)

    }
    for (auto e : CubeGaussObjects) {
    }
    for (auto e : CubeUniformObjects) {
    }
    for (auto e : SphereObjects) {
    }
    //search for smallest and biggest X,Y,Z coordinates and set BoxMin and BoxMax
}

const set<ContainerOption> &YamlParser::getContainerOptions() const { return containerOptions; }

const set<DataLayoutOption> &YamlParser::getDataLayoutOptions() const { return dataLayoutOptions; }

SelectorStrategyOption YamlParser::getSelectorStrategy() const { return selectorStrategy; }

const set<TraversalOption> &YamlParser::getTraversalOptions() const { return traversalOptions; }

TuningStrategyOption YamlParser::getTuningStrategyOption() const { return tuningStrategyOption; }

const set<Newton3Option> &YamlParser::getNewton3Options() const { return newton3Options; }

const autopas::NumberSet<double> &YamlParser::getCellSizeFactors() const { return *cellSizeFactors; }

double YamlParser::getCutoff() const { return cutoff; }

YamlParser::FunctorOption YamlParser::getFunctorOption() const { return functorOption; }

size_t YamlParser::getIterations() const { return iterations; }

spdlog::level::level_enum YamlParser::getLogLevel() const { return logLevel; }

bool YamlParser::getMeasureFlops() const { return measureFlops; }

unsigned int YamlParser::getTuningInterval() const { return tuningInterval; }

unsigned int YamlParser::getTuningSamples() const { return tuningSamples; }

unsigned int YamlParser::getTuningMaxEvidence() const { return tuningMaxEvidence; }

const string &YamlParser::getWriteVtk() const { return writeVTK; }

const string &YamlParser::getLogFileName() const { return logFileName; }

unsigned int YamlParser::getVerletRebuildFrequency() const { return verletRebuildFrequency; }

double YamlParser::getVerletSkinRadius() const { return verletSkinRadius; }

double YamlParser::getEpsilon() const { return epsilon; }

double YamlParser::getSigma() const { return sigma; }

double YamlParser::getDeltaT() const { return delta_t; }

double YamlParser::getMass() const { return mass; }

const vector<CubeGrid> &YamlParser::getCubeGrid() const { return CubeGridObjects; }

const vector<CubeGauss> &YamlParser::getCubeGauss() const { return CubeGaussObjects; }

const vector<CubeUniform> &YamlParser::getCubeUniform() const { return CubeUniformObjects; }

const vector<Sphere> &YamlParser::getSphere() const { return SphereObjects; }
