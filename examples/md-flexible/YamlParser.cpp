//
// Created by nicola on 15.07.19.
//

#include "YamlParser.h"

bool YamlParser::parseInput(int argc, char **argv) {
  using namespace std;
  bool displayHelp = false;
  int option, option_index;
  static struct option long_options[] = {{"filename", required_argument, nullptr, 'Y'},
                                         {"container", required_argument, nullptr, 'c'},
                                         {"selector-strategy", required_argument, nullptr, 'y'},
                                         {"cutoff", required_argument, nullptr, 'C'},
                                         {"cell-size-factor", required_argument, nullptr, 'a'},
                                         {"data-layout", required_argument, nullptr, 'd'},
                                         {"functor", required_argument, nullptr, 'f'},
                                         {"help", no_argument, nullptr, 'h'},
                                         {"iterations", required_argument, nullptr, 'i'},
                                         {"no-flops", no_argument, nullptr, 'F'},
                                         {"newton3", required_argument, nullptr, '3'},
                                         {"delta_t", required_argument, nullptr, 'D'},
                                         {"traversal", required_argument, nullptr, 't'},
                                         {"tuning-interval", required_argument, nullptr, 'I'},
                                         {"tuning-samples", required_argument, nullptr, 'S'},
                                         {"tuning-max-evidence", required_argument, nullptr, 'E'},
                                         {"tuning-strategy", required_argument, nullptr, 'T'},
                                         {"log-level", required_argument, nullptr, 'l'},
                                         {"log-file", required_argument, nullptr, 'L'},
                                         {"verlet-rebuild-frequency", required_argument, nullptr, 'v'},
                                         {"verlet-skin-radius", required_argument, nullptr, 'r'},
                                         {"vtk", required_argument, nullptr, 'w'},
                                         {nullptr, 0, nullptr, 0}};  // needed to signal the end of the array
  string strArg;
  // Yaml Parsing file parameter muss als erstes übergeben werden
  bool yamlparsed = false;
  option = getopt_long(argc, argv, "", long_options, &option_index);
  if (option == 'Y') {
    filename = optarg;
    try {
      YAML::LoadFile(filename);
      yamlparsed = true;
    } catch (const exception &) {
      cerr << "Error parsing Yaml File: " << filename << ", check filename or yaml Syntax" << endl;
    }
  } else {
    cout << "[INFO] No Yaml Parsing File specified" << endl << endl;
  }
  if (yamlparsed) {
    parseYamlFile();
  }
  optind = 1;
  while ((option = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    if (optarg != nullptr) strArg = optarg;
    transform(strArg.begin(), strArg.end(), strArg.begin(), ::tolower);
    switch (option) {
      case 'Y': {
        break;
      }
      case 'D': {
        try {
          delta_t = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing epsilon value: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case '3': {
        newton3Options = autopas::utils::StringUtils::parseNewton3Options(strArg, false);
        if (newton3Options.empty()) {
          cerr << "Unknown Newton3 option: " << strArg << endl;
          cerr << "Please use 'enabled' or 'disabled'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'c': {
        // overwrite default argument
        containerOptions = autopas::utils::StringUtils::parseContainerOptions(strArg, false);
        if (containerOptions.empty()) {
          cerr << "Unknown container option: " << strArg << endl;
          cerr << "Please use 'DirectSum', 'LinkedCells', 'VerletLists', 'VCells' or 'VCluster'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'C': {
        try {
          cutoff = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing cutoff Radius: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'a': {
        cellSizeFactors = autopas::utils::StringUtils::parseNumberSet(strArg);
        if (cellSizeFactors->isEmpty()) {
          cerr << "Error parsing cell size factors: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'd': {
        dataLayoutOptions = autopas::utils::StringUtils::parseDataLayout(strArg);
        if (dataLayoutOptions.empty()) {
          cerr << "Unknown data layouts: " << strArg << endl;
          cerr << "Please use 'AoS' or 'SoA'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'f': {
        if (strArg.find("avx") != string::npos) {
          functorOption = lj12_6_AVX;
        } else if (strArg.find("lj") != string::npos || strArg.find("lennard-jones") != string::npos) {
          functorOption = lj12_6;
        } else {
          cerr << "Unknown functor: " << strArg << endl;
          cerr << "Please use 'Lennard-Jones' or 'Lennard-Jones-AVX'" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'F': {
        measureFlops = false;
        break;
      }
      case 'h': {
        displayHelp = true;
        break;
      }
      case 'i': {
        try {
          iterations = stoul(strArg);
          if (iterations < 1) {
            cerr << "IterationNumber of iterations has to be a positive integer!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing number of iterations: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'I': {
        try {
          tuningInterval = (unsigned int)stoul(strArg);
          if (tuningInterval < 1) {
            cerr << "Tuning interval has to be a positive integer!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing tuning interval: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'y': {
        selectorStrategy = autopas::utils::StringUtils::parseSelectorStrategy(strArg);
        if (selectorStrategy == autopas::SelectorStrategyOption(-1)) {
          cerr << "Unknown Selector Strategy: " << strArg << endl;
          cerr << "Please use 'fastestAbs', 'fastestMean' or 'fastestMedian'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'l': {
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
          default: {
            cerr << "Unknown Log Level: " << strArg << endl;
            cerr << "Please use 'trace', 'debug', 'info', 'warning', 'error', 'critical' or 'off'." << endl;
            displayHelp = true;
          }
        }
        break;
      }
      case 'L': {
        logFileName = strArg;
        break;
      }
      case 'S': {
        try {
          tuningSamples = (unsigned int)stoul(strArg);
          if (tuningSamples < 1) {
            cerr << "Tuning samples has to be a positive integer!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing number of tuning samples: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'E': {
        try {
          tuningMaxEvidence = (unsigned int)stoul(strArg);
          if (tuningMaxEvidence < 1) {
            cerr << "Tuning max evidence has to be a positive integer!" << endl;
            displayHelp = true;
          }
        } catch (const exception &) {
          cerr << "Error parsing number of tuning max evidence: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 't': {
        traversalOptions = autopas::utils::StringUtils::parseTraversalOptions(strArg);
        if (traversalOptions.empty()) {
          cerr << "Unknown Traversal: " << strArg << endl;
          cerr << "Please use 'c08', 'c01', 'c18', 'sliced' or 'direct'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'T': {
        tuningStrategyOption = autopas::utils::StringUtils::parseTuningStrategyOption(strArg);
        if (tuningStrategyOption == autopas::TuningStrategyOption(-1)) {
          cerr << "Unknown Tuning Strategy: " << strArg << endl;
          cerr << "Please use 'full-search' or 'bayesian-search'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'v': {
        try {
          verletRebuildFrequency = (unsigned int)stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet-rebuild-frequency: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'w': {
        writeVTK = strArg;
        break;
      }
      case 'r': {
        try {
          verletSkinRadius = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet-skin-radius: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      default: {
        // error message handled by getopt
        displayHelp = true;
      }
    }
  }
  if (displayHelp) {
    cout << "Usage: " << argv[0] << endl;
    for (auto o : long_options) {
      if (o.name == nullptr) continue;
      cout << "    --" << setw(valueOffset + 2) << left << o.name;
      if (o.has_arg) cout << "option";
      cout << endl;
    }
    return false;
  }
  return true;
}

void YamlParser::parseYamlFile() {
  using namespace autopas;
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
    if (strArg.find("avx") != std::string::npos) {
      this->functorOption = lj12_6_AVX;
    } else if (strArg.find("lj") != std::string::npos || strArg.find("lennard-jones") != std::string::npos) {
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
    CubeGridObjects = {};
    for (YAML::const_iterator it = config["Objects"].begin(); it != config["Objects"].end(); ++it) {
      if (it->first.as<std::string>() == "CubeGrid") {
        CubeGrid C({it->second["particles-per-Dim"][0].as<unsigned long>(),
                    it->second["particles-per-Dim"][1].as<unsigned long>(),
                    it->second["particles-per-Dim"][2].as<unsigned long>()},
                   it->second["particleSpacing"].as<double>(),
                   {it->second["velocity"][0].as<double>(), it->second["velocity"][1].as<double>(),
                    it->second["velocity"][2].as<double>()},
                   {it->second["center"][0].as<double>(), it->second["center"][1].as<double>(),
                    it->second["center"][2].as<double>()});
        CubeGridObjects.emplace_back(C);
        continue;
      }
      if (it->first.as<std::string>() == "CubeGauss") {
        CubeGauss C(it->second["numberOfParticles"].as<size_t>(),
                    {it->second["box-length"][0].as<double>(), it->second["box-length"][1].as<double>(),
                     it->second["box-length"][2].as<double>()},
                    it->second["distribution-mean"].as<double>(), it->second["distribution-stddev"].as<double>(),
                    {it->second["velocity"][0].as<double>(), it->second["velocity"][1].as<double>(),
                     it->second["velocity"][2].as<double>()},
                    {it->second["center"][0].as<double>(), it->second["center"][1].as<double>(),
                     it->second["center"][2].as<double>()});
        CubeGaussObjects.emplace_back(C);
        continue;
      }
      if (it->first.as<std::string>() == "CubeUniform") {
        CubeUniform C(it->second["numberOfParticles"].as<size_t>(),
                      {it->second["box-length"][0].as<double>(), it->second["box-length"][1].as<double>(),
                       it->second["box-length"][2].as<double>()},
                      {it->second["velocity"][0].as<double>(), it->second["velocity"][1].as<double>(),
                       it->second["velocity"][2].as<double>()},
                      {it->second["center"][0].as<double>(), it->second["center"][1].as<double>(),
                       it->second["center"][2].as<double>()});
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
  this->calcAutopasBox();
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
  using namespace std;
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
  cout << setw(valueOffset) << left << "Allowed traversals"
       << ":  " << iterableToString(traversalOptions) << endl;
  cout << setw(valueOffset) << left << "Tuning Strategy"
       << ":  " << autopas::utils::StringUtils::to_string(tuningStrategyOption) << endl;
  cout << setw(valueOffset) << left << "Tuning Interval"
       << ":  " << tuningInterval << endl;
  cout << setw(valueOffset) << left << "Tuning Samples"
       << ":  " << tuningSamples << endl;
  cout << setw(valueOffset) << left << "Tuning Max evidence"
       << ":  " << tuningMaxEvidence << endl;
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
       << ":  " << iterations << endl
       << endl;

  cout << setw(valueOffset) << left << "Object Generation:" << endl;
  int i = 1;
  for (auto c : CubeGridObjects) {
    cout << "-Cube Grid Nr " << i << ":  " << endl;
    c.printConfig();
    i++;
  }
  i = 1;
  for (auto c : CubeGaussObjects) {
    cout << "-Cube Gauss Nr" << i << ":  " << endl;
    c.printConfig();
    i++;
  }
  i = 1;
  for (auto c : CubeUniformObjects) {
    cout << "-Cube Uniform Nr " << i << ":  " << endl;
    c.printConfig();
    i++;
  }
  i = 1;
  for (auto c : SphereObjects) {
    cout << "-Sphere Nr " << i << ":  " << endl;
    c.printConfig();
    i++;
  }
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
  std::vector<double> XcoordMin;
  std::vector<double> YcoordMin;
  std::vector<double> ZcoordMin;
  std::vector<double> XcoordMax;
  std::vector<double> YcoordMax;
  std::vector<double> ZcoordMax;
  for (auto c : CubeGridObjects) {
    XcoordMin.emplace_back(c.getBoxMin()[0]);
    YcoordMin.emplace_back(c.getBoxMin()[1]);
    ZcoordMin.emplace_back(c.getBoxMin()[2]);
    XcoordMax.emplace_back(c.getBoxMax()[0]);
    YcoordMax.emplace_back(c.getBoxMax()[1]);
    ZcoordMax.emplace_back(c.getBoxMax()[2]);
  }
  for (auto c : CubeGaussObjects) {
    XcoordMin.emplace_back(c.getBoxMin()[0]);
    YcoordMin.emplace_back(c.getBoxMin()[1]);
    ZcoordMin.emplace_back(c.getBoxMin()[2]);
    XcoordMax.emplace_back(c.getBoxMax()[0]);
    YcoordMax.emplace_back(c.getBoxMax()[1]);
    ZcoordMax.emplace_back(c.getBoxMax()[2]);
  }
  for (auto c : CubeUniformObjects) {
    XcoordMin.emplace_back(c.getBoxMin()[0]);
    YcoordMin.emplace_back(c.getBoxMin()[1]);
    ZcoordMin.emplace_back(c.getBoxMin()[2]);
    XcoordMax.emplace_back(c.getBoxMax()[0]);
    YcoordMax.emplace_back(c.getBoxMax()[1]);
    ZcoordMax.emplace_back(c.getBoxMax()[2]);
  }
  for (auto c : SphereObjects) {
    XcoordMin.emplace_back(c.getBoxMin()[0]);
    YcoordMin.emplace_back(c.getBoxMin()[1]);
    ZcoordMin.emplace_back(c.getBoxMin()[2]);
    XcoordMax.emplace_back(c.getBoxMax()[0]);
    YcoordMax.emplace_back(c.getBoxMax()[1]);
    ZcoordMax.emplace_back(c.getBoxMax()[2]);
  }
  if (not XcoordMin.empty()) {
    BoxMin = {*std::min_element(XcoordMin.begin(), XcoordMin.end()),
              *std::min_element(YcoordMin.begin(), YcoordMin.end()),
              *std::min_element(ZcoordMin.begin(), ZcoordMin.end())};
    BoxMax = {*std::max_element(XcoordMax.begin(), XcoordMax.end()),
              *std::max_element(YcoordMax.begin(), YcoordMax.end()),
              *std::max_element(ZcoordMax.begin(), ZcoordMax.end())};
  }
}

const std::set<autopas::ContainerOption> &YamlParser::getContainerOptions() const { return containerOptions; }

const std::set<autopas::DataLayoutOption> &YamlParser::getDataLayoutOptions() const { return dataLayoutOptions; }

autopas::SelectorStrategyOption YamlParser::getSelectorStrategy() const { return selectorStrategy; }

const std::set<autopas::TraversalOption> &YamlParser::getTraversalOptions() const { return traversalOptions; }

autopas::TuningStrategyOption YamlParser::getTuningStrategyOption() const { return tuningStrategyOption; }

const std::set<autopas::Newton3Option> &YamlParser::getNewton3Options() const { return newton3Options; }

const autopas::NumberSet<double> &YamlParser::getCellSizeFactors() const { return *cellSizeFactors; }

double YamlParser::getCutoff() const { return cutoff; }

YamlParser::FunctorOption YamlParser::getFunctorOption() const { return functorOption; }

size_t YamlParser::getIterations() const { return iterations; }

spdlog::level::level_enum YamlParser::getLogLevel() const { return logLevel; }

bool YamlParser::getMeasureFlops() const { return measureFlops; }

unsigned int YamlParser::getTuningInterval() const { return tuningInterval; }

unsigned int YamlParser::getTuningSamples() const { return tuningSamples; }

unsigned int YamlParser::getTuningMaxEvidence() const { return tuningMaxEvidence; }

const std::string &YamlParser::getWriteVtk() const { return writeVTK; }

const std::string &YamlParser::getLogFileName() const { return logFileName; }

unsigned int YamlParser::getVerletRebuildFrequency() const { return verletRebuildFrequency; }

double YamlParser::getVerletSkinRadius() const { return verletSkinRadius; }

double YamlParser::getEpsilon() const { return epsilon; }

double YamlParser::getSigma() const { return sigma; }

double YamlParser::getDeltaT() const { return delta_t; }

double YamlParser::getMass() const { return mass; }

const std::vector<CubeGrid> &YamlParser::getCubeGrid() const { return CubeGridObjects; }

const std::vector<CubeGauss> &YamlParser::getCubeGauss() const { return CubeGaussObjects; }

const std::vector<CubeUniform> &YamlParser::getCubeUniform() const { return CubeUniformObjects; }

const std::vector<Sphere> &YamlParser::getSphere() const { return SphereObjects; }

const std::array<double, 3> &YamlParser::getBoxMin() const { return BoxMin; }

const std::array<double, 3> &YamlParser::getBoxMax() const { return BoxMax; }

void YamlParser::setFilename(const std::string &filename) { YamlParser::filename = filename; }
