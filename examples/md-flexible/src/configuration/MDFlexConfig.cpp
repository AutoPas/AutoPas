/**
 * @file MDFlexConfig.cpp
 * @author F. Gratl
 * @date 18.10.2019
 */
#include "MDFlexConfig.h"

#include "MDFlexParser.h"
#include "src/ParticleSerializationTools.h"
#include "src/Thermostat.h"

namespace {
/**
 * Reads the next numberOfParticles lines from file. The lines are expected to contain xyz data.
 * @tparam dataType type of the data to be read. (Can be scalars or containers that support [])
 * @tparam size number of entries per dataType.
 * @param file
 * @param numberOfParticles
 * @return Vector of read data.
 */
template <class dataType, int size>
std::vector<dataType> readPayload(std::ifstream &file, size_t numberOfParticles) {
  std::vector<dataType> data(numberOfParticles);
  // loop over every line (=particle)
  for (size_t i = 0; i < numberOfParticles; ++i) {
    // loop over line (=coordinates)
    if constexpr (size == 1) {
      file >> data[i];
    } else {
      for (size_t j = 0; j < size; ++j) {
        file >> data[i][j];
      }
    }
  }
  return data;
}

/**
 * Searches the file word by word and sets the file accessor directly behind the first found position.
 * @param file
 * @param word
 */
void findWord(std::ifstream &file, const std::string &word) {
  std::string currentWord;
  while (currentWord != word) {
    file >> currentWord;
  }
}
}  // namespace

MDFlexConfig::MDFlexConfig(int argc, char **argv) {
  auto parserExitCode = MDFlexParser::parseInput(argc, argv, *this);
  if (parserExitCode != MDFlexParser::exitCodes::success) {
    if (parserExitCode == MDFlexParser::exitCodes::parsingError) {
      std::cout << "Error when parsing configuration file." << std::endl;
      exit(EXIT_FAILURE);
    }
    exit(EXIT_SUCCESS);
  }

  calcSimulationBox();

  if (tuningPhases.value > 0) {
    iterations.value = 0ul;
  }

  initializeParticlePropertiesLibrary();

  initializeObjects();
}

std::string MDFlexConfig::to_string() const {
  using namespace std;
  ostringstream os;

  auto passedContainerOptionsStr = autopas::utils::ArrayUtils::to_string(containerOptions.value);
  os << setw(valueOffset) << left << containerOptions.name << ":  " << passedContainerOptionsStr << endl;

  // since all containers are rebuilt only periodically print Verlet config always.
  os << setw(valueOffset) << left << verletRebuildFrequency.name << ":  " << verletRebuildFrequency.value << endl;
  os << setw(valueOffset) << left << verletSkinRadius.name << ":  " << verletSkinRadius.value << endl;
  if (passedContainerOptionsStr.find("luster") != std::string::npos) {
    os << setw(valueOffset) << left << verletClusterSize.name << ":  " << verletClusterSize.value << endl;
  }

  if (containerOptions.value.size() > 1 or traversalOptions.value.size() > 1 or dataLayoutOptions.value.size() > 1) {
    os << setw(valueOffset) << left << selectorStrategy.name << ":  " << selectorStrategy.value << endl;
  }

  os << setw(valueOffset) << left << dataLayoutOptions.name << ":  "
     << autopas::utils::ArrayUtils::to_string(dataLayoutOptions.value) << endl;
  os << setw(valueOffset) << left << traversalOptions.name << ":  "
     << autopas::utils::ArrayUtils::to_string(traversalOptions.value) << endl;
  os << setw(valueOffset) << left << tuningStrategyOption.name << ":  " << tuningStrategyOption.value << endl;
  if (tuningStrategyOption.value == autopas::TuningStrategyOption::bayesianSearch or
      tuningStrategyOption.value == autopas::TuningStrategyOption::bayesianClusterSearch) {
    os << setw(valueOffset) << left << acquisitionFunctionOption.name << ":  " << acquisitionFunctionOption.value
       << endl;
  }
  os << setw(valueOffset) << left << mpiStrategyOption.name << ":  " << mpiStrategyOption.value << endl;
  os << setw(valueOffset) << left << tuningInterval.name << ":  " << tuningInterval.value << endl;
  os << setw(valueOffset) << left << tuningSamples.name << ":  " << tuningSamples.value << endl;
  os << setw(valueOffset) << left << tuningMaxEvidence.name << ":  " << tuningMaxEvidence.value << endl;
  if (tuningStrategyOption.value == autopas::TuningStrategyOption::predictiveTuning) {
    os << setw(valueOffset) << left << relativeOptimumRange.name << ":  " << relativeOptimumRange.value << endl;
    os << setw(valueOffset) << left << maxTuningPhasesWithoutTest.name << ":  " << maxTuningPhasesWithoutTest.value
       << endl;
    os << setw(valueOffset) << left << relativeBlacklistRange.name << ":  " << relativeBlacklistRange.value << endl;
    os << setw(valueOffset) << left << evidenceFirstPrediction.name << ":  " << evidenceFirstPrediction.value << endl;
    os << setw(valueOffset) << left << extrapolationMethodOption.name << ":  " << extrapolationMethodOption.value
       << endl;
  }
  os << setw(valueOffset) << left << functorOption.name << ":  ";
  switch (functorOption.value) {
    case FunctorOption::lj12_6: {
      os << "Lennard-Jones (12-6)" << endl;
      break;
    }
    case FunctorOption::lj12_6_AVX: {
      os << "Lennard-Jones (12-6) AVX intrinsics" << endl;
      break;
    }
    case FunctorOption::lj12_6_Globals: {
      os << "Lennard-Jones (12-6) with globals" << endl;
      break;
    }
  }
  os << setw(valueOffset) << left << newton3Options.name << ":  "
     << autopas::utils::ArrayUtils::to_string(newton3Options.value) << endl;

  os << setw(valueOffset) << left << cutoff.name << ":  " << cutoff.value << endl;
  os << setw(valueOffset) << left << boxMin.name << ":  " << autopas::utils::ArrayUtils::to_string(boxMin.value)
     << endl;
  os << setw(valueOffset) << left << boxMax.name << ":  " << autopas::utils::ArrayUtils::to_string(boxMax.value)
     << endl;
  os << setw(valueOffset) << left << subdivideDimension.name << ":  "
     << autopas::utils::ArrayUtils::to_string(subdivideDimension.value) << endl;
  os << setw(valueOffset) << left << cellSizeFactors.name << ":  " << *cellSizeFactors.value << endl;
  os << setw(valueOffset) << left << deltaT.name << ":  " << deltaT.value << endl;
  // simulation length is either dictated by tuning phases or iterations
  if (tuningPhases.value > 0) {
    os << setw(valueOffset) << left << tuningPhases.name << ":  " << tuningPhases.value << endl;
  } else {
    os << setw(valueOffset) << left << iterations.name << ":  " << iterations.value << endl;
  }
  os << setw(valueOffset) << left << boolalpha << periodic.name << ":  " << periodic.value << endl;

  os << setw(valueOffset) << left << "Objects:" << endl;

  auto printObjectCollection = [](auto objectCollection, auto name, auto &os) {
    int objectId = 0;
    for (const auto &object : objectCollection) {
      os << "  " << name << ":" << endl;
      os << "    " << objectId << ":  " << endl;
      auto objectStr = object.to_string();
      // indent all lines of object
      objectStr = std::regex_replace(objectStr, std::regex("(^|\n)(.)"), "$1      $2");
      os << objectStr;  // no endl needed here because objectStr ends a line
      objectId++;
    }
  };

  printObjectCollection(cubeGridObjects, cubeGridObjectsStr, os);
  printObjectCollection(cubeGaussObjects, cubeGaussObjectsStr, os);
  printObjectCollection(cubeUniformObjects, cubeUniformObjectsStr, os);
  printObjectCollection(sphereObjects, sphereObjectsStr, os);

  if (not globalForceIsZero()) {
    os << setw(valueOffset) << left << globalForce.name << ":  "
       << autopas::utils::ArrayUtils::to_string(globalForce.value) << endl;
  }

  if (useThermostat.value) {
    os << useThermostat.name << ":" << endl;
    os << "  " << setw(valueOffset - 2) << left << initTemperature.name << ":  " << initTemperature.value << endl;
    os << "  " << setw(valueOffset - 2) << left << targetTemperature.name << ":  " << targetTemperature.value << endl;
    os << "  " << setw(valueOffset - 2) << left << deltaTemp.name << ":  " << deltaTemp.value << endl;
    os << "  " << setw(valueOffset - 2) << left << thermostatInterval.name << ":  " << thermostatInterval.value << endl;
    os << "  " << setw(valueOffset - 2) << left << addBrownianMotion.name << ":  " << addBrownianMotion.value << endl;
  }

  if (not vtkFileName.value.empty()) {
    os << setw(valueOffset) << left << vtkFileName.name << ":  " << vtkFileName.value << endl;
    os << setw(valueOffset) << left << vtkWriteFrequency.name << ":  " << vtkWriteFrequency.value << endl;
  }
  if (not checkpointfile.value.empty())
    os << setw(valueOffset) << left << checkpointfile.name << ":  " << checkpointfile.value << endl;

  os << setw(valueOffset) << logLevel.name << ":  " << (logLevel.value) << endl;

  os << setw(valueOffset) << dontMeasureFlops.name << ":  " << (not dontMeasureFlops.value) << endl;
  os << setw(valueOffset) << dontCreateEndConfig.name << ":  " << (not dontCreateEndConfig.value) << endl;
  os << setw(valueOffset) << dontShowProgressBar.name << ":  " << (dontShowProgressBar.value) << endl;
  return os.str();
}

void MDFlexConfig::calcSimulationBox() {
  const double interactionLength = cutoff.value + verletSkinRadius.value;

  // helper function so that we can do the same for every object collection
  // resizes the domain to the maximal extents of all objects
  auto resizeToObjectLimits = [&](const auto &objectCollection) {
    for (auto &object : objectCollection) {
      auto objectMin = object.getBoxMin();
      auto objectMax = object.getBoxMax();
      auto objectSpacing = object.getParticleSpacing();

      for (size_t i = 0; i < 3; ++i) {
        // pad domain such that periodic boundaries can work.
        // This is necessary if the given min/max is not at least half the spacing away of the farthest object.
        boxMin.value[i] = std::min(boxMin.value[i], objectMin[i] - objectSpacing / 2);
        boxMax.value[i] = std::max(boxMax.value[i], objectMax[i] + objectSpacing / 2);
      }
    }
  };

  resizeToObjectLimits(cubeGaussObjects);
  resizeToObjectLimits(cubeGridObjects);
  resizeToObjectLimits(cubeUniformObjects);
  resizeToObjectLimits(sphereObjects);
  resizeToObjectLimits(cubeClosestPackedObjects);

  // guarantee the box is at least of size interationLength
  for (int i = 0; i < 3; i++) {
    // needed for 2D Simulation, that BoxLength >= interactionLength for all Dimensions
    if (boxMax.value[i] - boxMin.value[i] < interactionLength) {
      std::cout << "WARNING: Simulation box in dimension " << i
                << " is shorter than interaction length and will be increased." << std::endl;
      boxMin.value[i] -= interactionLength / 2;
      boxMax.value[i] += interactionLength / 2;
    }
  }
}

void MDFlexConfig::addParticleType(unsigned long typeId, double epsilon, double sigma, double mass) {
  // check if type id is already existing and if there no error in input
  if (epsilonMap.value.count(typeId) == 1) {
    // check if type is already added
    if (epsilonMap.value.at(typeId) == epsilon and sigmaMap.value.at(typeId) == sigma and
        massMap.value.at(typeId) == mass) {
      return;
    } else {  // wrong initialization:
      throw std::runtime_error("Wrong Particle initialization: using same typeId for different properties");
    }
  } else {
    epsilonMap.value.emplace(typeId, epsilon);
    sigmaMap.value.emplace(typeId, sigma);
    massMap.value.emplace(typeId, mass);
  }
}

void MDFlexConfig::flushParticles() { _particles.clear(); }

void MDFlexConfig::initializeParticlePropertiesLibrary() {
  if (epsilonMap.value.empty()) {
    throw std::runtime_error("No properties found in particle properties library!");
  }

  if (epsilonMap.value.size() != sigmaMap.value.size() or epsilonMap.value.size() != massMap.value.size()) {
    throw std::runtime_error("Number of particle properties differ!");
  }

  _particlePropertiesLibrary = std::make_shared<ParticlePropertiesLibraryType>(cutoff.value);

  for (auto [type, epsilon] : epsilonMap.value) {
    _particlePropertiesLibrary->addType(type, epsilon, sigmaMap.value.at(type), massMap.value.at(type));
  }
  _particlePropertiesLibrary->calculateMixingCoefficients();
}

void MDFlexConfig::initializeObjects() {
  if (not checkpointfile.value.empty()) {
    loadParticlesFromCheckpoint();
  }
  for (const auto &object : cubeGridObjects) {
    object.generate(_particles);
  }
  for (const auto &object : cubeGaussObjects) {
    object.generate(_particles);
  }
  for (const auto &object : cubeUniformObjects) {
    object.generate(_particles);
  }
  for (const auto &object : sphereObjects) {
    object.generate(_particles);
  }
  for (const auto &object : cubeClosestPackedObjects) {
    object.generate(_particles);
  }
}

void MDFlexConfig::loadParticlesFromCheckpoint() {
  std::ifstream infile(checkpointfile.value);
  size_t numParticles;
  std::string dataType;

  if (not infile.is_open()) {
    std::cout << "Could not load checkpoint file " << checkpointfile.value << "." << std::endl;
    return;
  }

  findWord(infile, "POINTS");
  infile >> numParticles >> dataType;

  auto positions = readPayload<std::array<double, 3>, 3>(infile, numParticles);
  // next payload block is always preceded by the datatype
  findWord(infile, dataType);
  auto velocities = readPayload<std::array<double, 3>, 3>(infile, numParticles);
  findWord(infile, dataType);
  auto forces = readPayload<std::array<double, 3>, 3>(infile, numParticles);
  findWord(infile, "default");
  auto typeID = readPayload<size_t, 1>(infile, numParticles);

  // creating Particles from checkpoint:
  for (auto i = 0ul; i < numParticles; ++i) {
    ParticleType particle;

    particle.setR(positions[i]);
    particle.setV(velocities[i]);
    particle.setF(forces[i]);
    particle.setID(i);
    particle.setTypeId(typeID[i]);

    _particles.push_back(particle);
  }
}
