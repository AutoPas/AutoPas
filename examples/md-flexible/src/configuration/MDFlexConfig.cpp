/**
 * @file MDFlexConfig.cpp
 * @author F. Gratl
 * @date 18.10.2019
 */
#include "MDFlexConfig.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <string>
#include <vector>

#include "MDFlexParser.h"
#include "autopas/utils/WrapMPI.h"
#include "autopas/utils/isSmartPointer.h"
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
  std::vector<char> separators{' ', '/', '.', ',', '?', '!', '"', '\'', '<', '>', '=', ':', ';', '\n', '\t', '\r'};
  std::string currentWord;
  while (not file.eof() and currentWord != word) {
    char currentChar = file.get();
    if (std::find(separators.begin(), separators.end(), currentChar) != separators.end()) {
      currentWord = "";
    } else {
      currentWord += currentChar;
    }
  }
}

/**
 * Evaluates the number of pieces that the checkpoint contains. This is the number of vtu files belonging to the pvtu.
 * This is also the number of ranks that was used during the creation of the data file.
 *
 * @param filename The name of the pvtu checkpoint file.
 * @return Number of pieces.
 */
size_t getNumPiecesInCheckpoint(const std::string &filename) {
  size_t numPiecesInCheckpoint{0ul};
  std::ifstream inputStream(filename);
  findWord(inputStream, "Piece");
  while (not inputStream.eof()) {
    ++numPiecesInCheckpoint;
    findWord(inputStream, "Piece");
  }
  return numPiecesInCheckpoint;
}

/**
 * Loads the particles from a checkpoint written with a specific rank.
 * @param filename The name of the pvtu file.
 * @param rank The rank which created the respective vtu file.
 * @param particles Container for the particles recorded in the respective vts file.
 */
void loadParticlesFromRankRecord(std::string_view filename, const size_t &rank, std::vector<ParticleType> &particles) {
  const size_t endOfPath = filename.find_last_of('/');
  const auto filePath = filename.substr(0ul, endOfPath);
  const auto fileBasename = filename.substr(endOfPath + 1);

  // infer checkpoint data from the filename
  const size_t endOfScenarioName = fileBasename.find_last_of('_');
  const auto checkpointScenarioName = fileBasename.substr(0, endOfScenarioName);
  // +1 because we want to skip the separating '_'
  const auto checkpointIteration =
      fileBasename.substr(endOfScenarioName + 1, fileBasename.find_last_of('.') - endOfScenarioName - 1);

  std::string rankFilename{};
  rankFilename.append(filePath)
      .append("/data/")
      .append(checkpointScenarioName)
      .append("_")
      .append(std::to_string(rank))
      .append("_")
      .append(checkpointIteration)
      .append(".vtu");

  std::ifstream inputStream(rankFilename);
  if (not inputStream.is_open()) {
    throw std::runtime_error("Could not open rank-checkpoint file: " + rankFilename);
    return;
  }

  size_t numParticles{0ul};
  findWord(inputStream, "NumberOfPoints");
  inputStream.ignore(std::numeric_limits<std::streamsize>::max(), '"');
  inputStream >> numParticles;

  // sanity check
  if (numParticles == 0) {
    throw std::runtime_error("Could not determine the number of particles in the checkpoint file " + rankFilename);
    return;
  }

  findWord(inputStream, "velocities");
  inputStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  auto velocities = readPayload<std::array<double, 3>, 3>(inputStream, numParticles);

  findWord(inputStream, "forces");
  inputStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  auto forces = readPayload<std::array<double, 3>, 3>(inputStream, numParticles);

  findWord(inputStream, "typeIds");
  inputStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  auto typeIds = readPayload<size_t, 1>(inputStream, numParticles);

  findWord(inputStream, "ids");
  inputStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  auto ids = readPayload<size_t, 1>(inputStream, numParticles);

  findWord(inputStream, "positions");
  inputStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  auto positions = readPayload<std::array<double, 3>, 3>(inputStream, numParticles);

  // creating Particles from checkpoint:
  for (auto i = 0ul; i < numParticles; ++i) {
    ParticleType particle;

    particle.setR(positions[i]);
    particle.setV(velocities[i]);
    particle.setF(forces[i]);
    particle.setID(ids[i]);
    particle.setTypeId(typeIds[i]);

    particles.push_back(particle);
  }
}
}  // namespace

MDFlexConfig::MDFlexConfig(int argc, char **argv) {
  auto parserExitCode = MDFlexParser::parseInput(argc, argv, *this);
  if (parserExitCode != MDFlexParser::exitCodes::success) {
    if (parserExitCode == MDFlexParser::exitCodes::parsingError) {
      std::cerr << "Error when parsing configuration file." << std::endl;
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
  using autopas::utils::ArrayUtils::operator<<;
  ostringstream os;
  os << boolalpha;

  auto printOption = [&](const auto &option, int offsetFromOffset = 0) {
    os << setw(valueOffset + offsetFromOffset) << left << option.name << ":  ";
    // depending on whether option.value is a smart-pointer or not dereference it
    if constexpr (autopas::utils::is_smart_ptr<decltype(option.value)>::value) {
      os << *option.value;
    } else {
      os << option.value;
    }
    os << endl;
  };
  printOption(containerOptions);

  // since all containers are rebuilt only periodically print Verlet config always.
  printOption(verletRebuildFrequency);
  printOption(verletSkinRadiusPerTimestep);
  printOption(fastParticlesThrow);
  const auto passedContainerOptionsStr = autopas::utils::ArrayUtils::to_string(containerOptions.value);
  if (passedContainerOptionsStr.find("luster") != std::string::npos) {
    printOption(verletClusterSize);
  }

  if (containerOptions.value.size() > 1 or traversalOptions.value.size() > 1 or dataLayoutOptions.value.size() > 1) {
    printOption(selectorStrategy);
  }

  printOption(dataLayoutOptions);
  printOption(traversalOptions);
  printOption(tuningStrategyOption);
  if (tuningStrategyOption.value == autopas::TuningStrategyOption::bayesianSearch or
      tuningStrategyOption.value == autopas::TuningStrategyOption::bayesianClusterSearch) {
    printOption(acquisitionFunctionOption);
  }
  printOption(mpiStrategyOption);
  if (mpiStrategyOption.value == autopas::MPIStrategyOption::divideAndConquer) {
    printOption(MPITuningMaxDifferenceForBucket);
    printOption(MPITuningWeightForMaxDensity);
  }
  printOption(tuningInterval);
  printOption(tuningSamples);
  printOption(tuningMaxEvidence);
  if (tuningStrategyOption.value == autopas::TuningStrategyOption::predictiveTuning) {
    printOption(relativeOptimumRange);
    printOption(maxTuningPhasesWithoutTest);
    printOption(relativeBlacklistRange);
    printOption(evidenceFirstPrediction);
    printOption(extrapolationMethodOption);
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
    case FunctorOption::lj12_6_SVE: {
      os << "Lennard-Jones (12-6) SVE intrinsics" << endl;
      break;
    }
    case FunctorOption::lj12_6_Globals: {
      os << "Lennard-Jones (12-6) with globals" << endl;
      break;
    }
  }
  printOption(newton3Options);
  printOption(cutoff);
  printOption(boxMin);
  printOption(boxMax);
  printOption(cellSizeFactors);
  printOption(deltaT);
  // simulation length is either dictated by tuning phases or iterations
  if (tuningPhases.value > 0) {
    printOption(tuningPhases);
  } else {
    printOption(iterations);
  }
  printOption(boundaryOption);

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
  printObjectCollection(cubeClosestPackedObjects, cubeClosestPackedObjectsStr, os);
  printObjectCollection(sphereObjects, sphereObjectsStr, os);

  if (not globalForceIsZero()) {
    printOption(globalForce);
  }

  if (useThermostat.value) {
    os << useThermostat.name << ":" << endl;
    // since these parameters are
    constexpr int indentWidth = 2;
    const auto indent = std::string(indentWidth, ' ');
    os << indent;
    printOption(initTemperature, -indentWidth);
    os << indent;
    printOption(targetTemperature, -indentWidth);
    os << indent;
    printOption(deltaTemp, -indentWidth);
    os << indent;
    printOption(thermostatInterval, -indentWidth);
    os << indent;
    printOption(addBrownianMotion, -indentWidth);
  }

  if (not vtkFileName.value.empty()) {
    printOption(vtkFileName);
    printOption(vtkWriteFrequency);
  }
  if (not checkpointfile.value.empty()) {
    printOption(checkpointfile);
  }

  os << setw(valueOffset) << left << logLevel.name << ":  " << spdlog::level::to_string_view(logLevel.value).data()
     << endl;

  os << setw(valueOffset) << left << dontMeasureFlops.name << ":  " << (not dontMeasureFlops.value) << endl;
  os << setw(valueOffset) << left << dontCreateEndConfig.name << ":  " << (not dontCreateEndConfig.value) << endl;
  printOption(dontShowProgressBar);
  printOption(loadBalancer);
  printOption(loadBalancingInterval);
  printOption(subdivideDimension);
  return os.str();
}

void MDFlexConfig::calcSimulationBox() {
  const double interactionLength = cutoff.value + verletSkinRadiusPerTimestep.value * verletRebuildFrequency.value;

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

void MDFlexConfig::loadParticlesFromCheckpoint(const size_t &rank, const size_t &communicatorSize) {
  const std::string &filename = checkpointfile.value;

  std::ifstream inputStream(filename);
  if (not inputStream.is_open()) {
    throw std::runtime_error("Could not open checkpoint file: " + filename);
    return;
  }

  size_t checkpointCommunicatorSize{getNumPiecesInCheckpoint(filename)};
  if (communicatorSize == checkpointCommunicatorSize) {
    loadParticlesFromRankRecord(filename, rank, _particles);
  } else {
    for (size_t i = 0; i < checkpointCommunicatorSize; ++i) {
      loadParticlesFromRankRecord(filename, i, _particles);
    }
  }
}
