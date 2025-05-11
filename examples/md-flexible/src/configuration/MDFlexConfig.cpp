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
#include "autopas/utils/ArrayUtils.h"
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
 * @param particles Container for the particles recorded in the respective vtu file.
 */
void loadParticlesFromRankRecord(std::string_view filename, const size_t &rank, std::vector<ParticleType> &particles) {
  // Parse filename into path and filename: some/long/path/filename
  const size_t endOfPath = filename.find_last_of('/');
  const auto filePath = filename.substr(0ul, endOfPath);
  const auto fileBasename = filename.substr(endOfPath + 1);

  // infer checkpoint data from the filename: scenario_name_iteration.xxx
  const size_t endOfScenarioName = fileBasename.find_last_of('_');
  const auto checkpointScenarioName = fileBasename.substr(0, endOfScenarioName);
  // +1 because we want to skip the separating '_'
  const auto checkpointIteration =
      fileBasename.substr(endOfScenarioName + 1, fileBasename.find_last_of('.') - endOfScenarioName - 1);

  // Reconstruct the expected vtu (=piece) name from what we parsed.
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
  const auto velocities = readPayload<std::array<double, 3>, 3>(inputStream, numParticles);

  findWord(inputStream, "forces");
  inputStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  const auto forces = readPayload<std::array<double, 3>, 3>(inputStream, numParticles);

#if MD_FLEXIBLE_MODE == MULTISITE
  findWord(inputStream, "quaternions");
  inputStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  const auto quaternions = readPayload<std::array<double, 4>, 4>(inputStream, numParticles);

  findWord(inputStream, "angularVelocities");
  inputStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  const auto angularVelocities = readPayload<std::array<double, 3>, 3>(inputStream, numParticles);

  findWord(inputStream, "torques");
  inputStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  const auto torques = readPayload<std::array<double, 3>, 3>(inputStream, numParticles);
#endif

  findWord(inputStream, "typeIds");
  inputStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  const auto typeIds = readPayload<size_t, 1>(inputStream, numParticles);

  findWord(inputStream, "ids");
  inputStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  const auto ids = readPayload<size_t, 1>(inputStream, numParticles);

  findWord(inputStream, "positions");
  inputStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  const auto positions = readPayload<std::array<double, 3>, 3>(inputStream, numParticles);

  // creating Particles from checkpoint:
  particles.reserve(particles.size() + numParticles);
  for (auto i = 0ul; i < numParticles; ++i) {
    ParticleType particle;

    particle.setR(positions[i]);
    particle.setV(velocities[i]);
    particle.setF(forces[i]);
    particle.setID(ids[i]);
    particle.setTypeId(typeIds[i]);

#if MD_FLEXIBLE_MODE == MULTISITE
    particle.setQuaternion(quaternions[i]);
    particle.setAngularVel(angularVelocities[i]);
    particle.setTorque(torques[i]);
#endif

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

#if MD_FLEXIBLE_MODE == MULTISITE
  os << "# Running multi-site MD simulation.\n" << endl;
#else
  os << "# Running single-site MD simulation.\n" << endl;
#endif

  printOption(containerOptions);

  // since all containers are rebuilt only periodically print Verlet config always.
#ifndef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
  printOption(fastParticlesThrow);
#endif
  printOption(verletRebuildFrequency);
  printOption(verletSkinRadius);
  const auto passedContainerOptionsStr = autopas::utils::ArrayUtils::to_string(containerOptions.value);
  if (passedContainerOptionsStr.find("luster") != std::string::npos) {
    printOption(verletClusterSize);
  }

  if (containerOptions.value.size() > 1 or traversalOptions.value.size() > 1 or dataLayoutOptions.value.size() > 1) {
    printOption(selectorStrategy);
  }

  printOption(tuningStrategyOptions);

  // helper function to check if any options of a given list is in the tuningStrategyOptions.
  auto tuningStrategyOptionsContainAnyOf = [&](const std::vector<autopas::TuningStrategyOption> &needles) {
    return std::any_of(tuningStrategyOptions.value.begin(), tuningStrategyOptions.value.end(), [&](const auto &lhs) {
      return std::any_of(needles.begin(), needles.end(), [&](const auto &rhs) { return lhs == rhs; });
    });
  };

  if (tuningStrategyOptionsContainAnyOf({
          autopas::TuningStrategyOption::bayesianSearch,
          autopas::TuningStrategyOption::bayesianClusterSearch,
      })) {
    printOption(acquisitionFunctionOption);
  }
  if (tuningStrategyOptionsContainAnyOf({autopas::TuningStrategyOption::mpiDivideAndConquer})) {
    printOption(MPITuningMaxDifferenceForBucket);
    printOption(MPITuningWeightForMaxDensity);
  }
  printOption(tuningMetricOption);
  printOption(tuningInterval);
  printOption(tuningSamples);
  printOption(earlyStoppingFactor);
  printOption(useLOESSSmoothening);
  if (tuningStrategyOptionsContainAnyOf({
          autopas::TuningStrategyOption::randomSearch,
          autopas::TuningStrategyOption::bayesianSearch,
          autopas::TuningStrategyOption::bayesianClusterSearch,
      })) {
    printOption(tuningMaxEvidence);
  }
  if (tuningStrategyOptionsContainAnyOf({autopas::TuningStrategyOption::predictiveTuning})) {
    printOption(relativeOptimumRange);
    printOption(maxTuningPhasesWithoutTest);
    printOption(evidenceFirstPrediction);
    printOption(extrapolationMethodOption);
  }
  if (tuningStrategyOptionsContainAnyOf({autopas::TuningStrategyOption::slowConfigFilter})) {
    printOption(relativeBlacklistRange);
  }
  if (tuningStrategyOptionsContainAnyOf({autopas::TuningStrategyOption::ruleBasedTuning})) {
    printOption(ruleFilename);
  }
  if (tuningStrategyOptionsContainAnyOf({autopas::TuningStrategyOption::fuzzyTuning})) {
    printOption(fuzzyRuleFilename);
  }

  // TODO: C++20 Use contains instead of count
  if (getInteractionTypes().count(autopas::InteractionTypeOption::pairwise)) {
    os << setw(valueOffset) << left << "PairwiseInteraction:" << endl;
    constexpr int indentWidth = 2;
    const auto indent = std::string(indentWidth, ' ');
    os << indent;
    os << setw(valueOffset - indentWidth) << left << functorOption.name << ":  ";
    switch (functorOption.value) {
      case FunctorOption::none: {
        os << "None selected" << endl;
        break;
      }
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
      case FunctorOption::cellCollective: {
        os << "Cell Collective" << endl;
        break;
      }
    }
    os << indent;
    printOption(traversalOptions, -indentWidth);
    os << indent;
    printOption(dataLayoutOptions, -indentWidth);
    os << indent;
    printOption(newton3Options, -indentWidth);
  }

  // TODO c++20: use contains instead of count
  if (getInteractionTypes().count(autopas::InteractionTypeOption::triwise)) {
    os << setw(valueOffset) << left << "ThreeBodyInteraction:" << endl;
    constexpr int indentWidth = 2;
    const auto indent = std::string(indentWidth, ' ');
    os << indent;
    os << setw(valueOffset - indentWidth) << left << functorOption3B.name << ":  ";
    switch (functorOption3B.value) {
      case FunctorOption3B::none: {
        os << "None selected" << endl;
        break;
      }
      case FunctorOption3B::at: {
        os << "Axilrod-Teller" << endl;
        break;
      }
    }
    os << indent;
    printOption(traversalOptions3B, -indentWidth);
    os << indent;
    printOption(dataLayoutOptions3B, -indentWidth);
    os << indent;
    printOption(newton3Options3B, -indentWidth);
  }

  printOption(cutoff);
  printOption(boxMin);
  printOption(boxMax);
  printOption(cellSizeFactors);
  printOption(deltaT);
  printOption(pauseSimulationDuringTuning);
  printOption(sortingThreshold);
  // simulation length is either dictated by tuning phases or iterations
  if (tuningPhases.value > 0) {
    printOption(tuningPhases);
  } else {
    printOption(iterations);
  }
  printOption(boundaryOption);

  os << setw(valueOffset) << left << "Sites:" << endl;
  for (auto [siteId, epsilon] : epsilonMap.value) {
    os << "  " << siteId << ":" << endl;
    os << "    " << setw(valueOffset - 4) << left << epsilonMap.name << ":  " << epsilon << endl;
    os << "    " << setw(valueOffset - 4) << left << sigmaMap.name << ":  " << sigmaMap.value.at(siteId) << endl;
    os << "    " << setw(valueOffset - 4) << left << nuMap.name << ":  " << nuMap.value.at(siteId) << endl;
    os << "    " << setw(valueOffset - 4) << left << massMap.name << ":  " << massMap.value.at(siteId) << endl;
  }
#if MD_FLEXIBLE_MODE == MULTISITE
  os << setw(valueOffset) << left << "Molecules:" << endl;
  for (auto [molId, molToSiteId] : molToSiteIdMap) {
    os << "  " << molId << ":" << endl;
    os << "    " << setw(valueOffset - 4) << left << moleculeToSiteIdStr << ":  " << molToSiteId << endl;
    auto sitePosIter = molToSitePosMap.at(molId).begin();
    os << "    " << setw(valueOffset - 4) << left << moleculeToSitePosStr << ":  [" << *sitePosIter;
    sitePosIter++;
    for (; sitePosIter != molToSitePosMap.at(molId).end(); sitePosIter++) {
      os << ", " << *sitePosIter;
    }
    os << "]" << endl;
    os << "    " << setw(valueOffset - 4) << left << momentOfInertiaStr << ":  " << massMap.value.at(molId) << endl;
  }
#endif

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

  os << setw(valueOffset) << left << useTuningLogger.name << ":  " << useTuningLogger.value << endl;
  os << setw(valueOffset) << left << outputSuffix.name << ":  " << outputSuffix.value << endl;

  os << setw(valueOffset) << left << logLevel.name << ":  " << spdlog::level::to_string_view(logLevel.value).data()
     << endl;

  os << setw(valueOffset) << left << dontCreateEndConfig.name << ":  " << (not dontCreateEndConfig.value) << endl;
  printOption(dontShowProgressBar);
  printOption(loadBalancer);
  printOption(loadBalancingInterval);
  printOption(subdivideDimension);
  printOption(energySensorOption);
  return os.str();
}

void MDFlexConfig::calcSimulationBox() {
  const double interactionLength = cutoff.value + verletSkinRadius.value;
  const auto preBoxMin = boxMin.value;
  const auto preBoxMax = boxMax.value;

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

  if (boxMin.value != preBoxMin or boxMax.value != preBoxMax) {
    std::cout << "WARNING: Simulation box increased due to particles being too close to the boundaries." << std::endl;
  }

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

void MDFlexConfig::addSiteType(unsigned long siteId, double mass) {
  // check if siteId is already existing and if there is no error in input
  if (massMap.value.count(siteId) == 1) {
    // check if type is already added
    if (autopas::utils::Math::isNearRel(massMap.value.at(siteId), mass)) {
      return;
    } else {  // wrong initialization:
      throw std::runtime_error("Wrong Particle initialization: using same siteId for different properties");
    }
  } else {
    massMap.value.emplace(siteId, mass);
  }
}

void MDFlexConfig::addLJParametersToSite(unsigned long siteId, double epsilon, double sigma) {
  // check if siteId was already declared and mass was specified
  if (massMap.value.count(siteId) == 1) {
    if (epsilonMap.value.count(siteId) == 1) {
      if (autopas::utils::Math::isNearRel(epsilonMap.value.at(siteId), epsilon) and
          autopas::utils::Math::isNearRel(sigmaMap.value.at(siteId), sigma)) {
        return;
      } else {
        throw std::runtime_error("Wrong Particle initialization: using same siteId for different properties");
      }
    } else {
      epsilonMap.value.emplace(siteId, epsilon);
      sigmaMap.value.emplace(siteId, sigma);
    }
  } else {
    throw std::runtime_error("Initializing LJ parameters to a non existing site-ID");
  }
}

void MDFlexConfig::addATParametersToSite(unsigned long siteId, double nu) {
  // check if siteId was already declared and mass was specified
  if (massMap.value.count(siteId) == 1) {
    if (nuMap.value.count(siteId) == 1) {
      if (autopas::utils::Math::isNearRel(nuMap.value.at(siteId), nu)) {
        return;
      } else {
        throw std::runtime_error("Wrong Particle initialization: using same siteId for different properties");
      }
    } else {
      nuMap.value.emplace(siteId, nu);
    }
  } else {
    throw std::runtime_error("Initializing AT parameters to a non existing site-ID");
  }
}

void MDFlexConfig::addMolType(unsigned long molId, const std::vector<unsigned long> &siteIds,
                              const std::vector<std::array<double, 3>> &relSitePos,
                              std::array<double, 3> momentOfInertia) {
#if MD_FLEXIBLE_MODE == MULTISITE
  // check if siteId is already existing and if there no error in input
  if (molToSiteIdMap.count(molId) == 1) {
    // check if type is already added
    if (autopas::utils::ArrayMath::isEqual(molToSiteIdMap.at(molId), siteIds) and
        autopas::utils::ArrayMath::isNearRel(molToSitePosMap.at(molId), relSitePos) and
        autopas::utils::ArrayMath::isNearRel(momentOfInertiaMap.at(molId), momentOfInertia)) {
      return;
    } else {  // wrong initialization:
      throw std::runtime_error("Wrong Particle initialization: using same molId for different properties");
    }
  } else {
    molToSiteIdMap.emplace(molId, siteIds);
    molToSitePosMap.emplace(molId, relSitePos);
    momentOfInertiaMap.emplace(molId, momentOfInertia);
  }
#else
  throw std::runtime_error(
      "MDFlexConfig::addMolType was used without support for multi-site simulations being compiled");
#endif
}

void MDFlexConfig::flushParticles() { particles.clear(); }

void MDFlexConfig::initializeParticlePropertiesLibrary() {
  _particlePropertiesLibrary = std::make_shared<ParticlePropertiesLibraryType>(cutoff.value);

  // initialize site Ids with mandatory mass parameter
  for (auto [siteTypeId, mass] : massMap.value) {
    _particlePropertiesLibrary->addSiteType(siteTypeId, massMap.value.at(siteTypeId));
  }
  // check size of LJ site parameter vectors match
  if (epsilonMap.value.size() != sigmaMap.value.size()) {
    throw std::runtime_error(
        "MDFlexConfig::initializeParticlePropertiesLibrary(): Number of LJ site-level properties differ! Potentially "
        "missing epsilon or sigma for some LJ sites.");
  }
  // initialize LJ parameters
  for (auto [siteTypeId, epsilon] : epsilonMap.value) {
    _particlePropertiesLibrary->addLJParametersToSite(siteTypeId, epsilon, sigmaMap.value.at(siteTypeId));
  }
  // initialize AT parameters
  for (auto [siteTypeId, nu] : nuMap.value) {
    _particlePropertiesLibrary->addATParametersToSite(siteTypeId, nu);
  }

  // if doing Multi-site MD simulation, also check molecule level vectors match and initialize at molecular level
#if MD_FLEXIBLE_MODE == MULTISITE
  // check size of molecular level vectors match
  if (molToSiteIdMap.size() != molToSitePosMap.size()) {
    throw std::runtime_error("Number of molecular-level properties differ!");
  }

  // initialize at molecular level
  for (auto [molTypeId, siteTypeIds] : molToSiteIdMap) {
    _particlePropertiesLibrary->addMolType(molTypeId, siteTypeIds, molToSitePosMap.at(molTypeId),
                                           momentOfInertiaMap.at(molTypeId));
  }
#endif

  _particlePropertiesLibrary->calculateMixingCoefficients();
}

void MDFlexConfig::initializeObjects() {
  // @todo: the object generators should only generate particles relevant for the current rank's domain
  // All objects are generated on rank one. Particles are then distributed upon insertion.
  int myRank{};
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &myRank);
  if (myRank == 0) {
    for (const auto &object : cubeGridObjects) {
      object.generate(particles);
    }
    for (const auto &object : cubeGaussObjects) {
      object.generate(particles);
    }
    for (const auto &object : cubeUniformObjects) {
      object.generate(particles);
    }
    for (const auto &object : sphereObjects) {
      object.generate(particles);
    }
    for (const auto &object : cubeClosestPackedObjects) {
      object.generate(particles);
    }
  }
}

void MDFlexConfig::loadParticlesFromCheckpoint(const size_t &rank, const size_t &numRanks) {
  const std::string &filename = checkpointfile.value;

  std::ifstream inputStream(filename);
  if (not inputStream.is_open()) {
    throw std::runtime_error("Could not open checkpoint file: " + filename);
    return;
  }

  // If the checkpoint was written by a MPI parallel simulation the checkpoint is a pvtu file,
  // that points to several vtu files, called pieces.
  size_t numCheckpointPieces{getNumPiecesInCheckpoint(filename)};
  // If the number of currently used ranks matches the number of pieces in the checkpoint,
  // we assume that the domain layout matches and each rank only loads their own piece.
  // If it doesn't match the overlap will be redistributed when particles are inserted into the AutoPas containers.
  if (numRanks == numCheckpointPieces) {
    loadParticlesFromRankRecord(filename, rank, particles);
  } else if (rank == 0) {
    // Otherwise, rank 0 loads everything and distributes it later.
    for (size_t i = 0; i < numCheckpointPieces; ++i) {
      loadParticlesFromRankRecord(filename, i, particles);
    }
  }
}
