/**
 * @file MDFlexConfig.cpp
 * @author F. Gratl
 * @date 18.10.2019
 */

#include "MDFlexConfig.h"

#include "autopas/utils/StringUtils.h"

std::string MDFlexConfig::to_string() const {
  using namespace std;
  ostringstream os;

  auto passedContainerOptionsStr = autopas::utils::ArrayUtils::to_string(containerOptions);
  os << setw(valueOffset) << left << containerOptionsStr << ":  " << passedContainerOptionsStr << endl;

  // since all containers are rebuilt only periodically print Verlet config always.
  os << setw(valueOffset) << left << verletRebuildFrequencyStr << ":  " << verletRebuildFrequency << endl;
  os << setw(valueOffset) << left << verletSkinRadiusStr << ":  " << verletSkinRadius << endl;
  if (passedContainerOptionsStr.find("luster") != std::string::npos) {
    os << setw(valueOffset) << left << verletClusterSizeStr << ":  " << verletClusterSize << endl;
  }

  if (containerOptions.size() > 1 or traversalOptions.size() > 1 or dataLayoutOptions.size() > 1) {
    os << setw(valueOffset) << left << selectorStrategyStr << ":  " << selectorStrategy << endl;
  }

  os << setw(valueOffset) << left << dataLayoutOptionsStr << ":  "
     << autopas::utils::ArrayUtils::to_string(dataLayoutOptions) << endl;
  os << setw(valueOffset) << left << traversalOptionsStr << ":  "
     << autopas::utils::ArrayUtils::to_string(traversalOptions) << endl;
  os << setw(valueOffset) << left << tuningStrategyOptionsStr << ":  " << tuningStrategyOption << endl;
  if (tuningStrategyOption == autopas::TuningStrategyOption::bayesianSearch) {
    os << setw(valueOffset) << left << acquisitionFunctionOptionStr << ":  " << acquisitionFunctionOption << endl;
  }
  os << setw(valueOffset) << left << tuningIntervalStr << ":  " << tuningInterval << endl;
  os << setw(valueOffset) << left << tuningSamplesStr << ":  " << tuningSamples << endl;
  os << setw(valueOffset) << left << tuningMaxEvidenceStr << ":  " << tuningMaxEvidence << endl;
  os << setw(valueOffset) << left << relativeOptimumRangeStr << ":  " << relativeOptimumRange << endl;
  os << setw(valueOffset) << left << maxTuningPhasesWithoutTestStr << ":  " << maxTuningPhasesWithoutTest << endl;
  os << setw(valueOffset) << left << functorOptionStr << ":  ";
  switch (functorOption) {
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
  os << setw(valueOffset) << left << newton3OptionsStr << ":  " << autopas::utils::ArrayUtils::to_string(newton3Options)
     << endl;

  os << setw(valueOffset) << left << cutoffStr << ":  " << cutoff << endl;
  os << setw(valueOffset) << left << boxMinStr << ":  " << autopas::utils::ArrayUtils::to_string(boxMin) << endl;
  os << setw(valueOffset) << left << boxMaxStr << ":  " << autopas::utils::ArrayUtils::to_string(boxMax) << endl;
  os << setw(valueOffset) << left << cellSizeFactorsStr << ":  " << *cellSizeFactors << endl;
  os << setw(valueOffset) << left << deltaTStr << ":  " << deltaT << endl;
  os << setw(valueOffset) << left << iterationsStr << ":  " << iterations << endl;
  os << setw(valueOffset) << left << boolalpha << periodicStr << ":  " << periodic << endl;

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

  if (useThermostat) {
    os << thermostatStr << ":" << endl;
    os << "  " << setw(valueOffset - 2) << left << initTemperatureStr << ":  " << initTemperature << endl;
    os << "  " << setw(valueOffset - 2) << left << targetTemperatureStr << ":  " << targetTemperature << endl;
    os << "  " << setw(valueOffset - 2) << left << deltaTempStr << ":  " << deltaTemp << endl;
    os << "  " << setw(valueOffset - 2) << left << thermostatIntervalStr << ":  " << thermostatInterval << endl;
    os << "  " << setw(valueOffset - 2) << left << addBrownianMotionStr << ":  " << addBrownianMotion << endl;
  }

  if (not vtkFileName.empty()) {
    os << setw(valueOffset) << left << vtkFileNameStr << ":  " << vtkFileName << endl;
    os << setw(valueOffset) << left << vtkWriteFrequencyStr << ":  " << vtkWriteFrequency << endl;
  }
  if (not checkpointfile.empty())
    os << setw(valueOffset) << left << checkpointfileStr << ":  " << checkpointfile << endl;

  os << setw(valueOffset) << dontMeasureFlopsStr << ":  " << (not dontMeasureFlops) << endl;
  os << setw(valueOffset) << dontCreateEndConfigStr << ":  " << (not dontCreateEndConfig) << endl;
  return os.str();
}

void MDFlexConfig::calcSimulationBox() {
  const double interactionLength = cutoff + verletSkinRadius;
  std::array<std::vector<double>, 3> mins;
  std::array<std::vector<double>, 3> maxs;

  // helper function
  auto emplaceObjectLimits = [&](const auto &objectCollection) {
    for (auto &object : objectCollection) {
      for (size_t i = 0; i < 3; ++i) {
        mins[i].emplace_back(object.getBoxMin()[i]);
        maxs[i].emplace_back(object.getBoxMax()[i]);
      }
    }
  };

  emplaceObjectLimits(cubeGaussObjects);
  emplaceObjectLimits(cubeGridObjects);
  emplaceObjectLimits(cubeUniformObjects);
  emplaceObjectLimits(sphereObjects);

  for (int i = 0; i < 3; i++) {
    // pad domain such that periodic boundaries can work.
    // This is necessary if the given min/max is not at least half the spacing away of the farthest object.
    if (not mins[0].empty()) {
      auto objectMin = *std::min_element(mins[i].begin(), mins[i].end());
      auto objectMax = *std::max_element(maxs[i].begin(), maxs[i].end());
      boxMin[i] = std::min(boxMin[i], objectMin - particleSpacing / 2);
      boxMax[i] = std::max(boxMax[i], objectMax + particleSpacing / 2);
    }

    // needed for 2D Simulation, that BoxLength >= interactionLength for all Dimensions
    if (boxMax[i] - boxMin[i] < interactionLength) {
      std::cout << "WARNING: Simulation box in dimension " << i
                << " is shorter than interaction length and will be increased." << std::endl;
      boxMin[i] -= interactionLength / 2;
      boxMax[i] += interactionLength / 2;
    }
  }
}

void MDFlexConfig::addParticleType(unsigned long typeId, double epsilon, double sigma, double mass) {
  // check if type id is already existing and if there no error in input
  if (epsilonMap.count(typeId) == 1) {
    // check if type is already added
    if (epsilonMap.at(typeId) == epsilon and sigmaMap.at(typeId) == sigma and massMap.at(typeId) == mass) {
      return;
    } else {  // wrong initialization:
      throw std::runtime_error("Wrong Particle initialization: using same typeId for different properties");
    }
  } else {
    epsilonMap.emplace(typeId, epsilon);
    sigmaMap.emplace(typeId, sigma);
    massMap.emplace(typeId, mass);
  }
}
