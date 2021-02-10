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
