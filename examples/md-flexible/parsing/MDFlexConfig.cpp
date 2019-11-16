/**
 * @file MDFlexConfig.cpp
 * @author F. Gratl
 * @date 10/18/19
 */

#include "MDFlexConfig.h"
#include "autopas/utils/StringUtils.h"

std::string MDFlexConfig::to_string() const {
  using namespace std;
  ostringstream os;

  os << setw(valueOffset) << left << containerOptionsStr << ":  "
     << autopas::utils::ArrayUtils::to_string(containerOptions) << endl;

  // if verlet lists are in the container options print verlet config data
  if (autopas::utils::ArrayUtils::to_string(containerOptions).find("erlet") != std::string::npos) {
    os << setw(valueOffset) << left << verletRebuildFrequencyStr << ":  " << verletRebuildFrequency << endl;

    os << setw(valueOffset) << left << verletSkinRadiusStr << ":  " << verletSkinRadius << endl;
  }

  if (containerOptions.size() > 1 or traversalOptions.size() > 1 or dataLayoutOptions.size() > 1) {
    os << setw(valueOffset) << left << selectorStrategyStr << ":  " << selectorStrategy << endl;
  }

  os << setw(valueOffset) << left << dataLayoutOptionsStr << ":  "
     << autopas::utils::ArrayUtils::to_string(dataLayoutOptions) << endl;
  os << setw(valueOffset) << left << traversalOptionsStr << ":  "
     << autopas::utils::ArrayUtils::to_string(traversalOptions) << endl;
  os << setw(valueOffset) << left << tuningStrategyOptionsStr << ":  " << tuningStrategyOption << endl;
  os << setw(valueOffset) << left << tuningIntervalStr << ":  " << tuningInterval << endl;
  os << setw(valueOffset) << left << tuningSamplesStr << ":  " << tuningSamples << endl;
  os << setw(valueOffset) << left << tuningMaxEvidenceStr << ":  " << tuningMaxEvidence << endl;
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
  os << setw(valueOffset) << left << cellSizeFactorsStr << ":  "
     << "\"" << *cellSizeFactors << "\"" << endl;
  os << setw(valueOffset) << left << deltaTStr << ":  " << deltaT << endl;
  os << setw(valueOffset) << left << iterationsStr << ":  " << iterations << endl;
  os << setw(valueOffset) << left << boolalpha << periodicStr << ":  " << periodic << endl << endl;

  os << setw(valueOffset) << left << "Objects:" << endl;

  auto printObjectCollection = [](auto objectCollection, auto name, auto &os) {
    int objectId = 0;
    for (auto object : objectCollection) {
      os << "  " << name << ":" << endl;
      os << "    " << objectId << ":  " << endl;
      auto objectStr = object.to_string();
      // indent all lines of object
      objectStr = std::regex_replace(objectStr, std::regex("(^|\n)(.)"), "$1      $2");
      os << objectStr << endl;
      objectId++;
    }
  };

  printObjectCollection(cubeGridObjects, cubeGridObjectsStr, os);
  printObjectCollection(cubeGaussObjects, cubeGaussObjectsStr, os);
  printObjectCollection(cubeUniformObjects, cubeUniformObjectsStr, os);
  printObjectCollection(sphereObjects, sphereObjectsStr, os);

  if (useThermostat) {
    os << setw(valueOffset) << left << boolalpha << thermostatStr << endl;
    os << setw(valueOffset) << left << initTemperatureStr << ":  " << initTemperature << endl;
    os << setw(valueOffset) << left << targetTemperatureStr << ":  " << targetTemperature << endl;
    os << setw(valueOffset) << left << deltaTempStr << ":  " << deltaTemp << endl;
    os << setw(valueOffset) << left << thermostatIntervalStr << ":  " << thermostatInterval << endl;
    os << setw(valueOffset) << left << useCurrentTempForBrownianMotionStr << ":  " << useCurrentTempForBrownianMotion
       << endl;
  }

  return os.str();
}

void MDFlexConfig::calcSimulationBox() {
  const double interactionLength = cutoff + verletSkinRadius;
  // also account for explicitly set box sizes
  std::vector<double> xMins = {boxMin[0]};
  std::vector<double> yMins = {boxMin[1]};
  std::vector<double> zMins = {boxMin[2]};
  std::vector<double> xMaxs = {boxMax[0]};
  std::vector<double> yMaxs = {boxMax[1]};
  std::vector<double> zMaxs = {boxMax[2]};
  for (auto &c : cubeGridObjects) {
    xMins.emplace_back(c.getBoxMin()[0]);
    yMins.emplace_back(c.getBoxMin()[1]);
    zMins.emplace_back(c.getBoxMin()[2]);
    xMaxs.emplace_back(c.getBoxMax()[0]);
    yMaxs.emplace_back(c.getBoxMax()[1]);
    zMaxs.emplace_back(c.getBoxMax()[2]);
  }
  for (auto &c : cubeGaussObjects) {
    xMins.emplace_back(c.getBoxMin()[0]);
    yMins.emplace_back(c.getBoxMin()[1]);
    zMins.emplace_back(c.getBoxMin()[2]);
    xMaxs.emplace_back(c.getBoxMax()[0]);
    yMaxs.emplace_back(c.getBoxMax()[1]);
    zMaxs.emplace_back(c.getBoxMax()[2]);
  }
  for (auto &c : cubeUniformObjects) {
    xMins.emplace_back(c.getBoxMin()[0]);
    yMins.emplace_back(c.getBoxMin()[1]);
    zMins.emplace_back(c.getBoxMin()[2]);
    xMaxs.emplace_back(c.getBoxMax()[0]);
    yMaxs.emplace_back(c.getBoxMax()[1]);
    zMaxs.emplace_back(c.getBoxMax()[2]);
  }
  for (auto &c : sphereObjects) {
    xMins.emplace_back(c.getBoxMin()[0]);
    yMins.emplace_back(c.getBoxMin()[1]);
    zMins.emplace_back(c.getBoxMin()[2]);
    xMaxs.emplace_back(c.getBoxMax()[0]);
    yMaxs.emplace_back(c.getBoxMax()[1]);
    zMaxs.emplace_back(c.getBoxMax()[2]);
  }
  if (not xMins.empty()) {
    boxMin = {*std::min_element(xMins.begin(), xMins.end()), *std::min_element(yMins.begin(), yMins.end()),
              *std::min_element(zMins.begin(), zMins.end())};
    boxMax = {*std::max_element(xMaxs.begin(), xMaxs.end()), *std::max_element(yMaxs.begin(), yMaxs.end()),
              *std::max_element(zMaxs.begin(), zMaxs.end())};
  }
  // needed for 2D Simulation, that BoxLength >= interactionLength for all Dimensions
  for (int i = 0; i < 3; i++) {
    // pad domain such that periodic boundaries can work.
    boxMin[i] -= particleSpacing / 2;
    boxMax[i] += particleSpacing / 2;

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
      throw std::runtime_error("Wrong Particle initializaition: using same typeId for different properties");
    }
  } else {
    epsilonMap.emplace(typeId, epsilon);
    sigmaMap.emplace(typeId, sigma);
    massMap.emplace(typeId, mass);
  }
}