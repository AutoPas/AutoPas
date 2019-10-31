/**
 * @file MDFlexConfig.cpp
 * @author F. Gratl
 * @date 10/18/19
 */

#include "MDFlexConfig.h"
#include "autopas/utils/StringUtils.h"

void MDFlexConfig::print() {
  using namespace std;
  cout << setw(valueOffset) << left << "Container"
       << ":  " << autopas::utils::StringUtils::iterableToString(containerOptions) << endl;

  // if verlet lists are in the container options print verlet config data
  if (autopas::utils::StringUtils::iterableToString(containerOptions).find("erlet") != std::string::npos) {
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
       << ":  " << autopas::utils::StringUtils::iterableToString(dataLayoutOptions) << endl;
  cout << setw(valueOffset) << left << "Allowed traversals"
       << ":  " << autopas::utils::StringUtils::iterableToString(traversalOptions) << endl;
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
    case FunctorOption::lj12_6_Globals: {
      cout << "Lennard-Jones (12-6) with globals" << endl;
      break;
    }
  }
  cout << setw(valueOffset) << left << "Newton3"
       << ":  " << autopas::utils::StringUtils::iterableToString(newton3Options) << endl;

  cout << setw(valueOffset) << left << "Cutoff radius"
       << ":  " << cutoff << endl;
  cout << setw(valueOffset) << left << "boxMin"
       << ":  " << autopas::ArrayUtils::to_string(boxMin) << endl;
  cout << setw(valueOffset) << left << "boxMax"
       << ":  " << autopas::ArrayUtils::to_string(boxMax) << endl;
  cout << setw(valueOffset) << left << "Cell size factor"
       << ":  " << static_cast<std::string>(*cellSizeFactors) << endl;
  cout << setw(valueOffset) << left << "deltaT"
       << ":  " << deltaT << endl;
  cout << setw(valueOffset) << left << "Iterations"  // iterations * deltaT = time_end;
       << ":  " << iterations << endl;
  cout << setw(valueOffset) << left << "periodic boundaries"
       << ":  " << periodic << endl
       << endl;

  cout << setw(valueOffset) << left << "Object Generation:" << endl;
  int objectId = 1;
  for (auto c : cubeGridObjects) {
    cout << "-Cube Grid Nr " << objectId << ":  " << endl;
    cout << c << endl;
    objectId++;
  }
  objectId = 1;
  for (auto c : cubeGaussObjects) {
    cout << "-Cube Gauss Nr" << objectId << ":  " << endl;
    cout << c << endl;
    objectId++;
  }
  objectId = 1;
  for (auto c : cubeUniformObjects) {
    cout << "-Cube Uniform Nr " << objectId << ":  " << endl;
    cout << c << endl;
    objectId++;
  }
  objectId = 1;
  for (auto c : sphereObjects) {
    cout << "-Sphere Nr " << objectId << ":  " << endl;
    cout << c << endl;
    objectId++;
  }
  if (useThermostat) {
    cout << setw(valueOffset) << left << "Thermostat:" << endl;
    cout << setw(valueOffset) << left << "initializing velocites"
         << ":  " << initTemperature << endl;
    //@todo print usage of eather maxwellB or BM during initialization(after adding parsing options)
    cout << setw(valueOffset) << left << "initial Temperature"
         << ":  " << initTemperature << endl;
    cout << setw(valueOffset) << left << "number of TimeSteps"
         << ":  " << thermostatInterval << endl;
    cout << setw(valueOffset) << left << "target Temperature"
         << ":  " << targetTemperature << endl;
    cout << setw(valueOffset) << left << "deltaTemp"
         << ":  " << deltaTemp << endl;
  }
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