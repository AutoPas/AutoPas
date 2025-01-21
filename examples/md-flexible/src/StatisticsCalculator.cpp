/**
 * @file StatisticsCalculator.cpp
 * @author Joon Kim
 * @date 21.11.2024
 */

#include "StatisticsCalculator.h"

#include <cstddef>
#include <ios>
#include <limits>
#include <string>
#include <utility>

StatisticsCalculator::StatisticsCalculator(std::string sessionName, const std::string &outputFolder)
    : _sessionName(std::move(sessionName)) {
  tryCreateStatisticsFolders(_sessionName, outputFolder);

  std::vector<std::string> columnNames = {
      "Iteration",          "MeanPotentialEnergyZ",  "MeanKineticEnergyX",    "MeanKineticEnergyY",
      "MeanKineticEnergyZ", "MeanRotationalEnergyX", "MeanRotationalEnergyY", "MeanRotationalEnergyZ", "FlowVelocityY", "CrossSectionArea", "VolumetricFlowRate", "numConsideredGasParticles"};
  /**
  const std::vector<std::string> columnNames = {
      "Iteration", "TorqueIX",   "TorqueIY",   "TorqueIZ",     "AngularVelIX", "AngularVelIY", "AngularVelIZ",
      "TorqueJX",  "TorqueJY",   "TorqueJZ",   "AngularVelJX", "AngularVelJY", "AngularVelJZ", "ForceIX",
      "ForceIY",   "ForceIZ",    "VelocityIX", "VelocityIY",   "VelocityIZ",   "ForceJX",      "ForceJY",
      "ForceJZ",   "VelocityJX", "VelocityJY", "VelocityJZ", "Overlap", "DistanceBetweenCenters", "NotNeeded"};
      **/
  generateOutputFile(columnNames);
}

void StatisticsCalculator::recordStatistics(size_t currentIteration, const double globalForceZ,
                                            const autopas::AutoPas<ParticleType> &autoPasContainer,
                                            const ParticlePropertiesLibraryType &particlePropertiesLib) {
  const auto energyStatistics =
      calculateMeanPotentialKineticRotationalEnergy(autoPasContainer, globalForceZ, particlePropertiesLib);
  const auto flowRateStatistics = calculateVolumetricFlowRate(autoPasContainer, particlePropertiesLib);
  /**
  const auto statisticsI = calculateTorquesAndAngularVel(autoPasContainer, 1L);
  const auto statisticsJ = calculateTorquesAndAngularVel(autoPasContainer, 0L);
  const auto statisticsIForceVel = calculateForceAndVelocity(autoPasContainer, 1L);
  const auto statisticsJForceVel = calculateForceAndVelocity(autoPasContainer, 0L);
  const auto statisticsDistanceOverlap = calculateOverlapDistForceMagSum(autoPasContainer, particlePropertiesLib);

  auto combinedStatistics = std::tuple_cat(statisticsI, statisticsJ, statisticsIForceVel, statisticsJForceVel,
  statisticsDistanceOverlap);
   **/
  auto combinedStatistics = std::tuple_cat(energyStatistics, flowRateStatistics);
  StatisticsCalculator::writeRow(StatisticsCalculator::outputFile, currentIteration, combinedStatistics);
  if (currentIteration % 10000 == 0) {
    const std::vector<std::tuple<size_t, double>> binIndex_to_rdf = calculateRDF(autoPasContainer, particlePropertiesLib);
    for (const auto &pair : binIndex_to_rdf) {
      writeRow(outputFile_rdf, currentIteration, pair);
    }
  }
}

std::tuple<double, double, double, double, double, double> StatisticsCalculator::calculateTorquesAndAngularVel(
    const autopas::AutoPas<ParticleType> &autoPasContainer, const size_t typeId) {
  double torqueIX = 0.;
  double torqueIY = 0.;
  double torqueIZ = 0.;
  // double torqueJX = 0.;
  // double torqueJY = 0.;
  // double torqueJZ = 0.;

  double angularVelIX = 0.;
  double angularVelIY = 0.;
  double angularVelIZ = 0.;
  // double angularVelJX = 0.;
  // double angularVelJY = 0.;
  // double angularVelJZ = 0.;

  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    if (particle->getTypeId() == typeId) {  // particle i
      auto torqueI = particle->getTorque();
      torqueIX = torqueI[0];
      torqueIY = torqueI[1];
      torqueIZ = torqueI[2];
      auto angularVelI = particle->getAngularVel();
      angularVelIX = angularVelI[0];
      angularVelIY = angularVelI[1];
      angularVelIZ = angularVelI[2];
    }
    //} else {  // particle j
    //  auto torqueJ = particle->getTorque();
    //  torqueJX = torqueJ[0];
    //  torqueJY = torqueJ[1];
    //  torqueJZ = torqueJ[2];
    //  auto angularVelJ = particle->getAngularVel();
    //  angularVelJX = angularVelJ[0];
    //  angularVelJY = angularVelJ[1];
    //  angularVelJZ = angularVelJ[2];
    //}
  }

  return std::make_tuple(torqueIX, torqueIY, torqueIZ, angularVelIX, angularVelIY, angularVelIZ);
}

std::tuple<double, double, double, double, double, double> StatisticsCalculator::calculateForceAndVelocity(
    const autopas::AutoPas<ParticleType> &autoPasContainer, const size_t typeId) {
  double forceX = 0.;
  double forceY = 0.;
  double forceZ = 0.;
  double velocityX = 0.;
  double velocityY = 0.;
  double velocityZ = 0.;

  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    if (particle->getTypeId() == typeId) {  // particle i
      auto force = particle->getF();
      forceX = force[0];
      forceY = force[1];
      forceZ = force[2];
      auto velocity = particle->getV();
      velocityX = velocity[0];
      velocityY = velocity[1];
      velocityZ = velocity[2];
    }
  }
  return std::make_tuple(forceX, forceY, forceZ, velocityX, velocityY, velocityZ);
}

std::tuple<double, double, double, double, double, double, double>
StatisticsCalculator::calculateMeanPotentialKineticRotationalEnergy(
    const autopas::AutoPas<ParticleType> &autoPasContainer, const double globalForceZ,
    const ParticlePropertiesLibraryType &particlePropertiesLib) {
  using namespace autopas::utils::ArrayMath::literals;

  size_t particleCount = 0;
  double meanPotentialEnergy = 0.;
  std::array<double, 3> meanKineticEnergy = {0., 0., 0.};
  std::array<double, 3> meanRotationalEnergy = {0., 0., 0.};

  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    if (particle->getTypeId() != 0) continue;

    const double mass = particlePropertiesLib.getSiteMass(particle->getTypeId());
    const double radius = particlePropertiesLib.getRadius(particle->getTypeId());
    const std::array<double, 3> v = particle->getV();
    const std::array<double, 3> w = particle->getAngularVel();
    const std::array<double, 3> r = particle->getR();
    const double momentOfInertia = 0.4 * mass * radius * radius;

    meanPotentialEnergy += (mass * (-globalForceZ) * r[2]);
    meanKineticEnergy += (v * v * 0.5 * mass);
    meanRotationalEnergy += (w * w * 0.5 * momentOfInertia);

    ++particleCount;
  }

  meanPotentialEnergy = meanPotentialEnergy * (1. / particleCount);
  meanKineticEnergy = meanKineticEnergy * (1. / particleCount);
  meanRotationalEnergy = meanRotationalEnergy * (1. / particleCount);

  return std::make_tuple(meanPotentialEnergy, meanKineticEnergy[0], meanKineticEnergy[1], meanKineticEnergy[2],
                         meanRotationalEnergy[0], meanRotationalEnergy[1], meanRotationalEnergy[2]);
}

std::tuple<double, double, double> StatisticsCalculator::calculateOverlapDistForceMagSum(
    const autopas::AutoPas<ParticleType> &autoPasContainer,
    const ParticlePropertiesLibraryType &particlePropertiesLib) {
  using namespace autopas::utils::ArrayMath::literals;
  using namespace autopas::utils::ArrayMath;

  double overlapSum = 0.;
  double distSum = 0.;
  double forceMagSum = 0.;

  for (auto i = autoPasContainer.begin(autopas::IteratorBehavior::owned); i.isValid(); ++i) {
    for (auto j = autoPasContainer.begin(autopas::IteratorBehavior::owned); j.isValid(); ++j) {
      if (i->getID() == j->getID()) {
        continue;
      }

      const std::array<double, 3> x_i = i->getR();
      const std::array<double, 3> x_j = j->getR();
      const std::array<double, 3> displacement = x_i - x_j;
      const double dist = L2Norm(displacement);

      const double radius_i = particlePropertiesLib.getRadius(i->getTypeId());
      const double radius_j = particlePropertiesLib.getRadius(j->getTypeId());
      const double overlap = radius_i + radius_j - dist;

      const std::array<double, 3> force_i = i->getF();
      const std::array<double, 3> force_j = j->getF();
      const double forceMag = L2Norm(force_i) + L2Norm(force_j);

      overlapSum += (overlap > 0 ? overlap : 0);
      distSum += dist;
      forceMagSum += forceMag;
    }
  }

  overlapSum /= 2.;
  distSum /= 2.;
  forceMagSum /= 2.;

  return std::make_tuple(overlapSum, distSum, forceMagSum);
}

std::tuple<double, double, double, size_t> StatisticsCalculator::calculateVolumetricFlowRate(
    const autopas::AutoPas<ParticleType> &autoPasContainer,
    const ParticlePropertiesLibraryType &particlePropertiesLib) {
  using namespace autopas::utils::ArrayMath::literals;
  using namespace autopas::utils::ArrayMath;

  const std::array<double, 3> currentDomain = autoPasContainer.getBoxMax();
  const double crossSectionArea = currentDomain[0] * currentDomain[2];

  double flowVelYSum = 0.;
  size_t gasParticleCount = 0;

  for (auto i = autoPasContainer.begin(autopas::IteratorBehavior::owned); i.isValid(); ++i) {
    if (i->getTypeId() != 1) {  // Only consider gas particles
      continue;
    }
    if (i->getR()[1] < (currentDomain[1] / 5.)) {
      flowVelYSum += i->getV()[1];
      gasParticleCount++;
    }
  }

  const double preventDivisionByZero = 1e-6;
  const double flowVelocityY = flowVelYSum / (static_cast<double>(gasParticleCount) + preventDivisionByZero);
  const double volumetricFlowRate = flowVelocityY * crossSectionArea;

  return std::make_tuple(flowVelocityY, crossSectionArea, volumetricFlowRate, gasParticleCount);
}

//---------------------------------------------Helper Methods-----------------------------------------------------

void StatisticsCalculator::generateOutputFile(const std::vector<std::string> &columnNames) {
  std::ostringstream filename;
  filename << _statisticsFolderPath << _sessionName << "_statistics.csv";

  outputFile.open(filename.str(), std::ios::out);

  if (outputFile.is_open()) {
    for (size_t i = 0; i < columnNames.size(); ++i) {
      outputFile << columnNames[i];
      if (i < columnNames.size() - 1) {
        outputFile << ",";
      }
    }
    outputFile << "\n";
  } else {
    throw std::runtime_error("StatisticsCalculator::generateOutputFile(): Could not open file " + filename.str());
  }
  // ---------------------------------------------------------------------------------------------------------
  std::ostringstream filename_rdf;
  filename_rdf << _statisticsFolderPath << _sessionName << "_statistics_rdf.csv";

  outputFile_rdf.open(filename_rdf.str(), std::ios::out);

  if (outputFile_rdf.is_open()) {
    outputFile_rdf << "Iteration, Distance, RDF\n";
  } else {
    throw std::runtime_error("StatisticsCalculator::generateOutputFile(): Could not open file " + filename_rdf.str());
  }
}

void StatisticsCalculator::tryCreateStatisticsFolders(const std::string &name, const std::string &location) {
  if (not checkFileExists(location)) {
    tryCreateFolder(location, "./");
  }

  _sessionFolderPath = location + "/" + name + "/";
  tryCreateFolder(name, location);

  _statisticsFolderPath = _sessionFolderPath + "statistics/";
  tryCreateFolder("statistics", _sessionFolderPath);
}

void StatisticsCalculator::tryCreateFolder(const std::string &name, const std::string &location) {
  try {
    // took reference of ParallelVtkWriter.cpp
    const auto newDirectoryPath{location + "/" + name};
    mkdir(newDirectoryPath.c_str(), 0777);
  } catch (const std::exception &ex) {
    throw std::runtime_error("StatisticsCalculator::tryCreateFolder(): The output location " + location +
                             " passed to StatisticsCalculator is invalid: " + ex.what());
  }
}
std::vector<std::tuple<size_t, double>> StatisticsCalculator::calculateRDF(
    const autopas::AutoPas<ParticleType> &autoPasContainer,
    const ParticlePropertiesLibraryType &particlePropertiesLib) {
  using namespace autopas::utils::ArrayMath::literals;
  using namespace autopas::utils::ArrayMath;

  std::map<size_t, size_t> binIndex_to_counts;
  const double distance_unit = 1.;
  std::vector<std::tuple<size_t, double>> binIndex_to_rdf;

  for (auto i = autoPasContainer.begin(autopas::IteratorBehavior::owned); i.isValid(); ++i) {
    if (i->getTypeId() != 0) {  // Only consider solid particles
      continue;
    }
    auto j = i;
    ++j;  // Start `j` from the next particle after `i`.
    for (; j.isValid(); ++j) {
      if (j->getTypeId() != 0) {  // Only consider solid particles
        continue;
      }

      // No need to consider periodically imaged positions as for grain particles, only reflective BD are active.
      const std::array<double, 3> r_i = i->getR();
      const std::array<double, 3> r_j = j->getR();
      const double dist = L2Norm(r_i - r_j);
      const size_t binindex = std::floor(dist / distance_unit);
      binIndex_to_counts[binindex]++;
    }  // End of `j` loop
  }  // End of 'i' loop

  for (auto &pair : binIndex_to_counts) {
    const double dist_interval_start_casted = static_cast<double>(pair.first);
    const double counts_casted = static_cast<double>(pair.second);
    const double localDensity_Denominator =
        (4. / 3.) * M_PI *
        (std::pow(dist_interval_start_casted + distance_unit, 3) - std::pow(dist_interval_start_casted, 3));
    const double localDensity = counts_casted / localDensity_Denominator;

      binIndex_to_rdf.emplace_back(pair.first, localDensity);
  }

  return binIndex_to_rdf;
}
