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
      "MeanKineticEnergyZ", "MeanRotationalEnergyX", "MeanRotationalEnergyY", "MeanRotationalEnergyZ",
      "MeanNormalVectorXAbs",  "MeanNormalVectorYAbs",     "MeanNormalVectorZAbs",     "VarNormalVectorXAbs",
      "VarNormalVectorYAbs",   "VarNormalVectorZAbs"};
  /**
  const std::vector<std::string> columnNames = {
      "Iteration", "TorqueIX",   "TorqueIY",   "TorqueIZ",     "AngularVelIX", "AngularVelIY", "AngularVelIZ",
      "TorqueJX",  "TorqueJY",   "TorqueJZ",   "AngularVelJX", "AngularVelJY", "AngularVelJZ", "ForceIX",
      "ForceIY",   "ForceIZ",    "VelocityIX", "VelocityIY",   "VelocityIZ",   "ForceJX",      "ForceJY",
      "ForceJZ",   "VelocityJX", "VelocityJY", "VelocityJZ", "Overlap", "DistanceBetweenCenter  s", "NotNeeded"};
      **/
  generateOutputFile(columnNames);
}

void StatisticsCalculator::recordStatistics(size_t currentIteration, const double globalForceZ,
                                            const autopas::AutoPas<ParticleType> &autoPasContainer,
                                            const ParticlePropertiesLibraryType &particlePropertiesLib) {
  const auto energyStatistics =
      calculateMeanPotentialKineticRotationalEnergy(autoPasContainer, globalForceZ, particlePropertiesLib);
  const auto meanQ = calculateMeanAndVarNormalVectorAbs(autoPasContainer, particlePropertiesLib, 0L);
  /**
  const auto statisticsI = calculateTorquesAndAngularVel(autoPasContainer, 1L);
  const auto statisticsJ = calculateTorquesAndAngularVel(autoPasContainer, 0L);
  const auto statisticsIForceVel = calculateForceAndVelocity(autoPasContainer, 1L);
  const auto statisticsJForceVel = calculateForceAndVelocity(autoPasContainer, 0L);
  const auto statisticsDistanceOverlap = calculateOverlapDistForceMagSum(autoPasContainer, particlePropertiesLib);

  auto combinedStatistics = std::tuple_cat(statisticsI, statisticsJ, statisticsIForceVel, statisticsJForceVel,
  statisticsDistanceOverlap);
   **/
  auto combinedStatistics = std::tuple_cat(energyStatistics, meanQ);
  StatisticsCalculator::writeRow(currentIteration, combinedStatistics);
}

std::tuple<double, double, double, double, double, double> StatisticsCalculator::calculateMeanAndVarNormalVectorAbs(
    const autopas::AutoPas<ParticleType> &autoPasContainer, const ParticlePropertiesLibraryType &particlePropertiesLib,
    size_t typeId) {
  double nXAbsSum = 0.;
  double nYAbsSum = 0.;
  double nZAbsSum = 0.;
  size_t particleCount = 0;

  std::vector<std::array<double, 3>> normalVectorAbsList = {};

  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    if (particle->getTypeId() == typeId) {
      const std::array<double, 3> CoM = particle->getR();
      const std::vector<std::array<double, 3>> unrotatedSitePositions =
          particlePropertiesLib.getSitePositions(particle->getTypeId());
      const auto rotatedSitePositions =
          autopas::utils::quaternion::rotateVectorOfPositions(particle->getQuaternion(), unrotatedSitePositions);
      const std::array<double, 3> CoMToSite0 = autopas::utils::ArrayMath::sub(rotatedSitePositions[0], CoM);
      const std::array<double, 3> CoMToSite1 = autopas::utils::ArrayMath::sub(rotatedSitePositions[1], CoM);
      std::array<double, 3> normalVector = autopas::utils::ArrayMath::cross(CoMToSite0, CoMToSite1);
      const double normalVectorNorm = autopas::utils::ArrayMath::L2Norm(normalVector);
      normalVector = autopas::utils::ArrayMath::divScalar(normalVector, (normalVectorNorm + 1e-10));
      const std::array<double, 3> normalVectorAbs = std::array<double, 3>{std::abs(normalVector[0]),
                                                                          std::abs(normalVector[1]),
                                                                          std::abs(normalVector[2])};

      normalVectorAbsList.emplace_back(normalVectorAbs);

      nXAbsSum += normalVectorAbs[0];
      nYAbsSum += normalVectorAbs[1];
      nZAbsSum += normalVectorAbs[2];
      particleCount++;
    }
  }

  const double meanNXAbs = nXAbsSum / particleCount;
  const double meanNYAbs = nYAbsSum / particleCount;
  const double meanNZAbs = nZAbsSum / particleCount;

  if (particleCount == 0) {
    return std::make_tuple(0, 0, 0, 0, 0, 0);
  }

  double varNXAbs = 0.;
  double varNYAbs = 0.;
  double varNZAbs = 0.;

  for (const auto&normalVectorAbs : normalVectorAbsList) {
    varNXAbs += (normalVectorAbs[0] - meanNXAbs) * (normalVectorAbs[0] - meanNXAbs);
    varNYAbs += (normalVectorAbs[1] - meanNYAbs) * (normalVectorAbs[1] - meanNYAbs);
    varNZAbs += (normalVectorAbs[2] - meanNZAbs) * (normalVectorAbs[2] - meanNZAbs);
  }

  varNXAbs = varNXAbs / particleCount;
  varNYAbs = varNYAbs / particleCount;
  varNZAbs = varNZAbs / particleCount;

  return std::make_tuple(meanNXAbs, meanNYAbs, meanNZAbs, varNXAbs, varNYAbs, varNZAbs);
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
    if (particle->getTypeId() != 0L) {
      continue;
    }
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