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
  /**
  std::vector<std::string> columnNames = {
      "Iteration", "MeanPotentialEnergyZ", "MeanKineticEnergyX", "MeanKineticEnergyY",
      "MeanKineticEnergyZ", "MeanRotationalEnergyX", "MeanRotationalEnergyY", "MeanRotationalEnergyZ"
  };
   **/
  const std::vector<std::string> columnNames = {"Iteration", "Overlap", "Dist", "ForceIX", "MinusRelVelDotNormalUnit"};
  generateOutputFile(columnNames);
}

void StatisticsCalculator::recordStatistics(size_t currentIteration, const double globalForceZ,
                                            const autopas::AutoPas<ParticleType> &autoPasContainer,
                                            const ParticlePropertiesLibraryType &particlePropertiesLib) {
  // const auto statistics = calculateMeanPotentialKineticRotationalEnergy(autoPasContainer, globalForceZ,
  // particlePropertiesLib); TODO: change
  const auto statistics = calculateOverlapDistForceRelVelNormal(autoPasContainer, particlePropertiesLib);
  StatisticsCalculator::writeRow(currentIteration, statistics);
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

std::tuple<double, double, double, double> StatisticsCalculator::calculateOverlapDistForceRelVelNormal(
    const autopas::AutoPas<ParticleType> &autoPasContainer,
    const ParticlePropertiesLibraryType &particlePropertiesLib) {
  using namespace autopas::utils::ArrayMath::literals;
  using namespace autopas::utils::ArrayMath;

  // Force calculation for the linear normal contact force model in 2 particle setting.

  demLib::GranularDEM particle_i;
  demLib::GranularDEM particle_j;

  for (auto i = autoPasContainer.begin(autopas::IteratorBehavior::owned); i.isValid(); ++i) {
    if (i->getID() == 1) {
      particle_i = *i;
    } else if (i->getID() == 0) {
      particle_j = *i;
    }
  }

  const std::array<double, 3> x_i = particle_i.getR();
  const std::array<double, 3> x_j = particle_j.getR();
  const std::array<double, 3> displacement = x_i - x_j;
  const double dist = L2Norm(displacement);
  const std::array<double, 3> normalUnit = displacement / dist;

  const double radius_i = particlePropertiesLib.getRadius(particle_i.getTypeId());
  const double radius_j = particlePropertiesLib.getRadius(particle_j.getTypeId());
  double overlap = radius_i + radius_j - dist;
  overlap = overlap > 0 ? overlap : 0;

  const std::array<double, 3> force_i = particle_i.getF();
  const double forceIX = force_i[0];

  const std::array<double, 3> v_i = particle_i.getV();
  const std::array<double, 3> v_j = particle_j.getV();
  const std::array<double, 3> relVel = v_i - v_j;
  const double relVelDotNormalUnit = dot(relVel, normalUnit);

  return std::make_tuple(overlap, dist, forceIX, -relVelDotNormalUnit);
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