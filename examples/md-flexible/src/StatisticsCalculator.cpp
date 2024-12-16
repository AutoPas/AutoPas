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
  // const std::vector<std::string> columnNames = {"Iteration", "Overlap", "Dist", "ForceIX",
  // "MinusRelVelDotNormalUnit"};
  const std::vector<std::string> strainStressColumnNames = {"Iteration",  "DomainX",  "DomainY",      "DomainZ",
                                                            "Epsilon_yy", "Density",  "VolumeStrain", "Stress_xx",
                                                            "Stress_yy",  "Stress_zz"};
  generateOutputFile(strainStressColumnNames);
}

void StatisticsCalculator::recordStatistics(size_t currentIteration, const double globalForceZ,
                                            const autopas::AutoPas<ParticleType> &autoPasContainer,
                                            const ParticlePropertiesLibraryType &particlePropertiesLib,
                                            const double initialVolume, const double finalBoxMaxY,
                                            const double spring_stiffness) {
  // const auto statistics = calculateMeanPotentialKineticRotationalEnergy(autoPasContainer, globalForceZ,
  // particlePropertiesLib); TODO: change
  // const auto statistics = calculateOverlapDistForceRelVelNormal(autoPasContainer, particlePropertiesLib);
  const auto statistics = calculateStrainStressStatistics(autoPasContainer, particlePropertiesLib, initialVolume,
                                                          finalBoxMaxY, spring_stiffness);
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

std::array<double, 9> StatisticsCalculator::calculateStrainStressStatistics(
    const autopas::AutoPas<ParticleType> &autoPasContainer, const ParticlePropertiesLibraryType &particlePropertiesLib,
    const double initialVolume, const double finalBoxMaxY, const double spring_stiffness) {
  using namespace autopas::utils::ArrayMath::literals;
  using namespace autopas::utils::ArrayMath;
  // To calculate: DomainSize (x,y,z), Density, Volumetric Strain
  // DomainSize
  const std::array<double, 3> currentDomain = autoPasContainer.getBoxMax();

  // E_yy
  const double epsilon_yy = 1 - (currentDomain[1] / finalBoxMaxY);

  // Density
  const double currentVolume = currentDomain[0] * currentDomain[1] * currentDomain[2];
  double particleVolumeSum = 0.;
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    const double radius = particlePropertiesLib.getRadius(particle->getTypeId());
    particleVolumeSum += 4. / 3. * M_PI * radius * radius * radius;
  }
  const double density = particleVolumeSum / currentVolume;

  // Volumetric Strain
  const double volumetricStrain = (currentVolume - initialVolume) / initialVolume;

  // Stress
  std::array<double, 3> stress_tensor_diagonal = {0., 0., 0.};
  for (auto i = autoPasContainer.begin(autopas::IteratorBehavior::owned); i.isValid(); ++i) {
    // static stress due to velocity
    const double m_i = particlePropertiesLib.getSiteMass(i->getTypeId());
    const std::array<double, 3> v_i = i->getV();
    const std::array<double, 3> v_i_squared = v_i * v_i;
    stress_tensor_diagonal += (v_i_squared * m_i);

    auto j = i;
    ++j;  // Start `j` from the next particle after `i`.
    for (; j.isValid(); ++j) {
      // Process the pair (i, j)
      // Dynamic stress due to contacts
      const double radius_i = particlePropertiesLib.getRadius(i->getTypeId());  // Assume same radius for j
      const std::array<double, 3> x_i = i->getR();
      const std::array<double, 3> x_j = j->getR();
      const std::array<double, 3> displacement = x_i - x_j;
      const double dist = L2Norm(displacement);
      const double overlap = 2. * radius_i - dist;

      if (overlap <= 0) {
        continue;
      }  // No contact between particles `i` and `j`.

      const std::array<double, 3> normalUnit = displacement / dist;
      const std::array<double, 3> overlap_vectorized = normalUnit * (2. * radius_i) - displacement;

      const std::array<double, 3> contact_force = overlap_vectorized * -spring_stiffness;
      stress_tensor_diagonal += (contact_force * displacement);
    }  // End of `j` loop
  }  // End of `i` loop
  stress_tensor_diagonal = stress_tensor_diagonal / currentVolume;

  return {currentDomain[0],
          currentDomain[1],
          currentDomain[2],
          epsilon_yy,
          density,
          volumetricStrain,
          stress_tensor_diagonal[0],
          stress_tensor_diagonal[1],
          stress_tensor_diagonal[2]};
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
