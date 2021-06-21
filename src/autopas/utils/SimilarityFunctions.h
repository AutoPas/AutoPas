/**
 * @file SimilarityFunctions.h
 * @author J. Kroll
 * @date 01.06.2020
 */

#pragma once

#include "ThreeDimensionalMapping.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/molecularDynamics/MoleculeLJ.h"

namespace autopas::utils {

/**
 * calculates homogeneity and max density of given AutoPas simulation
 *
 * @tparam Particle
 * @param autopas
 * @return  [homogeneity, max_density]
 */
template<class Particle>
std::pair<double,double> calculateHomogeneityAndMaxDensity(const std::shared_ptr<autopas::ParticleContainerInterface<Particle>>& container) {
  size_t numberOfParticles = container->getNumParticles();
  // approximately the resolution we want to get.
  size_t numberOfCells = ceil(numberOfParticles / 10.);

  std::array<double, 3> startCorner = container->getBoxMin();
  std::array<double, 3> endCorner = container->getBoxMax();
  std::array<double, 3> domainSizePerDimension = {};
  for (int i = 0; i < 3; ++i) {
    domainSizePerDimension[i] = endCorner[i] - startCorner[i];
  }

  // get cellLength which is equal in each direction, derived from the domainsize and the requested number of cells
  double volume = domainSizePerDimension[0] * domainSizePerDimension[1] * domainSizePerDimension[2];
  double cellVolume = volume / numberOfCells;
  double cellLength = cbrt(cellVolume);

  // calculate the size of the boundary cells, which might be smaller then the other cells
  std::array<size_t, 3> cellsPerDimension = {};
  // size of the last cell layer per dimension. This cell might get truncated to fit in the domain.
  std::array<double, 3> outerCellSizePerDimension = {};
  for (int i = 0; i < 3; ++i) {
    outerCellSizePerDimension[i] =
        domainSizePerDimension[i] - (floor(domainSizePerDimension[i] / cellLength) * cellLength);
    cellsPerDimension[i] = ceil(domainSizePerDimension[i] / cellLength);
  }
  // Actual number of cells we end up with
  numberOfCells = cellsPerDimension[0] * cellsPerDimension[1] * cellsPerDimension[2];

  std::vector<size_t> particlesPerCell(numberOfCells, 0);
  std::vector<double> allVolumes(numberOfCells, 0);

  // add particles accordingly to their cell to get the amount of particles in each cell
  for (auto particleItr = container->begin(autopas::IteratorBehavior::owned); particleItr.isValid(); ++particleItr) {
    std::array<double, 3> particleLocation = particleItr->getR();
    std::array<size_t, 3> index = {};
    for (int i = 0; (size_t) i < particleLocation.size(); i++) {
      index[i] = particleLocation[i] / cellLength;
    }
    const size_t cellIndex = autopas::utils::ThreeDimensionalMapping::threeToOneD(index, cellsPerDimension);
    particlesPerCell[cellIndex] += 1;
    // calculate the size of the current cell
    allVolumes[cellIndex] = 1;
    for (int i = 0; (size_t) i < cellsPerDimension.size(); ++i) {
      // the last cell layer has a special size
      if (index[i] == cellsPerDimension[i] - 1) {
        allVolumes[cellIndex] *= outerCellSizePerDimension[i];
      } else {
        allVolumes[cellIndex] *= cellLength;
      }
    }
  }

  // calculate density for each cell
  double maxDensity{0.};
  std::vector<double> densityPerCell(numberOfCells, 0.0);
  for (int i = 0; (size_t) i < particlesPerCell.size(); i++) {
    densityPerCell[i] =
        (allVolumes[i] == 0) ? 0 : (particlesPerCell[i] / allVolumes[i]);  // make sure there is no division of zero
    if (densityPerCell[i] > maxDensity) {
      maxDensity = densityPerCell[i];
    }
  }

  // get mean and reserve variable for densityVariance
  double densityMean = numberOfParticles / volume;
  double densityVariance = 0.0;

  // calculate densityVariance
  for (int r = 0; (size_t) r < densityPerCell.size(); ++r) {
    double distance = densityPerCell[r] - densityMean;
    densityVariance += (distance * distance / densityPerCell.size());
  }

  // finally calculate standard deviation
  // return sqrt(densityVariance);
  double homogeneity = sqrt(densityVariance);
  return {homogeneity, maxDensity};
}

} // namespace autopas::utils
