/**
 * @file SimilarityFunctions.h
 * @author J. Kroll
 * @date 01.06.2020
 */

#pragma once

#include "ThreeDimensionalMapping.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/utils/WrapMPI.h"

namespace autopas::utils {

/**
 * Calculates homogeneity and max density of given AutoPas simulation.
 * Both values are computed at once to avoid iterating over the same space twice.
 * homogeneity > 0.0, normally < 1.0, but for extreme scenarios > 1.0
 * maxDensity > 0.0, normally < 3.0, but for extreme scenarios >> 3.0
 *
 * @tparam Particle
 * @param container container of current simulation
 * @return {homogeneity, maxDensity}
 */
template <class Container>
std::pair<double, double> calculateHomogeneityAndMaxDensity(const Container &container,
                                                            const std::array<double, 3> startCorner,
                                                            const std::array<double, 3> endCorner) {
  unsigned int numberOfParticles = static_cast<unsigned int>(container->getNumberOfParticles());
  autopas::AutoPas_MPI_Allreduce(&numberOfParticles, &numberOfParticles, 1, AUTOPAS_MPI_UNSIGNED_INT, AUTOPAS_MPI_SUM,
                                 AUTOPAS_MPI_COMM_WORLD);
  // approximately the resolution we want to get.
  size_t numberOfCells = ceil(numberOfParticles / 10.);
  std::array<double, 3> domainSizePerDimension = {};
  for (int i = 0; i < 3; ++i) {
    domainSizePerDimension[i] = endCorner[i] - startCorner[i];
  }

  // get cellLength which is equal in each direction, derived from the domainsize and the requested number of cells
  const double volume = domainSizePerDimension[0] * domainSizePerDimension[1] * domainSizePerDimension[2];
  const double cellVolume = volume / numberOfCells;
  const double cellLength = cbrt(cellVolume);

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
    const std::array<double, 3> particleLocation = particleItr->getR();
    std::array<size_t, 3> index = {};
    for (size_t i = 0; i < particleLocation.size(); i++) {
      index[i] = particleLocation[i] / cellLength;
    }
    const size_t cellIndex = autopas::utils::ThreeDimensionalMapping::threeToOneD(index, cellsPerDimension);
    particlesPerCell[cellIndex] += 1;
    // calculate the size of the current cell
    allVolumes[cellIndex] = 1;
    for (size_t i = 0; i < cellsPerDimension.size(); ++i) {
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
  for (size_t i = 0; i < particlesPerCell.size(); i++) {
    densityPerCell[i] =
        (allVolumes[i] == 0) ? 0 : (particlesPerCell[i] / allVolumes[i]);  // make sure there is no division of zero
    if (densityPerCell[i] > maxDensity) {
      maxDensity = densityPerCell[i];
    }
  }

  if (maxDensity < 0.0)
    throw std::runtime_error("maxDensity can never be smaller than 0.0, but is:" + std::to_string(maxDensity));

  // get mean and reserve variable for densityVariance
  const double densityMean = numberOfParticles / volume;
  double densityVariance = 0.0;

  // calculate densityVariance
  for (size_t r = 0; r < densityPerCell.size(); ++r) {
    double distance = densityPerCell[r] - densityMean;
    densityVariance += (distance * distance / densityPerCell.size());
  }

  // finally calculate standard deviation
  const double homogeneity = sqrt(densityVariance);
  // normally between 0.0 and 1.5
  if (homogeneity < 0.0)
    throw std::runtime_error("homogeneity can never be smaller than 0.0, but is:" + std::to_string(homogeneity));
  return {homogeneity, maxDensity};
}

}  // namespace autopas::utils
