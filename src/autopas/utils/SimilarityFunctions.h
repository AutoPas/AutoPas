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
 * maxDensity > 0.0, normally < 3.0, but for extreme scenarios > 3.0
 *
 * These values are calculated by dividing the domain into bins of perfectly equal cuboid shapes.
 *
 * @note The bin shapes between different AutoPas container will not match if the dimensions of their domains don't match.
 * Bins with greater surface area are probably more likely to feature fluctuations in the density.
 *
 * Not a rigorous proof, but should give an idea:
 * - Assume homogeneous distribution.
 * - For every particle assume it could move in any direction with equal probability (not a terrible assumption).
 * - So every particle near a bin face has a chance to more into the neighboring bin.
 * - And this probability is independent of the bin (i.e. the bin shape)
 * - With a larger surface area, there are more particles with a probability of moving into the neighboring bin.
 * - On average, same density, but each individual bin has a greater probability of losing or gaining a lot of particles.
 *
 * This is probably still fine for MPI Parallelized Tuning purposes, but be careful with comparing homogeneities and densities.
 *
 * @tparam Particle
 * @param container container of current simulation
 * @param startCorner lower left front corner of the box where the homogeneity shall be calculated.
 * @param endCorner upper right back corner of the box where the homogeneity shall be calculated.
 * @return {homogeneity, maxDensity}
 */
template <class Container>
std::pair<double, double> calculateHomogeneityAndMaxDensity(const Container &container,
                                                            const std::array<double, 3> &startCorner,
                                                            const std::array<double, 3> &endCorner) {
  using namespace autopas::utils::ArrayMath::literals;

  const auto numberOfParticles = container.getNumberOfParticles();
  const auto domainDimensions = endCorner - startCorner;
  const auto domainVolume = domainDimensions[0] * domainDimensions[1] * domainDimensions[2];

  // We scale the dimensions of the domain to bins with volumes which give approximately 10 particles per bin.
  const auto numberOfBins = std::ceil(numberOfParticles / 10.);
  const auto binVolume = domainVolume / (double)numberOfBins;
  const auto scalingFactor = binVolume / domainVolume;
  const auto binDimensions = domainDimensions * scalingFactor;


  // Calculate the actual number of bins per dimension. The rounding is needed in case the division slightly
  // underestimates the division.
  const int binsPerDimension = static_cast<int>(std::round(domainVolume / binVolume));

  std::vector<size_t> particlesPerBin(numberOfBins, 0);

  // add particles accordingly to their bin to get the amount of particles in each bin
  for (auto particleItr = container.begin(autopas::IteratorBehavior::owned); particleItr.isValid(); ++particleItr) {
    const auto &particleLocation = particleItr->getR();

    const auto binIndex3d = autopas::utils::ArrayMath::floorToInt((particleLocation-startCorner)/binDimensions);
    const auto binIndex1d = autopas::utils::ThreeDimensionalMapping::threeToOneD(binIndex3d, binsPerDimension);

    particlesPerBin[binIndex1d] += 1;
  }

  // calculate density for each bin and track max density
  double maxDensity{0.};
  std::vector<double> densityPerBin;
  densityPerBin.reserve(numberOfBins);
  for (size_t i = 0; i < numberOfBins; i++) {
    densityPerBin[i] = binVolume / (double)particlesPerBin[i];
    if (densityPerBin[i] > maxDensity) {
      maxDensity = densityPerBin[i];
    }
  }

  if (maxDensity < 0.0) {
    throw std::runtime_error("maxDensity can never be smaller than 0.0, but is:" + std::to_string(maxDensity));
  }
  // get mean and reserve variable for densityVariance
  const double densityMean = numberOfParticles / domainVolume;
  double densityVariance = 0.0;

  // calculate densityVariance
  for (size_t i = 0; i < numberOfBins; ++i) {
    const double densityDifference = densityPerBin[i] - densityMean;
    densityVariance += (densityDifference * densityDifference / densityPerBin[i]);
  }

  // finally calculate standard deviation
  const double homogeneity = std::sqrt(densityVariance);
  // normally between 0.0 and 1.5
  if (homogeneity < 0.0) {
    throw std::runtime_error("homogeneity can never be smaller than 0.0, but is:" + std::to_string(homogeneity));
  }
  return {homogeneity, maxDensity};
}

}  // namespace autopas::utils
