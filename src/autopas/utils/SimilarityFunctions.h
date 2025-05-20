/**
 * @file SimilarityFunctions.h
 * @author J. Kroll
 * @date 01.06.2020
 */

#pragma once

#include "ThreeDimensionalMapping.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapMPI.h"

namespace autopas::utils {

/**
 * Calculates homogeneity and max density of given AutoPas container.
 * Both values are computed at once to avoid iterating over the same space twice.
 * homogeneity > 0.0, normally < 1.0, but for extreme scenarios > 1.0
 * maxDensity > 0.0, normally < 3.0, but for extreme scenarios > 3.0
 *
 * @warning This function is somewhat questionable. Firstly, it does not include buffer particles. Secondly, it is
 * unclear if there are any advantages of this binning scheme over e.g. cell-bins in LiveInfo. Use with caution!
 *
 * A homogeneity of 0 describes a perfectly homogeneous scenario, and as the scenario becomes more heterogeneous, this
 * becomes greater. If there are no particles, homogeneity is forced as 0 as this scenario is perfectly homogeneous.
 *
 * These values are calculated by dividing the container's domain into bins of perfectly equal cuboid shapes.
 *
 * @note The bin shapes between different AutoPas containers will not match if the dimensions of their domains don't
 * match. Bins with greater surface area are probably more likely to feature fluctuations in the density.
 *
 * Not a rigorous proof, but should give an idea:
 * - Assume homogeneous distribution.
 * - For every particle assume it could move in any direction with equal probability (not a terrible assumption).
 * - So every particle near a bin face has a chance to move into the neighboring bin.
 * - And this probability is independent of the bin (i.e. the bin shape)
 * - With a larger surface area, there are more particles with a probability of moving into the neighboring bin.
 * - On average, same density, but each individual bin has a greater probability of losing or gaining a lot of
 * particles.
 *
 * This is probably still fine for MPI Parallelized Tuning purposes, but be careful with comparing homogeneities and
 * densities.
 *
 * @tparam Particle
 * @param container container of current simulation
 * @return {homogeneity, maxDensity}
 */
template <typename Particle>
std::pair<double, double> calculateHomogeneityAndMaxDensity(const ParticleContainerInterface<Particle> &container) {
  using namespace autopas::utils::ArrayMath::literals;

  const auto numberOfParticles = container.getNumberOfParticles();
  const auto domainDimensions = container.getBoxMax() - container.getBoxMin();
  const auto domainVolume = domainDimensions[0] * domainDimensions[1] * domainDimensions[2];

  // We scale the dimensions of the domain to bins with volumes which give approximately 10 particles per bin.
  // Todo The choice of 10 is arbitrary and probably can be optimized
  const auto targetNumberOfBins = std::ceil(numberOfParticles / 10.);
  const auto targetNumberOfBinsPerDim = std::cbrt(static_cast<double>(targetNumberOfBins));
  // This is probably not an integer, so we floor to get more than 10 particles per bin than too small bins
  const auto numberOfBinsPerDim = static_cast<int>(std::floor(targetNumberOfBinsPerDim));
  const auto binDimensions = domainDimensions / static_cast<double>(numberOfBinsPerDim);

  const auto numberOfBins = numberOfBinsPerDim * numberOfBinsPerDim * numberOfBinsPerDim;
  const auto binVolume = domainVolume / static_cast<double>(numberOfBins);

  std::vector<size_t> particlesPerBin(numberOfBins, 0);

  // add particles accordingly to their bin to get the amount of particles in each bin
  for (auto particleItr = container.begin(autopas::IteratorBehavior::owned); particleItr.isValid(); ++particleItr) {
    const auto &particleLocation = particleItr->getR();

    const auto binIndex3dUnsafe =
        autopas::utils::ArrayMath::floorAndCast((particleLocation - container.getBoxMin()) / binDimensions);

    // It is possible that floating point errors result in out of bounds indices.
    // e.g. if there are 7 bins in the x dimension, and that particle is close to the right domain boundary, the
    // division above might result in 7.0, which gets floored to 7 corresponding to a bin index that is out of bounds!

    // We therefore check for particle indices that are out of bounds and, if the particle is within the domain we
    // adjust the index. We don't care about floating point errors causing incorrect indices internally in the
    // domain as, for particles so close to a boundary, it is somewhat arbitrary which bin they fall into.

    const auto binIndex3DSafe = [&]() {
      auto newBinIndex3D = binIndex3dUnsafe;
      for (int dim = 0; dim < 3; ++dim) {
        if (particleLocation[dim] > container.getBoxMin()[dim] or particleLocation[dim] < container.getBoxMax()[dim]) {
          newBinIndex3D[dim] = std::clamp(binIndex3dUnsafe[dim], 0, numberOfBinsPerDim - 1);
        } else {
          autopas::utils::ExceptionHandler::exception(
              "calculateHomogeneityAndMaxDensity: Particle is outside the container!");
        }
      }
      return newBinIndex3D;
    }();

    const auto binIndex1d = autopas::utils::ThreeDimensionalMapping::threeToOneD(
        binIndex3DSafe, {numberOfBinsPerDim, numberOfBinsPerDim, numberOfBinsPerDim});

    particlesPerBin[binIndex1d] += 1;
  }

  // calculate density for each bin and track max density
  double maxDensity{0.};
  std::vector<double> densityPerBin{};
  densityPerBin.reserve(numberOfBins);
  for (auto i = 0; i < numberOfBins; i++) {
    densityPerBin.push_back(static_cast<double>(particlesPerBin[i]) / binVolume);
    maxDensity = std::max(maxDensity, densityPerBin[i]);
  }

  if (maxDensity < 0.0) {
    utils::ExceptionHandler::exception(
        "calculateHomogeneityAndMaxDensity(): maxDensity can never be smaller than 0.0, but is: {}", maxDensity);
  }
  // get mean and reserve variable for densityVariance
  const double densityMean = numberOfParticles / domainVolume;

  const double densityDifferenceSquaredSum = std::transform_reduce(densityPerBin.begin(), densityPerBin.end(), 0.0,
                                                                   std::plus<>(), [densityMean](double density) {
                                                                     double densityDifference = density - densityMean;
                                                                     return densityDifference * densityDifference;
                                                                   });

  const auto densityVariance = densityDifferenceSquaredSum / numberOfBins;

  // finally calculate standard deviation. If there are no particles, the above calculation leads to nan. This should
  // be forced to 0.
  const double homogeneity = numberOfParticles == 0 ? 0. : std::sqrt(densityVariance);
  // normally between 0.0 and 1.5
  if (homogeneity < 0.0) {
    utils::ExceptionHandler::exception(
        "calculateHomogeneityAndMaxDensity(): homogeneity can never be smaller than 0.0, but is: {}", homogeneity);
  }
  return {homogeneity, maxDensity};
}

}  // namespace autopas::utils
