/**
* @file ParticleBinStructure.cpp
* @author S. Newcome
* @date 13/12/2024
 */

//#include <array>
//#include <vector>

#include "ParticleBinStructure.h"

namespace autopas::utils {
ParticleBinStructure::ParticleBinStructure(std::array<size_t, 3> numBinsPerDim, std::array<double, 3> binLength, std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cutoff) {
  using namespace ArrayMath::literals;
  const auto numBins = numBinsPerDim[0] * numBinsPerDim[1] * numBinsPerDim[2];
  _particleCounts.resize(numBins);
  _numBinsPerDim = numBinsPerDim;
  _binLength = binLength;
  _binLengthReciprocal = 1. / binLength;
  _boxMin = boxMin;
  _boxMax = boxMax;
}

ParticleBinStructure::ParticleBinStructure(size_t numBinsPerDim, std::array<double, 3> binLength, std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cutoff) {
  using namespace ArrayMath::literals;
  const auto numBins = numBinsPerDim * numBinsPerDim * numBinsPerDim;
  _particleCounts.resize(numBins);
  _numBinsPerDim = {numBinsPerDim, numBinsPerDim, numBinsPerDim};
  _binLength = binLength;
  _binLengthReciprocal = 1. / binLength;
  _boxMin = boxMin;
  _boxMax = boxMax;
}

void ParticleBinStructure::countParticle(const std::array<double, 3> &particlePosition) {
  using namespace ArrayMath::literals;

  _statisticsCalculated = false;

  // Determine the 3D index of the bin the particle falls into
  const auto offsetIntoBox = particlePosition - _boxMin;
  const auto binIndex3DUnsafe = utils::ArrayMath::castedFloor<size_t>(offsetIntoBox * _binLengthReciprocal);

  // It is possible that floating point errors result in out of bounds indices.
  // e.g. if there are 7 bins in the x dimension, and that particle is close to the right domain boundary, the
  // division above might result in 7.0, which gets floored to 7 corresponding to a bin index that is out of bounds!

  // We therefore check for particle indices that are out of bounds and, if the particle is within the domain we
  // adjust the index. We don't care about floating point errors causing incorrect indices internally in the
  // domain as, for particles so close to a boundary, it is somewhat arbitrary which bin they fall into.

  const auto binIndex3DSafe = [&]() {
    auto newBinIndex3D = binIndex3DUnsafe;
    for (int dim = 0; dim < 3; ++dim) {
      if (particlePosition[dim] > _boxMin[dim] or particlePosition[dim] < _boxMax[dim]) {
        // Todo C++23 Use the size_t literal 0z
        newBinIndex3D[dim] = std::clamp(binIndex3DUnsafe[dim], 0UL, _numBinsPerDim[dim] - 1);
      } else {
        AutoPasLog(WARN, "Particle being counted is outside the box and will be ignored.");
      }
    }
    return newBinIndex3D;
  }();

  // Get the 1D index
  const auto binIndex1DSafe = utils::ThreeDimensionalMapping::threeToOneD<size_t>(binIndex3DSafe, _numBinsPerDim);

  // Add the particle to the bin and total
  _particleCounts[binIndex1DSafe] += 1;
  _totalParticleCount += 1;
}

void ParticleBinStructure::calculateStatistics() {
  // Sort particle counters in ascending order. Use a copy to maintain the structure of the particle counters.
  std::vector<std::size_t> sortedParticleCounts(_particleCounts);
  std::sort(sortedParticleCounts.begin(), sortedParticleCounts.end());

  // Get the minimum, maximum, median, and quartile particle counts
  _minimumNumberOfParticles = sortedParticleCounts.front();
  _maximumNumberOfParticles = sortedParticleCounts.back();
  _medianNumberOfParticles = sortedParticleCounts[sortedParticleCounts.size() / 2];
  _lowerQuartileNumberOfParticles = sortedParticleCounts[sortedParticleCounts.size() / 4];
  _upperQuartileNumberOfParticles = sortedParticleCounts[3 * sortedParticleCounts.size() / 4];

  // Determine the mean number of particles and density
  _meanNumberOfParticles = static_cast<double>(_totalParticleCount) / static_cast<double>(getNumberOfBins());
  _meanDensity = _meanNumberOfParticles / getBinVolume();


  // For the estimated number of neighbor interactions calculations, determine the estimated hit rate if the Linked
  // Cells method with cells the size of these bins was used, assuming a homogeneous distribution.
  // This is the ratio of cutoff sphere volume to the volume of 27 cells (number of cells within which one particle
  // could find neighbors given that the cell size factor is 1)
  const auto volumeOfCutoffSphere = 4. / 3. * M_PI * _cutoff * _cutoff * _cutoff;
  const auto potentialInteractionVolume = getBinVolume() * 27.;
  const auto estimatedHitRate = volumeOfCutoffSphere / potentialInteractionVolume;

  // Loop over all bins, determining the standard deviation in number of particles and densities, the number of empty
  // bins, and the estimated number of neighbor interactions
  double estimatedNumNeighborInteractionsSum = 0.;
  double numParticlesVarianceSum = 0.;
  double densityVarianceSum = 0.;
  size_t emptyBinCount = 0;
  double maximumDensity = 0;
  for (auto &particleCount : _particleCounts) {
    if (particleCount == 0) {
      ++emptyBinCount;
    } else {
      // Using the estimateHitRate, the number of particles in this cell, and the assumptions that particles are evenly
      // distributed within a bin and that this even distribution extends into neighboring bins and that newton3 isn't used, calculate the
      // number of neighbor interactions.
      // This is
      // - for every particle in this bin: [particleCount * ...]
      //   - it has distance checks with every particle in this bin and neighboring bin [... * (particleCount * 27)]
      //   - these checks have a hit rate of `estimatedHitRate`
      // - remove self interactions [... - particleCount]
      // In a very sparse situation, with a large potentialInteractionVolume and small particleCount, this could lead
      // to a negative estimate, so just take 0.
      estimatedNumNeighborInteractionsSum += std::max(
          static_cast<double>(particleCount * (particleCount * 27)) * estimatedHitRate  - particleCount, 0.);
    }

    const auto numParticlesDiffFromMean = particleCount - _meanNumberOfParticles;
    numParticlesVarianceSum += numParticlesDiffFromMean * numParticlesDiffFromMean;

    const auto density = particleCount / getBinVolume();
    const auto densityDiffFromMean = density - _meanDensity;
    densityVarianceSum += densityDiffFromMean * densityDiffFromMean;

    maximumDensity = std::max(maximumDensity, density);
  }

  _estimatedNumberOfNeighborInteractions = estimatedNumNeighborInteractionsSum;
  _standardDeviationNumberOfParticles = std::sqrt(numParticlesVarianceSum / static_cast<double>(getNumberOfBins()));
  _relativeStandardDeviationNumberOfParticles = _standardDeviationNumberOfParticles / _meanNumberOfParticles;
  _standardDeviationDensity = std::sqrt(densityVarianceSum / static_cast<double>(getNumberOfBins()));
  _numEmptyBins = emptyBinCount;
  _maxDensity = maximumDensity;

  _statisticsCalculated = true;
}

void ParticleBinStructure::resetCounters() {
  _totalParticleCount = 0;
  for (auto &particleCount : _particleCounts) {
    particleCount = 0;
  }
  _statisticsCalculated = false;
}

const size_t &ParticleBinStructure::getTotalParticleCount() const {
  return _totalParticleCount;
}

const std::array<std::size_t, 3> &ParticleBinStructure::getNumBinsPerDim() const {
  return _numBinsPerDim;
}

void ParticleBinStructure::setCellsPerDim(const std::array<std::size_t, 3> &numBinsPerDim) {
  _numBinsPerDim = numBinsPerDim;
}

const std::array<double, 3> &ParticleBinStructure::getBinLength() const {
  return _binLength;
}

const std::array<double, 3> &ParticleBinStructure::getBinLengthReciprocal() const {
  return _binLengthReciprocal;
}

void ParticleBinStructure::setBinLength(const std::array<double, 3> &binLength) {
  using namespace ArrayMath::literals;
  _binLength = binLength;
  _binLengthReciprocal = 1. / binLength;
}

const std::array<double, 3> &ParticleBinStructure::getBoxMin() const {
  return _boxMin;
}

void ParticleBinStructure::setBoxMin(const std::array<double, 3> &boxMin) {
  _boxMin = boxMin;
}

const std::array<double, 3> &ParticleBinStructure::getBoxMax() const {
  return _boxMax;
}

void ParticleBinStructure::setBoxMax(const std::array<double, 3> &boxMax) {
  _boxMax = boxMax;
}

std::size_t ParticleBinStructure::getNumberOfBins() { return _particleCounts.size(); }

void ParticleBinStructure::resize() {_particleCounts.resize(getNumberOfBins());}

double ParticleBinStructure::getBinVolume() { return _binLength[0] * _binLength[1] * _binLength[2]; }

double ParticleBinStructure::getMeanNumberOfParticles() {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _meanNumberOfParticles;
}

double ParticleBinStructure::getStandardDeviationNumberOfParticles() {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _standardDeviationNumberOfParticles;
}

double ParticleBinStructure::getRelativeStandardDeviationNumberOfParticles() {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _relativeStandardDeviationNumberOfParticles;
}

double ParticleBinStructure::getMeanDensity() {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _meanDensity;
}

double ParticleBinStructure::getStandardDeviationDensity() {
  if (!_statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _standardDeviationDensity;
}

double ParticleBinStructure::getMaxDensity() {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _maxDensity;
}

size_t ParticleBinStructure::getMaximumNumberOfParticles() {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _maximumNumberOfParticles;
}

size_t ParticleBinStructure::getMinimumNumberOfParticles() {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _minimumNumberOfParticles;
}

size_t ParticleBinStructure::getMedianNumberOfParticles() {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _medianNumberOfParticles;
}

size_t ParticleBinStructure::getLowerQuartileNumberOfParticles() {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _lowerQuartileNumberOfParticles;
}

size_t ParticleBinStructure::getUpperQuartileNumberOfParticles() {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _upperQuartileNumberOfParticles;
}

size_t ParticleBinStructure::getNumEmptyBins() {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet. Returning default/stored value for _emptyBinRatio.");
  }
  return _numEmptyBins;
}

double ParticleBinStructure::getEstimatedNumberOfNeighborInteractions() {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _estimatedNumberOfNeighborInteractions;
}


}