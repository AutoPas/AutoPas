/**
 * @file ParticleBinStructure.cpp
 * @author S. Newcome
 * @date 13/12/2024
 */

#include "ParticleBinStructure.h"

namespace autopas::utils {
ParticleBinStructure::ParticleBinStructure(std::array<size_t, 3> numBinsPerDim, std::array<double, 3> binLength,
                                           std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cutoff) {
  using namespace ArrayMath::literals;
  const auto numBins = numBinsPerDim[0] * numBinsPerDim[1] * numBinsPerDim[2];
  // Some error checking
  if (numBins < 1) {
    AutoPasLog(ERROR, "There must be at least one bin in the particle bin structure!");
  }
  if (std::any_of(binLength.begin(), binLength.end(), [](auto length) { return length <= 0; })) {
    AutoPasLog(ERROR, "Bin lengths must be positive!");
  }
  const auto boxDimensions = boxMax - boxMin;
  if (std::any_of(boxDimensions.begin(), boxDimensions.end(), [](auto length) { return length <= 0; })) {
    AutoPasLog(ERROR, "boxMax must be to the right/above boxMin!");
  }
  _particleCounts.resize(numBins);
  _numBinsPerDim = numBinsPerDim;
  _binLength = binLength;
  _binLengthReciprocal = 1. / binLength;
  _boxMin = boxMin;
  _boxMax = boxMax;
  _cutoff = cutoff;
}

template <typename T>
void ParticleBinStructure::countParticle(const std::array<T, 3> &particlePosition) {
  using namespace ArrayMath::literals;

  _statisticsCalculated = false;

  // Determine the 3D index of the bin the particle falls into
  const auto offsetIntoBox = particlePosition - _boxMin;
  const auto binIndex3DUnsafe = utils::ArrayMath::floorAndCast<size_t>(offsetIntoBox * _binLengthReciprocal);

  // It is possible that floating point errors result in out of bounds indices.
  // e.g. if there are 7 bins in the x dimension, and that particle is close to the right domain boundary, the
  // division above might result in 7.0, which gets floored to 7 corresponding to a bin index that is out of bounds!

  // We therefore check for particle indices that are out of bounds and, if the particle is within the domain we
  // adjust the index. We don't care about floating point errors causing incorrect indices internally in the
  // domain as, for particles so close to a boundary, it is somewhat arbitrary which bin they fall into.

  const auto binIndex3DSafe = [&]() {
    auto newBinIndex3D = binIndex3DUnsafe;
    for (int dim = 0; dim < 3; ++dim) {
      // Note: the following check also checks for "negative" indices (due to underflow.)
      if (binIndex3DUnsafe[dim] >= _numBinsPerDim[dim]) {
        // Add a little tolerance to what is considered inside the domain.
        if (particlePosition[dim] > _boxMin[dim] - 1e-12 and particlePosition[dim] < _boxMax[dim] + 1e-12) {
          // Todo C++23 Use the size_t literal 0z
          newBinIndex3D[dim] = std::clamp(binIndex3DUnsafe[dim], 0UL, _numBinsPerDim[dim] - 1);
        } else {
          AutoPasLog(WARN, "Particle being counted is outside the box and will be ignored.");
          // Assign index beyond domain to indicate invalidity
          newBinIndex3D[dim] = _numBinsPerDim[dim];
        }
      }
    }
    return newBinIndex3D;
  }();

  // If the particle is outside the domain, ignore it
  if (binIndex3DSafe[0] == _numBinsPerDim[0] or binIndex3DSafe[1] == _numBinsPerDim[1] or
      binIndex3DSafe[2] == _numBinsPerDim[2]) {
    return;
  }

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

  const bool noParticles =
      std::all_of(sortedParticleCounts.begin(), sortedParticleCounts.end(), [](auto count) { return count == 0; });

  // Get the minimum, maximum, median, and quartile particle counts
  _minimumParticlesPerBin = noParticles ? 0 : sortedParticleCounts.front();
  _maxParticlesPerBin = noParticles ? 0 : sortedParticleCounts.back();
  _medianParticlesPerBin = noParticles ? 0 : sortedParticleCounts[sortedParticleCounts.size() / 2];
  _lowerQuartileParticlesPerBin = noParticles ? 0 : sortedParticleCounts[sortedParticleCounts.size() / 4];
  _upperQuartileParticlesPerBin = noParticles ? 0 : sortedParticleCounts[3 * sortedParticleCounts.size() / 4];

  // Determine the mean number of particles and density
  _meanParticlesPerBin = static_cast<double>(_totalParticleCount) / static_cast<double>(getNumberOfBins());
  _meanDensity = _meanParticlesPerBin / getBinVolume();

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
      // distributed within a bin and that this even distribution extends into neighboring bins and that newton3 isn't
      // used, calculate the number of neighbor interactions. This is
      // - for every particle in this bin: [particleCount * ...]
      //   - it has distance checks with every particle in this bin and neighboring bin [... * (particleCount * 27)]
      //   - these checks have a hit rate of `estimatedHitRate`
      // - remove self interactions [... - particleCount]
      // In a very sparse situation, with a large potentialInteractionVolume and small particleCount, this could lead
      // to a negative estimate, so just take 0.
      estimatedNumNeighborInteractionsSum +=
          std::max(static_cast<double>(particleCount * (particleCount * 27)) * estimatedHitRate -
                       static_cast<double>(particleCount),
                   0.);
    }

    const auto numParticlesDiffFromMean = particleCount - _meanParticlesPerBin;
    numParticlesVarianceSum += numParticlesDiffFromMean * numParticlesDiffFromMean;

    const auto density = particleCount / getBinVolume();
    const auto densityDiffFromMean = density - _meanDensity;
    densityVarianceSum += densityDiffFromMean * densityDiffFromMean;

    maximumDensity = std::max(maximumDensity, density);
  }

  _estimatedNumberOfNeighborInteractions = estimatedNumNeighborInteractionsSum;
  _stdDevParticlesPerBin = std::sqrt(numParticlesVarianceSum / static_cast<double>(getNumberOfBins()));
  _relStdDevParticlesPerBin = noParticles ? 0 : _stdDevParticlesPerBin / _meanParticlesPerBin;
  _stdDevDensity = std::sqrt(densityVarianceSum / static_cast<double>(getNumberOfBins()));
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

const size_t &ParticleBinStructure::getTotalParticleCount() const { return _totalParticleCount; }

const std::vector<size_t> &ParticleBinStructure::getParticleCounts() const { return _particleCounts; }

const std::array<std::size_t, 3> &ParticleBinStructure::getNumBinsPerDim() const { return _numBinsPerDim; }

void ParticleBinStructure::setCellsPerDim(const std::array<std::size_t, 3> &numBinsPerDim) {
  _numBinsPerDim = numBinsPerDim;
}

const std::array<double, 3> &ParticleBinStructure::getBinLength() const { return _binLength; }

const std::array<double, 3> &ParticleBinStructure::getBinLengthReciprocal() const { return _binLengthReciprocal; }

void ParticleBinStructure::setBinLength(const std::array<double, 3> &binLength) {
  using namespace ArrayMath::literals;
  _binLength = binLength;
  _binLengthReciprocal = 1. / binLength;
}

const std::array<double, 3> &ParticleBinStructure::getBoxMin() const { return _boxMin; }

void ParticleBinStructure::setBoxMin(const std::array<double, 3> &boxMin) { _boxMin = boxMin; }

const std::array<double, 3> &ParticleBinStructure::getBoxMax() const { return _boxMax; }

void ParticleBinStructure::setBoxMax(const std::array<double, 3> &boxMax) { _boxMax = boxMax; }

std::size_t ParticleBinStructure::getNumberOfBins() const { return _particleCounts.size(); }

void ParticleBinStructure::resize() { _particleCounts.resize(getNumberOfBins()); }

double ParticleBinStructure::getBinVolume() const { return _binLength[0] * _binLength[1] * _binLength[2]; }

double ParticleBinStructure::getMeanParticlesPerBin() const {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _meanParticlesPerBin;
}

double ParticleBinStructure::getStdDevParticlesPerBin() const {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _stdDevParticlesPerBin;
}

double ParticleBinStructure::getRelStdDevParticlesPerBin() const {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _relStdDevParticlesPerBin;
}

double ParticleBinStructure::getMeanDensity() const {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _meanDensity;
}

double ParticleBinStructure::getStdDevDensity() const {
  if (!_statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _stdDevDensity;
}

double ParticleBinStructure::getMaxDensity() const {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _maxDensity;
}

size_t ParticleBinStructure::getMaxParticlesPerBin() const {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _maxParticlesPerBin;
}

size_t ParticleBinStructure::getMinParticlesPerBin() const {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _minimumParticlesPerBin;
}

size_t ParticleBinStructure::getMedianParticlesPerBin() const {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _medianParticlesPerBin;
}

size_t ParticleBinStructure::getLowerQuartileParticlesPerBin() const {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _lowerQuartileParticlesPerBin;
}

size_t ParticleBinStructure::getUpperQuartileParticlesPerBin() const {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _upperQuartileParticlesPerBin;
}

size_t ParticleBinStructure::getNumEmptyBins() const {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet. Returning default/stored value for _emptyBinRatio.");
  }
  return _numEmptyBins;
}

double ParticleBinStructure::getEstimatedNumberOfNeighborInteractions() const {
  if (not _statisticsCalculated) {
    AutoPasLog(WARN, "Statistics have not been calculated yet.");
  }
  return _estimatedNumberOfNeighborInteractions;
}

}  // namespace autopas::utils