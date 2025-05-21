/**
 * @file ParticleBinStructure.h
 * @author S. Newcome
 * @date 13/12/2024
 */

#pragma once

#include <array>
#include <vector>

#include "ArrayMath.h"
#include "ExceptionHandler.h"
#include "ThreeDimensionalMapping.h"

namespace autopas::utils {
/**
 * Particle-counting bin structure. Stores all relevant information about the bins and a counter per bin. Provides a
 * particle counting function and can return statistics based on the particle count. Assumes all bins are of the same
 * shape.
 *
 * Expected use order:
 * 1) Constructed
 * 2) Particles counted
 * 3) Statistics calculated
 * 4) Statistics can be extracted from getters
 *
 * A flag ensures (4) cannot be done before (3) and that if particles are added to the count after (3), (3) must be
 * called again before (4).
 *
 * The purpose of this is to calculate all statistics together because it is rather trivial and efficient but provides
 * getters for individual statistics in case some are not desired.
 */
class ParticleBinStructure {
 public:
  /**
   * Constructor, that takes 3D array of the number of bins in each dimension and the dimensions of each bin, and
   * resizes space in the particle count structure as well as setting other relevant statistics.
   * @param numBinsPerDim Number of bins per dimension.
   * @param binLength Dimensions of each bin.
   * @param boxMin Lower left corner of region considered for particle binning.
   * @param boxMax Upper right corner of region considered for particle binning.
   * @param cutoff Not needed for particle counting, but needed for estimatedNumberOfNeighborInteractions statistic.
   */
  ParticleBinStructure(std::array<size_t, 3> numBinsPerDim, std::array<double, 3> binLength,
                       std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cutoff);

  /**
   * Constructor, for if there are the same number of bins in each dimension. Takes a number of bins and the dimensions
   * of each bin, and resizes space in the particle count structure as well as setting other relevant statistics.
   * @param numBinsPerDim Number of bins per dimension
   * @param binLength Dimensions of each bin.
   * @param boxMin Lower left corner of region considered for particle binning.
   * @param boxMax Upper right corner of region considered for particle binning.
   * @param cutoff Not needed for particle counting, but needed for estimatedNumberOfNeighborInteractions statistic.
   */
  ParticleBinStructure(size_t numBinsPerDim, std::array<double, 3> binLength, std::array<double, 3> boxMin,
                       std::array<double, 3> boxMax, double cutoff)
      : ParticleBinStructure({numBinsPerDim, numBinsPerDim, numBinsPerDim}, binLength, boxMin, boxMax, cutoff){};

  /**
   * Determines the appropriate bin for the particle based on its position, and
   * increments the corresponding bin counter to reflect one more particle in that bin.
   * @param particlePosition
   */
  void countParticle(const std::array<double, 3> &particlePosition);

  /**
   * Calculates the following statistics:
   * - Mean number of particles per bin
   * - Standard deviation in number of particles per bin
   * - Standard deviation in number of particles per bin relative to the mean.
   * - Maximum number of particles per bin
   * - Minimum number of particles per bin
   * - Median number of particles per bin
   * - Lower quartile number of particles per bin
   * - Upper quartile number of particles per bin
   * - Number of bins that are empty
   * - Estimated number of neighbor interactions
   */
  void calculateStatistics();

  /**
   * Resets all counters to zero.
   */
  void resetCounters();

  /**
   * Getter for the total number of particles counted
   * @return
   */
  [[nodiscard]] const size_t &getTotalParticleCount() const;

  /**
   * Getter for the particle counts.
   * @return
   */
  [[nodiscard]] const std::vector<size_t> &getParticleCounts() const;

  /**
   * Getter for the number of bins per dimension
   * @return
   */
  [[nodiscard]] const std::array<std::size_t, 3> &getNumBinsPerDim() const;

  /**
   * Setter for the number of bins per dimension
   * @param numBinsPerDim
   */
  void setCellsPerDim(const std::array<std::size_t, 3> &numBinsPerDim);

  /**
   * Getter for the dimensions of each bin.
   * @return
   */
  [[nodiscard]] const std::array<double, 3> &getBinLength() const;

  /**
   * Getter for the reciprocal of the dimensions of each bin.
   * @return
   */
  [[nodiscard]] const std::array<double, 3> &getBinLengthReciprocal() const;

  /**
   * Setter for the dimension of each bin. Also sets the reciprocal of these lengths.
   * @param binLength
   */
  void setBinLength(const std::array<double, 3> &binLength);

  /**
   * Getter for the box min.
   * @return
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMin() const;

  /**
   * Setter for the box min.
   * @param boxMin
   */
  void setBoxMin(const std::array<double, 3> &boxMin);

  /**
   * Getter for the box max.
   * @return
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMax() const;

  /**
   * Setter for the box max.
   * @param boxMax
   */
  void setBoxMax(const std::array<double, 3> &boxMax);

  /**
   * Getter for the number of bins.
   * @return
   */
  [[nodiscard]] std::size_t getNumberOfBins() const;

  /**
   * Resizes the particle count structure.
   */
  void resize();

  /**
   * Getter for the volume of each bin.
   * @return
   */
  [[nodiscard]] double getBinVolume() const;

  /**
   * Returns the mean number of particles per bins.
   * Logs a warning if statistics have not been calculated yet.
   *
   * @return The mean number of particles per bin.
   */
  [[nodiscard]] double getMeanParticlesPerBin() const;

  /**
   * Returns the standard deviation of the number of particles per bins.
   * Logs a warning if statistics have not been calculated yet.
   *
   * @return The standard deviation of the number of particles per bin.
   */
  [[nodiscard]] double getStdDevParticlesPerBin() const;

  /**
   * Returns the standard deviation of the number of particles per bins relative to the mean.
   * Logs a warning if statistics have not been calculated yet.
   *
   * @return The relative standard deviation of the number of particles per bin.
   */
  [[nodiscard]] double getRelStdDevParticlesPerBin() const;

  /**
   * Returns the mean density of particles per bins.
   * Logs a warning if statistics have not been calculated yet.
   *
   * @return The mean density of particles.
   */
  [[nodiscard]] double getMeanDensity() const;

  /**
   * Returns the standard deviation of the density per bins.
   * Logs a warning if statistics have not been calculated yet.
   *
   * @return The standard deviation of densities.
   */
  [[nodiscard]] double getStdDevDensity() const;

  /**
   * Returns the maximum density of any bin.
   * @return
   */
  [[nodiscard]] double getMaxDensity() const;

  /**
   * Returns the maximum number of particles in a single bin.
   * Logs a warning if statistics have not been calculated yet.
   *
   * @return The maximum number of particles in any bin.
   */
  [[nodiscard]] size_t getMaxParticlesPerBin() const;

  /**
   * Returns the minimum number of particles in a single bin.
   * Logs a warning if statistics have not been calculated yet.
   *
   * @return The minimum number of particles in any bin.
   */
  [[nodiscard]] size_t getMinParticlesPerBin() const;

  /**
   * Returns the median number of particles across all bins.
   * Logs a warning if statistics have not been calculated yet.
   *
   * @return The median number of particles per bin.
   */
  [[nodiscard]] size_t getMedianParticlesPerBin() const;

  /**
   * Returns the lower quartile (25th percentile) number of particles across all bins.
   * Logs a warning if statistics have not been calculated yet.
   *
   * @return The lower quartile of the number of particles per bin.
   */
  [[nodiscard]] size_t getLowerQuartileParticlesPerBin() const;

  /**
   * Returns the upper quartile (75th percentile) number of particles across all bins.
   * Logs a warning if statistics have not been calculated yet.
   *
   * @return The upper quartile of the number of particles per bin.
   */
  [[nodiscard]] size_t getUpperQuartileParticlesPerBin() const;

  /**
   * Returns the number of empty bins.
   * Logs a warning if statistics have not been calculated yet.
   *
   * @return The empty bin ratio.
   */
  [[nodiscard]] size_t getNumEmptyBins() const;

  /**
   * Returns the estimated number of neighbor interactions across all bins. This takes the assumption that particles are
   * evenly distributed within a bin and that this even distribution extends into neighboring bins and that newton3
   * isn't used.
   * Logs a warning if statistics have not been calculated yet.
   *
   * @return The estimated number of neighbor interactions.
   */
  [[nodiscard]] double getEstimatedNumberOfNeighborInteractions() const;

 private:
  /**
   * The actual particle count structure.
   */
  std::vector<std::size_t> _particleCounts{0};

  /**
   * Total number of particles counted.
   */
  size_t _totalParticleCount{0};

  /**
   * The number of bins per dimension.
   */
  std::array<std::size_t, 3> _numBinsPerDim{};

  /**
   * The length of each bin in each dimension.
   */
  std::array<double, 3> _binLength{};

  /**
   * The reciprocal of the length of each bin in each dimension.
   */
  std::array<double, 3> _binLengthReciprocal{};

  /**
   * The lower left corner of the box which the bin structure partitions.
   */
  std::array<double, 3> _boxMin{};

  /**
   * The upper right corner of the box which the bin structure partitions.
   */
  std::array<double, 3> _boxMax{};

  /**
   * Cutoff. Only needed for the statistics.
   */
  double _cutoff{};

  /**
   * Flag for if statistics have been calculated. Is set to true after calculateStatistics is called. Is set to false
   * if a particle is added to the count after this.
   */
  bool _statisticsCalculated = false;

  /**
   * Statistics
   */

  /**
   * Mean number of particles per bin.
   */
  double _meanParticlesPerBin{};

  /**
   * Mean density (num particles in a bin / bin volume) per bin.
   */
  double _meanDensity{};

  /**
   * Standard deviation in the number of particles per bin.
   */
  double _stdDevParticlesPerBin{};

  /**
   * Standard deviation in the number of particles per bin relative to the number of particles.
   */
  double _relStdDevParticlesPerBin{};

  /**
   * Standard deviation in the bin density per bin.
   */
  double _stdDevDensity{};

  /**
   * Maximum density in any bin.
   */
  double _maxDensity{};

  /**
   * Maximum number of particles that any bin has.
   */
  size_t _maxParticlesPerBin{};

  /**
   * Minimum number of particles that any bin has.
   */
  size_t _minimumParticlesPerBin{};

  /**
   * Median number of particles that any bin has.
   */
  size_t _medianParticlesPerBin{};

  /**
   * Lower quartile number of particles that any bin has.
   */
  size_t _lowerQuartileParticlesPerBin{};

  /**
   * Upper quartile number of particles that any bin has.
   */
  size_t _upperQuartileParticlesPerBin{};

  /**
   * Number of bins which are empty.
   */
  size_t _numEmptyBins{};

  /**
   * Estimation in the number of particle pair interactions. This takes the assumption that particles are evenly
   * distributed within a bin and that this even distribution extends into neighboring bins and that newton3 isn't used.
   *
   * At a 3x3x3 bins level, with bin dimensions about the size of an interaction length, this is maybe okay but the
   * validity of this estimate with this assumption is not tested.
   */
  double _estimatedNumberOfNeighborInteractions{};
};
}  // namespace autopas::utils
