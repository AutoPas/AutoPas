/**
 * @file MPIParallelizedStrategy.h
 * @author W. Thieme
 * @date 27.05.2020
 */

#pragma once

#include <cstddef>
#include <random>
#include <type_traits>
#include <vector>

#include "TuningStrategyInterface.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/options/TuningStrategyOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "autopas/utils/AutoPasConfigurationCommunicator.h"
#include "autopas/utils/ConfigurationAndRankIteratorHandler.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/WrapMPI.h"
#include "options/ContainerOption.h"
#include "options/DataLayoutOption.h"
#include "options/LoadEstimatorOption.h"
#include "options/Newton3Option.h"
#include "options/TraversalOption.h"

namespace autopas {

/**
 * This strategy spreads the configuration queue in a round robin fashion over all ranks (with similar domain).
 *
 * The actual splitting of the search space and details of the communication logic is not currently handled by
 * this class, but by AutoPasConfigurationCommunicator.
 */
class MPIParallelizedStrategy : public TuningStrategyInterface {
 public:
  /**
   * Constructor.
   * @param fallbackConfiguration Generally applicable configuration that can be used if no local configuration are
   * applicable. It is suggested to use MPIParallelizedStrategy::createFallBackConfiguration() for this.
   * @param comm The communicator holding all ranks which participate in this tuning strategy.
   * @param mpiTuningMaxDifferenceForBucket
   * @param mpiTuningWeightForMaxDensity
   */
  MPIParallelizedStrategy(const Configuration &fallbackConfiguration, const AutoPas_MPI_Comm &comm,
                          double mpiTuningMaxDifferenceForBucket, double mpiTuningWeightForMaxDensity);

  void optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

  void receiveSmoothedHomogeneityAndMaxDensity(double homogeneity, double maxDensity) override;

  void reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

  [[nodiscard]] inline bool needsSmoothedHomogeneityAndMaxDensity() const override;

  void rejectConfiguration(const Configuration &configuration, bool indefinitely) override;

  /**
   * Create a resilient configuration that should always be applicable.
   *
   * The idea is to pick some configuration that is basically always applicable (lc_c08) and guess from the search space
   * which options Newton3 and Data Layout options are allowed.
   * @param searchSpace
   * @return The fall back configuration.
   */
  static Configuration createFallBackConfiguration(const std::set<Configuration> &searchSpace);

 private:
  /**
   * The smoothed homogeneity of the whole simulation.
   * See SimilarityFunctions::calculateHomogeneityAndMaxDensity().
   */
  double _smoothedHomogeneity{-1.};

  /**
   * The maximal density throughout the whole domain.
   * See SimilarityFunctions::calculateHomogeneityAndMaxDensity().
   */
  double _maxDensity{-1.};

  /**
   * The tuning strategy tuning locally
   */
  std::unique_ptr<TuningStrategyInterface> _tuningStrategy;

  /**
   * The full communicator of all AutoPas instances.
   */
  AutoPas_MPI_Comm _comm;
  /**
   * The communicator of AutoPas instances that this process shares its tuning with.
   */
  AutoPas_MPI_Comm _bucket{AUTOPAS_MPI_COMM_NULL};
  /**
   * Fallback configuration which should be always applicable, in case all local configurations are not applicable.
   */
  Configuration _fallbackConfiguration;
  /**
   * Indicator that we needed to use the fall back configuration.
   * If all ranks needed to do this we should abort.
   */
  bool _usingFallbackConfig{false};
  /**
   * Set of all configurations that were rejected in this tuning phase.
   */
  std::set<Configuration> _rejectedConfigurations{};
  /**
   * Maximum absolute difference in similarity metric for two ranks to fall in the same bucket.
   */
  double _mpiTuningMaxDifferenceForBucket;
  /**
   * Weight for maxDensity in the calculation for bucket distribution.
   */
  double _mpiTuningWeightForMaxDensity;
  /**
   * Random device used to determine the configuration if there are more ranks than configurations.
   */
  std::mt19937 _rng{std::random_device()()};
};
}  // namespace autopas
