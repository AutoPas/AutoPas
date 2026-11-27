/**
 * @file SimulatedAnnealing.h
 * @author Niklas Ladurner
 * @date 08.11.2026
 */

#pragma once

#include <optional>

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"

namespace autopas {

/// @brief Implements a stochastic tuning strategy based on Simulated Annealing. The selection process of which
/// configurations are proposed for sampling is guided by a temperature parameter. At high temperatures, the next
/// candidate configuration may vary significantly in its parameters from the current candidate. At lower temperatures,
/// the next candidate configuration is more likely to be similar in its parameters to the current candidate. The
/// temperature decreases geometrically over the tuning interval, so initially, the parameter space is explored and
/// later on the search settles down. The sampled proposal configuration is compared to the current configuration by
/// computation of  Δ = f(proposed)/f(current) - 1. Any configuration with Δ ≤ 0 is accepted. Any configuration with Δ >
/// 0 is accepted only with probability exp(-Δ / T_current).
class SimulatedAnnealing final : public TuningStrategyInterface {
 public:
  /// @brief Constructor for the Simulated Annealing tuning strategy.
  /// @param maxEvidence Maximum number of evidence to collect.
  /// @param seed Seed for the random number generator.
  /// @param initialTemperature Dimensionless T_0 > 0 used as the initial value for the temperature.
  /// @param endTemperature  Dimensionless T_end with T_0 > T_end > 0 used as the target value for the temperature in
  /// the last comparison of a tuning phase. The temperature decreases geometrically from T_0 to T_end over the whole
  /// tuning phase.
  /// @param exclusionFactor Untested candidate configurations that have been measured in a recent tuning phase (up to 2
  /// phases back) are not tested in the current tuning phase, if the recorded measurements predict a Δ >
  /// exclusionFactor * T_current. Can be disabled with exclusionFactor = inf.
  /// @param levelChangeVariance The tuning phase is divided up into overlapping normal distributions, for which this
  /// argument sets the variance. These distriutions indicate the probability of parameter change happening at a given
  /// "level" s.t. large changes such as swapping the container is more likely in the beginning of the tuning phase and
  /// smaller changes are more likely later on.
  explicit SimulatedAnnealing(size_t maxEvidence, unsigned long seed, double initialTemperature, double endTemperature,
                              double exclusionFactor, double levelChangeVariance);

  TuningStrategyOption getOptionType() const override;

  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

  bool optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

  bool reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

 private:
  /** Maximum number of evidence to collect. */
  size_t _maxEvidence;
  /** Current number of evidence collected. */
  size_t _numEvidenceCollected;

  /** Random number generator. */
  std::mt19937 _rng;

  /** Initial temperature of the annealing process. */
  double _initialTemperature;
  /** Current temperature of the annealing process. Decays geometrically from _initialTemperature to
  / _endTemperature. */
  double _currentTemperature;
  /** Target temperature of the annealing process. */
  double _endTemperature;

  // @todo _exclusionFactor and _levelChangeVariance should be derived empirically
  /** Factor at which to exclude candidates for which the predicted Δ is greater than exclusionFactor * T_current. */
  double _exclusionFactor;
  /** Variance parameter for the normal distribution handling the switching between levels over the course of the tuning
   * phase. */
  double _levelChangeVariance;

  /** Current state of the annealing process, contains a configuration and its sampled effective value. */
  struct AlgorithmState {
    Configuration config;
    long effectiveValue;
  };
  std::optional<AlgorithmState> _currentState;
  /** Current proposal configuration to sample next. */
  std::optional<Configuration> _pendingProposal;

  /** Factor by which the temperature cools down after each added evidence. */
  double _coolingFactor;
  /** Current tuning phase. */
  size_t _currentTuningPhase;
};

}  // namespace autopas
