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
/// later on the search settles down. Te sampled proposal configuration is compared to the current configuration with Δ
/// = f(proposed)/f(current) - 1. Any configuration with Δ ≤ 0 is accepted. Any configuration with Δ > 0 is accepted
/// only with probability exp(-Δ / T_current).
class SimulatedAnnealing final : public TuningStrategyInterface {
 public:
  explicit SimulatedAnnealing(size_t maxEvidence, unsigned long seed, double initialTemperature, double endTemperature,
                              double exclusionFactor, double levelChangeVariance);

  TuningStrategyOption getOptionType() const override;

  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

  bool optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

  bool reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

 private:
  /// @brief Maximum number of evidence to collect.
  size_t _maxEvidence;
  /// @brief Current number of evidence collected.
  size_t _numEvidenceCollected;

  /// @brief Random number generator.
  std::mt19937 _rng;

  /// @brief Initial temperature of the annealing process.
  double _initialTemperature;
  /// @brief Current temperature of the annealing process. Decays geometrically from _initialTemperature to
  /// _endTemperature.
  double _currentTemperature;
  /// @brief Target temperature of the annealing process.
  double _endTemperature;

  // @todo _exclusionFactor and _levelChangeVariance should be derived empirically
  /// @brief Factor at which to exclude candidates for which the predicted Δ is greater than exclusionFactor * T_current
  double _exclusionFactor;
  /// @brief Variance of the normal distribution for switching between the level to perform changes at over the course
  /// of the tuning phase.
  double _levelChangeVariance;

  /// @brief Current state of the annealing process, contains a configuration and its sampled effective value.
  struct AlgorithmState {
    Configuration config;
    long effectiveValue;
  };
  std::optional<AlgorithmState> _currentState;
  /// @brief Current proposal configuration to sample next.
  std::optional<Configuration> _pendingProposal;

  /// @brief Factor by which the temperature cools down after each added evidence.
  double _coolingFactor;
  /// @brief Current tuning phase.
  size_t _currentTuningPhase;
};

}  // namespace autopas
