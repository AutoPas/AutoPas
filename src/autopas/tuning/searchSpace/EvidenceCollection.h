/**
 * @file EvidenceCollection.h
 * @author F. Gratl
 * @date 28.06.23
 */

#pragma once

#include <spdlog/fmt/ostr.h>

#include <map>
#include <optional>
#include <tuple>
#include <vector>

#include "autopas/options/ContainerOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/searchSpace/Evidence.h"

namespace autopas {

/**
 * Class to manage all evidence.
 */
class EvidenceCollection {
 public:
  /**
   * Helper enum to switch between evidence measurements modes.
   */
  enum EvidenceMode {
    // Calculate the effective measurement for an average iteration
    EFFECTIVE,
    // Evidence of only the traversal measurement.
    TRAVERSAL,
    // Evidence of the rebuild + traversal measurement for a single iteration where a rebuild actually occurred.
    TOTAL
  };

  EvidenceCollection() = default;

  /**
   * Store a piece of evidence in the internal storage.
   * @param configuration
   * @param evidence
   */
  void addEvidence(const Configuration &configuration, const Evidence &evidence);

  /**
   * Returns all evidence collected for a given configuration.
   * @param configuration
   * @return If there is evidence a const pointer to the vector of collected evidence, otherwise nullptr.
   */
  [[nodiscard]] const std::vector<Evidence> *getEvidence(const Configuration &configuration) const;

  /**
   * Returns a modifiable reference to the last evidence of a given configuration.
   * @param configuration
   * @return Non-const reference to the last entry in the vector of a given configuration.
   */
  Evidence &modifyLastEvidence(const Configuration &configuration);

  /**
   * Retrieve the configuration with the best reduced evidence for a particular container plus cell size factor.
   * @param containerOption The container option to limit our configurations to.
   * @param cellSizeFactor The cellSizeFactor to limit our configurations to.
   * @return The optimal configuration and corresponding evidence
   */
  [[nodiscard]] std::tuple<Configuration, Evidence> getBestConfigForContainerAndCSF(ContainerOption containerOption,
                                                                                    double cellSizeFactor) const;

  /**
   * Retrieve the configuration with the best total evidence. Rebuild + Traversal.
   * @return The optimal configuration and corresponding evidence
   */
  [[nodiscard]] std::tuple<Configuration, Evidence> getBestConfigWithRebuild() const;

  /**
   * Retrieve the configuration with the lowest evidence value for the given tuning phase.
   * @param tuningPhase The tuning phase for which the optimum should be returned.
   * @param mode The Evidence criterion to look at (REDUCED, TRAVERSAL or TOTAL). Default = REDUCED
   * @param containerConstraint The container for which the optimum should be returned. Default = none
   * @param csfConstraint The cell size factor for which the optimum should be returned. Default = none
   * @return The optimal configuration.
   */
  [[nodiscard]] std::tuple<Configuration, Evidence> getOptimalConfiguration(
      size_t tuningPhase, EvidenceMode mode = EFFECTIVE,
      std::optional<ContainerOption> containerConstraint = std::nullopt,
      std::optional<double> csfConstraint = std::nullopt) const;

  /**
   * Retrieve the configuration with the lowest evidence value for the latest tuning phase.
   *
   * @return The optimal configuration.
   */
  [[nodiscard]] std::tuple<Configuration, Evidence> getLatestOptimalConfiguration() const;

  /**
   * Report if there is any evidence in the collection.
   * @param tuningPhase Optionally pass a tuning phase to check whether there is no evidence for this particular phase.
   * @return True if there is no evidence in the collection.
   */
  [[nodiscard]] bool empty(std::optional<size_t> tuningPhase = std::nullopt) const;

 private:
  /**
   * For each configuration the collection of all evidence (smoothed values) collected so far and in which iteration.
   * Configuration -> vector< iteration, time >
   */
  std::map<Configuration, std::vector<Evidence>> _evidenceMap{};

  /**
   * Number of the newest tuning phase for which evidence exists in this object.
   * This tuning phase is not necessarily completed yet.
   */
  size_t _latestTuningPhase{0};
};

/**
 * Enable readable logging for spdlog/fmt.
 * @param mode
 */
inline std::string format_as(EvidenceCollection::EvidenceMode mode) {
  switch (mode) {
    case EvidenceCollection::EFFECTIVE:
      return "Effective";
    case EvidenceCollection::TRAVERSAL:
      return "Traversal";
    case EvidenceCollection::TOTAL:
      return "Total";
    default:
      return "UnknownEvidenceMode";
  }
}

}  // namespace autopas