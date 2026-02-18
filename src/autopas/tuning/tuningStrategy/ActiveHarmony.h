/**
 * @file ActiveHarmony.h
 * @author Jakob Englhauser
 * @date 21.10.2019
 */

#pragma once

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>

#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/LoadEstimators.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "autopas/utils/AutoPasConfigurationCommunicator.h"
#include "autopas/utils/WrapMPI.h"

#ifdef AUTOPAS_ENABLE_HARMONY
#include "hclient.h"
#endif

namespace autopas {

/**
 * Interface to the Active Harmony (AH) tuning framework.
 * If a global AH server is provided and MPI is enabled, the server can be reached via the environment variables
 * HARMONY_HOST and HARMONY_PORT.
 * They have to be set to the host address and port of the server respectively.
 */
class ActiveHarmony : public TuningStrategyInterface {
 public:
  /**
   * Constructor. Note that ActiveHarmony assumes every traversal option is only applicable for one container.
   * @param interactionType
   * @param allowedContainerOptions
   * @param allowedCellSizeFactors
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedVecPatternOptions
   * @param mpiDivideAndConquer
   * @param comm
   */
  ActiveHarmony(const InteractionTypeOption &interactionType, const std::set<ContainerOption> &allowedContainerOptions,
                const NumberSet<double> &allowedCellSizeFactors,
                const std::set<TraversalOption> &allowedTraversalOptions,
                const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                const std::set<Newton3Option> &allowedNewton3Options,
                const std::set<VectorizationPatternOption> &allowedVecPatternOptions, bool mpiDivideAndConquer,
                AutoPas_MPI_Comm comm);

  ~ActiveHarmony() override;

  TuningStrategyOption getOptionType() const override;

  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

  bool optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

  /**
   * Indicate if the search space contains only one configuration.
   * @return
   */
  bool searchSpaceIsTrivial() const;

  /**
   * Indicate if the search space contains any configurations.
   * @return
   */
  bool searchSpaceIsEmpty() const;

  bool needsDomainSimilarityStatistics() const override;

  bool reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

 private:
#ifdef AUTOPAS_ENABLE_HARMONY
  /**
   * Pointer for the connection to the ActiveHarmony server.
   */
  hdesc_t *hdesc = nullptr;
  /**
   * Pointer to the ActiveHarmony tuning task defining the tuning parameters and tuning process.
   */
  htask_t *htask = nullptr;
#else
  using hdef_t = void *;
#endif

  const InteractionTypeOption _interactionType;

  std::set<ContainerOption> _allowedContainerOptions;
  std::unique_ptr<NumberSet<double>> _allowedCellSizeFactors;
  std::set<TraversalOption> _allowedTraversalOptions;
  std::set<LoadEstimatorOption> _allowedLoadEstimatorOptions;
  std::set<DataLayoutOption> _allowedDataLayoutOptions;
  std::set<Newton3Option> _allowedNewton3Options;
  std::set<VectorizationPatternOption> _allowedVecPatterns;

  size_t _tuningPhase{0};
  bool _mpiDivideAndConquer;
  AutoPas_MPI_Comm _comm;
  bool _nonLocalServer;

  /**
   * Resets the Harmony Server but keeps evidence.
   */
  void resetHarmony();

  /**
   * Fetch parameter-values from harmony server and update _currentConfig.
   */
  Configuration fetchConfiguration();

  /**
   * Invalidates the current configuration by reporting the worst possible performance to the harmony server.
   */
  void invalidateConfiguration();

  /**
   * Define the tuning parameter in the ActiveHarmony tuning definition as enum and set possible values.
   * @tparam OptionClass Option class.
   * @param hdef Pointer to ActiveHarmony tuning definition.
   * @param name Name of the tuning parameter.
   * @param options Set of possible values of the tuning parameter.
   */
  template <class OptionClass>
  void configureTuningParameter(hdef_t *hdef, const char *name, const std::set<OptionClass> &options);

  /**
   * Fetch value for enum-type tuning parameter.
   * @tparam OptionClass Option class.
   * @param name Name of the tuning parameter.
   * @param options Set of all allowed options.
   * @return Value for tuning parameter.
   */
  template <class OptionClass>
  OptionClass fetchTuningParameter(const char *name, const std::set<OptionClass> &options);

  void rejectConfiguration(const autopas::Configuration &configuration, bool indefinitely) override;

  void setupTuningParameters(int commSize, hdef_t *hdef);

  static constexpr int cellSizeSamples = 100;

  static constexpr const char *traversalOptionName = "traversalOption";
  static constexpr const char *loadEstimatorOptionName = "loadEstimatorOption";
  static constexpr const char *dataLayoutOptionName = "dataLayoutOption";
  static constexpr const char *cellSizeFactorsName = "cellSizeFactor";
  static constexpr const char *newton3OptionName = "newton3Option";
  static constexpr const char *vecPatternOptionName = "vectorPattern";
};

}  // namespace autopas
