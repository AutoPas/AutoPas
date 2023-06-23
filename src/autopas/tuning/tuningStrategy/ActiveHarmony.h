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
#include "autopas/options/MPIStrategyOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/utils/AutoPasConfigurationCommunicator.h"
#include "autopas/utils/WrapMPI.h"
#include "hclient.h"

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
   * @param allowedContainerOptions
   * @param allowedCellSizeFactors
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param mpiStrategyOption
   * @param comm
   */
  ActiveHarmony(const std::set<ContainerOption> &allowedContainerOptions,
                const NumberSet<double> &allowedCellSizeFactors,
                const std::set<TraversalOption> &allowedTraversalOptions,
                const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                const std::set<Newton3Option> &allowedNewton3Options, const MPIStrategyOption mpiStrategyOption,
                const AutoPas_MPI_Comm comm);

  ~ActiveHarmony() override;

  void addEvidence(long time, size_t iteration) override;

  long getEvidence(Configuration configuration) const override;

  bool tune(bool currentInvalid) override;

  void removeN3Option(Newton3Option option) override;

  bool searchSpaceIsTrivial() const override;

  bool searchSpaceIsEmpty() const override;

  bool smoothedHomogeneityAndMaxDensityNeeded() const override;

  const Configuration &getCurrentConfiguration() const override;

  void reset(size_t iteration) override;

  std::set<ContainerOption> getAllowedContainerOptions() const override;

 private:
  /**
   * Pointer for the connection to the ActiveHarmony server.
   */
  hdesc_t *hdesc = nullptr;
  /**
   * Pointer to the ActiveHarmony tuning task defining the tuning parameters and tuning process.
   */
  htask_t *htask = nullptr;

  std::set<ContainerOption> _allowedContainerOptions;
  std::unique_ptr<NumberSet<double>> _allowedCellSizeFactors;
  std::set<TraversalOption> _allowedTraversalOptions;
  std::set<LoadEstimatorOption> _allowedLoadEstimatorOptions;
  std::set<DataLayoutOption> _allowedDataLayoutOptions;
  std::set<Newton3Option> _allowedNewton3Options;

  Configuration _currentConfig;

  MPIStrategyOption _mpiStrategyOption;
  AutoPas_MPI_Comm _comm;
  bool _nonLocalServer;

  /**
   * Traversal times for configurations. Used to save evidence when resetting active-harmony server.
   */
  std::unordered_map<Configuration, size_t, ConfigHash> _traversalTimes;

  /**
   * Resets the Harmony Server but keeps evidence.
   */
  void resetHarmony();

  /**
   * Fetch parameter-values from harmony server and update _currentConfig.
   */
  void fetchConfiguration();

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

  void setupTuningParameters(int commSize, hdef_t *hdef);

  static constexpr int cellSizeSamples = 100;

  static constexpr const char *traversalOptionName = "traversalOption";
  static constexpr const char *loadEstimatorOptionName = "loadEstimatorOption";
  static constexpr const char *dataLayoutOptionName = "dataLayoutOption";
  static constexpr const char *cellSizeFactorsName = "cellSizeFactor";
  static constexpr const char *newton3OptionName = "newton3Option";
};

}  // namespace autopas
