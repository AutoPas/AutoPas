/**
 *@file ClusterBasedTuning.h
 *@author Frederik Maximilian Willemsen
 *@date 05.05.2025
 */

#pragma once

#include <Python.h>

#include "autopas/tuning/Configuration.h"
#include "tuning/tuningStrategy/TuningStrategyInterface.h"

namespace autopas {

/**
 * @class ClusterBasedTuning
 *
 * This class implements a k-means based tuning strategy. The k-means model makes a cluster prediction based on the
 * liveInfo during the simulation and a configuration mapping is then used to assign configurations to the clusters. The
 * k-means pipeline model is trained with Python's scikit learn library, and the prediction is also done in Python.
 */
class ClusterBasedTuning : public TuningStrategyInterface {
 public:
  /**
   * Constructor of the ClusterBasedTuning.
   * @param searchSpace set of configurations to be considered, not used, but required to be read in through the
   * constructor.
   * @param modelFileName name of the k-means prediction model path.
   * @param configuration_mapping name of the configuration mapping path.
   */
  explicit ClusterBasedTuning(const std::set<Configuration> &searchSpace, const std::string &modelFileName,
                              const std::string &configuration_mapping);
  /**
   * Destructor of the ClusterBasedTuning.
   */
  ~ClusterBasedTuning() override;

  /**
   * Method indicates whether the tuning strategy needs liveInfo during runtime.
   * @return
   */
  [[nodiscard]] bool needsLiveInfo() const override { return true; }

  /**
   * This method loads the Python module predict_kmeans from the examples/md-flexible/scripts directory, imports the
   * module and retrieves the main function to store it for later use.
   */
  void load_python_script();

  /**
   * This function saves the current simulation LiveInfo and stores it for later use.
   * @param info current simulation LiveInfo
   */
  void receiveLiveInfo(const LiveInfo &info) override;

  /**
   * Invokes the Python module and gets the prediction.
   * @return
   */
  std::string get_Python_prediction();

  /**
   * Update the configuration queue to the configuration prediction.
   * @param configurations the configuration predictions.
   * @param configQueue the configuration Queue.
   */
  void parse_configurations(std::string &configurations, std::vector<Configuration> &configQueue);

  bool optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

  bool reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const EvidenceCollection &evidenceCollection) override;

  [[nodiscard]] TuningStrategyOption getOptionType() const override;

 private:
  /**
   * Live Information for the current iteration used to make predictions.
   */
  std::map<std::string, double> _currentLiveInfo;

  /**
   * vector of configurations to be considered.
   */
  std::set<Configuration> _configurations;

  /**
   * Name of the file containing the k-means pipeline model.
   */
  std::string _modelFileName;

  /**
   * Name of the configuration mapping
   */
  std::string _configuration_mapping_FileName;

  /**
   * Pointer to the Python function main in the predict_kmeans.py script.
   */
  PyObject *_pFunc;
};

}  // namespace autopas
