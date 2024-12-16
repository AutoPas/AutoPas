#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"
#include <pybind11/embed.h>

namespace autopas {

/**
 * @class DecisionTreeTuning
 *
 * This class uses a decision tree model to predict the best configuration based on the current live
 * information. The decision tree model is trained using the Python scikit-learn library. The Python
 * script predict.py is used to load the model and make predictions.
 */
class DecisionTreeTuning : public TuningStrategyInterface {
 public:
  /**
   * Constructor of DecisionTreeTuning.
   * @param searchSpace Set of configurations to be considered.
   * @param modelFileName Name of the file containing the decision tree model.
   * @param confidenceThreshold Minimum confidence threshold for accepting predictions.
   */
  DecisionTreeTuning(const std::set<Configuration> &searchSpace, const std::string &modelFileName, double confidenceThreshold);

  ~DecisionTreeTuning() override;

  [[nodiscard]] bool needsLiveInfo() const override;
  void receiveLiveInfo(const LiveInfo &value) override;
  bool reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const EvidenceCollection &evidenceCollection) override;
  [[nodiscard]] TuningStrategyOption getOptionType() const override;

  bool optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

 private:
  /**
   * Get the prediction from the Python script.
   * @return Prediction from the Python script.
   */
  std::string getPredictionFromPython();
  /**
   * Update the configuration queue based on the prediction.
   * @param configQueue Configuration queue to be updated.
   * @param prediction Prediction from the Python script.
   */
  void updateConfigQueue(std::vector<Configuration> &configQueue, const std::string &prediction);

  /**
   * Live Information for the current iteration used to make predictions.
   */
  std::map<std::string, double> _currentLiveInfo;
  /**
   * Set of configurations to be considered.
   */
  std::set<Configuration> _configurations;
  /**
   * Name of the file containing the decision tree model.
   */
  std::string _modelFileName;
  /**
   * Confidence threshold for the prediction.
   */
  double _confidenceThreshold;
  /**
   * Pybind11 object representing the Python function `main`.
   */
  pybind11::object _pythonMainFunc;
};

}  // namespace autopas
