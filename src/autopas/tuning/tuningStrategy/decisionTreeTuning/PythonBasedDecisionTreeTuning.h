/**
 * @file PythonBasedDecisionTreeTuning.h
 * @author Abdulkadir Pazar
 * @date 20.06.24
 */

#pragma once

#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
#include <pybind11/embed.h>
#endif

#include <map>
#include <set>
#include <string>
#include <vector>

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"

namespace autopas {

/**
 * @class PythonBasedDecisionTreeTuning
 *
 * This class uses a decision tree model to predict the best configuration based on the current live
 * information. The decision tree model is trained using the Python scikit-learn library. The Python
 * script predict.py is used to load the model and make predictions.
 *
 * @note This tuning strategy requires the CMake variable AUTOPAS_ENABLE_PYTHON_BASED_TUNING set to ON
 */
class PythonBasedDecisionTreeTuning : public TuningStrategyInterface {
 public:
  /**
   * Constructor of PythonBasedDecisionTreeTuning.
   * @param searchSpace Set of configurations to be considered.
   * @param modelFileName Name of the file containing the random forest models.
   * @param confidenceThreshold Minimum confidence threshold for accepting predictions.
   * @param interactionType The interaction type (used to select the appropriate pairwise/triwise model).
   */
  PythonBasedDecisionTreeTuning(const std::set<Configuration> &searchSpace, const std::string &modelFileName,
                     double confidenceThreshold, InteractionTypeOption interactionType);

  ~PythonBasedDecisionTreeTuning() override;

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
   * The interaction type for which this tuning strategy predicts.
   */
  InteractionTypeOption _interactionType;

#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  pybind11::object _decisionTreeTuningPyObj;
#endif
};

}  // namespace autopas
