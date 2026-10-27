/**
 * @file PythonBasedTuning.h
 * @author Abdulkadir Pazar
 * @date 20.06.24
 */

#pragma once

#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
#include <pybind11/embed.h>
#endif

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"
#include "autopas/utils/PythonInterpreter.h"

namespace autopas {

/**
 * @class PythonBasedTuning
 *
 * This class uses a model written in Python to predict the best configuration based on the current live information.
 * The Python script randomForest/predict.py is used to load the model and make predictions. It provides a random forest
 * model, trained with randomForest/train.py using the Python scikit-learn library. For how to use a different model,
 * see docs/userdoc/PythonBasedMachineLearningTuning.md.
 *
 * @note This tuning strategy requires the CMake variable AUTOPAS_ENABLE_PYTHON_BASED_TUNING set to ON
 */
class PythonBasedTuning : public TuningStrategyInterface {
 public:
  /**
   * Constructor of PythonBasedTuning.
   * @param searchSpace Set of configurations to be considered.
   * @param modelFileName Name of the file containing the models.
   * @param confidenceThreshold Minimum confidence threshold for accepting predictions.
   * @param interactionType The interaction type (used to select the appropriate pairwise/triwise model).
   */
  PythonBasedTuning(const std::set<Configuration> &searchSpace, const std::string &modelFileName,
                    double confidenceThreshold, InteractionTypeOption interactionType);

  ~PythonBasedTuning() override;

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
   * Name of the file containing the model.
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

  /**
   * Keeps the embedded Python interpreter alive while this object exists. All PythonBasedTuning instances (e.g.
   * separate pairwise/triwise tuners) share it through PythonInterpreter::getSharedInterpreter(), since pybind11
   * permits only one interpreter per process.
   *
   * Keep this declaration before all `pybind11::object`s, so that they are destroyed before the interpreter is
   * released.
   */
  std::shared_ptr<utils::PythonInterpreter> _pythonInterpreter;

#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  /**
   * The actual python PythonBasedTuning object which we call .predict() on. Must be created and destroyed whilst
   * there is a live python interpreter.
   */
  pybind11::object _pythonBasedTuningPyObj;
#endif
};

}  // namespace autopas
