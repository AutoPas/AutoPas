/**
 * @file PythonBasedDecisionTreeTuning.cpp
 * @author Abdulkadir Pazar
 * @date 20.06.24
 */

#include "PythonBasedDecisionTreeTuning.h"

#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
#include <pybind11/embed.h>
#include <pybind11/stl.h>

#include <json.hpp>
#endif

#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
namespace py = pybind11;
#endif

PythonBasedDecisionTreeTuning::PythonBasedDecisionTreeTuning(const std::set<Configuration> &searchSpace, const std::string &modelFileName,
                                       double confidenceThreshold, InteractionTypeOption interactionType)
    : _configurations(searchSpace), _modelFileName(modelFileName), _confidenceThreshold(confidenceThreshold),
  _interactionType(interactionType) {
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  try {
    // Initialize the Python interpreter using scoped_interpreter
    static py::scoped_interpreter guard{};

    // Add the script directory to Python's path
    py::module::import("sys").attr("path").attr("append")(std::string(AUTOPAS_SOURCE_DIR) +
                                                          "/src/autopas/tuning/tuningStrategy/decisionTreeTuning");

    // Import the Python module and retrieve the 'load_model_and_encoder' and 'predict' functions
    py::module predictModule = py::module::import("predict");

    // Initialize the python object.
    _decisionTreeTuningPyObj = predictModule.attr("PythonBasedDecisionTreeTuning")(_modelFileName, _interactionType.to_string());

  } catch (const py::error_already_set &e) {
    utils::ExceptionHandler::exception("Failed to initialize Python environment: {}", e.what());
  }
#else
  utils::ExceptionHandler::exception("PythonBasedDecisionTreeTuning constructed but AUTOPAS_ENABLE_PYTHON_BASED_TUNING=OFF! "
"Set this CMake variable to ON to use this tuning strategy.");
#endif
}

PythonBasedDecisionTreeTuning::~PythonBasedDecisionTreeTuning() = default;

bool PythonBasedDecisionTreeTuning::needsLiveInfo() const { return true; }

void PythonBasedDecisionTreeTuning::receiveLiveInfo(const LiveInfo &info) {
  _currentLiveInfo.clear();
  for (const auto &infoEntry : info.get()) {
    const auto &name = infoEntry.first;
    const auto &entry = infoEntry.second;
    std::visit(
        [&](auto value) {
          if constexpr (std::is_arithmetic_v<decltype(value)>) {
            _currentLiveInfo[name] = value;
          }
        },
        entry);
  }
}

bool PythonBasedDecisionTreeTuning::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                               const EvidenceCollection &evidenceCollection) {
  std::string configPrediction = getPredictionFromPython();
  updateConfigQueue(configQueue, configPrediction);
  return true;
}

bool PythonBasedDecisionTreeTuning::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                             const EvidenceCollection &evidenceCollection) {
  return true;
}

TuningStrategyOption PythonBasedDecisionTreeTuning::getOptionType() const { return TuningStrategyOption::pythonBasedDecisionTreeTuning; }

std::string PythonBasedDecisionTreeTuning::getPredictionFromPython() {
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  try {
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    utils::Timer pythonPredictionTimer;
    pythonPredictionTimer.start();
#endif
    // Convert live info to JSON string
    const nlohmann::json liveInfoJson = _currentLiveInfo; // todo make this a reference
    // Call the Python function and get the result
    const py::object result = _decisionTreeTuningPyObj.attr("predict")(liveInfoJson.dump());
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    pythonPredictionTimer.stop();
    AutoPasLog(TRACE, "Python prediction took {} ms.", pythonPredictionTimer.getTotalTime());
#endif
    return result.cast<std::string>();

  } catch (const py::error_already_set &e) {
    utils::ExceptionHandler::exception("Error during Python function call: {}", e.what());
    return {};
  }
#endif
}

void PythonBasedDecisionTreeTuning::updateConfigQueue(std::vector<Configuration> &configQueue, const std::string &prediction) {
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  try {
    nlohmann::json predictionJson = nlohmann::json::parse(prediction);
    if (double confidence = predictionJson["confidence"]; confidence < _confidenceThreshold) {
      AutoPasLog(WARN, "Prediction confidence ({:.2f}) below threshold ({:.2f}), skipping update.", confidence,
                 _confidenceThreshold);
      return;
    }

    Configuration config;
    config.container = ContainerOption::parseOptionExact(predictionJson["Container"]);
    config.traversal = TraversalOption::parseOptionExact(predictionJson["Traversal"]);
    config.dataLayout = DataLayoutOption::parseOptionExact(predictionJson["Data Layout"]);
    config.newton3 = Newton3Option::parseOptionExact(predictionJson["Newton 3"]);
    // Sanity check that was stored as the expected type (i.e. a string).
    if (predictionJson["CellSizeFactor"].is_string()) {
      config.cellSizeFactor = std::stod(static_cast<std::string>(predictionJson["CellSizeFactor"]));
    } else {
      AutoPasLog(ERROR, "The Python predict.py script expected a string for CellSizeFactor, but a {} was returned.",
                 predictionJson["CellSizeFactor"].type_name());
    }
    config.loadEstimator = LoadEstimatorOption::parseOptionExact(predictionJson["Load Estimator"]);

    config.interactionType = _interactionType;

    configQueue.clear();
    configQueue.push_back(config);

  } catch (const std::exception &e) {
    AutoPasLog(ERROR, "Error parsing prediction: {}", e.what());
  }
#endif
}

}  // namespace autopas
