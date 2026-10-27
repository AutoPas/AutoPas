/**
 * @file DecisionTreeTuning.cpp
 * @author Abdulkadir Pazar
 * @date 20.06.24
 */

#include "DecisionTreeTuning.h"

#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
#include <pybind11/embed.h>
#include <pybind11/stl.h>

#include <nlohmann/json.hpp>
#include <sstream>
#endif

#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
namespace py = pybind11;
#endif

DecisionTreeTuning::DecisionTreeTuning(const std::set<Configuration> &searchSpace, const std::string &modelFileName,
                                       double confidenceThreshold, InteractionTypeOption interactionType)
    : _configurations(searchSpace),
      _modelFileName(modelFileName),
      _confidenceThreshold(confidenceThreshold),
      _interactionType(interactionType) {
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  _pythonInterpreter = utils::PythonInterpreter::getSharedInterpreter();
  try {
    // Add the script directory to Python's path
    py::module::import("sys").attr("path").attr("append")(std::string(AUTOPAS_SOURCE_DIR) +
                                                          "/src/autopas/tuning/tuningStrategy/decisionTreeTuning");

    // Import the Python module and retrieve the 'load_model_and_encoder' and 'predict' functions
    py::module predictModule = py::module::import("predict");

    // Initialize the python object.
    _decisionTreeTuningPyObj = predictModule.attr("DecisionTreeTuning")(_modelFileName, _interactionType.to_string());

    // The model must predict exactly the components that configurations consist of. Otherwise, it was trained before a
    // component was added, removed, or renamed.
    const auto modelComponentsList = _decisionTreeTuningPyObj.attr("config_columns").cast<std::vector<std::string>>();
    const std::set<std::string> modelComponents(modelComponentsList.begin(), modelComponentsList.end());
    std::set<std::string> configComponents;
    std::istringstream configHeader(Configuration().getCSVHeader());
    for (std::string component; std::getline(configHeader, component, ',');) {
      configComponents.insert(component);
    }
    if (modelComponents != configComponents) {
      utils::ExceptionHandler::exception(
          "The model {} predicts the configuration components {}, but AutoPas configurations consist of {}. Retrain "
          "the model with the current train.py.",
          _modelFileName, utils::ArrayUtils::to_string(modelComponents),
          utils::ArrayUtils::to_string(configComponents));
    }

  } catch (const py::error_already_set &e) {
    utils::ExceptionHandler::exception("Failed to initialize Python environment: {}", e.what());
  }
#else
  utils::ExceptionHandler::exception(
      "DecisionTreeTuning constructed but AUTOPAS_ENABLE_PYTHON_BASED_TUNING=OFF! "
      "Set this CMake variable to ON to use this tuning strategy.");
#endif
}

DecisionTreeTuning::~DecisionTreeTuning() = default;

bool DecisionTreeTuning::needsLiveInfo() const { return true; }

void DecisionTreeTuning::receiveLiveInfo(const LiveInfo &info) {
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

bool DecisionTreeTuning::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                               const EvidenceCollection &evidenceCollection) {
  std::string configPrediction = getPredictionFromPython();
  updateConfigQueue(configQueue, configPrediction);
  return true;
}

bool DecisionTreeTuning::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                             const EvidenceCollection &evidenceCollection) {
  return true;
}

TuningStrategyOption DecisionTreeTuning::getOptionType() const { return TuningStrategyOption::decisionTreeTuning; }

std::string DecisionTreeTuning::getPredictionFromPython() {
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  try {
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    utils::Timer pythonPredictionTimer;
    pythonPredictionTimer.start();
#endif
    // Convert live info to JSON string
    const nlohmann::json liveInfoJson = _currentLiveInfo;  // todo make this a reference
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
#else
  return {};
#endif
}

void DecisionTreeTuning::updateConfigQueue(std::vector<Configuration> &configQueue, const std::string &prediction) {
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  try {
    nlohmann::json predictionJson = nlohmann::json::parse(prediction);
    if (double confidence = predictionJson["confidence"]; confidence < _confidenceThreshold) {
      AutoPasLog(WARN, "Prediction confidence ({:.2f}) below threshold ({:.2f}), skipping update.", confidence,
                 _confidenceThreshold);
      return;
    }

    std::stringstream configStream;
    std::istringstream configHeader(Configuration().getCSVHeader());
    for (std::string component; std::getline(configHeader, component, ',');) {
      configStream << component << ": " << predictionJson.at(component).get<std::string>() << " , ";
    }
    Configuration config;
    configStream >> config;

    configQueue.clear();
    configQueue.push_back(config);

  } catch (const std::exception &e) {
    AutoPasLog(ERROR, "Error parsing prediction: {}", e.what());
  }
#endif
}

}  // namespace autopas
