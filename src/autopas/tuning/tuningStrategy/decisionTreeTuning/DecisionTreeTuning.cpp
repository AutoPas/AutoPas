/**
 * @file DecisionTreeTuning.cpp
 * @author Abdulkadir Pazar
 * @date 20.06.24
 */

#include "DecisionTreeTuning.h"

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

DecisionTreeTuning::DecisionTreeTuning(const std::set<Configuration> &searchSpace, const std::string &modelFileName,
                                       double confidenceThreshold)
    : _configurations(searchSpace), _modelFileName(modelFileName), _confidenceThreshold(confidenceThreshold) {
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  try {
    // Initialize the Python interpreter using scoped_interpreter
    static py::scoped_interpreter guard{};

    // Add the script directory to Python's path
    py::module::import("sys").attr("path").attr("append")(std::string(AUTOPAS_SOURCE_DIR) +
                                                          "/src/autopas/tuning/tuningStrategy/decisionTreeTuning");

    // Import the Python module and retrieve the 'main' function
    py::module predictModule = py::module::import("predict");
    _pythonMainFunc = predictModule.attr("main");
  } catch (const py::error_already_set &e) {
    utils::ExceptionHandler::exception("Failed to initialize Python environment: {}", e.what());
  }
#else
  utils::ExceptionHandler::exception("DecisionTreeTuning constructed but AUTOPAS_ENABLE_PYTHON_BASED_TUNING=OFF! "
"Set this CMake variable to ON to use this tuning strategy.");
#endif
}

DecisionTreeTuning::~DecisionTreeTuning() {}

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
    // Convert live info to JSON string
    nlohmann::json liveInfoJson = _currentLiveInfo;
    std::string modelPath = _modelFileName;
    // Call the Python function and get the result
    py::object result = _pythonMainFunc(modelPath, liveInfoJson.dump());
    return result.cast<std::string>();
  } catch (const py::error_already_set &e) {
    utils::ExceptionHandler::exception("Error during Python function call: {}", e.what());
    return {};
  }
#endif
}

void DecisionTreeTuning::updateConfigQueue(std::vector<Configuration> &configQueue, const std::string &prediction) {
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  nlohmann::json predictionJson;
  try {
    predictionJson = nlohmann::json::parse(prediction);
    double confidence = predictionJson["confidence"];
    if (confidence < _confidenceThreshold) {
      AutoPasLog(WARN, "Prediction confidence ({:.2f}) below threshold ({:.2f}), skipping update.", confidence,
                 _confidenceThreshold);
      return;
    }

    Configuration config;
    config.container = ContainerOption::parseOptionExact(predictionJson["Container"]);
    config.traversal = TraversalOption::parseOptionExact(predictionJson["Traversal"]);
    config.dataLayout = DataLayoutOption::parseOptionExact(predictionJson["Data Layout"]);
    config.newton3 = Newton3Option::parseOptionExact(predictionJson["Newton 3"]);
    config.cellSizeFactor = std::stod(getValue("CellSizeFactor"));
    config.loadEstimator = LoadEstimatorOption::parseOptionExact(getValue("Load Estimator"));

    config.interactionType = configQueue.front().interactionType;

    configQueue.clear();
    configQueue.push_back(config);

  } catch (const std::exception &e) {
    AutoPasLog(ERROR, "Error parsing prediction: {}", e.what());
  }
#endif
}

}  // namespace autopas
