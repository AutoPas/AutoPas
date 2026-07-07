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
#include <memory>
#include <mutex>
#endif

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
  try {
    // The code within the try attempts to create an embedded python interpreter, or attach to an already existing one
    // (there cannot be two!). This requires a relatively complicated juggling of ownership, to ensure multiple
    // instances of DecisionTreeTuning (notable pairwise/triwise but future-proofed for other possibilities) have access
    // to the same pointer to the interpreter without the pointer itself being static (which leads to problems with
    // finalizing/destructing the python interpreter - particularly in combination with some implementations of MPI).

    // We use a combination of the object owned shared pointer to the interpreter with a static weak_ptr to ensure each
    // object points to the same interpreter, if one exists. As weak_ptrs do not own that which they point to, we need
    // not worry about the static weak_ptr keeping the interpreter alive after this object's destruction, the same way
    // a static shared_ptr would.

    // We use a locked mutex to avoid multiple threads attempting this at once (not currently possible but provides us
    // with some future-proofing). This is static so all threads have the same mutex.
    {
      static std::mutex interpreterMutex;
      static std::weak_ptr<py::scoped_interpreter> interpreterWeak;
      // Lock the mutex to make this thread-safe.
      const std::scoped_lock lock(interpreterMutex);
      // Set the object owned shared_ptr to that pointed to by the static weak_ptr
      // weak_ptr::lock() returns a shared_ptr that is null if no interpreter has yet been created or non-null if it has
      // been created (set below).
      _pythonInterpreter = interpreterWeak.lock();
      if (not _pythonInterpreter) {
        // weak_ptr returned a null pointer => create a shared scoped_interpreter
        _pythonInterpreter = std::make_shared<py::scoped_interpreter>();
        // Point the static weak reference at the new interpreter so the next instance's .lock()
        // finds and shares it.
        interpreterWeak = _pythonInterpreter;
      }
    }

    // Add the script directory to Python's path
    py::module::import("sys").attr("path").attr("append")(std::string(AUTOPAS_SOURCE_DIR) +
                                                          "/src/autopas/tuning/tuningStrategy/decisionTreeTuning");

    // Import the Python module and retrieve the 'load_model_and_encoder' and 'predict' functions
    py::module predictModule = py::module::import("predict");

    // Initialize the python object.
    _decisionTreeTuningPyObj = predictModule.attr("DecisionTreeTuning")(_modelFileName, _interactionType.to_string());

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
    config.vecPattern = VectorizationPatternOption::parseOptionExact(predictionJson["VectorizationPattern"]);

    config.interactionType = _interactionType;

    configQueue.clear();
    configQueue.push_back(config);

  } catch (const std::exception &e) {
    AutoPasLog(ERROR, "Error parsing prediction: {}", e.what());
  }
#endif
}

}  // namespace autopas
