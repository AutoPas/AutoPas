/**
 * @file DecisionTreeTuning.cpp
 * @author Abdulkadir Pazar
 * @date 20.06.24
 */

#include "DecisionTreeTuning.h"

#include <Python.h>

#include <iostream>
#include <sstream>

#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

// Constructor
DecisionTreeTuning::DecisionTreeTuning(const std::set<Configuration> &searchSpace, const std::string &modelFileName)
    : _configurations(searchSpace), _modelFileName(modelFileName), _pFunc(nullptr) {
  // Initialize Python Interpreter
  Py_Initialize();
  loadScript();
}

// Destructor
DecisionTreeTuning::~DecisionTreeTuning() {
  // Finalize Python Interpreter
  if (_pFunc) {
    Py_DECREF(_pFunc);  // Clean up the reference to _pFunc
  }
  Py_Finalize();
}

void DecisionTreeTuning::loadScript() {
  // Add the current directory to the Python path
  PyObject *sysPath = PySys_GetObject((char *)"path");
  PyObject *cwd = PyUnicode_FromString(".");
  PyList_Append(sysPath, cwd);
  Py_DECREF(cwd);

  // Load the fixed Python script predict.py
  PyObject *pName = PyUnicode_DecodeFSDefault("predict");
  if (!pName) {
    utils::ExceptionHandler::exception("Failed to convert script name to Python string.");
  }

  PyObject *pModule = PyImport_Import(pName);
  Py_DECREF(pName);

  // Check if the module was loaded successfully
  // If not, print the error and throw an exception
  if (!pModule) {
    PyErr_Print();
    utils::ExceptionHandler::exception("Failed to load Python module: predict.py");
  }

  // Retrieve the module's dictionary object, which stores all its defined variables, functions, and classes.
  // pDict will be used to access functions or objects within the "predict" module.
  PyObject *pDict = PyModule_GetDict(pModule);
  Py_DECREF(pModule);

  // Check if the dictionary was retrieved successfully
  // If not, print the error and throw an exception
  if (!pDict) {
    PyErr_Print();
    utils::ExceptionHandler::exception("Failed to get dictionary from Python module: predict.py");
  }

  // Get the Python function 'main' from the module
  // If the function does not exist or is not callable, print the error and throw an exception
  PyObject *pFunc = PyDict_GetItemString(pDict, "main");
  if (!pFunc || !PyCallable_Check(pFunc)) {
    PyErr_Print();
    utils::ExceptionHandler::exception("Failed to get Python function: main");
  }

  // Increment the reference count to keep the Python function alive
  // as _pFunc is stored as a persistent member of the class.
  // Without this, Python may garbage-collect the function after this scope.
  _pFunc = pFunc;
  Py_INCREF(_pFunc);
}

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
  // Convert live info to a JSON string
  std::string liveInfoJson = "{";
  for (const auto &[key, value] : _currentLiveInfo) {
    liveInfoJson += "\"" + key + "\":" + std::to_string(value) + ",";
  }
  liveInfoJson.back() = '}';  // Replace the last comma with a closing bracket

  // Call the Python function 'main' with the model file and live info JSON
  PyObject *firstArg = PyUnicode_FromString(_modelFileName.c_str());
  PyObject *secondArg = PyUnicode_FromString(liveInfoJson.c_str());
  PyObject *pArgs = PyTuple_Pack(2, firstArg, secondArg);

  PyObject *pResult = PyObject_CallObject(_pFunc, pArgs);

  Py_DECREF(firstArg);
  Py_DECREF(secondArg);
  Py_DECREF(pArgs);

  if (!pResult) {
    PyErr_Print();
    utils::ExceptionHandler::exception("Error occurred during Python function call");
  }

  // Extract the result (which is a JSON string of predictions)
  const char *prediction = PyUnicode_AsUTF8(pResult);
  std::string configPrediction(prediction);
  Py_DECREF(pResult);

  return configPrediction;
}

void DecisionTreeTuning::updateConfigQueue(std::vector<Configuration> &configQueue, const std::string &prediction) {
  std::stringstream ss(prediction);
  std::string item;
  Configuration config;

  try {
    auto getValue = [&](const std::string &key) -> std::string {
      size_t startPos = prediction.find("\"" + key + "\":");
      if (startPos == std::string::npos) {
        throw std::runtime_error("Key '" + key + "' not found in prediction");
      }
      startPos = prediction.find("\"", startPos + key.size() + 3) + 1;
      size_t endPos = prediction.find("\"", startPos);
      return prediction.substr(startPos, endPos - startPos);
    };

    // Using parseOptionExact to parse values for each configuration option
    config.container = ContainerOption::parseOptionExact(getValue("Container"));
    config.traversal = TraversalOption::parseOptionExact(getValue("Traversal"));
    config.dataLayout = DataLayoutOption::parseOptionExact(getValue("Data Layout"));
    config.newton3 = Newton3Option::parseOptionExact(getValue("Newton 3"));
    // Load values from the config queue's first element for load estimator and cell size factor
    config.cellSizeFactor = configQueue.front().cellSizeFactor;
    config.loadEstimator = configQueue.front().loadEstimator;

  } catch (const std::exception &e) {
    AutoPasLog(ERROR, "Exception during parsing prediction: {}", e.what());
    return;  // In case of an error, return without modifying the configQueue
  }

  // Clear and update the configQueue with the new configuration
  configQueue.clear();
  configQueue.push_back(config);
}

}  // namespace autopas
