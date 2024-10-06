/**
 * @file DecisionTreeTuning.cpp
 * @author Abdulkadir Pazar
 * @date 20.06.24
 */

#include "DecisionTreeTuning.h"
#include <Python.h>
#include <iostream>
#include <sstream>

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
  Py_Finalize();
}

void DecisionTreeTuning::loadScript() {
  // Add the current directory to the Python path
  PyObject *sysPath = PySys_GetObject((char*)"path");
  PyObject *cwd = PyUnicode_FromString(".");
  PyList_Append(sysPath, cwd);
  Py_DECREF(cwd);

  // Load the fixed Python script predict.py
  PyObject *pName = PyUnicode_DecodeFSDefault("predict");
  if (!pName) {
    throw std::runtime_error("Failed to convert script name to Python string.");
  }

  PyObject *pModule = PyImport_Import(pName);
  Py_DECREF(pName);

  if (!pModule) {
    PyErr_Print();
    throw std::runtime_error("Failed to load Python module: predict.py");
  }

  PyObject *pDict = PyModule_GetDict(pModule);
  Py_DECREF(pModule);

  if (!pDict) {
    PyErr_Print();
    throw std::runtime_error("Failed to get dictionary from Python module: predict.py");
  }

  PyObject *pFunc = PyDict_GetItemString(pDict, "main");
  if (!pFunc || !PyCallable_Check(pFunc)) {
    PyErr_Print();
    throw std::runtime_error("Failed to find function 'main' in module: predict.py");
  }

  _pFunc = pFunc;
  Py_INCREF(_pFunc);
}

bool DecisionTreeTuning::needsLiveInfo() const {
  return true;
}

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

TuningStrategyOption DecisionTreeTuning::getOptionType() const {
  return TuningStrategyOption::decisionTreeTuning;
}

std::string DecisionTreeTuning::getPredictionFromPython() {
  // Convert live info to a JSON string
  std::string liveInfoJson = "{";
  for (const auto &[key, value] : _currentLiveInfo) {
    liveInfoJson += "\"" + key + "\":" + std::to_string(value) + ",";
  }
  liveInfoJson.back() = '}';  // Replace the last comma with a closing bracket

  // Call the Python function 'main' with the model file and live info JSON
  PyObject *pArgs = PyTuple_Pack(2, PyUnicode_FromString(_modelFileName.c_str()), PyUnicode_FromString(liveInfoJson.c_str()));
  PyObject *pResult = PyObject_CallObject(_pFunc, pArgs);
  Py_DECREF(pArgs);

  if (!pResult) {
    PyErr_Print();
    throw std::runtime_error("Error occurred during Python function call");
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
    config.loadEstimator = LoadEstimatorOption::parseOptionExact(getValue("Load Estimator"));
    config.dataLayout = DataLayoutOption::parseOptionExact(getValue("Data Layout"));
    config.newton3 = Newton3Option::parseOptionExact(getValue("Newton 3"));
    config.cellSizeFactor = std::stod(getValue("CellSizeFactor"));
  } catch (const std::exception &e) {
    AutoPasLog(ERROR, "Exception during parsing prediction: {}", e.what());
    return;  // In case of an error, return without modifying the configQueue
  }

  // Clear and update the configQueue with the new configuration
  configQueue.clear();
  configQueue.push_back(config);
}

}  // namespace autopas
