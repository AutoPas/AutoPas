//
// Created by frederik on 4/12/25.
//

#include "ClusterBasedTuning.h"

#include "nlohmann/json.hpp"
#include "spdlog/fmt/bundled/chrono.h"

namespace autopas {

ClusterBasedTuning::ClusterBasedTuning(const std::set<Configuration> &searchSpace, const std::string &modelFileName,
                                       const std::string &configuration_mapping)
    : _configurations(searchSpace),
      _modelFileName(modelFileName),
      _configuration_mapping_FileName(configuration_mapping),
      _pFunc(nullptr) {
  Py_Initialize();
  load_python_script();
}
ClusterBasedTuning::~ClusterBasedTuning() {
  if (_pFunc) {
    Py_DECREF(_pFunc);
  }
  Py_Finalize();
}

TuningStrategyOption ClusterBasedTuning::getOptionType() const { return TuningStrategyOption::clusterBasedTuning; }

bool ClusterBasedTuning::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                             const EvidenceCollection &evidenceCollection) {
  return true;
}

void ClusterBasedTuning::load_python_script() {
  // returns a pointer to the sys.path
  PyObject *sysPath = PySys_GetObject("path");

  // append the directory examples/md-flexible/scripts to the sys.path list where python is looking for the modules.
  PyList_Append(sysPath, PyUnicode_FromString("examples/md-flexible/scripts"));

  // import the module
  PyObject *pName = PyUnicode_FromString("predict_kmeans");
  PyObject *pModule = PyImport_Import(pName);

  Py_DECREF(pName);

  if (!pModule) {
    PyErr_Print();
    utils::ExceptionHandler::exception("Failed to load Python Module: predict_kmeans");
  }

  // get the function main from the module
  PyObject *pFunc = PyObject_GetAttrString(pModule, "main");
  Py_DECREF(pModule);

  if (!pFunc) {
    PyErr_Print();
    utils::ExceptionHandler::exception("Failed to load Python function: main");
  }
  if (!PyCallable_Check(pFunc)) {
    PyErr_Print();
    utils::ExceptionHandler::exception("Failed to call Python function: main");
  }

  // save the function so that it can be called later
  _pFunc = pFunc;
}

void ClusterBasedTuning::receiveLiveInfo(const LiveInfo &info) {
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

std::string ClusterBasedTuning::get_Python_prediction() {
  // Convert live info to a JSON string
  std::string liveInfoJson = "{";
  for (const auto &[key, value] : _currentLiveInfo) {
    liveInfoJson += "\"" + key + "\":" + std::to_string(value) + ",";
  }
  liveInfoJson.back() = '}';  // Replace the last comma with a closing bracket

  // modelFile name
  PyObject *firstArgument = PyUnicode_FromString(_modelFileName.c_str());

  // configuration mapping
  PyObject *secondArgument = PyUnicode_FromString(_configuration_mapping_FileName.c_str());

  // liveInfo
  PyObject *thirdArgument = PyUnicode_FromString(liveInfoJson.c_str());

  // pack args together
  PyObject *pArgs = PyTuple_Pack(3, firstArgument, secondArgument, thirdArgument);

  if (!PyCallable_Check(_pFunc)) {
    PyErr_Print();
    utils::ExceptionHandler::exception("Failed to call Python function: main");
  }

  // python call
  PyObject *result = PyObject_CallObject(_pFunc, pArgs);

  Py_DECREF(firstArgument);
  Py_DECREF(secondArgument);
  Py_DECREF(thirdArgument);
  Py_DECREF(pArgs);
  if (!result) {
    utils::ExceptionHandler::exception("Error occured in result of call to Python function: main");
    PyErr_Print();
  }

  const char *configuration_predictions = PyUnicode_AsUTF8(result);
  std::string configPrediction(configuration_predictions);

  Py_DECREF(result);

  Py_INCREF(_pFunc);

  return configPrediction;
}

void ClusterBasedTuning::parse_configurations(std::string &configurations, std::vector<Configuration> &configQueue) {
  try {
    // parse the configuration string to json
    auto json_configs = nlohmann::json::parse(configurations);

    // loops over the different configurations of that cluster
    for (const auto &elem : json_configs) {
      Configuration cfg;
      // read in from prediction
      cfg.container = ContainerOption::parseOptionExact(elem["Container"].get<std::string>());
      cfg.traversal = TraversalOption::parseOptionExact(elem["Traversal"].get<std::string>());
      cfg.newton3 = Newton3Option::parseOptionExact(elem["Newton 3"].get<std::string>());
      cfg.dataLayout = DataLayoutOption::parseOptionExact(elem["Data Layout"].get<std::string>());

      std::cout << "Configuration {" + cfg.container.to_string() + ", " + cfg.traversal.to_string() + ", " +
                       cfg.newton3.to_string() + "," + cfg.dataLayout.to_string() + "}"
                << std::endl;

      // set constant
      cfg.loadEstimator = LoadEstimatorOption::none;
      cfg.cellSizeFactor = 1.0;
      cfg.interactionType = InteractionTypeOption::pairwise;

      // adds that configuration to the configuration queue
      configQueue.push_back(cfg);
    }
  } catch (const std::exception &e) {
    AutoPasLog(ERROR, "Exception during parsing prediction: {}", e.what());
  }
}

bool ClusterBasedTuning::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                               const EvidenceCollection &evidenceCollection) {
  // clear the current configQueue
  configQueue.clear();

  // get the Python prediction
  std::string prediction = get_Python_prediction();

  // parses the prediction into the configQueue, configQueue is updated
  parse_configurations(prediction, configQueue);

  return true;
}

};  // namespace autopas
