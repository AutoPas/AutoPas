/**
 * @file TreeliteBasedDecisionTreeTuning.cpp
 * @author Elizaveta Polysaeva
 * @date 18.12.25
 */

#include "TreeliteBasedDecisionTreeTuning.h"

#ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
#include <algorithm>
#include <filesystem>
#include <json.hpp>
#include <string_view>
#endif
#include <type_traits>

#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

TreeliteBasedDecisionTreeTuning::TreeliteBasedDecisionTreeTuning(const std::set<Configuration> &searchSpace,
                                                                 const std::string &treeliteModelFileName,
                                                                 double confidenceThreshold)
    : _configurations(searchSpace), _modelFileName(treeliteModelFileName), _confidenceThreshold(confidenceThreshold) {
#ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  try {
    namespace fs = std::filesystem;

    if (_modelFileName.empty()) {
      utils::ExceptionHandler::exception("TreeliteBasedDecisionTreeTuning: model file path is empty: {}.",
                                         _modelFileName);
    }

    if (not _modelFileName.ends_with(".tl")) {
      utils::ExceptionHandler::exception(
          "TreeliteBasedDecisionTreeTuning expects a '.tl' Treelite model file, but it received '{}'.", _modelFileName);
    }

    // In next steps, files containing classes and features must be derived from the model file.

    // Files naming conventions:
    //   <model-prefix>_pairwise.tl / <model-prefix>_triwise.tl
    //   <model-prefix>_pairwise_classes.txt / <model-prefix>_triwise_classes.txt
    //   <model-prefix>_features.json

    // Resolve paths relative to current working directory.
    const fs::path modelPath(_modelFileName);           // full path to the .tl file
    const fs::path basePath = modelPath.parent_path();  // path to the parent directory

    const std::string modelStem = modelPath.stem().string();  // model filename without the .tl extension, e.g.
                                                              // my_model_pairwise part of my_model_pairwise.tl file

    // Get interaction type of the model.
    std::string_view interactionSuffix;
    if (modelStem.ends_with("_pairwise")) {
      interactionSuffix = "_pairwise";
    } else if (modelStem.ends_with("_triwise")) {
      interactionSuffix = "_triwise";
    } else {
      utils::ExceptionHandler::exception("Treelite model file '{}' must end with '_pairwise.tl' or '_triwise.tl'.",
                                         _modelFileName);
    }

    const std::string modelPrefix = modelStem.substr(
        0, modelStem.size() - interactionSuffix.size());  // e.g. model prefix of my_model_pairwise.tl is my_model

    const fs::path classesPath = basePath / (modelPrefix + std::string(interactionSuffix) +
                                             "_classes.txt");                   // e.g. my_model_pairwise_classes.txt
    const fs::path featuresPath = basePath / (modelPrefix + "_features.json");  // e.g. my_model_features.json

    // Accumulate missing file errors and throw once with all missing paths.
    std::vector<std::string> missingFiles;
    if (not fs::exists(modelPath)) {
      missingFiles.emplace_back(modelPath.string());
    }
    if (not fs::exists(classesPath)) {
      missingFiles.emplace_back(classesPath.string());
    }
    if (not fs::exists(featuresPath)) {
      missingFiles.emplace_back(featuresPath.string());
    }
    if (not missingFiles.empty()) {
      std::string missingFilesError{"Missing required Treelite files: "};
      for (size_t i = 0; i < missingFiles.size(); ++i) {
        missingFilesError += missingFiles[i];
        if (i + 1 < missingFiles.size()) {
          missingFilesError += ", ";
        }
      }
      utils::ExceptionHandler::exception("{}", missingFilesError);
    }

    // Construct Treelite predictor.
    _treeliteModel = std::make_unique<TreelitePredictor>(modelPath, classesPath, featuresPath);

  } catch (const std::exception &e) {
    utils::ExceptionHandler::exception("Failed to initialize Treelite model: {}", e.what());
  }
#else
  utils::ExceptionHandler::exception(
      "TreeliteBasedDecisionTreeTuning constructed but AUTOPAS_ENABLE_TREELITE_BASED_TUNING=OFF! "
      "Set this CMake variable to ON to use this tuning strategy.");
#endif
}

TreeliteBasedDecisionTreeTuning::~TreeliteBasedDecisionTreeTuning() = default;

bool TreeliteBasedDecisionTreeTuning::needsLiveInfo() const { return true; }

void TreeliteBasedDecisionTreeTuning::receiveLiveInfo(const LiveInfo &info) {
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

bool TreeliteBasedDecisionTreeTuning::reset(size_t iteration, size_t tuningPhase,
                                            std::vector<Configuration> &configQueue,
                                            const EvidenceCollection &evidenceCollection) {
  const std::string configPrediction = getPredictionFromTreelite();
  updateConfigQueue(configQueue, configPrediction);
  return true;
}

bool TreeliteBasedDecisionTreeTuning::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                                          const EvidenceCollection &evidenceCollection) {
  return true;
}

TuningStrategyOption TreeliteBasedDecisionTreeTuning::getOptionType() const {
  return TuningStrategyOption::treeliteBasedDecisionTreeTuning;
}

std::string TreeliteBasedDecisionTreeTuning::getPredictionFromTreelite() {
#ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  try {
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    utils::Timer timer;
    timer.start();
#endif

    if (not _treeliteModel) {
      utils::ExceptionHandler::exception("Treelite model not initialized.");
      return {};
    }

    const std::string result = _treeliteModel->predict(_currentLiveInfo);

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    timer.stop();
    AutoPasLog(TRACE, "Treelite prediction took {} ns.", timer.getTotalTime());
#endif
    return result;

  } catch (const std::exception &e) {
    utils::ExceptionHandler::exception("Error during Treelite prediction: {}", e.what());
    return {};
  }
#else
  utils::ExceptionHandler::exception(
      "getPredictionFromTreelite() called but AUTOPAS_ENABLE_TREELITE_BASED_TUNING=OFF.");
  return {};
#endif
}

void TreeliteBasedDecisionTreeTuning::updateConfigQueue(std::vector<Configuration> &configQueue,
                                                        const std::string &prediction) {
#ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
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
    if (predictionJson["CellSizeFactor"].is_string()) {
      config.cellSizeFactor = std::stod(static_cast<std::string>(predictionJson["CellSizeFactor"]));
    } else if (predictionJson["CellSizeFactor"].is_number()) {
      config.cellSizeFactor = predictionJson["CellSizeFactor"].get<double>();
    } else {
      AutoPasLog(ERROR, "Treelite prediction expected CellSizeFactor as string/number but got {}. Skipping update.",
                 predictionJson["CellSizeFactor"].type_name());
    }

    config.loadEstimator = LoadEstimatorOption::parseOptionExact(predictionJson["Load Estimator"]);

    configQueue.clear();
    configQueue.push_back(config);

  } catch (const std::exception &e) {
    AutoPasLog(ERROR, "Error parsing prediction: {}", e.what());
  }
#endif
}

}  // namespace autopas
