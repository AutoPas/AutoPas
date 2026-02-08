/**
 * @file TreeliteBasedDecisionTreeTuning.cpp
 * @author Elizaveta Polysaeva
 * @date 18.12.25
 */

#include "TreeliteBasedDecisionTreeTuning.h"

#ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
#include <filesystem>
#include <json.hpp>
#include <type_traits>
#endif

#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

TreeliteBasedDecisionTreeTuning::TreeliteBasedDecisionTreeTuning(const std::set<Configuration> &searchSpace,
                                                                 const std::string &modelPairwiseFileName,
                                                                 const std::string &modelTriwiseFileName,
                                                                 double confidenceThreshold,
                                                                 InteractionTypeOption interactionType)
    : _configurations(searchSpace),
      _modelPairwiseFileName(modelPairwiseFileName),
      _modelTriwiseFileName(modelTriwiseFileName),
      _confidenceThreshold(confidenceThreshold),
      _interactionType(interactionType) {
#ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  try {
    namespace fs = std::filesystem;

    // Select model based on interaction type.
    const std::string &modelFile =
        (_interactionType == InteractionTypeOption::pairwise) ? _modelPairwiseFileName : _modelTriwiseFileName;

    if (modelFile.empty()) {
      utils::ExceptionHandler::exception(
          "TreeliteBasedDecisionTreeTuning: model file path is empty for interaction type '{}'.",
          _interactionType.to_string());
    }

    if (!modelFile.ends_with(".tl")) {
      utils::ExceptionHandler::exception(
          "TreeliteBasedDecisionTreeTuning expects a '.tl' Treelite model file. Got '{}'.", modelFile);
    }

    // Resolve paths relative to current working directory.
    const fs::path modelPath(modelFile);
    const fs::path basePath = modelPath.parent_path();
    const std::string baseStem = modelPath.stem().string();

    const fs::path classesPath = basePath / (baseStem + "_classes.txt");
    const fs::path featuresPath = basePath / "features.json";

    // Fail fast if anything is missing.
    if (!fs::exists(modelPath)) {
      utils::ExceptionHandler::exception("Treelite model file not found: {}", modelPath.string());
    }
    if (!fs::exists(classesPath)) {
      utils::ExceptionHandler::exception("Treelite classes file not found: {}", classesPath.string());
    }
    if (!fs::exists(featuresPath)) {
      utils::ExceptionHandler::exception("Treelite features file not found: {}", featuresPath.string());
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

    if (!_treeliteModel) {
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
    config.interactionType = _interactionType;

    configQueue.clear();
    configQueue.push_back(config);

  } catch (const std::exception &e) {
    AutoPasLog(ERROR, "Error parsing prediction: {}", e.what());
  }
#endif
}

}  // namespace autopas
