/**
 * @file TreelitePredictor.cpp
 * @author Elizaveta Polysaeva
 * @date 18.12.25
 */

// #ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
#include "autopas/tuning/tuningStrategy/decisionTreeTuning/TreelitePredictor.h"

#include <cmath>
#include <fstream>
#include <json.hpp>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

namespace {

/**
 * Helper function that checks a Treelite return code and throws on error.
 * @param err Treelite return code.
 * @param msg Context message.
 */
void tlCheck(const int err, const std::string &msg) {
  if (err != 0) {
    autopas::utils::ExceptionHandler::exception("Treelite error: {} (code {}).", msg, err);
  }
}

/**
 * Helper function that removes a trailing carriage return from a line if present.
 * @param s Input string.
 * @return The string without a trailing carriage return.
 */
std::string stripCarriageReturn(std::string s) {
  if (!s.empty() && s.back() == '\r') {
    s.pop_back();
  }
  return s;
}

/**
 * Helper function that splits a string by a single-character delimiter.
 * @param s Input string.
 * @param delim Delimiter character.
 * @return The token list.
 */
std::vector<std::string> split(const std::string &s, const char delim) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream iss(s);
  while (std::getline(iss, token, delim)) {
    tokens.emplace_back(std::move(token));
  }
  return tokens;
}

}  // namespace

TreelitePredictor::TreelitePredictor(const std::string &modelPath, const std::string &classesPath,
                                     const std::string &featuresPath) {
  // Load features and classes before model to catch parsing error before model allocation.
  loadFeatures(featuresPath);
  loadClasses(classesPath);
  loadModel(modelPath);

  int numFeature = 0;
  tlCheck(TreeliteQueryNumFeature(_model.get(), &numFeature), "TreeliteQueryNumFeature failed.");

  if (numFeature != static_cast<int>(_features.size())) {
    autopas::utils::ExceptionHandler::exception("Treelite model expects {} features, but features.json provides {}. ",
                                                numFeature, _features.size());
  }

  initGtilConfig();

  _row.assign(_features.size(), 0.0f);
  _out.assign(_classes.size(), 0.0f);
}

TreelitePredictor::~TreelitePredictor() = default;

void TreelitePredictor::loadFeatures(const std::string &featuresPath) {
  std::ifstream file(featuresPath);
  if (!file.is_open()) {
    autopas::utils::ExceptionHandler::exception("Failed to open features.json: '{}'.", featuresPath);
  }

  nlohmann::json featuresJson;
  file >> featuresJson;

  if (!featuresJson.is_array()) {
    autopas::utils::ExceptionHandler::exception("features.json must be a JSON array: '{}'.", featuresPath);
  }

  _features.clear();
  _features.reserve(featuresJson.size());
  for (const auto &feature : featuresJson) {
    if (!feature.is_string()) {
      autopas::utils::ExceptionHandler::exception("features.json must contain only strings: '{}'.", featuresPath);
    }
    // Preserve order of features.
    _features.emplace_back(feature.get<std::string>());
  }

  if (_features.empty()) {
    autopas::utils::ExceptionHandler::exception("features.json is empty: '{}'.", featuresPath);
  }

  // Enforce exact feature list to prevent silent mismatch.
  if (_features != _expectedFeatures) {
    std::ostringstream oss;
    for (const auto &f : _expectedFeatures) oss << f << ' ';
    autopas::utils::ExceptionHandler::exception("features.json does not match expected list of features: {}",
                                                oss.str());
  }
}

void TreelitePredictor::loadClasses(const std::string &classesPath) {
  std::ifstream file(classesPath);
  if (!file.is_open()) {
    autopas::utils::ExceptionHandler::exception("Failed to open classes file: '{}'.", classesPath);
  }

  _classes.clear();
  std::string line;

  // Read classes line by line.
  while (std::getline(file, line)) {
    line = stripCarriageReturn(line);
    if (!line.empty()) {
      _classes.emplace_back(line);
    }
  }

  if (_classes.empty()) {
    autopas::utils::ExceptionHandler::exception("Classes file is empty: '{}'.", classesPath);
  }

  // Check that each class label can be parsed into a valid AutoPas option.
  for (const auto &cls : _classes) {
    const auto labels = split(cls, ';');

    if (labels.size() != _numLabels) {
      autopas::utils::ExceptionHandler::exception("TreelitePredictor expected {} labels, but got {}: '{}'.", _numLabels,
                                                  labels.size(), cls);
    }

    try {
      (void)ContainerOption::parseOptionExact(labels[0]);
      (void)TraversalOption::parseOptionExact(labels[1]);
      (void)LoadEstimatorOption::parseOptionExact(labels[2]);
      (void)DataLayoutOption::parseOptionExact(labels[3]);
      (void)Newton3Option::parseOptionExact(labels[4]);

      // CellSizeFactor must be numeric.
      size_t pos = 0;
      const double csf = std::stod(labels[5], &pos);
      if (pos != labels[5].size() || !std::isfinite(csf)) {
        autopas::utils::ExceptionHandler::exception("Invalid CellSizeFactor token '{}', class='{}'.", labels[5], cls);
      }
    } catch (const std::exception &e) {
      autopas::utils::ExceptionHandler::exception("Invalid label(s) in classes file. class='{}' error='{}'.", cls,
                                                  e.what());
    }
  }
}

void TreelitePredictor::loadModel(const std::string &modelPath) {
  TreeliteModelHandle modelHandle = nullptr;
  tlCheck(TreeliteDeserializeModelFromFile(modelPath.c_str(), &modelHandle),
          "TreeliteDeserializeModelFromFile failed for '" + modelPath + "'.");

  if (modelHandle == nullptr) {
    autopas::utils::ExceptionHandler::exception("Treelite returned a null model handle for '{}'.", modelPath);
  }

  // Take ownership of the new handle.
  _model.reset(modelHandle);
}

void TreelitePredictor::initGtilConfig() {
  // Minimal GTIL configuration:
  // - default prediction type
  // - single-threaded inference
  const char *cfgJson = R"({"predict_type":"default","nthread":1})";

  TreeliteGTILConfigHandle cfgHandle = nullptr;
  tlCheck(TreeliteGTILParseConfig(cfgJson, &cfgHandle), "TreeliteGTILParseConfig failed.");

  if (cfgHandle == nullptr) {
    autopas::utils::ExceptionHandler::exception("Treelite returned a null GTIL config handle.");
  }

  // Take ownership of the new handle.
  _cfg.reset(cfgHandle);
}

void TreelitePredictor::runInference(const std::map<std::string, double> &liveInfo) {
  if (_model == nullptr) {
    autopas::utils::ExceptionHandler::exception("Treelite model was not loaded.");
  }
  if (_cfg == nullptr) {
    autopas::utils::ExceptionHandler::exception("Treelite GTIL config was not initialized.");
  }

  // Build input row in the exact feature order used during training.
  // Missing features default to zero.
  for (size_t i = 0; i < _features.size(); ++i) {
    const auto it = liveInfo.find(_features[i]);
    const double value = (it != liveInfo.end()) ? it->second : 0.0;
    _row[i] = static_cast<float>(value);
  }

  // Run inference for a single row.
  tlCheck(TreeliteGTILPredict(_model.get(), static_cast<const void *>(_row.data()), "float32",
                              /*num_row=*/1, static_cast<void *>(_out.data()), _cfg.get()),
          "TreeliteGTILPredict failed.");

  for (float val : _out) {
    if (std::isnan(val) || std::isinf(val)) {
      utils::ExceptionHandler::exception("Treelite GTIL output contains NaN/Inf.");
    } else if (val < 0.0) {
      utils::ExceptionHandler::exception("Treelite GTIL output contains negative values.");
    }
  }
}

std::pair<int, double> TreelitePredictor::getConfidence() const {
  if (_out.empty()) {
    autopas::utils::ExceptionHandler::exception("Treelite predictor output buffer is empty.");
  }

  int bestIdx = 0;
  float bestVal = _out[0];
  double sum = 0.0;

  // Find argmax in output buffer and compute sum of all buffer values.
  for (int i = 0; i < static_cast<int>(_out.size()); ++i) {
    const float val = _out[i];
    sum += static_cast<double>(val);
    if (val > bestVal) {
      bestVal = val;
      bestIdx = i;
    }
  }

  // Normalize.
  if (sum <= 0.0) {
    autopas::utils::ExceptionHandler::exception("Treelite predictor output sum is non-positive.");
  }
  double confidence = static_cast<double>(bestVal) / sum;

  return {bestIdx, confidence};
}

std::pair<std::string, double> TreelitePredictor::getClassAndConfidence(const std::map<std::string, double> &liveInfo) {
  // Run inference on given live information.
  runInference(liveInfo);

  const auto [bestIdx, confidence] = getConfidence();

  if (bestIdx < 0 || bestIdx >= static_cast<int>(_classes.size())) {
    autopas::utils::ExceptionHandler::exception("Predicted class index {} is out of bounds.", bestIdx);
  }

  return {_classes[bestIdx], confidence};
}

std::string TreelitePredictor::predict(const std::map<std::string, double> &liveInfo) {
  const auto [combined, confidence] = getClassAndConfidence(liveInfo);

  // Get lables from a class string.
  const auto labels = split(combined, ';');

  nlohmann::json j;
  j["Container"] = labels[0];
  j["Traversal"] = labels[1];
  j["Load Estimator"] = labels[2];
  j["Data Layout"] = labels[3];
  j["Newton 3"] = labels[4];
  j["CellSizeFactor"] = labels[5];
  j["confidence"] = confidence;

  return j.dump();
}

}  // namespace autopas
// #endif
