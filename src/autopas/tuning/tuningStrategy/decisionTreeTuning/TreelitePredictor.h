/**
 * @file TreelitePredictor.h
 * @author Elizaveta Polysaeva
 * @date 18.12.25
 */

#pragma once

// #ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
#include <treelite/c_api.h>

#include <map>
#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "autopas/tuning/Configuration.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * A Treelite-based predictor for decision-tree tuning.
 *
 * The predictor loads:
 *  - A Treelite checkpoint (*.tl).
 *  - A classes file (*_classes.txt) mapping class index -> combined label string.
 *  - A features file (features.json) storing the feature order used during training.
 *
 * The predict() method returns a JSON string compatible with updateConfigQueue().
 *
 * @note This class requires AUTOPAS_ENABLE_TREELITE_BASED_TUNING to be set to ON.
 */
class TreelitePredictor {
 public:
  /**
   * Constructor of the TreelitePredictor.
   * @param modelPath Path to model checkpoint (*.tl).
   * @param classesPath Path to classes file (*_classes.txt).
   * @param featuresPath Path to features.json (JSON array of strings).
   */
  TreelitePredictor(const std::string &modelPath, const std::string &classesPath, const std::string &featuresPath);

  /**
   * Destructor releasing model and GTIL config.
   */
  ~TreelitePredictor();

  /**
   * Predicts and returns a JSON string compatible with updateConfigQueue().
   * The JSON contains:
   *  - Container, Traversal, Load Estimator, Data Layout, Newton 3, CellSizeFactor
   *  - confidence
   * @param liveInfo Map of live info feature name -> numeric value.
   * @return JSON string.
   */
  std::string predict(const std::map<std::string, double> &liveInfo);

 private:
  /**
   * Loads features from features.json.
   * @param featuresPath Path to features.json.
   */
  void loadFeatures(const std::string &featuresPath);

  /**
   * Loads classes from the classes file.
   * @param classesPath Path to the classes file.
   */
  void loadClasses(const std::string &classesPath);

  /**
   * Loads the Treelite model checkpoint.
   * @param modelPath Path to the model checkpoint (*.tl).
   */
  void loadModel(const std::string &modelPath);

  /**
   * Initializes the Treelite GTIL config for inference.
   */
  void initGtilConfig();

  /**
   * Runs Treelite inference and fills the output buffer.
   * @param liveInfo Map of LiveInfo feature name to numeric value.
   */
  void runInference(const std::map<std::string, double> &liveInfo);

  /**
   * Finds argmax in output buffer (confidence).
   * @return Pair of index and confidence.
   */
  std::pair<int, double> getConfidence() const;

  /**
   * Predicts a class with highest confidence.
   * @param liveInfo Map of live info feature name to numeric value.
   * @return Pair of class and its confidence.
   */
  std::pair<std::string, double> getClassAndConfidence(const std::map<std::string, double> &liveInfo);

  /**
   * Custom deleter for TreeliteModelHandle used by std::unique_ptr.
   * Defines a way to destroy the resources.
   */
  struct TreeliteModelDeleter {
    void operator()(TreeliteModelHandle handle) const {
      if (handle != nullptr) {
        TreeliteFreeModel(handle);
      }
    }
  };

  /**
   * Custom deleter for TreeliteGTILConfigHandle used by std::unique_ptr.
   * Defines a way to destroy the resources.
   */
  struct TreeliteConfigDeleter {
    void operator()(TreeliteGTILConfigHandle handle) const {
      if (handle != nullptr) {
        TreeliteGTILDeleteConfig(handle);
      }
    }
  };

  /**
   * RAII wrapper types for Treelite handles.
   * Treelite handles are pointers, so we store them in std::unique_ptr with custom deleters to free memory
   * appropriately.
   */
  using TreeliteModelPtr = std::unique_ptr<std::remove_pointer_t<TreeliteModelHandle>, TreeliteModelDeleter>;
  using TreeliteConfigPtr = std::unique_ptr<std::remove_pointer_t<TreeliteGTILConfigHandle>, TreeliteConfigDeleter>;

  /**
   * RAII handle for Treelite model.
   */
  TreeliteModelPtr _model{nullptr};

  /**
   * RAII handle for Treelite GTIL config.
   */
  TreeliteConfigPtr _cfg{nullptr};

  /**
   * Feature names in the training order.
   */
  std::vector<std::string> _features;

  /**
   * Combined class labels as read from the classes file.
   */
  std::vector<std::string> _classes;

  /**
   * Buffer for a single row of live information.
   * Input for TreeliteGTILPredict.
   */
  std::vector<float> _row;

  /**
   * Buffer for probabilities / confidence.
   * Output of TreeliteGTILPredict.
   */
  std::vector<float> _out;

  /**
   * The number of labels expected per class string.
   */
  inline static constexpr std::size_t _numLabels = 6;

  /**
   * The expected feature list in the training order.
   */
  inline static const std::vector<std::string> _expectedFeatures = {
      "meanParticlesPerCell", "medianParticlesPerCell",
      "maxParticlesPerCell",  "relativeParticlesPerCellStdDev",
      "threadCount",          "numCells",
      "numEmptyCells",        "skin",
  };
};

}  // namespace autopas
// #endif
