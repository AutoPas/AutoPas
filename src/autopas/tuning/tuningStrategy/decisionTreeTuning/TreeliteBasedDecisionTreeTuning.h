/**
 * @file TreeliteBasedDecisionTreeTuning.h
 * @author Elizaveta Polysaeva
 * @date 18.12.25
 */

#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"

#ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
#include <memory>

#include "autopas/tuning/tuningStrategy/decisionTreeTuning/TreelitePredictor.h"
#endif

namespace autopas {

/**
 * @class TreeliteBasedDecisionTreeTuning
 *
 * This class uses a decision tree model to predict the best configuration based on the current live
 * information. The decision tree model is trained using the Python scikit-learn library. Class
 * TreelitePredictor is used to load the model and make predictions.
 *
 * @note This tuning strategy requires the CMake variable AUTOPAS_ENABLE_TREELITE_BASED_TUNING set to ON
 */
class TreeliteBasedDecisionTreeTuning : public TuningStrategyInterface {
 public:
  /**
   * Constructor of TreeliteBasedDecisionTreeTuning.
   * @param searchSpace Set of configurations to be considered.
   * @param modelPairwiseFileName Name of the file containing the pairwise decision tree model.
   * @param modelTriwiseFileName Name of the file containing the triwise decision tree model.
   * @param confidenceThreshold Minimum confidence threshold for accepting predictions.
   * @param interactionType The interaction type (used to select the appropriate pairwise/triwise model).
   */
  TreeliteBasedDecisionTreeTuning(const std::set<Configuration> &searchSpace, const std::string &modelPairwiseFileName,
                                  const std::string &modelTriwiseFileName, double confidenceThreshold,
                                  InteractionTypeOption interactionType);

  ~TreeliteBasedDecisionTreeTuning() override;

  [[nodiscard]] bool needsLiveInfo() const override;
  void receiveLiveInfo(const LiveInfo &value) override;
  bool reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const EvidenceCollection &evidenceCollection) override;
  [[nodiscard]] TuningStrategyOption getOptionType() const override;

  bool optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

 private:
  /**
   * Runs a Treelite prediction and returns a JSON prediction string.
   * @return A JSON string containing the predicted configuration and a confidence value.
   */
  std::string getPredictionFromTreelite();

  /**
   * Updates the configuration queue based on a JSON prediction string.
   * @param configQueue Configuration queue to be updated.
   * @param prediction JSON prediction string.
   */
  void updateConfigQueue(std::vector<Configuration> &configQueue, const std::string &prediction);

  /**
   * Live Information for the current iteration used to make predictions.
   */
  std::map<std::string, double> _currentLiveInfo;

  /**
   * Set of configurations to be considered.
   */
  std::set<Configuration> _configurations;

  /**
   * Name of the file containing the pairwise decision tree model.
   */
  std::string _modelPairwiseFileName;

  /**
   * Name of the file containing the triwise decision tree model.
   */
  std::string _modelTriwiseFileName;

  /**
   * Confidence threshold for the prediction.
   */
  double _confidenceThreshold;

  /**
   * The interaction type for which this tuning strategy predicts.
   */
  InteractionTypeOption _interactionType;

#ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  /**
   * Treelite predictor used to run model inference.
   */
  std::unique_ptr<TreelitePredictor> _treeliteModel;
#endif
};

}  // namespace autopas
