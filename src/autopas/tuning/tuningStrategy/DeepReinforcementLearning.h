/**
 * @file DeepReinforcementLearning.h
 * @author P. Metscher
 * @date 18.05.2025
 */

#pragma once

#include <Eigen/Dense>
#include <map>
#include <random>

#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"
#include "autopas/utils/Random.h"

namespace autopas {

/**
 * Deep Reinforcement Learning tuning strategy.
 *
 * This class implements a deep reinforcement learning approach to tune the configurations. Therefore, the method learns
 * from previously collected evidence and will use this data to predict the next best configurations. This method does
 * not use information about the domain.
 *
 * This implementation is based on the thesis "Low-Data Algorithm Selection in AutoPas"
 * (https://mediatum.ub.tum.de/node?id=1797277). The constant parameters were found through the parameter tuning in the
 * thesis.
 */
class DeepReinforcementLearning final : public TuningStrategyInterface {
 public:
  /**
   * Enum class for the different exploration methods.
   */
  enum class ExplorationMethod {
    /**
     * The polynomial exploration method.
     *
     * Choose the next configuration based on the weighted age squared and the predicted evidence:
     *
     * @f[age^2 \cdot _phaseScale + predictedEvidence()@f]
     *
     * The predicted evidence is an interpolation of the last two tuning phases.
     */
    polynomial,
    /**
     * The random exploration method. Choose the next exploration samples randomly.
     */
    random,
    /**
     * The longest ago exploration method. Choose the exploration samples based on the longest ago search.
     */
    longestAgo,
  };

 private:
  /**
   * Enum class for the different tuning states.
   */
  enum class TuningState {
    /**
     * The first search phase, where the algorithm is initialized. In this phase a full search is performed.
     */
    firstSearch,
    /**
     * The exploration phase, where the algorithm explores the search space, to gather information.
     */
    exploration,
    /**
     * The exploitation phase, where the algorithm exploits the information gained in the exploration phase.
     */
    exploitation,
  };

  /**
   * Container for the evidence history of a single configuration.
   */
  class DRLHistoryData final {
   public:
    /**
     * Constructor.
     * @param inputLength The number of evidence points to keep.
     * @param value The initial value of the history.
     * @param expectedResampleInterval The expected interval after which an algorithm gets re-sampled.
     */
    DRLHistoryData(const size_t inputLength, const double value, const size_t expectedResampleInterval);

    /**
     * Destructor.
     */
    ~DRLHistoryData() = default;

    /**
     * Add the evidence to the history.
     * @param value The value that should be added to the history.
     */
    void addValue(const double value);

    /**
     * Get the vector storing the evidence history.
     *
     * The history is stored in the following format: [value1 * scale, value2 * scale, ..., valueN * scale,
     * passedTuningPhase1, passedTuningPhase2, ..., passedTuningPhaseN]. This can be used as the input for the neural
     * network.
     * @param scale The scale factor to apply to the values.
     * @return The value vector.
     */
    [[nodiscard]] Eigen::Matrix<double, Eigen::Dynamic, 1> getValueVector(const double scale) const;

    /**
     * Age the history data by incrementing the number of tuning phases passed.
     */
    void ageData();

    /**
     * Get how many tuning phases ago the algorithm was last searched.
     * @return The number of tuning phases ago the algorithm was last searched.
     */
    [[nodiscard]] size_t lastSearched() const;

    /**
     * Get the predicted evidence.
     *
     * The value is computed based on linear interpolation of the last two tuning results.
     *
     * @return The predicted evidence.
     */
    [[nodiscard]] double predictedEvidence() const;

   private:
    /**
     * Store a random number generator.
     */
    static inline thread_local Random _rand = Random(0x5CD85F7E3F2262B3);

    /**
     * Store the values (evidences) of the last searches.
     */
    std::vector<double> _values;

    /**
     * Store the number of tuning phases passed until the current search.
     */
    std::vector<size_t> _tuningPhases;

    /**
     * The number of input values to keep in the history.
     */
    size_t _inputLength;
  };

 public:
  // Forbid the default constructor.
  DeepReinforcementLearning() = delete;

  /**
   * Constructor for the DeepReinforcementLearning class.
   * @param retrain Whether to retrain the neural network or not. (I.e. make this simple a deep learning strategy)
   * @param numExplorationSamples The number of exploration samples to use.
   * @param explorationMethod The exploration method to use.
   */
  explicit DeepReinforcementLearning(const bool retrain, const size_t numExplorationSamples = 3,
                                     const ExplorationMethod explorationMethod = ExplorationMethod::polynomial);

  /**
   * Destructor.
   */
  ~DeepReinforcementLearning() = default;

  /**
   * Get the type of the tuning strategy.
   * @return The deepReinforcementLearning type.
   */
  [[nodiscard]] TuningStrategyOption getOptionType() const override;

  /**
   * Add evidence to the history.
   * @param configuration The configuration for which the evidence is added.
   * @param evidence The evidence to be added.
   */
  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

  /**
   * Optimize the suggestions for the next iteration.
   *
   * This method will reorder and shrink the config queue based on the phase (add exploration or exploitation samples).
   *
   * @param configQueue The queue of configurations to be optimized.
   * @param evidenceCollection The collection of evidence to be used for the optimization.
   * @return Whether the queue was wiped intentionally.
   */
  bool optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

  /**
   * Reset the tuning process at the beginning of a tuning phase.
   * @param iteration The current iteration.
   * @param tuningPhase The current tuning phase.
   * @param configQueue The queue of configurations to be optimized.
   * @param evidenceCollection The collection of evidence to be used for the optimization.
   * @return Whether the queue was wiped intentionally.
   */
  bool reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const EvidenceCollection &evidenceCollection) override;

 private:
  /**
   * Get the exploration samples based on the _explorationMethod.
   * @return A set of configurations to be used for exploration.
   */
  [[nodiscard]] std::set<Configuration> getExplorationSamples() const;

  /**
   * Define the activation function to be used.
   *
   * For deep reinforcement learning, the ELU activation function is being used: @f$elu(x) = \max(0, x) + \min(0,
   * \alpha(\exp(x) - 1))@f$
   *
   * @param x The input to the activation function.
   * @return The output of the activation function.
   */
  [[nodiscard]] static inline constexpr double activationFunction(const double x) {
    // Use the ELU activation function
    return (x > 0) ? x : std::exp(x) - 1;
  }

  /**
   * Define the derivative of the activation function.
   *
   * This is the derivative of the ELU activation function: @f$elu'(x) = \max(0, 1) + \min(0, \exp(x))@f$
   *
   * @param x The input to the derivative.
   * @return The output of the derivative.
   */
  [[nodiscard]] static inline constexpr double activationFunctionDerivative(const double x) {
    // Derivative of the ELU activation function
    return (x > 0) ? 1 : std::exp(x);
  }

  /**
   * Evaluate the sample using the neural network.
   *
   * This value will return the expected runtime of the configuration.
   *
   * @param sample The sample to be evaluated.
   * @return The expected runtime of the configuration.
   */
  [[nodiscard]] double evaluate(const DRLHistoryData &sample) const;

  /**
   * Get the exploitation samples based on the current neural network.
   * @return A set of configurations to be used for exploitation.
   */
  [[nodiscard]] std::set<Configuration> getExploitationSamples() const;

  /**
   * Retrain the neural network. This is done based on the _retrainingInputData and _retrainingOutputData.
   */
  void retrainNeuralNetwork();

  /**
   * The random number generator.
   */
  static inline thread_local Random _rand = Random(0xDB91B2DB10E659CF);

  /**
   * The exploration method to be used.
   */
  ExplorationMethod _explorationMethod = ExplorationMethod::polynomial;

  /**
   * Store if the neural network should be retrained with reinforcement learning.
   */
  bool _retrain = false;

  /**
   * The configurations that have been explored.
   */
  std::set<Configuration> _exploredConfigurations;

  /**
   * The history of the search.
   */
  std::unordered_map<Configuration, DRLHistoryData, ConfigHash> _history;

  /**
   * The input length of the neural network. This is equivalent to the number of evidences stored in the history for a
   * particular configuration.
   */
  size_t _inputLength = 4;

  /**
   * The size of the hidden layer.
   */
  size_t _hiddenLayerSize = 16;

  /**
   * The hidden layer weights. This is a _inputLength * 2 x _hiddenLayerSize matrix for the hidden layer.
   */
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _hiddenLayer;

  /**
   * The hidden layer bias vector. This is a _hiddenLayerSize x 1 vector for the hidden layer biases.
   */
  Eigen::Matrix<double, Eigen::Dynamic, 1> _hiddenLayerBias;

  /**
   * The output layer weights. This is a _hiddenLayerSize x 1 vector for the output layer.
   */
  Eigen::Matrix<double, Eigen::Dynamic, 1> _outputLayer;

  /**
   * The output layer bias. This is a scalar value for the output layer bias.
   */
  double _outputLayerBias = 2.715298880101034151e-01;

  /**
   * The retraining input data. This is a vector of _inputLength x 1 matrices for the retraining input data.
   */
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> _retrainingInputData;

  /**
   * The retraining output data. This is a vector of scalar values for the retraining output data.
   */
  std::vector<double> _retrainingOutputData;

  /**
   * The number of retraining iterations.
   */
  size_t _retrainingIterations = 34;

  /**
   * The learning rate.
   */
  double _learningRate = 1.771665733218049020e-09;

  /**
   * The current state of the tuning process.
   */
  TuningState _state = TuningState::firstSearch;

  /**
   * The minimum evidence found during the full search.
   */
  double _minEvidence = std::numeric_limits<double>::infinity();

  /**
   * The number of samples to be used for exploration.
   */
  size_t _numExplorationSamples = 4;

  /**
   * The number of samples to be used for exploitation.
   */
  size_t _numExploitationSamples = 1;

  /**
   * The scaling factor used for the priority of the age of a data in polynomial exploration.
   *
   * The phase scale factor defines how much the priority depends on the squared age of the data. It is only used if the
   * polynomial method is being used for selecting the exploration samples.
   */
  double _phaseScale = 1.162197398643888948e-05;
};
}  // namespace autopas
