/**
 * @file DeepReinforcementLearning.cpp
 * @author P. Metscher
 * @date 18.05.2025
 */

#include "DeepReinforcementLearning.h"

#include "autopas/utils/logging/Logger.h"

autopas::DeepReinforcementLearning::DRLHistoryData::DRLHistoryData(const size_t inputLength, const double value,
                                                                   const size_t expectedResampleInterval) {
  // Initialize the vectors with the given input length
  _inputLength = inputLength;

  // Initialize the values
  _values.resize(_inputLength);
  _values.assign(_inputLength, value);

  // Initialize the tuning phases
  _tuningPhases.resize(_inputLength);
  std::normal_distribution<> dist(expectedResampleInterval, expectedResampleInterval / 4);
  _tuningPhases[0] = 0;

  for (size_t i = 1; i < _inputLength; i++) {
    _tuningPhases[i] = _tuningPhases[i - 1] + dist(_rand);
  }
}

void autopas::DeepReinforcementLearning::DRLHistoryData::addValue(const double value) {
  // Rotate everything back
  for (size_t i = _inputLength - 1; i > 0; i--) {
    _values[i] = _values[i - 1];
    _tuningPhases[i] = _tuningPhases[i - 1];
  }

  // Add new value
  _values[0] = value;
  _tuningPhases[0] = 0;
}

void autopas::DeepReinforcementLearning::DRLHistoryData::ageData() {
  // Increment tuning phases to age the data
  for (size_t i = 0; i < _inputLength; i++) {
    _tuningPhases[i]++;
  }
}

size_t autopas::DeepReinforcementLearning::DRLHistoryData::lastSearched() const { return _tuningPhases[0]; }

double autopas::DeepReinforcementLearning::DRLHistoryData::predictedEvidence() const {
  return _values[0] - _tuningPhases[0] * (_values[1] - _values[0]) / (_tuningPhases[1] - _tuningPhases[0]);
}

Eigen::Matrix<double, Eigen::Dynamic, 1> autopas::DeepReinforcementLearning::DRLHistoryData::getValueVector(
    const double scale) const {
  std::vector<double> ret;

  // Construct the data vector as follows:
  // [value1, value2, ..., valueN, tuningPhase1, tuningPhase2, ..., tuningPhaseN]
  // where N is the inputLength.
  for (size_t i = 0; i < _inputLength; i++) {
    ret.emplace_back(_values[i] * scale);
  }

  for (size_t i = 0; i < _inputLength; i++) {
    ret.emplace_back(_tuningPhases[i]);
  }

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec = Eigen::Map<Eigen::VectorXd>(ret.data(), ret.size());
  return vec;
}

autopas::DeepReinforcementLearning::DeepReinforcementLearning(const std::set<Configuration> &searchSpace,
                                                              const bool train, const size_t explorationSamples,
                                                              const ExplorationMethod explorationMethod)
    : _explorationMethod(explorationMethod),
      _train(train),
      _searchSpace(searchSpace),
      _explorationSamples(explorationSamples) {
  // Check if the search space is empty.
  if (_searchSpace.empty()) {
    utils::ExceptionHandler::exception("DeepReinforcementLearning: The searchspace may not be empty.");
  }

  // Check if the exploration samples are not too many
  if (_explorationSamples + _exploitationSamples >= _searchSpace.size()) {
    utils::ExceptionHandler::exception(
        "DeepReinforcementLearning: The exploration and exploitation samples may not be greater than the search space "
        "size.");
  }

  // Test that there is at least one exploration sample
  if (_explorationSamples < 1) {
    utils::ExceptionHandler::exception("DeepReinforcementLearning: The exploration samples may not be less than 1.");
  }

  // clang-format off
  _hiddenLayer = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(_inputLength * 2, _hiddenLayerSize);
  _hiddenLayer <<
      7.168054529013406784e-01, 7.168040585285946967e-03, 2.704026286165521209e-01, 4.630724881088673500e-03, 1.687074397039885210e-01, 1.585413257664988751e-01, 1.760893397538867422e-01, 5.590097621677001527e-01, 1.050048910109433742e-01, 8.255051644726370610e-03, -7.082428308776984004e-03, -6.861423195762151456e-03, -7.394085678574022213e-02, -1.014808086514065039e-02, -6.885103686948274362e-03, 4.554077184914036465e-03,
      1.255160466506496063e-02, 4.570199778548201218e-03, 5.195389527127776234e-03, 8.720365495863871552e-03, -1.061224377335888815e-03, -2.245135679628698904e-03, -2.402129639579072672e-03, 1.013535188850097669e-02, 1.541235646191406765e-03, -8.175324462705345699e-03, -6.268562710917340713e-03, -7.532738693801125686e-03, 3.261654187958517789e-03, -1.904582227782710161e-03, 5.067618028809222665e-03, 5.320846747132212604e-03,
      -1.933812015150767824e-02, -4.072587613029335536e-03, -8.593698148292543484e-03, -3.896110854711215803e-03, 1.009982038397928307e-02, 1.171835333860572238e-02, 4.760041874346415702e-04, -2.995534291799160573e-02, 1.535694586485886019e-02, 6.687425050013228248e-03, 2.209119530600127477e-03, -9.215260312158383852e-04, -6.646528482441664814e-04, 5.537976887816015159e-03, 4.684117594951018801e-03, 3.248904098057020344e-03,
      -1.939117861948497917e-02, 1.852853689618171794e-05, 1.042208088381155716e-02, -2.464591078294206539e-03, 2.521341489924397786e-02, 1.990424454435281085e-02, 1.822782233966340273e-02, -1.637183458252445406e-02, 1.061672227881722309e-02, 1.870169973569269775e-03, -6.764724490092286803e-03, 2.749449175137412115e-03, -8.309988717452769860e-03, 4.452278143525482783e-03, -3.264791007548906773e-03, 2.941309560030008936e-03,
      -9.113462877831761400e-04, -7.774784026942291142e-03, -1.080761851494165720e-03, 5.983924625579903908e-04, -1.953424542214928818e-03, 6.643152323881822661e-03, -9.064919407502492962e-03, 1.047234792270523894e-03, 9.544400149830324989e-03, 3.650441089737995175e-03, -6.367477597637214749e-03, 8.213278670087788280e-03, -3.674876647716281749e-03, -2.872067547698399503e-03, 8.489885224207271638e-03, 6.371575835637979801e-03,
      6.228553160500610803e-04, -4.831682220226504905e-03, 1.673969917514481849e-04, -4.090232047752387075e-04, 1.740639394486662563e-05, -2.836681311814974633e-04, 1.943580839178212253e-03, 7.010354482840374128e-04, 3.416293078557566636e-03, 6.765815878390201397e-03, -4.829020878320208071e-03, -7.935035298489160477e-03, 1.514479651106706570e-02, 9.398901511109950674e-03, 5.082480381349336988e-03, -6.741953471496572869e-03,
      -4.458449655766699802e-04, 2.367484330904916310e-03, 8.110222568302746056e-04, 1.690516230758963567e-04, 1.120366953742010438e-02, -4.158466489582061075e-03, 2.589206703273883286e-03, -2.417770889894326272e-03, 4.024996616782129145e-03, -8.486769132638146904e-03, 9.074710982512930110e-03, 3.958934413415912505e-03, 1.039974959705559053e-02, -8.982923583128720954e-03, 1.004297451722017129e-03, 6.408255850043616079e-04,
      -1.531084020497015374e-03, -8.919777177402516549e-03, -1.012026627843938152e-03, -9.103637899890491478e-03, -4.269797253534059721e-03, 4.771564812481082253e-03, 8.355063093490230561e-03, -2.330336218393410633e-04, 2.361900839276299610e-03, -6.936814157230741795e-03, -9.013219541378972299e-03, -4.830965752382156142e-03, 1.169971873925122578e-02, 7.931049502331625426e-04, -2.196480952793811517e-03, -9.612086438292812354e-03;
  // clang-format on

  _hiddenLayerBias = Eigen::Matrix<double, Eigen::Dynamic, 1>(_hiddenLayerSize);
  _hiddenLayerBias <<             //
      -3.135705915255077806e-03,  //
      -6.243448213205700023e-03,  //
      -1.717444889700283110e-03,  //
      -8.942428778361341188e-03,  //
      -1.017224253555813751e-02,  //
      -8.893363991131367566e-03,  //
      -7.355436019233416861e-03,  //
      -1.455493202264515964e-02,  //
      4.594051174181832732e-04,   //
      -8.844623220133039643e-03,  //
      7.804188427476882982e-03,   //
      8.943758451508279986e-03,   //
      -7.804225827731178441e-03,  //
      8.528590591742253244e-03,   //
      -5.446880809276771979e-04,  //
      -4.879841284888828887e-03;  //

  _outputLayer = Eigen::Matrix<double, Eigen::Dynamic, 1>(_hiddenLayerSize);
  _outputLayer <<                 //
      7.145312041908945533e-01,   //
      -3.894996838784765131e-03,  //
      2.695768165356492996e-01,   //
      4.033694103345464451e-04,   //
      1.706423712976079532e-01,   //
      1.601416487127459154e-01,   //
      1.764293396649790624e-01,   //
      5.577547497622016293e-01,   //
      1.060730897150632401e-01,   //
      -6.040194951323567263e-03,  //
      -6.991749510121741951e-03,  //
      -3.582789992043560866e-03,  //
      -7.673545289944561087e-02,  //
      -1.237576000350001058e-02,  //
      -4.503395418219970255e-03,  //
      -3.863324328016031774e-04;  //
}

autopas::TuningStrategyOption autopas::DeepReinforcementLearning::getOptionType() const {
  return TuningStrategyOption::deepReinforcementLearning;
}

void autopas::DeepReinforcementLearning::addEvidence(const Configuration &configuration, const Evidence &evidence) {
  switch (_state) {
    case TuningState::firstSearch: {
      if (_searchSpace.find(configuration) == _searchSpace.end()) {
        utils::ExceptionHandler::exception("DeepReinforcementLearning: Configuration not found in search space.");
      }

      // In the first search state, we simply collect evidence without any further processing.
      _history.emplace(configuration, DRLHistoryData(_inputLength, evidence.value,
                                                     static_cast<double>(_searchSpace.size()) / _explorationSamples));
      _minEvidence = std::min(_minEvidence, static_cast<double>(evidence.value));
      break;
    }

    case TuningState::exploration:
    case TuningState::exploitation: {
      // In the exploration and exploitation state, we collect evidence and update the exploration samples.
      _trainingInputData.push_back(_history.at(configuration).getValueVector(_minEvidence));
      _trainingOutputData.push_back(static_cast<double>(evidence.value));
      _history.at(configuration).addValue(evidence.value);

      break;
    }

    default: {
      utils::ExceptionHandler::exception("DeepReinforcementLearning: Unknown tuning state.");
      break;
    }
  }
}

bool autopas::DeepReinforcementLearning::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                                             const EvidenceCollection &evidenceCollection) {
  switch (_state) {
    case TuningState::firstSearch: {
      // During the first search, we simply perform a full search.
      if (configQueue.empty()) {
        if (_searchSpace.size() != _history.size()) {
          utils::ExceptionHandler::exception(
              "DeepReinforcementLearning: The search space size does not match the history size. Not every "
              "configuration was searched.");
        }
      }

      return true;
    }

    case TuningState::exploration: {
      // Code for optimizing suggestions in the exploration state

      // Select the exploration samples if the config queue is larger than the exploration samples
      if (configQueue.size() > _explorationSamples) {
        std::set<Configuration> samples = getExplorationSamples();
        configQueue.clear();

        for (const auto &sample : samples) {
          // If the sample is not in the history, add it to the config queue.
          configQueue.push_back(sample);
        }
      }

      // Search the exploration samples
      if (!configQueue.empty()) {
        return false;
      }

      // Code for optimizing suggestions in the training state
      if (_train) {
        trainNeuralNetwork();
      }

      // Clear the training data
      _trainingInputData.clear();
      _trainingOutputData.clear();

      // Select the exploitation samples
      std::set<Configuration> samples = getExploitationSamples();
      configQueue.clear();

      for (const auto &sample : samples) {
        // If the sample is not in the history, add it to the config queue.
        configQueue.push_back(sample);
      }

      // Transition to the exploitation if the training has finished
      _state = TuningState::exploitation;
      return true;
    }

    case TuningState::exploitation: {
      // Search all exploitation samples and finish the search
      return true;
    }

    default: {
      utils::ExceptionHandler::exception("DeepReinforcementLearning: Unknown tuning state.");
      return false;
    }
  }
}

bool autopas::DeepReinforcementLearning::reset(size_t iteration, size_t tuningPhase,
                                               std::vector<Configuration> &configQueue,
                                               const EvidenceCollection &evidenceCollection) {
  // Validate the input
  if (configQueue.size() != _searchSpace.size()) {
    utils::ExceptionHandler::exception(
        "DeepReinforcementLearning: The searchspace size does not match the config queue size.");
    return false;
  }

  // Reset some internal parameters
  if (tuningPhase == 0) {
    _state = TuningState::firstSearch;
  } else {
    _state = TuningState::exploration;
  }

  // Age the history table
  for (auto &[config, history] : _history) {
    history.ageData();
  }

  return false;
}

std::set<autopas::Configuration> autopas::DeepReinforcementLearning::getExplorationSamples() const {
  if (_explorationSamples > _searchSpace.size()) {
    utils::ExceptionHandler::exception("DeepReinforcementLearning: Exploration samples exceed search space size.");
    return _searchSpace;  // Return the entire search space if samples exceed size
  }

  switch (_explorationMethod) {
    case ExplorationMethod::polynomial: {
      // Implementation for polynomial selector exploration
      std::vector<std::pair<Configuration, double>> explorationPriority;
      double timeScale = 1.0 / _minEvidence;
      double phaseScale = 1.0 / _searchSpace.size();

      for (const Configuration &config : _searchSpace) {
        double priority = 0.0;

        // Calculate the priority based on the last tuning phase
        double lastTuningPhase = _history.at(config).lastSearched() * phaseScale;

        priority += lastTuningPhase * lastTuningPhase * _phaseScale;
        priority -= _history.at(config).predictedEvidence() * timeScale;

        explorationPriority.emplace_back(config, priority);
      }

      // Sort the configuration queue
      std::sort(explorationPriority.begin(), explorationPriority.end(),
                [](const auto &a, const auto &b) { return a.second > b.second; });

      // Generate the subset
      std::set<autopas::Configuration> subSet;

      for (size_t i = 0; i < _explorationSamples; i++) {
        subSet.emplace(explorationPriority[i].first);
      }

      return subSet;
    }

    case ExplorationMethod::random: {
      // Implementation for random exploration
      return _rand.randomSubset(_searchSpace, _explorationSamples);
    }

    case ExplorationMethod::longestAgo: {
      // Implementation for longest ago exploration
      std::vector<std::pair<Configuration, size_t>> sortedHistory;
      for (const auto &[config, history] : _history) {
        sortedHistory.emplace_back(config, history.lastSearched());
      }

      // Shuffle the array to prevent order based preferences
      std::shuffle(sortedHistory.begin(), sortedHistory.end(), _rand);
      std::sort(sortedHistory.begin(), sortedHistory.end(),
                [](const auto &a, const auto &b) { return a.second < b.second; });

      std::set<Configuration> ret;
      for (size_t i = 0; i < _explorationSamples; i++) {
        ret.emplace(sortedHistory[i].first);
      }

      return ret;
    }

    default: {
      utils::ExceptionHandler::exception("DeepReinforcementLearning: Unknown exploration method.");
      return _rand.randomSubset(_searchSpace, _explorationSamples);  // Fallback to random exploration
    }
  }
}

double autopas::DeepReinforcementLearning::evaluate(const DRLHistoryData &sample) const {
  // Prevent theoretical errors, if the execution is extremely fast and the timing sensor has a low resolution
  double scale = (_minEvidence > 0.0) ? (1.0 / _minEvidence) : 1.0;

  // Evaluate the sample using the neural network
  Eigen::Matrix<double, 1, Eigen::Dynamic> input = sample.getValueVector(scale).transpose();
  auto hiddenOut = input * _hiddenLayer + _hiddenLayerBias.transpose();
  auto activatedHidden = hiddenOut.unaryExpr([](const double x) { return activationFunction(x); });
  double output = _outputLayer.dot(activatedHidden) + _outputLayerBias;

  if (output <= 0.0) [[unlikely]] {
    AutoPasLog(DEBUG, "DeepReinforcementLearning: Neural Network output is {} (this value should always be > 0.0).",
               output);
  }

  return output;
}

std::set<autopas::Configuration> autopas::DeepReinforcementLearning::getExploitationSamples() const {
  // Handle invalid input sizes
  if (_exploitationSamples > _searchSpace.size()) {
    utils::ExceptionHandler::exception("DeepReinforcementLearning: Exploitation samples exceed search space size.");
    return _searchSpace;  // Return the entire search space if samples exceed size
  }

  if (_history.size() != _searchSpace.size()) {
    utils::ExceptionHandler::exception("DeepReinforcementLearning: History size does not match search space size.");
    return {};  // Return an empty set if history and search space size do not match
  }

  // Collect the exploitation samples based on the history
  std::set<Configuration> exploitationSamples;
  std::vector<std::pair<Configuration, double>> sortedConfigs;

  // Choose the best configurations based on the minimum evaluate(history[config])
  for (const auto &[config, DRLHistoryData] : _history) {
    double score = evaluate(DRLHistoryData);
    sortedConfigs.emplace_back(config, score);
  }

  // Sort the configurations by their scores
  std::sort(sortedConfigs.begin(), sortedConfigs.end(),
            [](const auto &a, const auto &b) { return a.second < b.second; });

  // Select the top configurations (no overflow possible due to the input validation)
  for (size_t i = 0; i < _exploitationSamples; i++) {
    const Configuration config = sortedConfigs[i].first;

    // Prevent duplicate searches
    if (_history.at(config).lastSearched() == 0) {
      continue;
    }

    exploitationSamples.insert(config);
  }

  return exploitationSamples;
}

void autopas::DeepReinforcementLearning::trainNeuralNetwork() {
  for (size_t i = 0; i < _trainingIterations; i++) {
    for (size_t j = 0; j < _trainingInputData.size(); j++) {
      // Forward propagate
      const auto &input = _trainingInputData[j];
      double target = _trainingOutputData[j];

      // Hidden layer: input * weights + bias
      Eigen::Matrix<double, 1, Eigen::Dynamic> hiddenInput = input.transpose() * _hiddenLayer;
      hiddenInput += _hiddenLayerBias.transpose();

      // Activation
      Eigen::Matrix<double, 1, Eigen::Dynamic> hiddenOutput =
          hiddenInput.unaryExpr([](const double x) { return activationFunction(x); });

      // Output layer: hidden * weights + bias
      double output = hiddenOutput * _outputLayer + _outputLayerBias;

      // Compute loss and gradients (mean squared error)
      double error = output - target;
      double dLoss_dOutput = 2.0 * error;

      // Gradients for output layer
      Eigen::Matrix<double, Eigen::Dynamic, 1> dOutputLayer = hiddenOutput.transpose() * dLoss_dOutput;
      double dOutputLayerBias = dLoss_dOutput;

      // Gradients for hidden layer
      Eigen::Matrix<double, Eigen::Dynamic, 1> dHidden = _outputLayer * dLoss_dOutput;
      Eigen::Matrix<double, Eigen::Dynamic, 1> dHiddenActivation =
          hiddenInput.unaryExpr([](const double x) { return activationFunctionDerivative(x); })
              .transpose()
              .cwiseProduct(dHidden);

      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dHiddenLayer = input * dHiddenActivation.transpose();
      Eigen::Matrix<double, Eigen::Dynamic, 1> dHiddenLayerBias = dHiddenActivation;

      // Update weights and biases
      _outputLayer -= _learningRate * dOutputLayer;
      _outputLayerBias -= _learningRate * dOutputLayerBias;
      _hiddenLayer -= _learningRate * dHiddenLayer;
      _hiddenLayerBias -= _learningRate * dHiddenLayerBias;
    }
  }
}
