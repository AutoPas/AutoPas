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
  std::normal_distribution<> dist(static_cast<double>(expectedResampleInterval), static_cast<double>(expectedResampleInterval) / 4);
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
  for (auto &tuningPhase : _tuningPhases) {
    tuningPhase++;
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

autopas::DeepReinforcementLearning::DeepReinforcementLearning(const bool train, const size_t explorationSamples,
                                                              const ExplorationMethod explorationMethod)
    : _explorationMethod(explorationMethod), _train(train), _explorationSamples(explorationSamples) {
  // Test that there is at least two exploration sample
  if (_explorationSamples <= 1) {
    utils::ExceptionHandler::exception("DeepReinforcementLearning: The exploration samples may not be less than 2.");
  }

  // clang-format off
  _hiddenLayer = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(_inputLength * 2, _hiddenLayerSize);
  _hiddenLayer <<
      -4.445636816225461557e-03, -1.665181529429264964e-02,  3.948844164874587004e-01, -3.328680133656577306e-02,  4.761291317370201781e-02,  4.812055042532794973e-03, -1.158019336760089495e-02,  9.575951123575705826e-02,  8.010436488579358860e-01,  8.848617998091938552e-03,  1.278131294927795369e-01,  5.510907005489978222e-02, -1.242152981509650131e-02, -6.686573178232044223e-02,  3.804136883699282649e-03,  3.915208173370351985e-01, //
      -7.633845370950443320e-03,  3.978045070409219619e-03,  6.296137947190323832e-03, -1.663389489684525754e-02,  4.423323635143023951e-03, -1.102122305391365145e-02, -3.103517540944439063e-03, -6.869874467793332488e-03,  2.179243814483783001e-02,  4.113896755473781322e-03, -9.563278360982070334e-03, -3.263499921348785677e-04, -1.038390640175560278e-03,  7.929403181207262535e-03,  1.365657500471845328e-03,  1.028319828419055511e-02, //
      -1.292413095196153633e-02, -1.263909335094434878e-02,  1.286477690600031781e-02, -1.929819712039636953e-02, -8.874919772069408616e-03,  2.917594502946532862e-03,  4.366148983805265803e-03,  1.956774969432472332e-03, -2.230559460514659753e-02,  7.221623609717045013e-03,  3.137178062620358133e-03, -3.766609012297297918e-03, -1.894852545415773564e-03, -2.508157605015383070e-04,  6.603336717711409470e-05, -3.665467449024812068e-03, //
      -1.576019160555879386e-02, -1.710650772634293975e-03,  2.028192063958074096e-02, -3.100664541782366390e-02, -4.698000194671637236e-03, -1.268033390795530665e-02,  2.542887179199748108e-03,  8.831480046930724123e-03, -1.226518402303907254e-02,  8.907365855020788967e-03, -1.028455799770529484e-03,  8.551982609435414967e-03, -1.428357308665839133e-02, -1.033772512507194138e-02, -3.729005237113039786e-03, -3.746236648285931024e-03, //
      -1.929572162441826888e-02, -2.427980440197391107e-02,  1.675352909002334302e-03, -5.680269960641709526e-02,  9.396822531984412238e-03, -1.914754190642114356e-02, -9.650801378095440891e-03, -5.323402991063825272e-03, -1.095233639560363147e-02,  3.687400089823864024e-03,  2.708654535713080243e-03, -7.035538904511221166e-03, -1.857961562538220160e-02, -1.979691965928451810e-03,  5.933508868737029894e-03,  2.451727045728428755e-03, //
      -7.150324393665755812e-03,  1.019052950994035929e-02, -4.570639988442704894e-03, -1.889589086578598353e-03, -6.743227042709448033e-03,  8.528905169081408416e-03,  1.686385683996903051e-03,  4.560863548836005359e-03,  4.082837679598005226e-03, -9.626381286583923520e-04, -7.774081937877609302e-03, -1.217637883107271267e-03, -1.348896594022675323e-04,  6.643555883739759964e-03, -6.957367419675847785e-03,  3.523583180287542002e-04, //
       3.223803392012643612e-03, -1.297416814299122726e-03,  3.008832125304946020e-03,  1.687387502165784217e-02,  1.171719431598491033e-02,  5.042003770513094063e-03,  3.039948669642006124e-03, -1.525123583403633753e-03,  1.166275754304003411e-03, -7.522049706793341138e-03,  6.792187393978685490e-03, -2.422812678573114248e-03,  9.688317331526804382e-03,  1.218945308803522608e-02,  4.512620676407223808e-03, -1.011106520064284937e-03, //
       1.686545305362980099e-02,  9.298976521838031584e-03,  6.157913372749283798e-03,  2.873757871656030735e-02,  1.744541763427530021e-03,  5.301813342739749397e-04,  1.666203836892147594e-03,  3.489948399435068883e-03,  6.389972098075923958e-04,  4.396764635958176799e-03,  4.157415765306895247e-04,  1.030690619294660076e-02,  3.891572171006626221e-03,  7.176661029600370743e-03, -2.603046562777876737e-03, -1.463776074465084455e-03; //
  // clang-format on

  _hiddenLayerBias = Eigen::Matrix<double, Eigen::Dynamic, 1>(_hiddenLayerSize);
  _hiddenLayerBias <<             //
      5.264647513546491693e-03,   //
      -2.940947794651221740e-03,  //
      -3.608824054843444082e-03,  //
      -4.532957276108481350e-03,  //
      7.385659791226325925e-03,   //
      3.993377577912109871e-03,   //
      3.368151622040083722e-03,   //
      -6.890303195688882165e-03,  //
      -1.643420056105823063e-02,  //
      -8.219386138704772027e-03,  //
      -7.537791033706692234e-03,  //
      6.454158585728376162e-03,   //
      4.194664315279335483e-03,   //
      -6.637979031233082541e-03,  //
      -7.640891325529360761e-03,  //
      -1.112110321001389367e-02;  //

  _outputLayer = Eigen::Matrix<double, Eigen::Dynamic, 1>(_hiddenLayerSize);
  _outputLayer <<                 //
      -3.295714186738026019e-02,  //
      -3.050614890021472944e-02,  //
      3.948492470871525772e-01,   //
      -8.338086645962415611e-02,  //
      4.755454493101033547e-02,   //
      -2.303950505755288092e-02,  //
      -1.941024036569514652e-02,  //
      9.551460402974001895e-02,   //
      7.998975727746606701e-01,   //
      -1.189319968520019359e-02,  //
      1.278703742170215485e-01,   //
      5.501077709947461286e-02,   //
      -2.952029626891259501e-02,  //
      -6.819006630639656841e-02,  //
      -6.963486161686841501e-03,  //
      3.907018762648295507e-01;   //
}

autopas::TuningStrategyOption autopas::DeepReinforcementLearning::getOptionType() const {
  return TuningStrategyOption::deepReinforcementLearning;
}

void autopas::DeepReinforcementLearning::addEvidence(const Configuration &configuration, const Evidence &evidence) {
  switch (_state) {
    case TuningState::firstSearch: {
      // In the first search state, we simply collect evidence without any further processing.
      _searchSpace.insert(configuration);
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
      _history.at(configuration).addValue(static_cast<double>(evidence.value));

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
      return true;
    }

    case TuningState::exploration: {
      // Code for optimizing suggestions in the exploration state

      // Select the exploration samples if the config queue is larger than the exploration samples
      if (configQueue.size() > _explorationSamples) {
        const std::set<Configuration> samples = getExplorationSamples();
        configQueue.clear();

        for (const auto &sample : samples) {
          // If the sample is not in the history, add it to the config queue.
          configQueue.push_back(sample);
        }
      }

      // Search the exploration samples
      if (not configQueue.empty()) {
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
      const std::set<Configuration> samples = getExploitationSamples();
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
  // Reset some internal parameters
  if (tuningPhase == 0) {
    _state = TuningState::firstSearch;
  } else {
    _state = TuningState::exploration;
  }

  // Age the history table
  for (auto &history: _history | std::views::values) {
    history.ageData();
  }

  // The first search may never contain an empty config queue.
  optimizeSuggestions(configQueue, evidenceCollection);

  return false;
}

std::set<autopas::Configuration> autopas::DeepReinforcementLearning::getExplorationSamples() const {
  switch (_explorationMethod) {
    case ExplorationMethod::polynomial: {
      // Implementation for polynomial selector exploration
      std::vector<std::pair<Configuration, double>> explorationPriority;
      double timeScale = 1.0 / _minEvidence;

      for (const Configuration &config : _searchSpace) {
        // Calculate the priority based on when the configuration was last tested.
        double lastTuningPhase = _history.at(config).lastSearched();
        
        const auto priority = lastTuningPhase * lastTuningPhase * _phaseScale - _history.at(config).predictedEvidence() * timeScale;

        explorationPriority.emplace_back(config, priority);
      }

      // Sort the configuration queue
      std::ranges::sort(explorationPriority,
            [](const auto &a, const auto &b) { return a.second > b.second; });

      // Generate the subset
      std::set<Configuration> subset;

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

      // Shuffle the array to prevent order-based preferences in case of multiple algorithms with the same age.
      std::ranges::shuffle(sortedHistory, _rand);
      std::ranges::sort(sortedHistory, [](const auto &a, const auto &b) { return a.second < b.second; });

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
