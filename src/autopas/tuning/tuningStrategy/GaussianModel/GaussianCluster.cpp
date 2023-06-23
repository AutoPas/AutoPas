/**
 * @file GaussianCluster.cpp
 * @author F. Gratl
 * @date 17.11.22
 */

#include "GaussianCluster.h"

autopas::GaussianCluster::GaussianCluster(const std::vector<int> &dimRestriction, size_t continuousDims,
                                          autopas::GaussianCluster::WeightFunction weightFun, double sigma,
                                          Random &rngRef, const GaussianModelTypes::VectorToStringFun &vectorToString,
                                          const std::string &outputSuffix)
    : _dimRestriction(dimRestriction),
      _continuousDims(continuousDims),
      _clusters(),
      _weightFun(weightFun),
      _evidenceMinValue(0),
      _evidenceMaxValue(0),
      _numEvidence(0),
      _sigma(sigma),
      _rng(rngRef),
      _logger(std::make_unique<GaussianClusterLogger>(vectorToString, outputSuffix)) {
  initClusters();
}

autopas::GaussianCluster::~GaussianCluster() = default;

const std::vector<int> &autopas::GaussianCluster::getDimensions() const { return _dimRestriction; }

void autopas::GaussianCluster::setDimensions(const std::vector<int> &newValue) {
  _dimRestriction = newValue;
  initClusters();
}

const autopas::GaussianProcess &autopas::GaussianCluster::getCluster(size_t index1D) const {
  return _clusters[index1D];
}

void autopas::GaussianCluster::clear() {
  for (auto cluster : _clusters) {
    cluster.clear();
  }
  _numEvidence = 0;
}

size_t autopas::GaussianCluster::numEvidence() const { return _numEvidence; }

void autopas::GaussianCluster::addEvidence(const autopas::GaussianModelTypes::VectorDiscrete &inputDiscrete,
                                           const autopas::GaussianModelTypes::VectorContinuous &inputContinuous,
                                           double output) {
  size_t clusterIdx = getIndex(inputDiscrete);
  if (static_cast<size_t>(inputContinuous.size()) != _continuousDims) {
    utils::ExceptionHandler::exception(
        "GaussianCluster: size of continuous input {} does not match specified dimensions {}", inputContinuous.size(),
        _continuousDims);
  }

  if (_numEvidence == 0) {
    // first evidence
    _evidenceMinValue = _evidenceMaxValue = output;
    _evidenceMaxVector = std::make_pair(inputDiscrete, inputContinuous);
  } else if (output < _evidenceMinValue) {
    _evidenceMinValue = output;
  } else if (output > _evidenceMaxValue) {
    _evidenceMaxValue = output;
    _evidenceMaxVector = std::make_pair(inputDiscrete, inputContinuous);
  }

  _clusters[clusterIdx].addEvidence(inputContinuous, output, false);
  ++_numEvidence;

  auto [sample_means, sample_thetas, sample_dimScales] = GaussianProcess::generateHyperparameterSamples(
      hp_sample_size, _rng, _continuousDims, _sigma, _evidenceMinValue, _evidenceMaxValue);

  // set hyperparameter of each cluster
  for (auto &cluster : _clusters) {
    cluster.setHyperparameters(sample_means, sample_thetas, sample_dimScales);
  }

  // combine score of each cluster
  for (size_t i = 0; i < hp_sample_size; ++i) {
    // set the score of a hyperparameter sample to the product of all scores
    double combinedScore = 0.;
    for (auto &cluster : _clusters) {
      combinedScore += cluster.getHyperparameters()[i].score;
    }

    if (std::isnan(combinedScore) or std::isinf(combinedScore)) {
      utils::ExceptionHandler::exception("GaussianCluster: Score of hyperparameter is {}", combinedScore);
    }

    for (auto &cluster : _clusters) {
      cluster.getHyperparameters()[i].score = combinedScore;
    }
  }

  // normalize all hyperparameters
  for (auto &cluster : _clusters) {
    cluster.normalizeHyperparameters();
  }
}

void autopas::GaussianCluster::addEvidence(const autopas::GaussianModelTypes::VectorPairDiscreteContinuous &input,
                                           double output) {
  addEvidence(input.first, input.second, output);
}

autopas::GaussianModelTypes::VectorPairDiscreteContinuous autopas::GaussianCluster::getEvidenceMax() const {
  if (_numEvidence == 0) {
    utils::ExceptionHandler::exception("GaussianCluster has no evidence");
  }

  return _evidenceMaxVector;
}

std::vector<autopas::GaussianModelTypes::VectorAcquisition> autopas::GaussianCluster::sampleAcquisition(
    autopas::AcquisitionFunctionOption af, const autopas::GaussianModelTypes::NeighbourFunction &neighbourFun,
    const std::vector<GaussianModelTypes::VectorContinuous> &continuousSamples) const {
  std::vector<GaussianModelTypes::VectorAcquisition> acquisitions;
  // pair up all discrete with all continuous samples
  acquisitions.reserve(_clusters.size() * continuousSamples.size());

  auto neighbourWeights = initNeighbourWeights(neighbourFun);

  for (const auto &continuousSample : continuousSamples) {
    auto [means, vars, stddevs] = precalculateDistributions(continuousSample);
    updateNeighbourWeights(neighbourWeights, means, vars, stddevs);
    auto currentAcquisisitions = precalculateAcquisitions(af, means, vars);

    _logger->add(_clusters, _discreteVectorMap, continuousSample, means, vars, neighbourWeights);

    // get acquisition considering neighbours
    for (size_t i = 0; i < _clusters.size(); ++i) {
      // target cluster gets weight 1.
      double mixAcquisition = currentAcquisisitions[i];
      double weightSum = 1.;

      // weighted sum over neighbours
      for (const auto &[n, _, weight] : neighbourWeights[i]) {
        mixAcquisition += currentAcquisisitions[n] * weight;
        weightSum += weight;
      }

      // normalize
      mixAcquisition /= weightSum;

      acquisitions.emplace_back(std::make_pair(_discreteVectorMap[i], continuousSample), mixAcquisition);
    }
  }

  _logger->flush();

  return acquisitions;
}

void autopas::GaussianCluster::logDebugGraph(
    const autopas::GaussianModelTypes::NeighbourFunction &neighbourFun,
    const std::vector<GaussianModelTypes::VectorContinuous> &continuousSamples) const {
  if (autopas::Logger::get()->level() > autopas::Logger::LogLevel::trace) {
    return;
  }

  auto neighbourWeights = initNeighbourWeights(neighbourFun);
  for (const auto &continuousSample : continuousSamples) {
    auto [means, vars, stddevs] = precalculateDistributions(continuousSample);
    updateNeighbourWeights(neighbourWeights, means, vars, stddevs);

    _logger->add(_clusters, _discreteVectorMap, continuousSample, means, vars, neighbourWeights);
  }

  _logger->flush();
}

void autopas::GaussianCluster::setVectorToStringFun(const autopas::GaussianModelTypes::VectorToStringFun &fun) {
  _logger->setVectorToStringFun(fun);
}

autopas::GaussianModelTypes::VectorAcquisition autopas::GaussianCluster::sampleAcquisitionMax(
    autopas::AcquisitionFunctionOption af, const autopas::GaussianModelTypes::NeighbourFunction &neighbourFun,
    const std::vector<GaussianModelTypes::VectorContinuous> &continuousSamples) const {
  // generate all combinations and acquisitions
  auto acquisitions = sampleAcquisition(af, neighbourFun, continuousSamples);

  // return maximum
  return *std::max_element(
      acquisitions.begin(), acquisitions.end(),
      [](const GaussianModelTypes::VectorAcquisition &va1, const GaussianModelTypes::VectorAcquisition &va2) {
        auto acquisition1 = va1.second;
        auto acquisition2 = va2.second;
        return acquisition1 < acquisition2;
      });
}

std::vector<autopas::GaussianModelTypes::VectorPairDiscreteContinuous>
autopas::GaussianCluster::sampleOrderedByAcquisition(
    autopas::AcquisitionFunctionOption af, const autopas::GaussianModelTypes::NeighbourFunction &neighbourFun,
    const std::vector<GaussianModelTypes::VectorContinuous> &continuousSamples) const {
  // generate all combinations and acquisitions
  auto acquisitions = sampleAcquisition(af, neighbourFun, continuousSamples);

  // order by acquisition
  std::sort(acquisitions.begin(), acquisitions.end(),
            [](const GaussianModelTypes::VectorAcquisition &va1, const GaussianModelTypes::VectorAcquisition &va2) {
              auto acquisition1 = va1.second;
              auto acquisition2 = va2.second;
              return acquisition1 < acquisition2;
            });

  // project to vectors, removing acquisitions
  std::vector<GaussianModelTypes::VectorPairDiscreteContinuous> orderedVectors;
  orderedVectors.reserve(acquisitions.size());
  for (const auto &[vec, acq] : acquisitions) {
    orderedVectors.push_back(vec);
  }

  return orderedVectors;
}

std::string autopas::GaussianCluster::defaultVecToString(
    const autopas::GaussianModelTypes::VectorPairDiscreteContinuous &vec) {
  std::stringstream result;
  const auto &[discreteVec, continuousVec] = vec;

  result << "(" << discreteVec[0];
  for (long d = 1; d < discreteVec.size(); ++d) {
    result << "," << discreteVec[d];
  }

  for (long c = 0; c < continuousVec.size(); ++c) {
    result << "," << continuousVec[c];
  }

  result << ")";

  return result.str();
}

void autopas::GaussianCluster::initClusters() {
  size_t numClusters = 1;
  for (auto restriction : _dimRestriction) {
    if (restriction <= 0) {
      utils::ExceptionHandler::exception("GaussianCluster: dimension-restriction is {} but has to be positive",
                                         restriction);
    }

    numClusters *= restriction;
  }

  _clusters.clear();
  _discreteVectorMap.clear();
  _numEvidence = 0;

  _clusters.reserve(numClusters);
  _discreteVectorMap.reserve(numClusters);

  GaussianModelTypes::VectorDiscrete currentDiscrete = Eigen::VectorXi::Zero(_dimRestriction.size());
  for (size_t i = 0; i < numClusters; ++i, discreteIncrement(currentDiscrete)) {
    _clusters.emplace_back(_continuousDims, _sigma, _rng);
    _discreteVectorMap.emplace_back(currentDiscrete);
  }
}

size_t autopas::GaussianCluster::getIndex(const autopas::GaussianModelTypes::VectorDiscrete &x) const {
  if (static_cast<size_t>(x.size()) != _dimRestriction.size()) {
    utils::ExceptionHandler::exception(
        "GaussianCluster: size of discrete input {} does not match specified dimensions {}", x.size(),
        _dimRestriction.size());
  }

  size_t result = 0;
  for (long i = x.size() - 1; i >= 0; --i) {
    if (x[i] < 0 or x[i] >= _dimRestriction[i]) {
      utils::ExceptionHandler::exception("GaussianCluster: The {}th dimension is {} but is restricted to [0,{})", i,
                                         x[i], _dimRestriction[i]);
    }

    result = result * _dimRestriction[i] + x[i];
  }

  return result;
}

void autopas::GaussianCluster::discreteIncrement(autopas::GaussianModelTypes::VectorDiscrete &x) const {
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    if (++x[i] < _dimRestriction[i]) {
      break;
    }
    x[i] = 0;
  }
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
autopas::GaussianCluster::precalculateDistributions(
    const autopas::GaussianModelTypes::VectorContinuous &continuousTuple) const {
  std::vector<double> means;
  std::vector<double> vars;
  std::vector<double> stddevs;
  means.reserve(_clusters.size());
  vars.reserve(_clusters.size());
  stddevs.reserve(_clusters.size());
  for (const auto &cluster : _clusters) {
    means.push_back(cluster.predictMean(continuousTuple));

    double var = cluster.predictVar(continuousTuple);
    vars.push_back(var);
    stddevs.push_back(std::sqrt(var));
  }

  return std::make_tuple(means, vars, stddevs);
}

std::vector<double> autopas::GaussianCluster::precalculateAcquisitions(autopas::AcquisitionFunctionOption af,
                                                                       const std::vector<double> &means,
                                                                       const std::vector<double> &vars) const {
  std::vector<double> result;
  result.reserve(_clusters.size());

  for (size_t i = 0; i < _clusters.size(); ++i) {
    result.push_back(AcquisitionFunction::calcAcquisition(af, means[i], vars[i], _evidenceMaxValue));
  }

  return result;
}

autopas::GaussianModelTypes::NeighboursWeights autopas::GaussianCluster::initNeighbourWeights(
    const autopas::GaussianModelTypes::NeighbourFunction &neighbourFun) const {
  GaussianModelTypes::NeighboursWeights result;
  result.reserve(_clusters.size());

  // for each cluster create a neighbour list
  for (size_t i = 0; i < _clusters.size(); ++i) {
    auto neighbours = neighbourFun(_discreteVectorMap[i]);

    // calculate initial weight for each neighbour
    std::vector<std::tuple<size_t, double, double>> neighbourWeights;
    neighbourWeights.reserve(neighbours.size());
    for (const auto &[n, priorWeight] : neighbours) {
      size_t n_index = getIndex(n);
      // ignore clusters without evidence
      if (_clusters[n_index].numEvidence() > 0) {
        // calculate initial value for given weight function
        double weight = 0.;
        switch (_weightFun) {
          case evidenceMatchingProbabilityGM: {
            const auto &[inputs, outputs] = _clusters[n_index].getEvidence();
            weight = 1.;
            // product of probability densitiy over all evidence in neighbouring cluster if provided to the target
            // cluster
            for (size_t e = 0; e < inputs.size(); ++e) {
              weight *= _clusters[i].predictOutputPDF(inputs[e], outputs[e]);
            }
            // geometric mean
            weight = std::pow(weight, 1.0 / inputs.size());
            break;
          }
          case evidenceMatchingScaledProbabilityGM: {
            const auto &[inputs, outputs] = _clusters[n_index].getEvidence();
            weight = 1.;
            // product of probability densitiy over all evidence in neighbouring cluster if provided to the target
            // cluster
            for (size_t e = 0; e < inputs.size(); ++e) {
              weight *= _clusters[i].predictOutputScaledPDF(inputs[e], outputs[e]);
            }
            // geometric mean
            weight = std::pow(weight, 1.0 / inputs.size());
            break;
          }
          case wasserstein2:
            break;
        }
        neighbourWeights.emplace_back(n_index, priorWeight, priorWeight * weight);
      }
    }
    result.push_back(std::move(neighbourWeights));
  }

  return result;
}

void autopas::GaussianCluster::updateNeighbourWeights(autopas::GaussianModelTypes::NeighboursWeights &neighbourWeights,
                                                      const std::vector<double> &means, const std::vector<double> &vars,
                                                      const std::vector<double> &stddevs) const {
  // for each cluster update neighbour-weight list
  GaussianModelTypes::VectorDiscrete currentDiscrete = Eigen::VectorXi::Zero(_dimRestriction.size());
  for (size_t i = 0; i < _clusters.size(); ++i) {
    // for each neighbour update weight
    for (auto &[n, priorWeight, weight] : neighbourWeights[i]) {
      switch (_weightFun) {
        case evidenceMatchingProbabilityGM:
        case evidenceMatchingScaledProbabilityGM:
          // keep inital value
          break;
        case wasserstein2:
          weight = vars[i] / (vars[i] + std::pow(means[n] - means[i], 2) + std::pow(stddevs[n] - stddevs[i], 2));
          weight *= priorWeight;
          break;
      }
    }
  }
}
