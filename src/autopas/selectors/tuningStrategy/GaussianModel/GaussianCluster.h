/**
 * @file GaussianCluster.h
 * @author Jan Nguyen
 * @date 05.04.20
 */

#pragma once

#include <Eigen/Core>
#include <utility>

#include "GaussianClusterLogger.h"
#include "GaussianProcess.h"
#include "autopas/options/AcquisitionFunctionOption.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Math.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/Random.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * Model to predicts the output of a blackbox function f(x) for given input x.
 * The model separates discrete and continuous dimensions of x. For each possible
 * discrete tuple a Gaussian process is assigned to estimate f(x) if the tuple is fixed.
 * Some sample input-output pairs (x,f(x)) should be provided as evidence.
 */
class GaussianCluster {
  // number of samples to find optimal hyperparameters
  static constexpr size_t hp_sample_size = 500;
  // number of hyperparameters
  static constexpr size_t hp_size = 25;

 public:
  /**
   * Different weight functions between clusters
   */
  enum WeightFunction {
    /**
     * geometric mean of probability densitiy over all evidence in neighbouring cluster if provided to the target
     * cluster
     */
    evidenceMatchingProbabilityGM,
    /**
     * geometric mean of scaled probability densitiy over all evidence in neighbouring cluster if provided to the target
     * cluster Probability densities are scaled such that the maximum is 1.
     */
    evidenceMatchingScaledProbabilityGM,
    /**
     * Wasserstein-2 distance of normal distributions of cluster given a continuous cluster
     */
    wasserstein2
  };

  /**
   * Constructor
   * @param dimRestriction restrict the i-th dimension to a integer between 0 and dimRestriction[i]-1
   * @param continuousDims additional unrestricted dimensions
   * @param weightFun function to calculate weight between clusters
   * @param sigma fixed noise
   * @param rngRef reference to random number generator
   * @param vectorToString function to convert vectors to a readable string
   */
  GaussianCluster(const std::vector<int> &dimRestriction, size_t continuousDims, WeightFunction weightFun, double sigma,
                  Random &rngRef, const GaussianModelTypes::VectorToStringFun &vectorToString = defaultVecToString)
      : _dimRestriction(dimRestriction),
        _continuousDims(continuousDims),
        _clusters(),
        _weightFun(weightFun),
        _evidenceMinValue(0),
        _evidenceMaxValue(0),
        _numEvidence(0),
        _sigma(sigma),
        _rng(rngRef),
        _vecToStringFun(vectorToString) {
    initClusters();
  }

  /**
   * Get the number of clusters in each dimension.
   * @return
   */
  [[nodiscard]] const std::vector<int> &getDimensions() const { return _dimRestriction; }

  /**
   * Change the number of cluster in all dimension.
   * This will discard all evidence.
   * @param newValue new number of clusters in each dimension
   */
  void setDimensions(const std::vector<int> &newValue) {
    _dimRestriction = newValue;
    initClusters();
  }

  /**
   * Get the underlying GaussianProcess of a cluster.
   * This function should only be used for testing purposes.
   * @param index1D
   * @return
   */
  [[nodiscard]] const GaussianProcess &getCluster(size_t index1D) const { return _clusters[index1D]; }

  /**
   * Discard all evidence.
   */
  void clear() {
    for (auto cluster : _clusters) {
      cluster.clear();
    }
    _numEvidence = 0;
  }

  /**
   * Get the number of evidence provided.
   * @return
   */
  [[nodiscard]] size_t numEvidence() const { return _numEvidence; }

  /**
   * Provide a input-output pair as evidence.
   * Each evidence improve the quality of future predictions.
   * @param inputDiscrete x
   * @param inputContinuous y
   * @param output f((x,y))
   */
  void addEvidence(const GaussianModelTypes::VectorDiscrete &inputDiscrete,
                   const GaussianModelTypes::VectorContinuous &inputContinuous, double output) {
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

  /**
   * Provide a input-output pair as evidence.
   * Each evidence improves the quality of future predictions.
   * @param input (x,y)
   * @param output f((x,y))
   */
  inline void addEvidence(const GaussianModelTypes::VectorPairDiscreteContinuous &input, double output) {
    addEvidence(input.first, input.second, output);
  }

  /**
   * Get the evidence with the highest output value
   * @return input of max
   */
  [[nodiscard]] GaussianModelTypes::VectorPairDiscreteContinuous getEvidenceMax() const {
    if (_numEvidence == 0) {
      utils::ExceptionHandler::exception("GaussianCluster has no evidence");
    }

    return _evidenceMaxVector;
  }

  /**
   * Generate all possible combinations of discrete tuples and continuous tuples in samples
   * and calculate their corresponding acquisition.
   * @param af function to optimize
   * @param neighbourFun function which generates neighbours of given discrete tuple
   * @param continuousSamples continuous tuples
   * @return all discrete-continuous tuples paired with their corresponding acquisition
   */
  [[nodiscard]] std::vector<GaussianModelTypes::VectorAcquisition> sampleAcquisition(
      AcquisitionFunctionOption af, const GaussianModelTypes::NeighbourFunction &neighbourFun,
      const std::vector<GaussianModelTypes::VectorContinuous> &continuousSamples) const {
    std::vector<GaussianModelTypes::VectorAcquisition> acquisitions;
    // pair up all discrete with all continuous samples
    acquisitions.reserve(_clusters.size() * continuousSamples.size());

    GaussianClusterLogger graphLogger(_vecToStringFun);

    auto neighbourWeights = initNeighbourWeights(neighbourFun);

    for (const auto &continuousSample : continuousSamples) {
      auto [means, vars, stddevs] = precalculateDistributions(continuousSample);
      updateNeighbourWeights(neighbourWeights, means, vars, stddevs);

      graphLogger.add(_clusters, _discreteVectorMap, continuousSample, means, vars, neighbourWeights);

      // get acquisition considering neighbours
      for (size_t i = 0; i < _clusters.size(); ++i) {
        // target cluster gets weight 1.
        double mixMean = means[i];
        double mixVar = vars[i];
        double weightSum = 1.;

        // sum over neighbours
        for (const auto &[n, weight] : neighbourWeights[i]) {
          mixMean += means[n] * weight;
          mixVar += vars[n] * weight * weight;
          weightSum += weight;
        }

        // normalize
        mixMean /= weightSum;
        mixVar /= weightSum * weightSum;

        // calc acquisition
        double currentValue = AcquisitionFunction::calcAcquisition(af, mixMean, mixVar, _evidenceMaxValue);
        acquisitions.emplace_back(std::make_pair(_discreteVectorMap[i], continuousSample), currentValue);
      }
    }

    graphLogger.end();

    return acquisitions;
  }

  /**
   * Generate all possible combinations of discrete tuples and continuous tuples in samples
   * and calculate their weight to each other. Output graph in AutoPasLog Trace.
   * @param neighbourFun function which generates neighbours of given discrete tuple
   * @param continuousSamples continuous tuples
   */
  void logDebugGraph(
      const std::function<std::vector<GaussianModelTypes::VectorDiscrete>(const GaussianModelTypes::VectorDiscrete &)>
          &neighbourFun,
      const std::vector<GaussianModelTypes::VectorContinuous> &continuousSamples) const {
    if (autopas::Logger::get()->level() > autopas::Logger::LogLevel::trace) {
      return;
    }

    GaussianClusterLogger graphLogger(_vecToStringFun);

    auto neighbourWeights = initNeighbourWeights(neighbourFun);
    for (const auto &continuousSample : continuousSamples) {
      auto [means, vars, stddevs] = precalculateDistributions(continuousSample);
      updateNeighbourWeights(neighbourWeights, means, vars, stddevs);

      graphLogger.add(_clusters, _discreteVectorMap, continuousSample, means, vars, neighbourWeights);
    }

    graphLogger.end();
  }
  /**
   * Change the used function to convert from vector to string.
   * @param fun new converter
   */
  void setVectorToStringFun(const GaussianModelTypes::VectorToStringFun &fun) { _vecToStringFun = fun; }

  /**
   * Generate all possible combinations of discrete tuples and continuous tuples in samples and
   * returns the vector with maximum acquisition.
   * @param af function to optimize
   * @param neighbourFun function which generates neighbours of given discrete tuple
   * @param continuousSamples continuous tuples
   * @return pair of vector and corresponding maximum acquistion
   */
  [[nodiscard]] GaussianModelTypes::VectorAcquisition sampleAcquisitionMax(
      AcquisitionFunctionOption af, const GaussianModelTypes::NeighbourFunction &neighbourFun,
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

  /**
   * Generate all possible combinations of discrete tuples and continuous tuples in samples and
   * order them by their acquisition.
   * @param af function to optimize
   * @param neighbourFun function which generates neighbours of given discrete tuple
   * @param continuousSamples continuous tuples
   * @return all discrete-continuous tuples ordered by their corresponding acquisition
   */
  [[nodiscard]] std::vector<GaussianModelTypes::VectorPairDiscreteContinuous> sampleOrderedByAcquisition(
      AcquisitionFunctionOption af, const GaussianModelTypes::NeighbourFunction &neighbourFun,
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

 private:
  /**
   * Create a GaussianProcess for each cluster and precalculate DiscreteVector for each cluster index.
   */
  void initClusters() {
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

  /**
   * Get the cluster index of a discrete tuple.
   * @param x the discrete tuple
   * @return
   */
  [[nodiscard]] size_t getIndex(const GaussianModelTypes::VectorDiscrete &x) const {
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

  /**
   * Increment the given discrete tuple x, resulting
   * in the tuple following x in the cluster list.
   * The increment works in a date-like fashion:
   * One increment adds one to the first vector entry.
   * If this entry overflows also the next entry is incremented and so on.
   * @param x discrete tuple
   */
  void discreteIncrement(GaussianModelTypes::VectorDiscrete &x) const {
    for (Eigen::Index i = 0; i < x.size(); ++i) {
      if (++x[i] < _dimRestriction[i]) {
        break;
      }
      x[i] = 0;
    }
  }

  /**
   * Calculate mean, variance and stddev for all clusters for given continous tuple.
   * @param continuousTuple
   * @return Vectors means, vars, stddevs containing corresponding values for each cluster
   */
  [[nodiscard]] std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> precalculateDistributions(
      const GaussianModelTypes::VectorContinuous &continuousTuple) const {
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

  /**
   * Initalize the neighbour-weight list for each cluster.
   * @param neighbourFun function which generates neighbours of given discrete tuple
   * @return
   */
  [[nodiscard]] GaussianModelTypes::NeighboursWeights initNeighbourWeights(
      const GaussianModelTypes::NeighbourFunction &neighbourFun) const {
    GaussianModelTypes::NeighboursWeights result;
    result.reserve(_clusters.size());

    // for each cluster create a neighbour list
    for (size_t i = 0; i < _clusters.size(); ++i) {
      auto neighbours = neighbourFun(_discreteVectorMap[i]);

      // calculate initial weight for each neighbour
      std::vector<std::pair<size_t, double>> neighbourWeights;
      neighbourWeights.reserve(neighbours.size());
      for (const auto &n : neighbours) {
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
          neighbourWeights.emplace_back(n_index, weight);
        }
      }
      result.push_back(std::move(neighbourWeights));
    }

    return result;
  }

  /**
   * Update the neighbour-weight list for each cluster. This function is called for each continuous tuple.
   * Given mean, var and stddev should be evaluated for the current continuous tuple.
   * @param neighbourWeights list to update
   * @param means mean for each cluster
   * @param vars variance for each cluster
   * @param stddev standard deviation for each cluster
   * @return
   */
  void updateNeighbourWeights(GaussianModelTypes::NeighboursWeights &neighbourWeights, const std::vector<double> &means,
                              const std::vector<double> &vars, const std::vector<double> &stddevs) const {
    // for each cluster update neighbour-weight list
    GaussianModelTypes::VectorDiscrete currentDiscrete = Eigen::VectorXi::Zero(_dimRestriction.size());
    for (size_t i = 0; i < _clusters.size(); ++i) {
      // for each neighbour update weight
      for (auto &[n, weight] : neighbourWeights[i]) {
        switch (_weightFun) {
          case evidenceMatchingProbabilityGM:
          case evidenceMatchingScaledProbabilityGM:
            // keep inital value
            break;
          case wasserstein2:
            weight = vars[i] / (vars[i] + std::pow(means[n] - means[i], 2) + std::pow(stddevs[n] - stddevs[i], 2));
            break;
        }
      }
    }
  }

  /**
   * Default function used to convert vectors to readable strings.
   * @note begin() and end() currently not available for Eigen::Vector, so AutoPas ArrayUtils cannot be used.
   * @param vec
   * @return string with format (a,b,...,n) beginning with discrete values.
   */
  static std::string defaultVecToString(const GaussianModelTypes::VectorPairDiscreteContinuous &vec) {
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

  /**
   * Number of clusters per discrete dimension.
   */
  std::vector<int> _dimRestriction;
  /**
   * Number of additional unrestricted continuous dimensions.
   */
  const size_t _continuousDims;

  /**
   * Gaussian process for each discrete tuple.
   */
  std::vector<GaussianProcess> _clusters;

  /**
   * Vector to easiliy map from index to DiscreteVector
   */
  std::vector<GaussianModelTypes::VectorDiscrete> _discreteVectorMap;

  /**
   * Function used to calculate the weight between two clusters.
   */
  WeightFunction _weightFun;

  /**
   * Current smallest evidence output.
   */
  double _evidenceMinValue;
  /**
   * Current greatest evidence output.
   */
  double _evidenceMaxValue;
  /**
   * Current greatest evidence input
   */
  GaussianModelTypes::VectorPairDiscreteContinuous _evidenceMaxVector;
  /**
   * Current number of evidence.
   */
  size_t _numEvidence;

  /**
   * Fixed noise assumed.
   */
  const double _sigma;

  Random &_rng;

  /**
   * Function to convert vectors to strings.
   */
  GaussianModelTypes::VectorToStringFun _vecToStringFun;
};
}  // namespace autopas
