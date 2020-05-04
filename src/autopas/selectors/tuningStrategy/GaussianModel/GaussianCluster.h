/**
 * @file GaussianCluster.h
 * @author Jan Nguyen
 * @date 05.04.20
 */

#pragma once

#include <Eigen/Core>
#include <utility>

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
  using VectorDiscrete = Eigen::VectorXi;
  using VectorContinuous = Eigen::VectorXd;

 public:
  /**
   * Pair of input tuples and corresponding expected acquisition
   */
  using VectorAcquisition = std::pair<std::pair<VectorDiscrete, VectorContinuous>, double>;

  /**
   * Constructor
   * @param dimRestriction restrict the i-th dimension to a integer between 0 and dimRestriction[i]-1
   * @param continuousDims additional unrestricted dimensions
   * @param sigma fixed noise
   * @param rngRef reference to random number generator
   */
  GaussianCluster(const std::vector<int> &dimRestriction, size_t continuousDims, double sigma, Random &rngRef)
      : _dimRestriction(dimRestriction),
        _continuousDims(continuousDims),
        _clusters(),
        _evidenceMinValue(0),
        _evidenceMaxValue(0),
        _numEvidence(0),
        _sigma(sigma),
        _rng(rngRef) {
    initClusters();
  }

  /**
   * Get the number of clusters in each dimension.
   * @return
   */
  const std::vector<int> &getDimensions() const { return _dimRestriction; }

  /**
   * Change the number of cluster in a dimension
   * @param dim the dimension to change
   * @param newValue new number of clusters
   */
  void setDimension(size_t dim, int newValue) {
    _dimRestriction[dim] = newValue;
    initClusters();
  }

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
  void addEvidence(const VectorDiscrete &inputDiscrete, const VectorContinuous &inputContinuous, double output) {
    size_t clusterIdx = getIndex(inputDiscrete);
    if (static_cast<size_t>(inputContinuous.size()) != _continuousDims) {
      utils::ExceptionHandler::exception(
          "GaussianCluster: size of continuous input {} does not match specified dimensions {}", inputContinuous.size(),
          _continuousDims);
    }

    if (_numEvidence == 0) {
      // first evidence
      _evidenceMinValue = _evidenceMaxValue = output;
      _evidenceMinVector = _evidenceMaxVector = std::make_pair(inputDiscrete, inputContinuous);
    } else if (output < _evidenceMinValue) {
      _evidenceMinValue = output;
      _evidenceMinVector = std::make_pair(inputDiscrete, inputContinuous);
    } else if (output > _evidenceMaxValue) {
      _evidenceMaxValue = output;
      _evidenceMaxVector = std::make_pair(inputDiscrete, inputContinuous);
    }

    _clusters[clusterIdx].addEvidence(inputContinuous, output);
    ++_numEvidence;
  }

  /**
   * Get the evidence with the smallest output value
   * @return input of min
   */
  [[nodiscard]] std::pair<VectorDiscrete, VectorContinuous> getEvidenceMin() const {
    if (_numEvidence == 0) {
      utils::ExceptionHandler::exception("GaussianCluster has no evidence");
    }

    return _evidenceMinVector;
  }

  /**
   * Get the evidence with the highest output value
   * @return input of max
   */
  [[nodiscard]] std::pair<VectorDiscrete, VectorContinuous> getEvidenceMax() const {
    if (_numEvidence == 0) {
      utils::ExceptionHandler::exception("GaussianCluster has no evidence");
    }

    return _evidenceMaxVector;
  }

  /**
   * Calculate the acquisition for all possible combinations of discrete tuples and continuous tuples in samples
   * @param af function to optimize
   * @param neighbourFun function which generates neighbours of given discrete tuple
   * @param samples continuous tuples
   * @return vector of pairs - vectors and corresponding acquisition
   */
  std::vector<VectorAcquisition> sampleAcquisition(
      AcquisitionFunctionOption af, const std::function<std::vector<Eigen::VectorXi>(Eigen::VectorXi)> &neighbourFun,
      const std::vector<VectorContinuous> &samples) {
    std::vector<VectorAcquisition> acquisitions;
    // pair up all discrete with all continuous samples
    acquisitions.reserve(_clusters.size() * samples.size());

    for (const auto &currentContinuous : samples) {
      // calc means of given continuous tuple for each cluster
      std::vector<double> means;
      std::vector<double> vars;
      std::vector<double> stddevs;
      means.reserve(_clusters.size());
      vars.reserve(_clusters.size());
      stddevs.reserve(_clusters.size());
      for (const auto &cluster : _clusters) {
        means.push_back(cluster.predictMean(currentContinuous));

        double var = cluster.predictVar(currentContinuous);
        vars.push_back(var);
        stddevs.push_back(std::sqrt(var));
      }

      // get maximum acquisition considering neighbours
      VectorDiscrete currentDiscrete = Eigen::VectorXi::Zero(_dimRestriction.size());
      for (size_t i = 0; i < _clusters.size(); ++i, discreteIncrement(currentDiscrete)) {
        auto neighbours = neighbourFun(currentDiscrete);

        double mixMean = means[i];
        double mixVar = vars[i];
        double weightSum = 1.;

        // sum over neighbours
        for (const auto &n : neighbours) {
          size_t n_index = getIndex(n);
          // ignore clusters without evidence
          if (_clusters[i].numEvidence() > 0) {
            // weight based on 2-wasserstein distance
            double distance = std::pow(means[n_index] - means[i], 2) + std::pow(stddevs[n_index] - stddevs[i], 2);
            double n_weight = vars[i] / (vars[i] + distance);
            mixMean += means[n_index] * n_weight;
            mixVar += vars[n_index] * n_weight * n_weight;
            weightSum += n_weight;
          }
        }

        // normalize
        mixMean /= weightSum;
        mixVar /= weightSum;

        // calc acquisition
        double currentValue;
        switch (af) {
          case AcquisitionFunctionOption::probabilityOfDecrease:
          case AcquisitionFunctionOption::expectedDecrease:
            currentValue = AcquisitionFunction::calcAcquisition(af, mixMean, mixVar, _evidenceMinValue);
          case AcquisitionFunctionOption::mean:
          case AcquisitionFunctionOption::variance:
          case AcquisitionFunctionOption::lowerConfidenceBound:
          case AcquisitionFunctionOption::upperConfidenceBound:
            currentValue = AcquisitionFunction::calcAcquisition(af, mixMean, mixVar);
        }

        acquisitions.emplace_back(std::make_pair(currentDiscrete, currentContinuous), currentValue);
      }
    }

    return acquisitions;
  }

 private:
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
    _numEvidence = 0;

    _clusters.reserve(numClusters);
    for (size_t i = 0; i < numClusters; ++i) {
      _clusters.emplace_back(_continuousDims, _sigma, _rng);
    }
  }

  /**
   * Get the cluster index of a discrete tuple
   * @param x the discrete tuple
   * @return
   */
  [[nodiscard]] size_t getIndex(const VectorDiscrete &x) const {
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
  void discreteIncrement(VectorDiscrete &x) const {
    for (Eigen::Index i = 0; i < x.size(); ++i) {
      if (++x[i] < _dimRestriction[i]) {
        break;
      }
      x[i] = 0;
    }
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
   * Current smallest evidence output.
   */
  double _evidenceMinValue;
  /**
   * Current smallest evidence input.
   */
  std::pair<VectorDiscrete, VectorContinuous> _evidenceMinVector;
  /**
   * Current greatest evidence output.
   */
  double _evidenceMaxValue;
  /**
   * Current greatest evidence input
   */
  std::pair<VectorDiscrete, VectorContinuous> _evidenceMaxVector;
  /**
   * Current number of evidence.
   */
  size_t _numEvidence;

  /**
   * Fixed noise assumed.
   */
  const double _sigma;

  Random &_rng;
};
}  // namespace autopas
