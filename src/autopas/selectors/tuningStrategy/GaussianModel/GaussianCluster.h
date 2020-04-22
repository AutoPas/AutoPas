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
        _numEvidence(0),
        _evidenceMinValue(0),
        _evidenceMaxValue(0) {
    size_t numClusters = 1;
    for (auto restriction : _dimRestriction) {
      if (restriction <= 0) {
        utils::ExceptionHandler::exception("GaussianCluster: dimension-restriction is {} but has to be positive",
                                           restriction);
      }

      numClusters *= restriction;
    }

    _clusters.reserve(numClusters);
    for (size_t i = 0; i < numClusters; ++i) {
      _clusters.emplace_back(continuousDims, sigma, rngRef);
    }
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
  size_t numEvidence() const { return _numEvidence; }

  /**
   * Provide a input-output pair as evidence.
   * Each evidence improve the quality of future predictions.
   * @param inputDiscrete x
   * @param inputContinuous y
   * @param output f((x,y))
   */
  void addEvidence(VectorDiscrete inputDiscrete, VectorContinuous inputContinuous, double output) {
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
      ;
    } else if (output > _evidenceMaxValue) {
      _evidenceMaxValue = output;
      _evidenceMaxVector = std::make_pair(inputDiscrete, inputContinuous);
      ;
    }

    _clusters[clusterIdx].addEvidence(inputContinuous, output);
    ++_numEvidence;
  }

  /**
   * Get the evidence with the smallest output value
   * @return input of min
   */
  std::pair<VectorDiscrete, VectorContinuous> getEvidenceMin() {
    if (_numEvidence == 0) {
      utils::ExceptionHandler::exception("GaussianCluster has no evidence");
    }

    return _evidenceMinVector;
  }

  /**
   * Get the evidence with the highest output value
   * @return input of max
   */
  std::pair<VectorDiscrete, VectorContinuous> getEvidenceMax() {
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

    for (auto currentContinuous : samples) {
      // calc means of given continuous tuple for each cluster
      std::vector<double> means;
      std::vector<double> vars;
      std::vector<double> stddevs;
      means.reserve(_clusters.size());
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
            double n_weight = 1. / (1 + distance / vars[i]);
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
          default:
            currentValue = AcquisitionFunction::calcAcquisition(af, mixMean, mixVar);
        }

        acquisitions.emplace_back(std::make_pair(currentDiscrete, currentContinuous), currentValue);
      }
    }

    return acquisitions;
  }

 private:
  /**
   * Get the cluster index of a discrete tuple
   * @param x the discrete tuple
   * @return
   */
  size_t getIndex(const VectorDiscrete &x) const {
    if (static_cast<size_t>(x.size()) != _dimRestriction.size()) {
      utils::ExceptionHandler::exception(
          "GaussianCluster: size of discrete input {} does not match specified dimensions {}", x.size(),
          _dimRestriction.size());
    }

    auto i = x.size() - 1;

    if (x[i] < 0 || x[i] >= _dimRestriction[i]) {
      utils::ExceptionHandler::exception("GaussianCluster: The {}th dimension is {} but is restricted to [0,{})", i,
                                         x[i], _dimRestriction[i]);
    }

    size_t result = x[i];
    for (--i; i >= 0; --i) {
      if (x[i] < 0 || x[i] >= _dimRestriction[i]) {
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
   * restrict the i-th dimension to a integer between 0 and _dimRestriction[i]-1
   */
  const std::vector<int> _dimRestriction;
  /**
   * number of additional unrestricted dimensions
   */
  const size_t _continuousDims;

  /**
   * Gaussian process for each discrete tuple
   */
  std::vector<GaussianProcess> _clusters;

  /**
   * Current smallest evidence output
   */
  double _evidenceMinValue;
  /**
   * Current smallest evidence input
   */
  std::pair<VectorDiscrete, VectorContinuous> _evidenceMinVector;
  /**
   * Current greatest evidence output
   */
  double _evidenceMaxValue;
  /**
   * Current greatest evidence input
   */
  std::pair<VectorDiscrete, VectorContinuous> _evidenceMaxVector;
  /**
   * Current number of evidence
   */
  size_t _numEvidence;
};
}  // namespace autopas
