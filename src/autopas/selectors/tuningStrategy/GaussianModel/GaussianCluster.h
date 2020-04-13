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
   * Constructor
   * @param dimRestriction restrict the i-th dimension to a integer between 0 and dimRestriction[i]-1
   * @param continuousDims additional unrestricted dimensions
   * @param sigma fixed noise
   * @param rngRef reference to random number generator
   */
  GaussianCluster(std::vector<int> dimRestriction, size_t continuousDims, double sigma, Random &rngRef)
      : _dimRestriction(dimRestriction), _continuousDims(continuousDims), _clusters(), _numEvidence(0) {
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
   * @param input x
   * @param output f(x)
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
   * Calculates the acquisition function for given input.
   * @param af acquisition function a:input->double
   * @param input i
   * @return a(i). This value can be compared with values a(x) of other inputs x to weigh which input would give the
   * most gain if its evidence were provided.
   */
  inline double calcAcquisition(AcquisitionFunctionOption af, const VectorDiscrete &inputDiscrete,
                                const VectorContinuous &inputContinuous) const {
    return _clusters[getIndex(inputDiscrete)].calcAcquisition(af, inputContinuous);
  }

  /**
   * Find the discrete tuple and the continuous tuple in samples which maximizes given aquisition function.
   * @param af function to maximize
   * @param samples
   * @return
   */
  std::pair<VectorDiscrete, VectorContinuous> sampleAquisitionMax(AcquisitionFunctionOption af,
                                                                  const std::vector<VectorContinuous> &samples) const {
    VectorDiscrete v = Eigen::VectorXi::Zero(_dimRestriction.size());

    // find the maximum of the first cluster
    auto maxLocal = _clusters[0].sampleAquisitionMaxValue(af, samples);
    VectorDiscrete maxVectorDiscrete = v;
    VectorContinuous maxVectorContinuous = maxLocal.first;
    double maxValue = maxLocal.second;

    // search for other clusters for a greater maximum
    for (size_t i = 1; i < _clusters.size(); ++i, discreteIncrement(v)) {
      maxLocal = _clusters[i].sampleAquisitionMaxValue(af, samples);
      if (maxLocal.second > maxValue) {
        maxVectorDiscrete = v;
        maxVectorContinuous = maxLocal.first;
        maxValue = maxLocal.second;
      }
    }

    return std::make_pair(maxVectorDiscrete, maxVectorContinuous);
  }

  /**
   * Find the discrete tuple and the continuous tuple in samples which minimizes given aquisition function.
   * @param af function to minimize
   * @param samples
   * @return
   */
  std::pair<VectorDiscrete, VectorContinuous> sampleAquisitionMin(AcquisitionFunctionOption af,
                                                                  const std::vector<VectorContinuous> &samples) const {
    VectorDiscrete v = Eigen::VectorXi::Zero(_dimRestriction.size());

    // find the minimum of the first cluster
    auto minLocal = _clusters[0].sampleAquisitionMinValue(af, samples);
    VectorDiscrete minVectorDiscrete = v;
    VectorContinuous minVectorContinuous = minLocal.first;
    double minValue = minLocal.second;

    // search for other clusters for a lesser minimum
    for (size_t i = 1; i < _clusters.size(); ++i, discreteIncrement(v)) {
      minLocal = _clusters[i].sampleAquisitionMinValue(af, samples);
      if (minLocal.second < minValue) {
        minVectorDiscrete = v;
        minVectorContinuous = minLocal.first;
        minValue = minLocal.second;
      }
    }

    return std::make_pair(minVectorDiscrete, minVectorContinuous);
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
      } else {
        x[i] = 0;
      }
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
