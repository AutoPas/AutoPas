/**
 * @file GaussianProcess.h
 * @author Jan Nguyen
 * @date 17.05.19
 */

#pragma once

#include <Eigen/Dense>
#include "autopas/options/AcquisitionFunctionOption.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Random.h"

namespace autopas {

/**
 * Gaussian process is a stochastical model. It predicts the
 * output of a blackbox function f(x) for given input x. To do so, some sample
 * input-output pairs (x,f(x)) should be provided as evidence.
 *
 * Currently the default mean is 0 and squared exponential kernel is used.
 * TODO: maybe offer some options.
 * @tparam Vector class should be subtractable and convertible to Eigen::VectorXd
 */
template <class Vector>
class GaussianProcess {
 public:
  /**
   * Constructor
   * @param dims number of input dimensions
   * @param theta default variance
   * @param sigma fixed noise
   * @param rngRef reference to rng
   */
  GaussianProcess(size_t dims, double theta, double sigma, Random &rngRef)
      : _inputs(),
        _outputs(),
        _dims(dims),
        _defaultTheta(theta),
        _sigma(sigma),
        _covMat(),
        _covMatInv(),
        _weights(),
        _rng(rngRef) {
    updateHyperparameters();
  }

  /**
   * Discard all evidence.
   */
  void clear() {
    _inputs.clear();
    updateHyperparameters();
  }

  /**
   * Get the number of evidence provided.
   * @return
   */
  size_t numEvidence() const { return _inputs.size(); }

  /**
   * Provide a input-output pair as evidence.
   * Each evidence improve the quality of future predictions.
   * @param input x
   * @param output f(x)
   */
  void addEvidence(Vector input, double output) {
    auto inputVec = static_cast<Eigen::VectorXd>(input);
    if (static_cast<size_t>(inputVec.size()) != _dims) {
      utils::ExceptionHandler::exception("GaussianProcess: size of input {} does not match specified dimensions {}",
                                         inputVec.size(), _dims);
    }

    _inputs.push_back(input);
    long newSize = _inputs.size();

    // extend output vector
    _outputs.conservativeResize(newSize, Eigen::NoChange_t());
    _outputs(newSize - 1) = output;

    updateHyperparameters();
  }

  /**
   * Get the evidence with the smallest output value
   * @return input of min
   */
  Vector getEvidenceMin() {
    if (_inputs.empty()) {
      utils::ExceptionHandler::exception("GaussianProcess has no evidence");
    }

    double min = _outputs[0];
    Vector minVec = _inputs[1];

    for (size_t i = 1; i < _inputs.size(); ++i) {
      if (_outputs[i] < min) {
        min = _outputs[i];
        minVec = _inputs[i];
      }
    }

    return minVec;
  }

  /**
   * Get the evidence with the highest output value
   * @return input of max
   */
  Vector getEvidenceMax() {
    if (_inputs.empty()) {
      utils::ExceptionHandler::exception("GaussianProcess has no evidence");
    }

    double max = _outputs[0];
    Vector maxVec = _inputs[1];

    for (size_t i = 1; i < _inputs.size(); ++i) {
      if (_outputs[i] > max) {
        max = _outputs[i];
        maxVec = _inputs[i];
      }
    }

    return maxVec;
  }

  /**
   * Try to predict f(x) using the evidence
   * provided so far.
   * @param input x
   * @return expected output of f(x)
   */
  double predictMean(const Vector &input) const {
    auto inputVec = static_cast<Eigen::VectorXd>(input);
    if (static_cast<size_t>(inputVec.size()) != _dims) {
      utils::ExceptionHandler::exception("GaussianProcess: size of input {} does not match specified dimensions {}",
                                         inputVec.size(), _dims);
    }

    if (_inputs.size() == 0) return 0.;

    return kernelVector(input).dot(_weights);
  }

  /**
   * The variance of the predicted f(x) from predictMean().
   * @param input x
   * @return variance
   */
  double predictVar(const Vector &input) const {
    auto inputVec = static_cast<Eigen::VectorXd>(input);
    if (static_cast<size_t>(inputVec.size()) != _dims) {
      utils::ExceptionHandler::exception("GaussianProcess: size of input {} does not match specified dimensions {}",
                                         inputVec.size(), _dims);
    }

    if (_inputs.size() == 0) return kernel(input, input);

    Eigen::VectorXd kVec = kernelVector(input);
    return kernel(input, input) - kVec.dot(_covMatInv * kVec);
  }

  /**
   * Calculates the acquisition function for given input.
   * @param af acquisition function a:input->double
   * @param input i
   * @return a(i). This value can be compared with values a(x) of other inputs x to weigh which input would give the
   * most gain if its evidence were provided.
   */
  inline double calcAcquisition(AcquisitionFunctionOption af, const Vector &input) const {
    switch (af) {
      case ucb: {
        return predictMean(input) + std::sqrt(predictVar(input));
      }
      case lcb: {
        return predictMean(input) - std::sqrt(predictVar(input));
      }
      case mean: {
        return predictMean(input);
      }
    }

    autopas::utils::ExceptionHandler::exception("GaussianProcess.calcAcquisition: Unknown acquisition function {}.",
                                                af);
    return 0;
  }

  /**
   * Find the input in samples which maximizes given aquisition function.
   * TODO: maybe add parameters for hyperparameters of aquisition functions
   * @param af function to maximize
   * @param samples
   * @return
   */
  Vector sampleAquisitionMax(AcquisitionFunctionOption af, const std::vector<Vector> &samples) const {
    int maxIdx = -1;
    double maxVal = 0.;

    // find maximum from samples
    for (unsigned i = 0; i < samples.size(); ++i) {
      double val = calcAcquisition(af, samples[i]);

      if (maxIdx == -1 || val > maxVal) {
        maxIdx = i;
        maxVal = val;
      }
    }

    return samples[maxIdx];
  }

  /**
   * Find the input in samples which minimizes given aquisition function.
   * TODO: maybe add parameters for hyperparameters of aquisition functions
   * @param af function to minimize
   * @param samples
   * @return
   */
  Vector sampleAquisitionMin(AcquisitionFunctionOption af, const std::vector<Vector> &samples) const {
    int minIdx = -1;
    double minVal = 0.;

    // find minimum from samples
    for (unsigned i = 0; i < samples.size(); ++i) {
      double val = calcAcquisition(af, samples[i]);

      if (minIdx == -1 || val < minVal) {
        minIdx = i;
        minVal = val;
      }
    }

    return samples[minIdx];
  }

 private:
  /**
   * Update the hyperparameters: theta, dimScale.
   * To do so, hyperparameter-samples are randomly generated.
   * The samples are combined using a weighted average. The weight of a sample
   * equals to the probability that given evidence and hyperparameter-sample
   * generates given output.
   */
  inline void updateHyperparameters() {
    // distribution of theta: gamma distribution with expectation _defaultTheta.
    std::gamma_distribution<double> thetaDistribution(3., _defaultTheta / 3.);
    // distribution of dimScale: gamma distribution with expectation 1.
    std::gamma_distribution<double> dimScaleDistribution(3., 1. / 3.);

    // size of monte carlo simulation
    const size_t mt_size = 1000;

    // number of evidence
    size_t newSize = _inputs.size();

    // if no evidence
    if (newSize == 0) {
      // use default values
      _theta = _defaultTheta;
      _dimScale = Eigen::VectorXd::Ones(_dims);
      return;
    }

    _covMat.resize(newSize, newSize);

    // initialize sums to 0
    double scoreSum = 0;
    double thetaSum = 0;
    Eigen::VectorXd dimScaleSum = Eigen::VectorXd::Zero(_dims);

    // @TODO parallelize?
    for (size_t t = 0; t < mt_size; ++t) {
      // generate theta
      double theta = thetaDistribution(_rng);

      // generate dimScale
      std::vector<double> dimScaleData;
      for (size_t d = 0; d < _dims; ++d) {
        dimScaleData.push_back(dimScaleDistribution(_rng));
      }
      Eigen::VectorXd dimScale = Eigen::Map<Eigen::VectorXd>(dimScaleData.data(), dimScaleData.size());

      // calculate covariance matrix
      for (size_t i = 0; i < newSize; ++i) {
        _covMat(i, i) = kernel(_inputs[i], _inputs[i], theta, dimScale) + _sigma;
        for (size_t j = i + 1; j < newSize; ++j) {
          _covMat(i, j) = _covMat(j, i) = kernel(_inputs[i], _inputs[j], theta, dimScale);
        }
      }

      // weight tested hyperparameter by probability to fit evidence
      double score = std::exp(-0.5 * (_outputs.dot(_covMat.llt().solve(_outputs)) + std::log(_covMat.determinant())));
      thetaSum += theta * score;
      dimScaleSum += dimScale * score;
      scoreSum += score;
    }

    // get weighted average;
    _theta = thetaSum / scoreSum;
    _dimScale = dimScaleSum / scoreSum;

    // calculate needed matrix and vector for predictions
    for (size_t i = 0; i < newSize; ++i) {
      _covMat(i, i) = kernel(_inputs[i], _inputs[i]) + _sigma;
      for (size_t j = i + 1; j < newSize; ++j) {
        _covMat(i, j) = _covMat(j, i) = kernel(_inputs[i], _inputs[j]);
      }
    }
    _covMatInv = _covMat.llt().solve(Eigen::MatrixXd::Identity(newSize, newSize));
    _weights = _covMatInv * _outputs;
  }

  /**
   * Kernel function to describe similarity between two inputs
   * using given hyperparameters.
   * @param input1
   * @param input2
   * @param theta
   * @param dimScale
   * @return
   */
  inline double kernel(const Vector &input1, const Vector &input2, double theta, Eigen::VectorXd dimScale) const {
    Eigen::VectorXd r = static_cast<Eigen::VectorXd>(input1 - input2);
    Eigen::VectorXd rSquared = r.array().square();
    return theta * exp(-rSquared.dot(dimScale));
  }

  /**
   * Kernel function to describe similarity between two inputs.
   * @param input1
   * @param input2
   * @return
   */
  inline double kernel(const Vector &input1, const Vector &input2) const {
    return kernel(input1, input2, _theta, _dimScale);
  }

  /**
   * Calculates the kernel between input and all evidence.
   * @param input
   * @return Vector of covariances
   */
  Eigen::VectorXd kernelVector(const Vector &input) const {
    std::vector<double> k;
    for (auto &d : _inputs) {
      k.push_back(kernel(input, d));
    }
    return Eigen::Map<Eigen::VectorXd>(k.data(), k.size());
  }

  std::vector<Vector> _inputs;
  Eigen::VectorXd _outputs;

  /**
   * input dimensions
   */
  const size_t _dims;

  /**
   * default prior variance
   */
  const double _defaultTheta;
  /**
   * prior variance
   */
  double _theta;
  /**
   * Scale for each input dimension
   */
  Eigen::VectorXd _dimScale;
  /**
   * fixed noise assumed
   */
  const double _sigma;

  Eigen::MatrixXd _covMat;
  Eigen::MatrixXd _covMatInv;
  Eigen::VectorXd _weights;

  Random &_rng;
};
}  // namespace autopas
