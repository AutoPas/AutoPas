/**
 * @file GaussianProcess.h
 * @author Jan Nguyen
 * @date 17.05.19
 */

#pragma once

#include <Eigen/Dense>
#include "autopas/options/AcquisitionFunctionOption.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberSet.h"
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
  // number of samples to find optimal hyperparameters
  static constexpr size_t hp_sample_size = 10000;

 public:
  /**
   * Constructor
   * @param dims number of input dimensions
   * @param sigma fixed noise
   * @param rngRef reference to rng
   */
  GaussianProcess(size_t dims, double sigma, Random &rngRef)
      : _inputs(), _outputs(), _dims(dims), _sigma(sigma), _covMat(), _covMatInv(), _weights(), _rng(rngRef) {
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
    Vector minVec = _inputs[0];

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
    Vector maxVec = _inputs[0];

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

    if (_inputs.size() == 0) return _mean;

    return _mean + kernelVector(input).dot(_weights);
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
    // number of evidence
    size_t newSize = _inputs.size();
    _covMat.resize(newSize, newSize);

    // if no evidence
    if (newSize == 0) {
      // use default values
      _mean = 0;
      _theta = 1.;
      _dimScale = Eigen::VectorXd::Ones(_dims);
      return;
    }

    if (newSize == 1) {
      // default values for one evidence
      _mean = _outputs[0];
      _theta = _mean * _mean;
      _dimScale = Eigen::VectorXd::Ones(_dims);
    } else {
      // output bounds
      double sampleMin = _outputs.minCoeff();
      double sampleMax = _outputs.maxCoeff();

      // range of mean
      // inside bounds of evidence outputs
      NumberInterval<double> meanRange(sampleMin, sampleMax);
      // range of theta
      // max sample stddev: (max - min)
      // max stddev from zero: abs(min) & abs(max)
      double thetaMax = std::pow(std::max({sampleMax - sampleMin, std::abs(sampleMin), std::abs(sampleMax)}), 2);
      NumberInterval<double> thetaRange(_sigma, thetaMax);
      // range of dimScale
      // Assuming most distances are greater equal 1.
      // For a dimScale d > 5 + ln(thetaMax): theta * exp(-d r) < 1%. So choosing
      // a greater dimScale may lead to many kernels close to zero.
      // But if needed the upper bound can be increased.
      NumberInterval<double> dimScaleRange(0., 5. + std::max(0., std::log(thetaMax)));

      // generate mean
      auto sample_means = meanRange.uniformSample(hp_sample_size, _rng);

      // generate theta
      auto sample_thetas = thetaRange.uniformSample(hp_sample_size, _rng);

      // generate dimScale
      std::vector<std::vector<double>> sample_dimScaleData;
      for (size_t d = 0; d < _dims; ++d) {
        sample_dimScaleData.push_back(dimScaleRange.uniformSample(hp_sample_size, _rng));
      }
      // convert dimScales to Vectors
      std::vector<Eigen::VectorXd> sample_dimScales;
      for (size_t t = 0; t < hp_sample_size; ++t) {
        std::vector<double> dimScaleData;
        for (size_t d = 0; d < _dims; ++d) {
          dimScaleData.push_back(sample_dimScaleData[d][t]);
        }
        Eigen::VectorXd dimScale = Eigen::Map<Eigen::VectorXd>(dimScaleData.data(), dimScaleData.size());
        sample_dimScales.push_back(std::move(dimScale));
      }

      double scoreMax = std::numeric_limits<double>::lowest();
      // @TODO parallelize?
      // find max likelyhood in samples
      for (size_t t = 0; t < hp_sample_size; ++t) {
        double mean = sample_means[t];
        double theta = sample_thetas[t];
        Eigen::VectorXd dimScale = sample_dimScales[t];

        // mean of output shifted to zero
        Eigen::VectorXd outputCentered = _outputs - mean * Eigen::VectorXd::Ones(newSize);

        // calculate covariance matrix
        for (size_t i = 0; i < newSize; ++i) {
          _covMat(i, i) = kernel(_inputs[i], _inputs[i], theta, dimScale) + _sigma;
          for (size_t j = i + 1; j < newSize; ++j) {
            _covMat(i, j) = _covMat(j, i) = kernel(_inputs[i], _inputs[j], theta, dimScale);
          }
        }

        // cholesky decomposition
        Eigen::LLT<Eigen::MatrixXd> llt = _covMat.llt();
        Eigen::MatrixXd l = llt.matrixL();

        // log likelihood of evidence given parameters
        double score = -0.5 * outputCentered.dot(llt.solve(outputCentered)) - std::log(l.diagonal().prod());

        if (std::isnan(score)) {
          // error score calculation failed
          utils::ExceptionHandler::exception("GaussianProcess: invalid score ", score);
        }

        if (score > scoreMax) {
          scoreMax = score;
          _mean = mean;
          _theta = theta;
          _dimScale = dimScale;
        }
      }
    }

    // calculate needed matrix and vector for predictions
    for (size_t i = 0; i < newSize; ++i) {
      _covMat(i, i) = kernel(_inputs[i], _inputs[i]) + _sigma;
      for (size_t j = i + 1; j < newSize; ++j) {
        _covMat(i, j) = _covMat(j, i) = kernel(_inputs[i], _inputs[j]);
      }
    }
    // mean of output shifted to zero
    _covMatInv = _covMat.llt().solve(Eigen::MatrixXd::Identity(newSize, newSize));
    _weights = _covMatInv * (_outputs - _mean * Eigen::VectorXd::Ones(newSize));
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
   * prior mean
   */
  double _mean;
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
