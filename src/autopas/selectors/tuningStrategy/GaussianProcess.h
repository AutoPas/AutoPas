/**
 * @file GaussianProcess.h
 * @author Jan Nguyen
 * @date 17.05.2019
 */

#pragma once

#include <Eigen/Core>
#include <utility>

#include "autopas/options/AcquisitionFunctionOption.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Math.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/Random.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * Gaussian process is a stochastical model. It predicts the
 * output of a blackbox function f(x) for given input x. To do so, some sample
 * input-output pairs (x,f(x)) should be provided as evidence.
 *
 * Currently the squared exponential kernel is used.
 * TODO: maybe offer some options.
 */
class GaussianProcess {
  using Vector = Eigen::VectorXd;

  // number of samples to find optimal hyperparameters
  static constexpr size_t hp_sample_size = 10000;
  // number of hyperparameters
  static constexpr size_t hp_size = 100;

  /**
   * Hyperparameters and derived matrices used for prediction
   */
  struct Hyperparameters {
    /**
     * score used to weight hyperparameters
     */
    double score;
    /**
     * prior mean
     */
    double mean;
    /**
     * prior variance
     */
    double theta;
    /**
     * Scale for each input dimension
     */
    Eigen::VectorXd dimScales;

    /**
     * Covariance Matrix Inverse
     */
    Eigen::MatrixXd covMatInv;
    /**
     * Weights used for predictions
     */
    Eigen::VectorXd weights;

    /**
     * Default Constructor
     */
    Hyperparameters()
        : mean(0.),
          theta(1.),
          dimScales(Eigen::VectorXd::Ones(1)),
          covMatInv(Eigen::MatrixXd::Ones(1, 1)),
          weights(Eigen::VectorXd::Ones(1)) {}

    /**
     * Constructor
     * @param mean prior mean of Gaussian process
     * @param theta prior variance
     * @param dimScales scale for each input dimension
     */
    Hyperparameters(double mean, double theta, Eigen::VectorXd dimScales)
        : mean(mean), theta(theta), dimScales(std::move(dimScales)) {}

    /**
     * Precalculate matrices needed for predictions
     * @param sigma assumed noise
     * @param inputs evidence input
     * @param outputs evidence output
     */
    void precalculate(double sigma, const std::vector<Eigen::VectorXd> &inputs, const Eigen::VectorXd &outputs);
  };

 public:
  /**
   * Constructor
   * @param dims number of input dimensions
   * @param sigma fixed noise
   * @param rngRef reference to rng
   */
  GaussianProcess(size_t dims, double sigma, Random &rngRef)
      : _inputs(), _outputs(), _dims(dims), _sigma(sigma), _hypers(), _rng(rngRef) {}

  /**
   * Discard all evidence.
   */
  void clear() { _inputs.clear(); }

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

    if (_inputs.empty()) {
      // first evidence
      _evidenceMinValue = _evidenceMaxValue = output;
      _evidenceMinVector = _evidenceMaxVector = input;
    } else if (output < _evidenceMinValue) {
      _evidenceMinValue = output;
      _evidenceMinVector = input;
    } else if (output > _evidenceMaxValue) {
      _evidenceMaxValue = output;
      _evidenceMaxVector = input;
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

    return _evidenceMinVector;
  }

  /**
   * Get the evidence with the highest output value
   * @return input of max
   */
  Vector getEvidenceMax() {
    if (_inputs.empty()) {
      utils::ExceptionHandler::exception("GaussianProcess has no evidence");
    }

    return _evidenceMaxVector;
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

    // default mean 0.
    if (_inputs.size() == 0) return 0.;

    double result = 0.;
    for (auto &hyper : _hypers) {
      result += hyper.score * (hyper.mean + kernelVector(input, hyper.theta, hyper.dimScales).dot(hyper.weights));
    }

    return result;
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

    // default variance 1.
    if (_inputs.size() == 0) return 1.;

    double result = 0.;
    for (auto &hyper : _hypers) {
      Eigen::VectorXd kVec = kernelVector(input, hyper.theta, hyper.dimScales);
      result += hyper.score * (kernel(input, input, hyper.theta, hyper.dimScales) - kVec.dot(hyper.covMatInv * kVec));
    }
    return result;
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
      case AcquisitionFunctionOption::upperConfidenceBound: {
        return predictMean(input) + 2 * std::sqrt(predictVar(input));
      }
      case AcquisitionFunctionOption::lowerConfidenceBound: {
        return predictMean(input) - 2 * std::sqrt(predictVar(input));
      }
      case AcquisitionFunctionOption::mean: {
        return predictMean(input);
      }
      case AcquisitionFunctionOption::variance: {
        return predictVar(input);
      }
      case AcquisitionFunctionOption::probabilityOfDecrease: {
        return utils::Math::normalCDF((_evidenceMinValue - predictMean(input)) / std::sqrt(predictVar(input)));
      }
      case AcquisitionFunctionOption::expectedDecrease: {
        double mean = predictMean(input);
        double stddev = std::sqrt(predictVar(input));
        double minNormed = (_evidenceMinValue - mean) / stddev;
        return (_evidenceMinValue - mean) * utils::Math::normalCDF(minNormed) +
               stddev * utils::Math::normalPDF(minNormed);
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
    _hypers.clear();

    // if no evidence
    if (newSize == 0) {
      // use default values
      return;
    }

    if (newSize == 1) {
      // default values for one evidence
      _hypers.emplace_back(_outputs[0], _outputs[0] * _outputs[0], Eigen::VectorXd::Ones(_dims));
      _hypers[0].precalculate(_sigma, _inputs, _outputs);
    } else {
      // range of mean
      // inside bounds of evidence outputs
      NumberInterval<double> meanRange(_evidenceMinValue, _evidenceMaxValue);
      // range of theta
      // max sample stddev: (max - min)
      // max stddev from zero: abs(min) & abs(max)
      double thetaMax = std::pow(
          std::max({_evidenceMaxValue - _evidenceMinValue, std::abs(_evidenceMinValue), std::abs(_evidenceMaxValue)}),
          2);
      // at least sigma
      thetaMax = std::max(thetaMax, _sigma);
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

      // initialize hyperparameter samples
      _hypers.reserve(hp_sample_size);
      for (size_t t = 0; t < hp_sample_size; ++t) {
        _hypers.emplace_back(sample_means[t], sample_thetas[t], sample_dimScales[t]);
      }

      // precalculate matrices for all hyperparameters
      // @TODO find sensible chunkSize
#ifdef AUTOPAS_OPENMP
      const size_t chunkSize = std::max(hp_sample_size / (autopas_get_num_threads() * 10), 1ul);
#pragma omp parallel for schedule(dynamic, chunkSize)
#endif
      for (size_t t = 0; t < hp_sample_size; ++t) {
        _hypers[t].precalculate(_sigma, _inputs, _outputs);
      }

      // sort by score
      std::sort(_hypers.begin(), _hypers.end(),
                [](const Hyperparameters &h1, const Hyperparameters &h2) { return h1.score > h2.score; });

      // only keep top 'hp_size' hyperparameters
      _hypers.resize(hp_size);
    }

    // normalize scores
    double scoreSum = 0.;
    for (auto &hyper : _hypers) {
      scoreSum += hyper.score;
    }
    for (auto &hyper : _hypers) {
      hyper.score /= scoreSum;
    }
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
  static inline double kernel(const Vector &input1, const Vector &input2, double theta,
                              const Eigen::VectorXd &dimScale) {
    double dot = 0;
    for (int i = 0; i < input1.size(); ++i) {
      double dist = input1[i] - input2[i];
      dist *= dist * dimScale[i];
      dot += dist;
    }
    return theta * std::exp(-dot);
  }

  /**
   * Calculates the kernel between input and all evidence.
   * @param input
   * @return Vector of covariances
   */
  Eigen::VectorXd kernelVector(const Vector &input, double theta, const Eigen::VectorXd &dimScale) const {
    std::vector<double> k(_inputs.size());
    for (size_t i = 0; i < k.size(); ++i) {
      k[i] = kernel(input, _inputs[i], theta, dimScale);
    }
    return Eigen::Map<Eigen::VectorXd>(k.data(), k.size());
  }

  std::vector<Vector> _inputs;
  Eigen::VectorXd _outputs;

  /**
   * Current smallest evidence output
   */
  double _evidenceMinValue;
  /**
   * Current smallest evidence input
   */
  Vector _evidenceMinVector;
  /**
   * Current greatest evidence output
   */
  double _evidenceMaxValue;
  /**
   * Current greatest evidence input
   */
  Vector _evidenceMaxVector;

  /**
   * input dimensions
   */
  const size_t _dims;

  /**
   * fixed noise assumed
   */
  const double _sigma;

  /**
   * hyperparameters
   */
  std::vector<Hyperparameters> _hypers;

  Random &_rng;
};
}  // namespace autopas
