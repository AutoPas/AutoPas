/**
 * @file GaussianProcess.h
 * @author Jan Nguyen
 * @date 17.05.19
 */

#pragma once

#include <Eigen/Core>
#include <utility>

#include "AcquisitionFunction.h"
#include "GaussianHyperparameters.h"
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
  friend GaussianHyperparameters;

  using Vector = Eigen::VectorXd;

  // number of samples to find optimal hyperparameters
  static constexpr size_t hp_sample_size = 10000;
  // number of hyperparameters
  static constexpr size_t hp_size = 100;

 public:
  /**
   * Constructor
   * @param dims number of input dimensions
   * @param sigma fixed noise
   * @param rngRef reference to rng
   */
  GaussianProcess(size_t dims, double sigma, Random &rngRef)
      : _inputs(),
        _outputs(),
        _dims(dims),
        _evidenceMinValue(0),
        _evidenceMaxValue(0),
        _sigma(sigma),
        _hypers(),
        _rng(rngRef) {
    tuneHyperparameters();
  }

  /**
   * Discard all evidence.
   */
  void clear() {
    _inputs.clear();
    _outputs = Eigen::VectorXd::Zero(0);
    tuneHyperparameters();
  }

  /**
   * Get the number of evidence provided.
   * @return
   */
  [[nodiscard]] size_t numEvidence() const { return _inputs.size(); }

  /**
   * Provide a input-output pair as evidence.
   * Each evidence improve the quality of future predictions.
   * @param input x
   * @param output f(x)
   * @param tuneHypers if false hyperparemeters need to be set manually
   */
  void addEvidence(const Vector &input, double output, bool tuneHypers) {
    if (static_cast<size_t>(input.size()) != _dims) {
      utils::ExceptionHandler::exception("GaussianProcess: size of input {} does not match specified dimensions {}",
                                         input.size(), _dims);
    }

    if (_inputs.empty()) {
      // first evidence
      _evidenceMinValue = _evidenceMaxValue = output;
      _evidenceMaxVector = input;
    } else if (output < _evidenceMinValue) {
      _evidenceMinValue = output;
    } else if (output > _evidenceMaxValue) {
      _evidenceMaxValue = output;
      _evidenceMaxVector = input;
    }

    _inputs.push_back(input);
    long newSize = _inputs.size();

    // extend output vector
    _outputs.conservativeResize(newSize, Eigen::NoChange_t());
    _outputs(newSize - 1) = output;

    if (tuneHypers) {
      tuneHyperparameters();
    } else {
      // hyperparameters should be recalculated
      _hypers.clear();
    }
  }

  /**
   * Get the evidence with the highest output value
   * @return input of max
   */
  [[nodiscard]] Vector getEvidenceMax() {
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
  [[nodiscard]] double predictMean(const Vector &input) const {
    if (static_cast<size_t>(input.size()) != _dims) {
      utils::ExceptionHandler::exception("GaussianProcess: size of input {} does not match specified dimensions {}",
                                         input.size(), _dims);
    }

    double result = 0.;
    if (_inputs.empty()) {
      // no evidence
      for (const auto &hyper : _hypers) {
        result += hyper.score * hyper.mean;
      }
    } else {
      for (const auto &hyper : _hypers) {
        result += hyper.score * (hyper.mean + kernelVector(input, hyper.theta, hyper.dimScales).dot(hyper.weights));
      }
    }

    return result;
  }

  /**
   * The variance of the predicted f(x) from predictMean().
   * @param input x
   * @return variance
   */
  [[nodiscard]] double predictVar(const Vector &input) const {
    if (static_cast<size_t>(input.size()) != _dims) {
      utils::ExceptionHandler::exception("GaussianProcess: size of input {} does not match specified dimensions {}",
                                         input.size(), _dims);
    }

    double result = 0.;
    if (_inputs.empty()) {
      // no evidence
      for (const auto &hyper : _hypers) {
        result += hyper.score * hyper.theta;
      }
    } else {
      for (const auto &hyper : _hypers) {
        Eigen::VectorXd kVec = kernelVector(input, hyper.theta, hyper.dimScales);
        result += hyper.score * (kernel(input, input, hyper.theta, hyper.dimScales) - kVec.dot(hyper.covMatInv * kVec));
      }
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
  [[nodiscard]] inline double calcAcquisition(AcquisitionFunctionOption af, const Vector &input) const {
    return AcquisitionFunction::calcAcquisition(af, predictMean(input), predictVar(input), _evidenceMaxValue);
  }

  /**
   * Find the input in samples which maximizes given aquisition function.
   * TODO: maybe add parameters for hyperparameters of aquisition functions
   * @param af function to maximize
   * @param samples
   * @return
   */
  [[nodiscard]] Vector sampleAquisitionMax(AcquisitionFunctionOption af, const std::vector<Vector> &samples) const {
    size_t bestIdx = 0;
    double bestVal = calcAcquisition(af, samples[0]);

    // find optimmum from samples
    for (size_t i = 1; i < samples.size(); ++i) {
      double val = calcAcquisition(af, samples[i]);

      if (val > bestVal) {
        bestIdx = i;
        bestVal = val;
      }
    }

    return samples[bestIdx];
  }

  /**
   * Generate hyperparameter samples.
   * @param sampleSize size
   * @param rng random number generator
   * @param dims number of input dimension
   * @param sigma fixed noise
   * @param evidenceMinValue current lowest evidence output
   * @param evidenceMaxValue current highest evidence output
   * @return tuple [means, thetas, dimScales]
   */
  [[nodiscard]] static std::tuple<std::vector<double>, std::vector<double>, std::vector<Eigen::VectorXd>>
  generateHyperparameterSamples(size_t sampleSize, Random &rng, size_t dims, double sigma, double evidenceMinValue,
                                double evidenceMaxValue) {
    // range of mean
    // inside bounds of evidence outputs
    NumberInterval<double> meanRange(evidenceMinValue, evidenceMaxValue);
    // range of theta
    // max sample stddev: (max - min)
    // max stddev from zero: abs(min) & abs(max)
    double thetaMax = std::pow(
        std::max({evidenceMaxValue - evidenceMinValue, std::abs(evidenceMinValue), std::abs(evidenceMaxValue)}), 2);
    // at least sigma
    thetaMax = std::max(thetaMax, sigma);
    NumberInterval<double> thetaRange(sigma, thetaMax);
    // range of dimScale
    // Assuming most distances are greater equal 1.
    // For a dimScale d > 5 + ln(thetaMax): theta * exp(-d r) < 1%.
    // So choosing a greater dimScale may lead to many kernels close to zero.
    // But if needed the upper bound can be increased.
    NumberInterval<double> dimScaleRange(0., 5. + std::max(0., std::log(thetaMax)));

    // generate mean
    auto sample_means = meanRange.uniformSample(sampleSize, rng);

    // generate theta
    auto sample_thetas = thetaRange.uniformSample(sampleSize, rng);

    // generate dimScale
    std::vector<std::vector<double>> sample_dimScaleData;
    sample_dimScaleData.reserve(dims);
    for (size_t d = 0; d < dims; ++d) {
      sample_dimScaleData.emplace_back(dimScaleRange.uniformSample(sampleSize, rng));
    }
    // convert dimScales to Vectors
    std::vector<Eigen::VectorXd> sample_dimScales;
    sample_dimScales.reserve(sampleSize);
    for (size_t t = 0; t < sampleSize; ++t) {
      std::vector<double> dimScaleData;
      dimScaleData.reserve(dims);
      for (size_t d = 0; d < dims; ++d) {
        dimScaleData.push_back(sample_dimScaleData[d][t]);
      }
      sample_dimScales.emplace_back(Eigen::Map<Eigen::VectorXd>(dimScaleData.data(), dimScaleData.size()));
    }

    return std::make_tuple(sample_means, sample_thetas, sample_dimScales);
  }

  /**
   * Get current hyperparameters.
   * @return
   */
  [[nodiscard]] std::vector<GaussianHyperparameters> &getHyperparameters() { return _hypers; }

  /**
   * Set the hyperparameters: means, theta, dimScale.
   * The samples are scored equal to the probability that given evidence and hyperparameter-sample
   * generates given output. Hyperparameters weights should be normalized.
   * @param sample_means
   * @param sample_thetas
   * @param sample_dimScales
   */
  void setHyperparameters(const std::vector<double> &sample_means, const std::vector<double> &sample_thetas,
                          const std::vector<Eigen::VectorXd> &sample_dimScales) {
    size_t hyperSize = sample_means.size();
    _hypers.clear();

    // initialize hyperparameter samples
    _hypers.reserve(hyperSize);
    for (size_t t = 0; t < hyperSize; ++t) {
      _hypers.emplace_back(sample_means[t], sample_thetas[t], sample_dimScales[t]);
    }

    // precalculate matrices for all hyperparameters
    // @TODO find sensible chunkSize
#ifdef AUTOPAS_OPENMP
    const size_t chunkSize = std::max(hyperSize / (autopas_get_num_threads() * 10), 1ul);
#pragma omp parallel for schedule(dynamic, chunkSize)
#endif
    for (size_t t = 0; t < hyperSize; ++t) {
      _hypers[t].precalculate(_sigma, _inputs, _outputs);
    }
  }

  /**
   * Normalize weights of hyperparameters and truncate lowest weights.
   */
  void normalizeHyperparameters() {
    // sort by score
    std::sort(_hypers.begin(), _hypers.end(),
              [](const GaussianHyperparameters &h1, const GaussianHyperparameters &h2) { return h1.score > h2.score; });

    // only keep hp_size highest scores
    if (_hypers.size() > hp_size) {
      _hypers.erase(_hypers.begin() + hp_size, _hypers.end());
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

 private:
  /**
   * Hyperparameter means, theta and dimScale are randomly generated.
   * The samples are combined using a weighted average. The weight of a sample
   * equals to the probability that given evidence and hyperparameter-sample
   * generates given output. The lowest weights are truncated.
   */
  void tuneHyperparameters() {
    // number of evidence
    size_t newSize = _inputs.size();
    _hypers.clear();

    // if no evidence
    if (newSize == 0) {
      // use default values
      _hypers.emplace_back(0., 1., Eigen::VectorXd::Ones(_dims));
      _hypers[0].precalculate(_sigma, _inputs, _outputs);
      _hypers[0].score = 1.;
      return;
    }

    auto [sample_means, sample_thetas, sample_dimScales] =
        generateHyperparameterSamples(hp_sample_size, _rng, _dims, _sigma, _evidenceMinValue, _evidenceMaxValue);
    setHyperparameters(sample_means, sample_thetas, sample_dimScales);
    normalizeHyperparameters();
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
  [[nodiscard]] static inline double kernel(const Vector &input1, const Vector &input2, double theta,
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
  [[nodiscard]] Eigen::VectorXd kernelVector(const Vector &input, double theta, const Eigen::VectorXd &dimScale) const {
    std::vector<double> k(_inputs.size());
    for (size_t i = 0; i < k.size(); ++i) {
      k[i] = kernel(input, _inputs[i], theta, dimScale);
    }
    return Eigen::Map<Eigen::VectorXd>(k.data(), k.size());
  }

  std::vector<Vector> _inputs;
  Eigen::VectorXd _outputs;

  /**
   * Number of input dimensions.
   */
  const size_t _dims;

  /**
   * Current smallest evidence output.
   */
  double _evidenceMinValue;
  /**
   * Current greatest evidence output.
   */
  double _evidenceMaxValue;
  /**
   * Current greatest evidence input.
   */
  Vector _evidenceMaxVector;

  /**
   * Fixed noise assumed.
   */
  const double _sigma;

  /**
   * Sampled hyperparameters including precalculated matrices and score
   */
  std::vector<GaussianHyperparameters> _hypers;

  Random &_rng;
};
}  // namespace autopas
