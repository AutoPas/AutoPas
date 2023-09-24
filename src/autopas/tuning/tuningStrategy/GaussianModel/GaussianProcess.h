/**
 * @file GaussianProcess.h
 * @author Jan Nguyen
 * @date 17.05.2019
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
  GaussianProcess(size_t dims, double sigma, Random &rngRef);

  virtual ~GaussianProcess();

  /**
   * Change input dimension. Current evidence will be discarded.
   * @param dims
   */
  void setDimension(size_t dims);

  /**
   * Discard all evidence.
   */
  void clear();

  /**
   * Get the number of evidence provided.
   * @return
   */
  [[nodiscard]] size_t numEvidence() const;

  /**
   * Get all currently stored evidence
   * @return pair of inputs and outputs
   */
  [[nodiscard]] std::pair<const std::vector<Vector> &, const Vector &> getEvidence() const;

  /**
   * Provide a input-output pair as evidence.
   * Each evidence improve the quality of future predictions.
   * @param input x
   * @param output f(x)
   * @param tuneHypers if false hyperparemeters need to be set manually
   */
  void addEvidence(const Vector &input, double output, bool tuneHypers);

  /**
   * Get the evidence with the highest output value
   * @return input of max
   */
  [[nodiscard]] const Vector &getEvidenceMax() const;

  /**
   * Try to predict f(x) using the evidence
   * provided so far.
   * @param input x
   * @return expected output of f(x)
   */
  [[nodiscard]] double predictMean(const Vector &input) const;

  /**
   * Get the variance if evidence are ignored.
   * @return
   */
  [[nodiscard]] double getDefaultVar() const;

  /**
   * The variance of the predicted f(x) from predictMean().
   * @param input x
   * @return variance
   */
  [[nodiscard]] double predictVar(const Vector &input) const;

  /**
   * Calculate the probability density of provided output given provided input.
   * @param input
   * @param output
   * @return
   */
  [[nodiscard]] double predictOutputPDF(const Vector &input, double output) const;

  /**
   * Calculate the scaled probability density of provided output given provided input.
   * The probability density is scaled such that the maximum is 1.
   * @param input
   * @param output
   * @return
   */
  [[nodiscard]] double predictOutputScaledPDF(const Vector &input, double output) const;

  /**
   * Calculates the acquisition function for given input.
   * @param af acquisition function a:input->double
   * @param input i
   * @return a(i). This value can be compared with values a(x) of other inputs x to weigh which input would give the
   * most gain if its evidence were provided.
   */
  [[nodiscard]] double calcAcquisition(AcquisitionFunctionOption af, const Vector &input) const;

  /**
   * Find the input in samples which maximizes given aquisition function.
   * TODO: maybe add parameters for hyperparameters of aquisition functions
   * @param af function to maximize
   * @param samples
   * @return
   */
  [[nodiscard]] Vector sampleAquisitionMax(AcquisitionFunctionOption af, const std::vector<Vector> &samples) const;

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
  [[nodiscard]] static std::tuple<std::vector<double>, std::vector<double>,
                                  std::vector<autopas::GaussianProcess::Vector>>
  generateHyperparameterSamples(size_t sampleSize, Random &rng, size_t dims, double sigma, double evidenceMinValue,
                                double evidenceMaxValue);

  /**
   * Get current hyperparameters.
   * @return
   */
  [[nodiscard]] std::vector<GaussianHyperparameters> &getHyperparameters();

  /**
   * Set the hyperparameters: means, theta, dimScale.
   * The samples are scored equal to the probability that given evidence and hyperparameter-sample
   * generates given output. Hyperparameters weights should be normalized.
   * @param sample_means
   * @param sample_thetas
   * @param sample_dimScales
   */
  void setHyperparameters(const std::vector<double> &sample_means, const std::vector<double> &sample_thetas,
                          const std::vector<autopas::GaussianProcess::Vector> &sample_dimScales);

  /**
   * Normalize weights of hyperparameters and truncate lowest weights.
   */
  void normalizeHyperparameters();

 private:
  /**
   * Hyperparameter means, theta and dimScale are randomly generated.
   * The samples are combined using a weighted average. The weight of a sample
   * equals to the probability that given evidence and hyperparameter-sample
   * generates given output. The lowest weights are truncated.
   */
  void tuneHyperparameters();

  /**
   * Kernel function to describe similarity between two inputs
   * using given hyperparameters.
   * @param input1
   * @param input2
   * @param theta
   * @param dimScale
   * @return
   */
  [[nodiscard]] static double kernel(const Vector &input1, const Vector &input2, double theta,
                                     const autopas::GaussianProcess::Vector &dimScale);

  /**
   * Calculates the kernel between input and all evidence.
   * @param input
   * @param theta
   * @param dimScale
   * @return Vector of covariances
   */
  [[nodiscard]] autopas::GaussianProcess::Vector kernelVector(const Vector &input, double theta,
                                                              const autopas::GaussianProcess::Vector &dimScale) const;

  std::vector<Vector> _inputs;
  Vector _outputs;

  /**
   * Number of input dimensions.
   */
  size_t _dims;

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
