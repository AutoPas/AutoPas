/**
 * @file GaussianHyperparameters.h
 * @author Jan Nguyen
 * @date 11.05.20
 */

#pragma once

#include <Eigen/Core>

namespace autopas {
/**
 * Hyperparameters of Gaussian models and derived matrices used for prediction
 */
class GaussianHyperparameters {
 public:
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
   * Constructor
   * @param mean prior mean of Gaussian process
   * @param theta prior variance
   * @param dimScales scale for each input dimension
   */
  GaussianHyperparameters(double mean, double theta, Eigen::VectorXd dimScales);

  GaussianHyperparameters(const GaussianHyperparameters &gaussianHyperparameters);

  GaussianHyperparameters(GaussianHyperparameters &&gaussianHyperparameters) noexcept;

  GaussianHyperparameters &operator=(const GaussianHyperparameters &gaussianHyperparameters);

  GaussianHyperparameters &operator=(GaussianHyperparameters &&gaussianHyperparameters) noexcept;

  ~GaussianHyperparameters();

  /**
   * Precalculate matrices needed for predictions
   * @param sigma assumed noise
   * @param inputs evidence input
   * @param outputs evidence output
   */
  void precalculate(double sigma, const std::vector<Eigen::VectorXd> &inputs, const Eigen::VectorXd &outputs);
};
}  // namespace autopas
