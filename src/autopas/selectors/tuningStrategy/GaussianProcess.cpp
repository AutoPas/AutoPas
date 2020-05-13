/**
 * @file GaussianProcess.cpp
 * @author seckler
 * @date 07.02.2020
 */

#include "GaussianProcess.h"

#include <Eigen/Cholesky>

void autopas::GaussianProcess::Hyperparameters::precalculate(double sigma, const std::vector<Eigen::VectorXd> &inputs,
                                                             const Eigen::VectorXd &outputs) {
  size_t size = outputs.size();
  // mean of output shifted to zero
  Eigen::VectorXd outputCentered = outputs - mean * Eigen::VectorXd::Ones(size);

  Eigen::MatrixXd covMat(size, size);
  // calculate covariance matrix
  for (size_t i = 0; i < size; ++i) {
    covMat(i, i) = kernel(inputs[i], inputs[i], theta, dimScales) + sigma;
    for (size_t j = i + 1; j < size; ++j) {
      covMat(i, j) = covMat(j, i) = kernel(inputs[i], inputs[j], theta, dimScales);
    }
  }

  // cholesky decomposition
  Eigen::LLT<Eigen::MatrixXd> llt = covMat.llt();
  Eigen::MatrixXd l = llt.matrixL();

  // precalculate inverse of covMat and weights for predictions
  covMatInv = llt.solve(Eigen::MatrixXd::Identity(size, size));
  weights = covMatInv * outputCentered;

  // likelihood of evidence given parameters
  score = std::exp(-0.5 * outputCentered.dot(weights)) / l.diagonal().prod();

  if (std::isnan(score)) {
    // error score calculation failed
    utils::ExceptionHandler::exception("GaussianProcess: invalid score ", score);
  }
}
