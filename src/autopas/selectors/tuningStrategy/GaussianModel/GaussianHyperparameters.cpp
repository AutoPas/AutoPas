/**
 * @file GaussianHyperparameters.cpp
 * @author Jan Nguyen
 * @date 11.05.20
 */

#include "GaussianHyperparameters.h"

#include <Eigen/Cholesky>

#include "GaussianProcess.h"

autopas::GaussianHyperparameters::GaussianHyperparameters(double mean, double theta, Eigen::VectorXd dimScales)
    : score(0), mean(mean), theta(theta), dimScales(std::move(dimScales)) {}

autopas::GaussianHyperparameters::~GaussianHyperparameters() = default;

void autopas::GaussianHyperparameters::precalculate(double sigma, const std::vector<Eigen::VectorXd> &inputs,
                                                    const Eigen::VectorXd &outputs) {
  size_t size = outputs.size();
  // mean of output shifted to zero
  Eigen::VectorXd outputCentered = outputs - mean * Eigen::VectorXd::Ones(size);

  Eigen::MatrixXd covMat(size, size);
  // calculate covariance matrix
  for (size_t i = 0; i < size; ++i) {
    covMat(i, i) = GaussianProcess::kernel(inputs[i], inputs[i], theta, dimScales) + sigma;
    for (size_t j = i + 1; j < size; ++j) {
      covMat(i, j) = covMat(j, i) = GaussianProcess::kernel(inputs[i], inputs[j], theta, dimScales);
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
