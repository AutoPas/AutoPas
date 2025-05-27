/**
 * @file PairwiseKernel.h
 *
 * @date 29.04.2025
 * @author Luis Gall
 */

#pragma once

namespace mdLib {

template <class CRTP_T>
class PairwiseKernel {
 public:
  explicit PairwiseKernel(){};

  virtual double calculateDerivative(double dr) = 0;

  virtual double calculate(double dr) = 0;
};

}  // namespace mdLib