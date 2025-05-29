/**
 * @file Kernel.h
 *
 * @date 29.04.2025
 * @author Luis Gall
 */

#pragma once

namespace mdLib {

template <class CRTP_T>
class Kernel {
 public:
  explicit Kernel(){};

  virtual double calculatePairDerivative(double dr) = 0;

  virtual double calculatePair(double dr) = 0;

  virtual double calculateTriplet(double dr1, double dr2, double dr3) = 0;
};

}  // namespace mdLib