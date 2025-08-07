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

  virtual double calculatePairForce(double dr) = 0;

  virtual double calculatePairPotential(double dr) = 0;

  virtual double calculateTripletPotential(double dr1, double dr2, double dr3) = 0;

  virtual std::array<double, 3> calculateTripletForce(double dr1, double dr2, double dr3) = 0;
};

}  // namespace mdLib