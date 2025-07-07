//
// Created by sliep on 05.05.2025.
//

#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>

#include "ParentLUT.h"

namespace mdLib {
class LUT2B : public ParentLUT {
 public:
  void logToFile(const std::string &message) {
    std::ofstream outFile("iteration_10000_NN_small_fallingDrop_distances.txt", std::ios::app);  // Open in append mode
    if (outFile.is_open()) {
      outFile << message << std::endl;
    } else {
      std::cerr << "Unable to open file for writing." << std::endl;
    }
  }

  void logIndexToFile(const std::string &message) {
    std::ofstream outFile("iteration_10000_NN_small_fallingDrop_res_100_in_lut.txt",
                          std::ios::app);  // Open in append mode
    if (outFile.is_open()) {
      outFile << message << std::endl;
    } else {
      std::cerr << "Unable to open file for writing." << std::endl;
    }
  }

  LUT2B(int resolution, double cutoffSquared, double sigsquared = 0, double eps = 0, double shift = 0)
      : ParentLUT(resolution, cutoffSquared) {
    //    _lut2B.reserve(resolution);

    epsilon24 = eps;
    sigmaSquared = sigsquared * sigsquared;
    sqrtSigma = std::sqrt(sigmaSquared);
    shift6 = shift;
    // logToFile(std::to_string(cutoffSquared));
    // logToFile(std::to_string(_pointDistance));
  }

  LUT2B(int resolution, double cutoffSquared, bool global, bool delay = false, double sigsquared = 0, double eps = 0,
        double shift = 0)
      : ParentLUT(resolution, cutoffSquared) {
    //    _lut2B.reserve(resolution);

    epsilon24 = eps;
    sigmaSquared = sigsquared;
    shift6 = shift;
    // logToFile(std::to_string(cutoffSquared));
    // logToFile(std::to_string(_pointDistance));

    if (resolution == 0) {
      return;
    }

    fill(cutoffSquared, global, delay);
  }
  template <class Functor>
  float retrieveValues(const Functor &functor, float distanceSquared) {
    // if smaller than the smallest distance saved in lut return manually calculated value
    if ((!delay && distanceSquared < (_pointDistance / 2)) || (delay && (distanceSquared < (sqrtSigma * 0.9)))) {
      return getLUTValues(distanceSquared);
    }

    // return getNextNeighbor(functor, distanceSquared);

       return getLinear(functor,distanceSquared);
  }

  template <class Functor>
  void fill(const Functor &functor, double cutoffSquared, bool globsLUT) {
    if (!globsLUT) {
      fill_plain(functor, cutoffSquared);
    }
    if (globsLUT) {
      fill_global(functor, cutoffSquared);
    }
  }

  void fill(double cutoffSquared, bool globsLUT, bool _delay) {
    delay = _delay;

    if (!globsLUT) {
      fill_plain(delay);
    }
    if (globsLUT) {
      fill_global(cutoffSquared);
    }
  }

  template <class Functor>
  float getLinear(const Functor &functor, float distance) {

    // std::cout<< "distance  " << distance << std::endl;
    if (distance >= _cutoffSquared - _pointDistance / 2) {

      return _lut2B.at(_resolution -1);
    }
    if (distance <= (_pointDistance / 2)) {
      return functor.getLUTValues(distance);
    }
    if (std::fmod(distance, _pointDistance) == 0) {
      //        logIndexToFile(std::to_string(distance/_pointDistance));
      return _lut2B.at(distance / _pointDistance);
    }

    // index
    float lowerX = std::floor((distance - _pointDistance / 2) / _pointDistance);
    // std::cout<< lowerX << std::endl;

    if (lowerX >= _resolution) {

      return _lut2B.at(_resolution -1);
    }

    // index
    float upperX = std::ceil((distance - _pointDistance / 2) / _pointDistance);
    // std::cout<< upperX<< std::endl;

    //      const float lowerY = _lut2B.at(lowerX);

    float lowerY;
    if (lowerX >= _numberOfPoints) {
      //        logIndexToFile(std::to_string(_numberOfPoints-1));
      lowerY = _lut2B.at(_numberOfPoints - 1);
    } else {
      //        logIndexToFile(std::to_string(lowerX));
      lowerY = _lut2B.at(lowerX);
    }

    float upperY;
    if (upperX >= _numberOfPoints) {
      //        logIndexToFile(std::to_string(_numberOfPoints-1));
      upperY = _lut2B.at(_numberOfPoints - 1);
    } else {
      //        logIndexToFile(std::to_string(upperX));
      upperY = _lut2B.at(upperX);
    }
    // std::cout<< lowerY<< std::endl;
    // std::cout<< upperY<< std::endl;
    lowerX = lowerX * _pointDistance + _pointDistance / 2;
    upperX = upperX * _pointDistance + _pointDistance / 2;


    // std::cout<<"to dist lower " << lowerX << std::endl;
    // std::cout<<"to dist upper " << upperX << std::endl;
    auto res = lowerY + (distance - lowerX) * (upperY - lowerY) / (upperX - lowerX);
    // std::cout << "--------------------------------------------------------res   " <<res <<  std::endl;
    return lowerY + (distance - lowerX) * (upperY - lowerY) / (upperX - lowerX);
  }

  template <class Functor>
  float getNextNeighbor(const Functor &functor, float dr2) {
    if (dr2 >= _cutoffSquared) {
      //        logIndexToFile(std::to_string(_numberOfPoints - 1));
      return _lut2B.at(_numberOfPoints - 1);
    }
    // Compute perfect value for lowest dr2 where error is greatest
    if (dr2 < _pointDistance / 2) {
      return functor.getLUTValues(dr2);
    }

    //                   std::cout << "in next Neighbour \n";

    // TODO: Try multiplication
    auto res = std::floor(dr2 / _pointDistance);
    ////                   std::cout << res;
    if (res >= _resolution) {
      //                     logIndexToFile(std::to_string(_resolution -1 ));
      return _lut2B.at(_resolution - 1);
    }
    //                   logIndexToFile( std::to_string(std::floor(dr2/ _pointDistance)));
    return _lut2B.at(std::floor(dr2 / _pointDistance));
    ;
  }

  // alternative LUT2B implementations

  double calculateGlobalFactor(float distance2) {
    double invdr2 = 1. / distance2;
    double lj6 = sigmaSquared * invdr2;
    lj6 = lj6 * lj6 * lj6;
    double lj12 = lj6 * lj6;
    double lj12m6 = lj12 - lj6;
    return epsilon24 * lj12m6 + shift6;
  }

  // TODO maybe just delete cutoff squared from the parameters cause you dont use it?
  template <class Functor>
  void fill_plain(const Functor &functor, bool delay) {
    if (!delay) {
      for (auto i = 0; i < _resolution; i++) {
        auto x = (_pointDistance / 2) + (i * _pointDistance);
        // std::cout << x << std::endl;
        _lut2B.push_back(functor.getLUTValues((_pointDistance / 2) + (i * _pointDistance)));
      }
    } else {
      fill_delayed_start();
    }
  }

  template <class Functor>
  void fill_delayed_start(const Functor &functor, double cutoffSquared) {
    double delay = std::sqrt(sigmaSquared) * 0.9;
    for (auto i = 0; i < _resolution; i++) {
      auto x = delay + (i * _pointDistance);
      // std::cout << x << std::endl;
      _lut2B.push_back(functor.getLUTValues((_pointDistance / 2) + (i * _pointDistance)));
    }
  }

  template <class Functor>
  void fill_global(const Functor &functor, double cutoffSquared) {
    _lut2B_globals.reserve(_resolution);
    for (auto i = 0; i < _resolution; i++) {
      double distance = (_pointDistance / 2) + (i * _pointDistance);
      _lut2B_globals.push_back({functor.getLUTValues(distance), calculateGlobalFactor(distance)});

    }
  }

  void fill_plain(bool delay) {
    if (!delay) {
      for (auto i = 0; i < _resolution; i++) {
        auto x = (_pointDistance / 2) + (i * _pointDistance);
        // std::cout << x << std::endl;
        _lut2B.push_back(getLUTValues((_pointDistance / 2) + (i * _pointDistance)));
      }
    } else {
      fill_delayed_start();
    }
  }

  void fill_delayed_start() {
    double delay = std::sqrt(sigmaSquared) * 0.9;
    for (auto i = 0; i < _resolution; i++) {
      auto x = delay + (i * _pointDistance);
      // std::cout << x << std::endl;
      _lut2B.push_back(getLUTValues((_pointDistance / 2) + (i * _pointDistance)));
    }
  }

  void fill_global(double cutoffSquared) {
    _lut2B_globals.reserve(_resolution);
    for (auto i = 0; i < _resolution; i++) {
      double distance = (_pointDistance / 2) + (i * _pointDistance);
      _lut2B_globals.push_back({getLUTValues(distance), calculateGlobalFactor(distance)});

    }
  }

  template <class Functor>
  std::array<double, 2> getNextNeighbor_global(const Functor &functor, float dr2) {
    if (dr2 >= _cutoffSquared) {
      //        logIndexToFile(std::to_string(_numberOfPoints - 1));
      return _lut2B_globals.at(_numberOfPoints - 1);
    }
    // Compute perfect value for lowest dr2 where error is greatest
    if (dr2 < _pointDistance / 2) {
      return {functor.getLUTValues(dr2), calculateGlobalFactor(dr2)};
    }

    // TODO: Try multiplication
    auto res = std::floor(dr2 / _pointDistance);
    ////                   std::cout << res;
    if (res >= _resolution) {
      //                     logIndexToFile(std::to_string(_resolution -1 ));
      return _lut2B_globals.at(_resolution - 1);
    }
    //                   logIndexToFile( std::to_string(std::floor(dr2/ _pointDistance)));
    return _lut2B_globals.at(std::floor(dr2 / _pointDistance));
    ;
  }

  template <class Functor>
  std::array<double, 2> retrieveValues_global(const Functor &functor, float distanceSquared) {
    if ((!delay && distanceSquared < (_pointDistance / 2)) || (delay && (distanceSquared < (sqrtSigma * 0.9)))) {
      return {getLUTValues(distanceSquared), calculateGlobalFactor(distanceSquared)};
    }

    return getNextNeighbor_global(functor, distanceSquared);

    //    return getLinear(functor,distanceSquared);
  }

  // alternative spacing LUT

  template <class Functor>
  void fill_unevenly(const Functor &functor, double cutoffSquared) {
    double increment = _cutoff / _resolution;

    for (auto i = 0; i < _resolution; i++) {
      auto distance = (increment * i) * (increment * i);
      _lut2B.push_back(functor.getLUTValues(distance));
    }
  }

  void setSigmaSquared(double sigmaSquared) { LUT2B::sigmaSquared = sigmaSquared; }

  void setEpsilon24(double epsilon24) { LUT2B::epsilon24 = epsilon24; }

  void setShift6(double shift6) { LUT2B::shift6 = shift6; }

  float getLUTValues(double distance2) const {
    // TOOD why is this hardcoded get the actual one from the mdconfig. save it in the parentLUT
    // this is force magnitude for true r not r^2 even though distance= r^2
    //    auto sigmaSquared = 1;
    //      auto sigmaSquared = sigmaSquared;
    //      auto epsilon24 = epsilon24;
    //    auto shift6 = _shift6AoS;

    double invdr2 = 1. / distance2;
    double lj6 = sigmaSquared * invdr2;
    lj6 = lj6 * lj6 * lj6;
    double lj12 = lj6 * lj6;
    double lj12m6 = lj12 - lj6;
    return epsilon24 * (lj12 + lj12m6) * invdr2;
  }

 private:
  std::vector<double> _lut2B;
  std::vector<std::array<double, 2>> _lut2B_globals;

  double sigmaSquared;
  double epsilon24;
  double shift6;
  bool delay = false;
  double sqrtSigma;
};
}  // namespace mdLib