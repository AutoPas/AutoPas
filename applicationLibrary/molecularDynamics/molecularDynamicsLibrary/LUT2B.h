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


  LUT2B(int resolution, double cutoffSquared, double sigsquared = 0, double eps = 0, double shift = 0)
      : ParentLUT(resolution, cutoffSquared) {


    epsilon24 = eps;
    sigmaSquared = sigsquared;
    sqrtSigma = std::sqrt(sigmaSquared);
    shift6 = shift;

  }

  LUT2B(int resolution, double cutoffSquared, bool global, bool delay = false, double sigsquared = 0, double eps = 0,
        double shift = 0)
      : ParentLUT(resolution, cutoffSquared) {


    epsilon24 = eps;
    sigmaSquared = sigsquared;
    shift6 = shift;
    sqrtSigma = std::sqrt(sigmaSquared);


    if (resolution == 0) {
      return;
    }

    fill(cutoffSquared, global, delay);
  }
  template <class Functor>
  float retrieveValues(const Functor &functor, float distanceSquared) {
    // if smaller than the smallest distance saved in lut return manually calculated value

if ((!delay && distanceSquared < ((_pointDistance /2 ) + _pointDistance) )|| (delay && (distanceSquared <  (delayValue +_pointDistance)))) {
      return functor.getLUTValues(distanceSquared);
    }
        return getNextNeighbor( distanceSquared);
//    return getLinear(functor,distanceSquared);
  }


  template <class Functor>
  std::array<double, 2> retrieveValues_global(const Functor &functor, float distanceSquared) {

    if (((!delay) && distanceSquared < ((_pointDistance /2 ) + _pointDistance)) || (delay && (distanceSquared < (delayValue + _pointDistance) ))) {
      return {getLUTValues(distanceSquared), calculateGlobalFactor(distanceSquared)};
    }

    //   return getNextNeighbor_global( distanceSquared);
    return getLinear_global(functor,distanceSquared);
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
      fill_global(delay);
    }
  }

  template <class Functor>
  float getLinear(const Functor &functor, float distance) {

    if (std::fmod(distance, _pointDistance) == 0) {

      return _lut2B.at(distance / _pointDistance);
    }

    int lutindexfloor = 0;
    int lutindexceil = 0;
    // index
    if(!delay) {

    lutindexfloor = std::floor((distance - ((_pointDistance /2)+ _pointDistance))/ _pointDistance);
      lutindexceil = std::ceil((distance - ((_pointDistance /2)+ _pointDistance))/ _pointDistance);

    }else{
      lutindexfloor = std::floor((distance - (delayValue + _pointDistance)) / _pointDistance);
      lutindexceil = std::ceil((distance - (delayValue + _pointDistance)) / _pointDistance);
    }

    float Y1=0;

    if (lutindexfloor >= _resolution) {

      Y1 = _lut2B.at(_resolution - 1);
    } else {
      Y1 = _lut2B.at(lutindexfloor);
    }

    float Y2=0;
    if (lutindexceil >= _resolution) {

      Y2 = _lut2B.at(_resolution - 1);
    } else {

      Y2 = _lut2B.at(lutindexceil);
    }

    double X1 = 0;
    double X2 = 0;


    if(!delay) {


      X1 = (lutindexfloor * _pointDistance) +(_pointDistance /2) + _pointDistance;
      X2 = (lutindexceil * _pointDistance) + (_pointDistance /2)+ _pointDistance;
    }else {
      X1 = (lutindexfloor * _pointDistance) + (delayValue + _pointDistance);
      X2 = (lutindexceil * _pointDistance) + (delayValue + _pointDistance);
    }

    if (fabs(Y2 - Y1) < 1e-12) {
      return Y1;
    }

    if (fabs(X2 - X1 ) < 1e-12) {
        return Y1;
    }

    if (fabs(distance - X1) < 1e-12) {
        return Y1;
    }
    return  Y1 + (distance - X1) * ((Y2 - Y1) / (X2 - X1 ));
  }


  float getNextNeighbor( float distance) {
    if (distance >= _cutoffSquared) {
      return _lut2B.at(_numberOfPoints - 1);
    }

    int index = 0;
    if(!delay){

           index = std::floor((distance - ((_pointDistance /2)+ _pointDistance ))/ _pointDistance);
    }else{
       index = std::floor((distance - (_lutCutoff + _pointDistance)) / _pointDistance);
    }


    if (index >= _resolution) {
      return _lut2B.at(_resolution - 1);
    }

    return _lut2B.at(index);
  }



  std::array<double, 2> getNextNeighbor_global( float distance) {
    if (distance >= _cutoffSquared) {

      return _lut2B_globals.at(_numberOfPoints - 1);
    }

    int index = 0;
    if(!delay){
      index = std::floor((distance - ((_pointDistance /2)+ _pointDistance ))/ _pointDistance);
    }else{
       index = std::floor((distance - (_lutCutoff + _pointDistance)) / _pointDistance);
    }
    if (index >= _resolution) {

      return _lut2B_globals.at(_resolution - 1);
    }

    return _lut2B_globals.at(index);

  }

  double calculateGlobalFactor(float distance2) {
    double invdr2 = 1. / distance2;
    double lj6 = sigmaSquared * invdr2;
    lj6 = lj6 * lj6 * lj6;
    double lj12 = lj6 * lj6;
    double lj12m6 = lj12 - lj6;
    return epsilon24 * lj12m6 + shift6;
  }




  void fill_plain(bool delay) {
    if (!delay) {
      for (auto i = 1; i <= _resolution; i++) {
       auto x =(_pointDistance /2) + (i * _pointDistance);

        _lut2B.push_back(getLUTValues((x)));
//    std::cout << "Index  " << i-1 << "  distance :   " << x << std::endl;

      }
    } else {
      fill_delayed_start();
    }
  }

  void fill_delayed_start() {
    delayValue = sqrtSigma * 0.9;

    //delayed start changes the lutcutoff
    _lutCutoff = delayValue;
    _lutDistance = _cutoffSquared - _lutCutoff;
    _lutFactor = _resolution / _lutDistance;
    _pointDistance = _lutDistance / static_cast<double>(_resolution);

    for (auto i = 0; i < _resolution; i++) {
      auto x = (delayValue + _pointDistance) + (i * _pointDistance);
      _lut2B.push_back(getLUTValues(x));
//      std::cout << "Index  " << i << "  distance :   " << x << std::endl;
    }
  }

  void fill_global( bool delay) {

    if(!delay) {
      _lut2B_globals.reserve(_resolution);
      for (auto i = 1; i <= _resolution; i++) {
        auto x =(_pointDistance /2) + (i * _pointDistance);
       _lut2B_globals.push_back({getLUTValues(x), calculateGlobalFactor(x)});
//        std::cout << "Index  " << i-1 << "  distance :   " << x << std::endl;
      }
    }if(delay){
      fill_global_delay();


    }
  }

  void fill_global_delay() {

    delayValue = sqrtSigma * 0.9;

    //delayed start changes the lutcutoff
    _lutCutoff = delayValue;
    _lutDistance = _cutoffSquared - _lutCutoff;
    _lutFactor = _resolution / _lutDistance;
    _pointDistance = _lutDistance / static_cast<double>(_resolution);

    _lut2B_globals.reserve(_resolution);
    for (auto i = 0; i < _resolution; i++) {
      double distance = (delayValue +_pointDistance ) + (i * _pointDistance);
      _lut2B_globals.push_back({getLUTValues(distance), calculateGlobalFactor(distance)});
//      std::cout << "Index  " << i << "  distance :   " << distance << std::endl;

    }
  }





  template <class Functor>
  std::array<double , 2>  getLinear_global(const Functor &functor, float distance) {

    if (std::fmod(distance, _pointDistance) == 0) {

      return _lut2B_globals.at(distance / _pointDistance);
    }

    int lutindexfloor = 0;
    int lutindexceil = 0;
    // index
    if(!delay) {


      lutindexfloor = std::floor((distance - ((_pointDistance /2)+ _pointDistance) )/ _pointDistance);
      lutindexceil = std::ceil((distance - ((_pointDistance /2)+ _pointDistance) )/ _pointDistance);
    }else{
      lutindexfloor = std::floor((distance - (delayValue + _pointDistance)) / _pointDistance);
      lutindexceil = std::ceil((distance -(delayValue + _pointDistance)) / _pointDistance);
    }


    std::array<double ,2 >  Y1={0,0};

    if (lutindexfloor >= _resolution) {

      Y1 = _lut2B_globals.at(_resolution - 1);
    } else {

      Y1 = _lut2B_globals.at(lutindexfloor);
    }

    std::array<double ,2 >  Y2={0,0};
    if (lutindexceil >= _resolution) {

      Y2 = _lut2B_globals.at(_resolution - 1);
    } else {

      Y2 = _lut2B_globals.at(lutindexceil);
    }

    double X1 = 0;
    double X2 = 0;


    if(!delay) {


            X1 = (lutindexfloor * _pointDistance) +(_pointDistance /2) + _pointDistance;
      X2 = (lutindexceil * _pointDistance) +(_pointDistance /2) + _pointDistance;
    }else {
      X1 = (lutindexfloor * _pointDistance) + (delayValue + _pointDistance);
      X2 = (lutindexceil * _pointDistance) + (delayValue + _pointDistance);
    }

    if (fabs(Y2[0] - Y1[0]) < 1e-12) {

      return Y1;
    }

    if (fabs(X2 - X1 ) < 1e-12) {
        return Y1;
    }

    if (fabs(distance - X1) < 1e-12) {

        return Y1;
    }

    auto force = Y1[0] + (distance - X1) * ((Y2[0] - Y1[0]) / (X2 - X1 ));
    auto globs = Y1[1] + (distance - X1) * ((Y2[1] - Y1[1]) / (X2 - X1 ));
    return  {force,globs};
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

 private:
  std::vector<double> _lut2B;
  std::vector<std::array<double, 2>> _lut2B_globals;

  double sigmaSquared;
  double epsilon24;
  double shift6;
  bool delay = false;
  double sqrtSigma;
  double delayValue = 0;
};
}  // namespace mdLib