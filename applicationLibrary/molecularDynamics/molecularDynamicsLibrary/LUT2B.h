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


  void logToFile(const std::string& message) {
    std::ofstream outFile("iteration_10000_NN_small_fallingDrop_distances.txt", std::ios::app); // Open in append mode
    if (outFile.is_open()) {
      outFile << message << std::endl;
    } else {
      std::cerr << "Unable to open file for writing." << std::endl;
    }
  }

  void logIndexToFile(const std::string& message) {
    std::ofstream outFile("iteration_10000_NN_small_fallingDrop_res_100_in_lut.txt", std::ios::app); // Open in append mode
    if (outFile.is_open()) {
      outFile << message << std::endl;
    } else {
      std::cerr << "Unable to open file for writing." << std::endl;
    }
  }


  LUT2B( int resolution, double cutoffSquared) : ParentLUT(resolution, cutoffSquared) {
//    _lut2B.reserve(resolution);

logToFile(std::to_string(cutoffSquared));
logToFile(std::to_string(_pointDistance));



  }

  template<class Functor>
  float retrieveValues(const Functor &functor,float distanceSquared) {

    AutoPasLog(DEBUG, "Retrieving from LUT2B {}", distanceSquared);
logToFile(std::to_string(distanceSquared));
return getNextNeighbor(functor,distanceSquared);

//    return getLinear(functor,distanceSquared);

  }

  template<class Functor>
  void fill(const Functor &functor, double cutoffSquared) {
  for (auto i = 0; i < _resolution; i++) {
//    auto x = (_pointDistance / 2) + (i * _pointDistance);
//    std::cout<< x<< std::endl;
    _lut2B.push_back(functor.getLUTValues((_pointDistance / 2) + (i * _pointDistance)));
  }}

  template<class Functor>
  float getLinear(const Functor &functor, float distance){

      if (distance >= _cutoffSquared - _pointDistance / 2) {
        logIndexToFile(std::to_string(_numberOfPoints-1));
        return _lut2B.at(_numberOfPoints - 1);
      }
      if (distance <= (_pointDistance / 2)) {

        return functor.getLUTValues(distance);
      }
      if (std::fmod(distance, _pointDistance) == 0) {
        logIndexToFile(std::to_string(distance/_pointDistance));
        return _lut2B.at(distance / _pointDistance);
      }

      //index
      float lowerX = std::floor((distance - _pointDistance / 2) / _pointDistance);
//      std::cout<< lowerX << std::endl;

      //index
      float upperX = std::ceil((distance - _pointDistance / 2) / _pointDistance);
//      std::cout<< upperX<< std::endl;


//      const float lowerY = _lut2B.at(lowerX);

      float lowerY;
      if(lowerX >= _numberOfPoints){
        logIndexToFile(std::to_string(_numberOfPoints-1));
        lowerY = _lut2B.at(_numberOfPoints -1 );
      }else {
        logIndexToFile(std::to_string(lowerX));
        lowerY = _lut2B.at(lowerX);
      }





      float upperY;
      if(upperX >= _numberOfPoints){
        logIndexToFile(std::to_string(_numberOfPoints-1));
         upperY = _lut2B.at(_numberOfPoints -1 );
      }else {
        logIndexToFile(std::to_string(upperX));
        upperY = _lut2B.at(upperX);
      }
//      std::cout<< lowerY<< std::endl;
//      std::cout<< upperY<< std::endl;
      lowerX = lowerX * _pointDistance + _pointDistance / 2;
      upperX = upperX * _pointDistance + _pointDistance / 2;


//           std::cout << "in 2b Get/ Lineear \n";
      return lowerY + (distance - lowerX) * (upperY - lowerY) / (upperX - lowerX);

  }

  template<class Functor>
  float getNextNeighbor(const Functor &functor, float dr2) {

      if (dr2 >= _cutoffSquared) {
        logIndexToFile(std::to_string(_numberOfPoints - 1));
        return _lut2B.at(_numberOfPoints - 1);
      }
      // Compute perfect value for lowest dr2 where error is greatest
      if (dr2 < _pointDistance / 2) {

        return functor.getLUTValues(dr2);
      }

//                   std::cout << "in next Neighbour \n";

      //TODO: Try multiplication
                   auto res = std::floor(dr2 / _pointDistance);
////                   std::cout << res;
                   if(res >= _resolution){
                     logIndexToFile(std::to_string(_resolution -1 ));
                     return _lut2B.at(_resolution-1);
                   }
                   logIndexToFile( std::to_string(std::floor(dr2/ _pointDistance)));
      return _lut2B.at(std::floor(dr2 / _pointDistance));;
    }


    //alternative LUT2B implementations


    double calculateGlobalFactor( float distance2){

      double invdr2 = 1. / distance2;
      double lj6 = sigmaSquared * invdr2;
      lj6 = lj6 * lj6 * lj6;
      double lj12 = lj6 * lj6;
      double lj12m6 = lj12 - lj6;
      return  epsilon24 * lj12m6 + shift6;

    }

    template<class Functor>
    void fill_global(const Functor &functor, double cutoffSquared) {
      for (auto i = 0; i < _resolution; i++) {

        double distance = (_pointDistance / 2) + (i * _pointDistance);
        _lut2B_globals[i][0] = functor.getLUTValues(distance);
        _lut2B_globals[i][1] = calculateGlobalFactor(distance);
      }}


    template<class Functor>
    std::array<double, 2> getNextNeighbor_global(const Functor &functor, float dr2) {

      if (dr2 >= _cutoffSquared) {
        //        logIndexToFile(std::to_string(_numberOfPoints - 1));
        return _lut2B_globals.at(_numberOfPoints - 1);
      }
      // Compute perfect value for lowest dr2 where error is greatest
      if (dr2 < _pointDistance / 2) {

        return {functor.getLUTValues(dr2) , calculateGlobalFactor(dr2)};
      }



      //TODO: Try multiplication
      auto res = std::floor(dr2 / _pointDistance);
      ////                   std::cout << res;
      if(res >= _resolution){
        //                     logIndexToFile(std::to_string(_resolution -1 ));
        return _lut2B_globals.at(_resolution-1);
      }
      //                   logIndexToFile( std::to_string(std::floor(dr2/ _pointDistance)));
      return _lut2B_globals.at(std::floor(dr2 / _pointDistance));;
    }


    template<class Functor>
    std::array<double, 2> retrieveValues_global(const Functor &functor,float distanceSquared) {

      AutoPasLog(DEBUG, "Retrieving from LUT2B {}", distanceSquared);
      //logToFile(std::to_string(distanceSquared));
      return getNextNeighbor(functor,distanceSquared);

      //    return getLinear(functor,distanceSquared);

    }




    // alternative spacing LUT

    template<class Functor>
    void fill_unevenly(const Functor &functor, double cutoffSquared) {

      double increment=  _cutoff/ _resolution;

      for (auto i = 0; i < _resolution; i++) {

        auto distance = (increment * i) *(increment * i);
        _lut2B.push_back(functor.getLUTValues(distance));
      }}




 private:
  std::vector<double > _lut2B;
  std::vector<std::array<double, 2>> _lut2B_globals;


  double sigmaSquared;
  double epsilon24;
  double shift6 ;




};}