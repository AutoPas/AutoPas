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


  LUT2B( int resolution, double cutoffSquared) : ParentLUT(resolution, cutoffSquared) {
//    _lut2B.reserve(resolution);



// Open the file
std::filesystem::path cwd = std::filesystem::current_path();
std::cout << "Current working directory: " << cwd << std::endl;

std::ofstream file("logging_distances.csv");
// Check if the file is opened successfully
if (!file.is_open()) {
  std::cout << "Error opening file!" << std::endl;

}
std::cout << "cutoffSquared   " <<cutoffSquared << std::endl;
std::cout << "pointDistance   " <<_pointDistance << std::endl;
std::cout << "distanceSquared" << std::endl;

  }

  template<class Functor>
  float retrieveValues(const Functor &functor,float distanceSquared) {
//    return getNextNeighbor(functor,distanceSquared);


std::cout << std::fixed << std::setprecision(15) << distanceSquared << std::endl;

    return getLinear(functor,distanceSquared);

  }

  template<class Functor>
  void fill(const Functor &functor, double cutoffSquared) {
  for (auto i = 0; i < _resolution; i++) {
    _lut2B.push_back(functor.getLUTValues((_pointDistance / 2) + (i * _pointDistance)));
  }}

  template<class Functor>
  float getLinear(const Functor &functor, float distance){

      if (distance >= _cutoffSquared - _pointDistance / 2) {

        return _lut2B.at(_numberOfPoints - 1);
      }
      if (distance <= (_pointDistance / 2)) {

        return functor.getLUTValues(distance);
      }
      if (std::fmod(distance, _pointDistance) == 0) {

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
        lowerY = _lut2B.at(_numberOfPoints -1 );
      }else {
        lowerY = _lut2B.at(lowerX);
      }





      float upperY;
      if(upperX >= _numberOfPoints){
         upperY = _lut2B.at(_numberOfPoints -1 );
      }else {
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

        return _lut2B.at(_numberOfPoints - 1);
      }
      // Compute perfect value for lowest dr2 where error is greatest
      if (dr2 < _pointDistance / 2) {

        return functor.getLUTValues(dr2);
      }


//                   std::cout << "in next Neighbour \n";

      //TODO: Try multiplication
                   auto res = std::floor(dr2 / _pointDistance);
//                   std::cout << res;
                   if(res >= _resolution){
                     return _lut2B.at(_resolution-1);
                   }
      return _lut2B.at(std::floor(dr2 / _pointDistance));;
    }



 private:
  std::vector<double > _lut2B;






};}