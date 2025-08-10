//
// Created by sliep on 05.05.2025.
//

#pragma once


#include <array>
#include <vector>

#include "autopas/utils/ArrayMath.h"

namespace mdLib {

enum InterpolationMethod{
  NN = 1,
  LIN = 2,
  MP =3
};


 class ParentLUT{
  public:
      ParentLUT(int resolution, double cutoffSquared){
        _resolution = resolution;
        _cutoffSquared = cutoffSquared;
        _cutoff = std::sqrt(cutoffSquared);
      _lutCutoff = cutoffSquared / 10.;
//      _lutCutoff = 2.5;
      _lutDistance = _cutoffSquared - _lutCutoff;
      _lutFactor = _resolution / _lutDistance;
      _pointDistance = _lutDistance / static_cast<double>(_resolution);



      //
      _numberOfPoints= resolution;
    }

    ParentLUT(){}



    virtual void fill() =0;

    protected:


     double _cutoff{};
  double _cutoffSquared{};
  double _lutCutoff{};
  double _lutFactor{};
  double _lutDistance{};
  int _resolution{};
  double _numberOfPoints{};
  double _pointDistance{};

 public:
  void setResolution(int resolution) {
      _resolution = resolution;

      _lutCutoff = _cutoffSquared / 10.;

      _lutDistance = _cutoffSquared - _lutCutoff;
      _lutFactor = _resolution / _lutDistance;
      _pointDistance = _lutDistance / static_cast<double>(_resolution);
      _numberOfPoints= resolution;

  }

  int calculateNumberOfPoints3B(int resolution) {

    return (((resolution+ 1) * (resolution +2) * (resolution+3))/6);
  }
 };}