//
// Created by sliep on 05.05.2025.
//

#pragma once


#include <array>
#include <vector>

#include "autopas/utils/ArrayMath.h"

namespace mdLib {


 class ParentLUT{
  public:
      ParentLUT(int resolution, double cutoffSquared){
        _resolution = resolution;
        _cutoffSquared = cutoffSquared;
        _cutoff = std::sqrt(cutoffSquared);
      _lutCutoff = cutoffSquared / 10.;
      _lutDistance = _cutoffSquared - _lutCutoff;
      _lutFactor = _resolution / _lutDistance;
      _pointDistance = _lutDistance / static_cast<double>(_resolution);



      //
      _numberOfPoints= resolution;
    }

    ParentLUT(){}


//     template<class Functor>
//     void fill(const Functor &functor, double cutoffSquared) {}

    protected:

  //Different vector for 2b and 3b functor
//  std::vector<std::vector<std::vector<std::array<double, 3>>>> _lut3B;
//  std::vector<double> _lut2B;

     double _cutoff{};
  double _cutoffSquared{};
  double _lutCutoff{};
  double _lutFactor{};
  double _lutDistance{};
  int _resolution{};
  bool _isFilled = false;
  double _numberOfPoints{};
  double _pointDistance{};





};}