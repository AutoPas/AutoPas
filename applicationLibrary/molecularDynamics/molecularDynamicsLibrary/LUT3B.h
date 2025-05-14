//
// Created by sliep on 05.05.2025.
//
#pragma once

#include "ParentLUT.h"


namespace mdLib {
class LUT3B : public ParentLUT{

public:
 LUT3B(int resolution, double cutoffSquared) : ParentLUT(resolution,cutoffSquared) {

    }

      template <class Functor>
  [[nodiscard]] std::pair<const std::array<double, 3>, std::array<u_int8_t, 3>> retrieveValues(const Functor &functor, double dist1, double dist2, double dist3) const {
//    using namespace autopas::utils::ArrayMath::literals;
    if (dist1 < _lutCutoff or dist2 < _lutCutoff or dist3 < _lutCutoff) {
      return std::make_pair(functor.getLUTValues(dist1, dist2, dist3), std::array<u_int8_t, 3>({0, 1, 2}));
    }

    size_t index1, index2, index3;
    std::array<u_int8_t,3> order{};  // Array to track the original indices of the distances

    if (dist1 >= dist2 && dist2 >= dist3) {
      index1 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor );
      index2 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor );
      index3 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor );
      order[0] = 0; order[1] = 1; order[2] = 2;
    } else if (dist1 >= dist3 && dist3 >= dist2) {
      index1 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor );
      index2 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor );
      index3 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor );
      order[0] = 0; order[1] = 2; order[2] = 1;
    } else if (dist2 >= dist1 && dist1 >= dist3) {
      index1 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor );
      index2 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor );
      index3 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor );
      order[0] = 1; order[1] = 0; order[2] = 2;
    } else if (dist2 >= dist3 && dist3 >= dist1) {
      index1 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor );
      index2 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor );
      index3 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor);
      order[0] = 1; order[1] = 2; order[2] = 0;
    } else if (dist3 >= dist1 && dist1 >= dist2) {
      index1 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor);
      index2 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor);
      index3 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor);
      order[0] = 2; order[1] = 0; order[2] = 1;
    } else {
      index1 = static_cast<size_t>((dist3 - _lutCutoff) * _lutFactor);
      index2 = static_cast<size_t>((dist2 - _lutCutoff) * _lutFactor);
      index3 = static_cast<size_t>((dist1 - _lutCutoff) * _lutFactor);
      order[0] = 2; order[1] = 1; order[2] = 0;
    }



//if Linear Interpolation
    auto result = getLinear(index1, index2, index3);



 // space for other interpolation
//auto result = getNextNeighbour(index1, index2, index3);

    auto resultPair = std::make_pair(result, order);

       std::cout << "RETRIEVE VALUE ";

    return resultPair;
  }



   std::array<double, 3> getLinear( size_t index1, size_t index2, size_t index3) const {
    using namespace autopas::utils::ArrayMath::literals;
    // Retrieve values from the LUT by linear interpolation

    auto resultLeft = _lut3B[index1][index2][index3];

    std::array<double, 3> resultRight{} ;

    auto z_max = _lut3B[index1][index2].size();
    auto y_max = _lut3B[index1].size();
    auto x_max = _lut3B.size();

    if(index3 < z_max -1 ){
      resultRight = _lut3B[index1][index2][index3+1];
      return  (resultRight + resultLeft) * 0.5;
    }else if( index2 < y_max-1){
      resultRight = _lut3B[index1][index2+1][index3];
      return  (resultRight + resultLeft) * 0.5;
    }
    resultRight = _lut3B[index1+1][index2][index3];

           std::cout << "in LINEAR ";

    return  (resultRight + resultLeft) * 0.5;
  }

     template<class Functor>
     void fill(const Functor &functor, double cutoffSquared) {


       _lut3B.resize(_resolution + 1);
       for (auto a = 0; a <= _resolution; a++) {
         double distA = _lutCutoff + a * _pointDistance;
         _lut3B[a].resize(a + 1);
         for (auto b = 0; b <= a; b++) {
           double distB = _lutCutoff + b * _pointDistance;
           _lut3B[a][b].reserve(b + 1);
           for (auto c = 0; c <= b; c++) {
             double distC = _lutCutoff + c * _pointDistance;
             _lut3B[a][b].push_back(functor.getLUTValues(distA, distB, distC));

           }
         }
       }

       std::cout << "filled ";
    }

    std::array<double, 3> getNextNeighbour( size_t index1, size_t index2, size_t index3) const {

      //just a function for consistancy, floored values are basically next neighbours
      return _lut3B[index1][index2][index3];

    }







   private:
    std::vector<std::vector<std::vector<std::array<double, 3>>>> _lut3B;


};
}