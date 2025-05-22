//
// Created by sliep on 05.05.2025.
//
#pragma once

#include "ParentLUT.h"
#include "autopas/utils/Math.h"


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
//    auto result = getLinear(index1, index2, index3);



 // space for other interpolation
auto result = getNextNeighbour(index1, index2, index3);

    auto resultPair = std::make_pair(result, order);

       std::cout << "RETRIEVE VALUE ";

    return resultPair;
  }



  std::array<double, 3> getWeightedAverage( size_t index1, size_t index2, size_t index3, size_t actual_dist0, size_t actual_dist1, size_t actual_dist2) const {

    auto V_Middle = _lut3B[index1][index2][index3];
    auto distances_middle = getDistances(index1, index2, index3);



    std::array<double, 3> V_Above{};
    std::array<double, 3> distances_above;
    auto z_max = _lut3B[index1][index2].size();
    auto y_max = _lut3B[index1].size();
    auto x_max = _lut3B.size();





    //index is enough doesnt need specific distance since the can be directly mapped to each other and bigger index always == bigger distance
    //no use the actual distances cause need to comapre to transported distances
    if(index3 < z_max -1 ){
      V_Above = _lut3B[index1][index2][index3+1];
       distances_above = getDistances(index1, index2, index3+1);

    }else if( index2 < y_max-1){
      V_Above = _lut3B[index1][index2+1][index3];
       distances_above = getDistances(index1, index2+1, index3);

    }else {

      V_Above = _lut3B[index1 + 1][index2][index3];
       distances_above = getDistances(index1+1, index2, index3);
    }






    std::array<double, 3> V_Below{};
    std::array<double, 3> distances_below;



    std::array<double, 3> resultBelow {};
    if(index3 > 0 ){
      resultBelow = _lut3B[index1][index2][index3-1];
      distances_below = getDistances(index1, index2, index3-1);

    }else if( index2 >0){
      resultBelow = _lut3B[index1][index2-1][index3];
      distances_below = getDistances(index1, index2-1, index3);

    }else {
      if (index1> 0){
      resultBelow = _lut3B[index1 - 1][index2][index3];
      distances_below = getDistances(index1-1, index2, index3);


      }
      else{resultBelow = _lut3B[0][0][0];
        distances_below = getDistances(0, 0, 0);}
    }





    //weights
    auto a = pow(actual_dist0 -distances_middle[0] , 2) + pow(actual_dist1-distances_middle[1] , 2) + pow(actual_dist2-distances_middle[2] , 2) ;
    auto b = pow(actual_dist0 -distances_below[0] , 2) + pow(actual_dist1-distances_below[1] , 2) + pow(actual_dist2-distances_below[2] , 2) ;
    auto c = pow(actual_dist0 -distances_above[0] , 2) + pow(actual_dist1-distances_above[1] , 2) + pow(actual_dist2-distances_above[2] , 2) ;


    auto w1 = 1 /a;
    auto w2 = 1 /b;
    auto w3 = 1 /c;

    auto w_total = w1 +w2 +w3;

    w1= w1 / w_total;
    w2= w2 / w_total;
    w3= w3 / w_total;


    auto v0 = (V_Middle[0]* w1) +(V_Below[0] * w2)+ (V_Above[0]* w3);
    auto v1 = (V_Middle[1]* w1) +(V_Below[1] * w2)+ (V_Above[1]* w3);
    auto v2= (V_Middle[2]* w1) +(V_Below[2] * w2)+ (V_Above[2]* w3);

    return std::array<double, 3> {v0, v1, v2};



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


    std::array<double, 3> getHybrid( size_t index1, size_t index2, size_t index3) const {

      auto resultMiddle = _lut3B[index1][index2][index3];



      auto resultAbove = getNextNeighbour(index1, index2, index3);



      std::array<double, 3> resultBelow {};
      if(index3 > 0 ){
        resultBelow = _lut3B[index1][index2][index3-1];

      }else if( index2 >0){
        resultBelow = _lut3B[index1][index2-1][index3];

      }else {
        resultBelow = _lut3B[index1 - 1][index2][index3];
      }



      // V = alpha *
      double  alpha, beta,  gamma;
      //resultAbove = alpha * distance[1] + beta * distance[2] + gamma * distance[3] + delta
      //note : distance has to be recalulated with the lutfactor like in retrieve value for the other ones +/ - pointdistance
      //the eigen lib here might already have all the stuff for solving the LGS
      Eigen::Matrix<double, 3,3> matrix;
      Eigen::Vector3d v;
//      matrix <<
      //wait does this even make sense cause V in lut is vector of 3





    }






   private:
    std::vector<std::vector<std::vector<std::array<double, 3>>>> _lut3B;


    //helper functions

    std::array<double, 3> getDistances(size_t i1,size_t i2, size_t i3) const{

      double d1 = (i1 / _lutFactor) + _lutCutoff;
      double d2 = (i2 / _lutFactor) + _lutCutoff;
      double d3 = (i3 / _lutFactor) + _lutCutoff;


      return  std::array<double, 3>{d1,d2,d3};

    }

};
}