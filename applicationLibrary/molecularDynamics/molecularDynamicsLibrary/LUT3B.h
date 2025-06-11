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

//       std::cout << "RETRIEVE VALUE ";

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
     void fill(const Functor &functor, double cutoffSquared, bool withGlobals) {

if(!withGlobals) {
  fill_plain(functor, cutoffSquared);
}if(withGlobals){
//  fill_gobal(functor, cutoffSquared);
 fill_gobal_all(functor, cutoffSquared);

}

       std::cout << "filled ";
    }

    std::array<double, 3> getNextNeighbour( size_t index1, size_t index2, size_t index3) const {

      //just a function for consistency, floored values are basically next neighbours
      return _lut3B[index1][index2][index3];

    }





//alternative global LUT

    std::array<double, 2> calculateGlobalFac(  double distSquaredIJ, double distSquaredJK, double distSquaredKI) const {


      const double allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
      const double allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
//      const double factor = 3.0 * nu / allDistsTo5;
      const double factor =   3.0 * 1.6152500E-3  / allDistsTo5;

      return{allDistsSquared, factor};

    }

    std::array<double, 2> calculateFullGlobal(  double distSquaredIJ, double distSquaredJK, double distSquaredKI) const {


      const double allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
      const double allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
      //      const double factor = 3.0 * nu / allDistsTo5;
      const double factor =   3.0 * 1.6152500E-3  / allDistsTo5;

      //        // Dot products of both distance vectors going from one particle  -> over law of cosine
      const double IJDotKI = 0.5 * (distSquaredIJ + distSquaredKI - distSquaredJK);
      const double IJDotJK = 0.5 * (distSquaredIJ + distSquaredJK - distSquaredKI);
      const double JKDotKI = 0.5 * (distSquaredJK + distSquaredKI - distSquaredIJ);
      const double allDotProducts = IJDotKI * IJDotJK * JKDotKI;
      const double potentialEnergy3 = factor * (allDistsSquared - 3.0 * allDotProducts);

      return{potentialEnergy3, 0};

    }


    std::array<double, 2> calculatePotEnergyTest( double distSquaredIJ, double distSquaredJK, double distSquaredKI) const {


      double distIJ = std::sqrt(distSquaredIJ);
      double distJK = std::sqrt(distSquaredJK);
      double distKI = std::sqrt(distSquaredKI);

      double distIJ3 = distSquaredIJ * distIJ;
      double distJK3 = distSquaredJK * distJK;
      double distKI3 = distSquaredKI * distKI;
      double KIcosIJ = (distSquaredIJ + distSquaredKI - distSquaredJK) / (2  * distIJ * distKI);
      double IJcosJK= (distSquaredIJ + distSquaredJK - distSquaredKI) / (2  * distIJ * distJK);
      double JKcosKI = (distSquaredJK + distSquaredKI - distSquaredIJ) / (2  * distJK * distKI);
      double res = 3* (1.6152500E-3  * (1+ (3 * KIcosIJ * IJcosJK * JKcosKI) ) )/(distIJ3* distKI3 *distJK3 );

      return {res, 0};
    }


    template<class Functor>
    void fill_plain(const Functor &functor, double cutoffSquared) {
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
    }

    template<class Functor>
    void fill_gobal(const Functor &functor, double cutoffSquared) {


      _lut3B_global.resize(_resolution + 1);
      for (auto a = 0; a <= _resolution; a++) {
        double distA = _lutCutoff + a * _pointDistance;
        _lut3B_global[a].resize(a + 1);
        for (auto b = 0; b <= a; b++) {
          double distB = _lutCutoff + b * _pointDistance;
          _lut3B_global[a][b].reserve(b + 1);
          for (auto c = 0; c <= b; c++) {
            double distC = _lutCutoff + c * _pointDistance;
            auto global = calculateGlobalFac(distA, distB, distC);
            auto lutVal = functor.getLUTValues(distA, distB, distC);
            std::array<double, 5> res = {lutVal[0],lutVal[1],lutVal[2], global[0], global[1]  };
            _lut3B_global[a][b].push_back(res);
          }
        }
      }

      std::cout << "filled ";
    }

    template<class Functor>
    void fill_gobal_all(const Functor &functor, double cutoffSquared) {


      _lut3B_global.resize(_resolution + 1);
      for (auto a = 0; a <= _resolution; a++) {
        double distA = _lutCutoff + a * _pointDistance;
        _lut3B_global[a].resize(a + 1);
        for (auto b = 0; b <= a; b++) {
          double distB = _lutCutoff + b * _pointDistance;
          _lut3B_global[a][b].reserve(b + 1);
          for (auto c = 0; c <= b; c++) {
            double distC = _lutCutoff + c * _pointDistance;
//            auto global = calculateFullGlobal(distA, distB, distC);
            //testing
            auto global = calculatePotEnergyTest(distA, distB, distC);
            auto lutVal = functor.getLUTValues(distA, distB, distC);
            std::array<double, 5> res = {lutVal[0],lutVal[1],lutVal[2], global[0], 0  };
            _lut3B_global[a][b].push_back(res);
          }
        }
      }

      std::cout << "filled ";
    }


    template <class Functor>
//    [[nodiscard]] std::pair<const std::array<double, 3>, std::pair<std::array<u_int8_t, 3>, std::array<double,2>> retrieveValues_global(const Functor &functor, double dist1, double dist2, double dist3) const {
    [[nodiscard]] std::pair<const std::array<double, 5>, std::array<u_int8_t, 3>> retrieveValues_global(const Functor &functor, double dist1, double dist2, double dist3) const {

    //    using namespace autopas::utils::ArrayMath::literals;
      if (dist1 < _lutCutoff or dist2 < _lutCutoff or dist3 < _lutCutoff) {
          auto globs = calculateGlobalFac(dist1, dist2, dist3);
          auto vals = functor.getLUTValues(dist1, dist2, dist3);
          std:std::array<double, 5> res = {vals[0],vals[1],vals[2],globs[0], globs[1] };
        return std::make_pair(res, std::array<u_int8_t, 3>({0, 1, 2}));
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
      auto result = getNextNeighbour_global(index1, index2, index3);
//      auto val ={result[0], result[1], result[2]};
//      auto fac = {result[3], result[4]};
//
//      auto pait = std::make_pair(order, fac);
//      auto resultPair = std::make_pair(val, pait);
      auto resultPair = std::make_pair(result, order);

//      std::cout << "RETRIEVE VALUE ";

      return resultPair;
    }


    std::array<double, 5> getNextNeighbour_global( size_t index1, size_t index2, size_t index3) const {

      //just a function for consistency, floored values are basically next neighbours
      return _lut3B_global[index1][index2][index3];

    }

    template<class Functor>
    void setNu(double nu ){
      this->nu = nu;
    }

   private:
    std::vector<std::vector<std::vector<std::array<double, 3>>>> _lut3B;
    std::vector<std::vector<std::vector<std::array<double, 5>>>> _lut3B_global;
    double nu;



    //helper functions

    std::array<double, 3> getDistances(size_t i1,size_t i2, size_t i3) const{

      double d1 = (i1 / _lutFactor) + _lutCutoff;
      double d2 = (i2 / _lutFactor) + _lutCutoff;
      double d3 = (i3 / _lutFactor) + _lutCutoff;


      return  std::array<double, 3>{d1,d2,d3};

    }

};
}