/**
* @file TriwiseLUT.h
 * @author muehlhaeusser
 * @date 03.09.2024
 */

#pragma once

#include <array>
#include <vector>
#include <cmath>
#include <algorithm>

namespace mdLib {

template<class TriwiseFunctor>
class TriwiseLUT {
  public:
  TriwiseLUT(double cutoffSquared, int resolution, TriwiseFunctor &functor) : _lutCutoff(cutoffSquared / 10.), _resolution(resolution) {

    double lutDistance = cutoffSquared - _lutCutoff;
    double pointDistance = lutDistance / static_cast<double>(_resolution);

    _lut.reserve(_resolution + 1);
    for (auto a = 0; a <= _resolution; a++) {
      double distA = _lutCutoff + a * pointDistance;
      _lut[a].reserve(a + 1);
      for (auto b = 0; b <= a; b++) {
        double distB = _lutCutoff + b * pointDistance;
        _lut[a][b].reserve(b + 1);
        for (auto c = 0; c <= b; c++) {
          double distC = _lutCutoff + c * pointDistance;
          _lut[a][b].push_back(functor.getLUTValues(distA, distB, distC));
        }
      }
    }
  }

  [[nodiscard]] std::array<double, 3> retrieveValues(double dist1, double dist2, double dist3) const {

    if (std::any_of({dist1, dist2, dist3}, [&](auto r) { return r < _lutCutoff; })) {
      return TriwiseFunctor::getLUTValues(dist1, dist2, dist3);
    }

    auto index1 = static_cast<size_t>((dist1 - _lutCutoff) / _resolution);
    auto index2 = static_cast<size_t>((dist2 - _lutCutoff) / _resolution);
    auto index3 = static_cast<size_t>((dist3 - _lutCutoff) / _resolution);

    return _lut[index1][index2][index3];
  }


  private:

  std::vector<std::vector<std::vector<std::array<double, 3>>>> _lut;

  double _lutCutoff{};
  int _resolution{};
};
}