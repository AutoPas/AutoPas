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

class TriwiseLUT {
  public:
  TriwiseLUT(int resolution) : _resolution(resolution) {}

  template <class Functor>
  [[nodiscard]] std::pair<const std::array<double, 3>, std::array<u_int8_t, 3>> retrieveValues(const Functor &functor, double dist1, double dist2, double dist3) const {

    if (dist1 < _lutCutoff or dist2 < _lutCutoff or dist3 < _lutCutoff) {
      return std::make_pair(functor.getLUTValues(dist1, dist2, dist3), std::array<u_int8_t, 3>({0, 1, 2}));
    }

    size_t index1, index2, index3;
    std::array<u_int8_t,3> order{};  // Array to track the original indices of the distances

    if (dist1 >= dist2 && dist2 >= dist3) {
      index1 = static_cast<size_t>((dist1 - _lutCutoff) / _resolution);
      index2 = static_cast<size_t>((dist2 - _lutCutoff) / _resolution);
      index3 = static_cast<size_t>((dist3 - _lutCutoff) / _resolution);
      order[0] = 0; order[1] = 1; order[2] = 2;
    } else if (dist1 >= dist3 && dist3 >= dist2) {
      index1 = static_cast<size_t>((dist1 - _lutCutoff) / _resolution);
      index2 = static_cast<size_t>((dist3 - _lutCutoff) / _resolution);
      index3 = static_cast<size_t>((dist2 - _lutCutoff) / _resolution);
      order[0] = 0; order[1] = 2; order[2] = 1;
    } else if (dist2 >= dist1 && dist1 >= dist3) {
      index1 = static_cast<size_t>((dist2 - _lutCutoff) / _resolution);
      index2 = static_cast<size_t>((dist1 - _lutCutoff) / _resolution);
      index3 = static_cast<size_t>((dist3 - _lutCutoff) / _resolution);
      order[0] = 1; order[1] = 0; order[2] = 2;
    } else if (dist2 >= dist3 && dist3 >= dist1) {
      index1 = static_cast<size_t>((dist2 - _lutCutoff) / _resolution);
      index2 = static_cast<size_t>((dist3 - _lutCutoff) / _resolution);
      index3 = static_cast<size_t>((dist1 - _lutCutoff) / _resolution);
      order[0] = 1; order[1] = 2; order[2] = 0;
    } else if (dist3 >= dist1 && dist1 >= dist2) {
      index1 = static_cast<size_t>((dist3 - _lutCutoff) / _resolution);
      index2 = static_cast<size_t>((dist1 - _lutCutoff) / _resolution);
      index3 = static_cast<size_t>((dist2 - _lutCutoff) / _resolution);
      order[0] = 2; order[1] = 0; order[2] = 1;
    } else {
      index1 = static_cast<size_t>((dist3 - _lutCutoff) / _resolution);
      index2 = static_cast<size_t>((dist2 - _lutCutoff) / _resolution);
      index3 = static_cast<size_t>((dist1 - _lutCutoff) / _resolution);
      order[0] = 2; order[1] = 1; order[2] = 0;
    }

    // Retrieve values from the LUT
    std::pair<const std::array<double, 3>, std::array<u_int8_t, 3>> result = std::make_pair(_lut[index1][index2][index3], order);

    return result;
  }

  template<class Functor>
  void fill(const Functor &functor, double cutoffSquared, bool forceRefill = false) {
    if (not _isFilled or forceRefill) {
      _cutoffSquared = cutoffSquared;
      _lutCutoff = cutoffSquared / 10.;
      double lutDistance = _cutoffSquared - _lutCutoff;
      double pointDistance = lutDistance / static_cast<double>(_resolution);

      _lut.resize(_resolution + 1);
      for (auto a = 0; a <= _resolution; a++) {
        double distA = _lutCutoff + a * pointDistance;
        _lut[a].resize(a + 1);
        for (auto b = 0; b <= a; b++) {
          double distB = _lutCutoff + b * pointDistance;
          _lut[a][b].reserve(b + 1);
          for (auto c = 0; c <= b; c++) {
            double distC = _lutCutoff + c * pointDistance;
            _lut[a][b].push_back(functor.getLUTValues(distA, distB, distC));
          }
        }
      }
      _isFilled = true;
    }
  }

  private:

  std::vector<std::vector<std::vector<std::array<double, 3>>>> _lut;

  double _cutoffSquared{};
  double _lutCutoff{};
  int _resolution{};
  bool _isFilled = false;
};
}