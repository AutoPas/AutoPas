/**
 * @file RotationalAnalysis.h
 * @author D. Martin
 * @date 19.05.2025
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "autopas/utils/Quaternion.h"

struct Quaternion {
  double w, x, y, z;

  Quaternion inverse() const { return {w, -x, -y, -z}; }

  Quaternion operator*(const Quaternion &q) const {
    return {w * q.w - x * q.x - y * q.y - z * q.z, w * q.x + x * q.w + y * q.z - z * q.y,
            w * q.y - x * q.z + y * q.w + z * q.x, w * q.z + x * q.y - y * q.x + z * q.w};
  }

  double dot(const Quaternion &q) const { return w * q.w + x * q.x + y * q.y + z * q.z; }
};

class RotationalAnalysis {
 public:
  void setValues(size_t maxLagSteps, size_t stepInterval) {
    _numLags = maxLagSteps / stepInterval + 1;
    _maxLagSteps = maxLagSteps;
    _stepInterval = stepInterval;
    _msdSums.resize(_numLags, 0.0);
    _oacfSums.resize(_numLags, 0.0);
    _counts.resize(_numLags, 0);
  }

  void recordOrientation(size_t particleId, size_t timestep, const Quaternion &q) {
    auto &buffer = _orientationBuffers[particleId];
    buffer.push_back({timestep, q});
    if (buffer.size() > _numLags) buffer.pop_front();

    for (size_t lagIdx = 0; lagIdx < _numLags; ++lagIdx) {
      size_t lag = lagIdx * _stepInterval;
      if (buffer.size() <= lag) continue;

      const auto &[tRef, qRef] = buffer[buffer.size() - 1 - lag];
      Quaternion qRel = qRef.inverse() * q;
      double angle = 2.0 * std::acos(std::clamp(qRel.w, -1.0, 1.0));
      double angleSq = angle * angle;

      double dot = qRef.dot(q);
      double oacfVal = dot * dot;

#pragma omp atomic
      _msdSums[lagIdx] += angleSq;

#pragma omp atomic
      _oacfSums[lagIdx] += oacfVal;

#pragma omp atomic
      _counts[lagIdx] += 1;
    }
  }

  void writeResults(std::string outputFolder, std::string filename) const {
    // Ensure the output folder exists
    std::filesystem::create_directories(outputFolder);

    std::ofstream msdOut(outputFolder + "/" + filename + "_msd_rot.csv");
    std::ofstream oacfOut(outputFolder + "/" + filename + "_oacf.csv");
    msdOut << "# lagTime,msd_rot\n";
    oacfOut << "# lagTime,oacf\n";
    for (size_t i = 0; i < _numLags; ++i) {
      double time = i * _stepInterval;
      size_t count = _counts[i];
      double msdAvg = count > 0 ? _msdSums[i] / count : 0.0;
      double oacfAvg = count > 0 ? _oacfSums[i] / count : 0.0;
      msdOut << time << "," << msdAvg << "\n";
      oacfOut << time << "," << oacfAvg << "\n";
    }
  }

 private:
  size_t _maxLagSteps;
  size_t _stepInterval;
  size_t _numLags;

  std::unordered_map<size_t, std::deque<std::pair<size_t, Quaternion>>> _orientationBuffers;

  std::vector<double> _msdSums;
  std::vector<double> _oacfSums;
  std::vector<size_t> _counts;
};
