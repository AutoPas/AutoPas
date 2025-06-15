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

  // Returns the inverse of the quaternion (conjugate, assuming unit quaternion)
  Quaternion inverse() const { return {w, -x, -y, -z}; }

  // Quaternion multiplication (Hamilton product)
  Quaternion operator*(const Quaternion &q) const {
    return {w * q.w - x * q.x - y * q.y - z * q.z, w * q.x + x * q.w + y * q.z - z * q.y,
            w * q.y - x * q.z + y * q.w + z * q.x, w * q.z + x * q.y - y * q.x + z * q.w};
  }

  // Dot product between two quaternions
  double dot(const Quaternion &q) const { return w * q.w + x * q.x + y * q.y + z * q.z; }
};

class RotationalAnalysis {
 public:
  /**
   * Set the parameters for the rotational correlation analysis.
   * @param maxLagSteps Maximum number of lag timesteps to consider.
   * @param stepInterval Interval between evaluated lag timesteps.
   */
  void setValues(const std::shared_ptr<autopas::AutoPas<ParticleType>> autoPasContainer, size_t maxLagSteps,
                 size_t stepInterval) {
    _numLags = maxLagSteps / stepInterval;
    _maxLagSteps = maxLagSteps;
    _stepInterval = stepInterval;
    _msdSums.resize(_numLags, 0.0);
    _oacfSums.resize(_numLags, 0.0);
    _counts.resize(_numLags, 0);
    _autoPasContainer = autoPasContainer;
    for (size_t id = 0; id < _autoPasContainer->getNumberOfParticles(); ++id) {
      _orientationBuffers[id] = {};
    }
  }

  /**
   * Record the orientation of all particles at a specific timestep and update
   * rotational MSD and OACF for all lag times.
   *
   * @param timestep Current timestep of the simulation.
   */
  void recordOrientations(size_t timestep) {
    if (!_autoPasContainer) return;

    int numThreads = autopas::autopas_get_max_threads();
    std::vector<std::vector<double>> msdSumsLocal(numThreads, std::vector<double>(_numLags, 0.0));
    std::vector<std::vector<double>> oacfSumsLocal(numThreads, std::vector<double>(_numLags, 0.0));
    std::vector<std::vector<size_t>> countsLocal(numThreads, std::vector<size_t>(_numLags, 0));

    AUTOPAS_OPENMP(parallel) {
      int threadID = autopas::autopas_get_thread_num();

      for (auto iter = _autoPasContainer->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
        const auto particleId = iter->getID();
        const auto [w, x, y, z] = iter->getQuaternion();
        const auto q = Quaternion{w, x, y, z};

        auto &buffer = _orientationBuffers[particleId];
        buffer.push_back({timestep, q});

        while (!buffer.empty() and
               static_cast<int>(buffer.front().first) < (static_cast<int>(timestep) - static_cast<int>(_maxLagSteps))) {
          // std::cout << "pop " << ", buffer.front().first: " << buffer.front().first
          //           << ", timestep - _maxLagSteps: " << (static_cast<int>(timestep) - static_cast<int>(_maxLagSteps))
          //           << std::endl;
          buffer.pop_front();
        }

        for (size_t lagIdx = 0; lagIdx < _numLags; ++lagIdx) {
          size_t lag = lagIdx * _stepInterval;
          if (timestep < lag) {
            // std::cout << "continue in timestep < lag with timestep " << timestep << " and lag " << lag << " and
            // lagIdx "
            //           << lagIdx << std::endl;
            continue;
          }

          size_t targetTime = timestep - lag;

          auto it = std::find_if(buffer.begin(), buffer.end(),
                                 [targetTime](const auto &entry) { return entry.first == targetTime; });

          if (it == buffer.end()) {
            // std::cout << "continue in it == buffer.end() with timestep " << timestep << " and lagIdx " << lagIdx
            //           << " and lag " << lag << std::endl;
            continue;
          }

          const auto &[tRef, qRef] = *it;

          Quaternion qRel = qRef.inverse() * q;

          // Angular displacement (modulo 2π), derived from quaternion w-component
          double angle = 2.0 * std::acos(std::clamp(qRel.w, -1.0, 1.0));
          double angleSq = angle * angle;

          // Orientation auto-correlation function: square of the dot product
          double dot = qRef.dot(q);
          double oacfVal = dot * dot;

          msdSumsLocal[threadID][lagIdx] += angleSq;
          oacfSumsLocal[threadID][lagIdx] += oacfVal;
          countsLocal[threadID][lagIdx] += 1;
        }
      }
    }

    for (size_t threadID = 0; threadID < numThreads; ++threadID) {
      for (size_t lagIdx = 0; lagIdx < _numLags; ++lagIdx) {
        _msdSums[lagIdx] += msdSumsLocal[threadID][lagIdx];
        _oacfSums[lagIdx] += oacfSumsLocal[threadID][lagIdx];
        _counts[lagIdx] += countsLocal[threadID][lagIdx];
      }
    }
  }

  /**
   * Write the time-lagged rotational MSD and OACF results to CSV files.
   * @param outputFolder Folder to write the results to.
   * @param filename Base name for the output files.
   */
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

  // Maps each particle to its orientation history (quaternion, timestep)
  std::unordered_map<size_t, std::deque<std::pair<size_t, Quaternion>>> _orientationBuffers;

  // Rotational mean square displacement (MSD) sums
  std::vector<double> _msdSums;

  // Orientation autocorrelation function (OACF) sums
  std::vector<double> _oacfSums;

  // Number of contributing samples per lag
  std::vector<size_t> _counts;

  /**
   * The the nodes' AutoPas container used for simulation.
   */
  std::shared_ptr<autopas::AutoPas<ParticleType>> _autoPasContainer;
};
