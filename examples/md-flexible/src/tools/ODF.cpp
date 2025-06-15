/**
 * @file ODF.cpp
 * @author D. Martin
 * @date 12.06.2025
 */
#include "ODF.h"

ODF::ODF(const std::shared_ptr<autopas::AutoPas<ParticleType>> autoPasContainer, double radiusMin, double radiusMax,
         size_t numBins, double guardArea, bool periodicBoundaries)
    : _autoPasContainer{autoPasContainer},
      _radiusMin{radiusMin},
      _radiusMax{radiusMax},
      _numBins{numBins},
      _guardArea{guardArea},
      _periodicBoundaries{periodicBoundaries} {
  _finalOdf.resize(_numBins);
  for (size_t i = 0; i < _numBins; ++i) {
    _finalOdf[i] =
        std::pair<double, double>{static_cast<double>(i + 1) / static_cast<double>(_numBins) * _radiusMax, 0};
  }
}

ODF::ODF(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << " for reading." << std::endl;
    return;
  }

  _finalOdf.clear();
  std::string line;
  std::getline(file, line);  // Skip header

  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string rStr, gStr;

    if (!std::getline(ss, rStr, ',') || !std::getline(ss, gStr, ',')) {
      std::cerr << "Warning: Malformed line in CSV: " << line << std::endl;
      continue;
    }

    try {
      double r = std::stod(rStr);
      double g = std::stod(gStr);
      _finalOdf.emplace_back(r, g);
    } catch (const std::exception &e) {
      std::cerr << "Error: Could not convert values in line: " << line << "\n";
    }
  }
  file.close();

  _odfFinished = true;
}

void ODF::reset() {
  _odfFinished = false;
  _finalOdf.clear();
  _finalOdf.resize(_numBins);
  _numODFs = 0;
  for (size_t i = 0; i < _numBins; ++i) {
    _finalOdf[i] =
        std::pair<double, double>{static_cast<double>(i + 1) / static_cast<double>(_numBins) * _radiusMax, 0};
  }
}

void ODF::captureODF() {
  using namespace autopas::utils::ArrayMath::literals;

  if (_odfFinished) {
    throw std::runtime_error("Final ODF is already computed!");
  }

  // Parameters for ODF calculation
  const double binWidth = (_radiusMax - _radiusMin) / _numBins;

  const auto &boxMin = _autoPasContainer->getBoxMin();
  const auto &boxMax = _autoPasContainer->getBoxMax();

  std::vector<double> odfBins(_numBins, 0);
  std::vector<double> odfBinsCount(_numBins, 0);
  std::vector<double> odf(_numBins, 0);
  size_t totalNumParticlesInInnerArea = 0;

  auto isInGuardArea = [&](const std::array<double, 3> &pos) {
    for (int dim = 0; dim < 3; dim++) {
      if (pos[dim] < boxMin[dim] + _guardArea || pos[dim] > boxMax[dim] - _guardArea) {
        return true;
      }
    }
    return false;
  };

  auto calcMinImageDistSquared = [&](const std::array<double, 3> &pos1, const std::array<double, 3> &pos2) {
    std::array<double, 3> delta;
    std::array<double, 3> boxLen = {boxMax[0] - boxMin[0], boxMax[1] - boxMin[1], boxMax[2] - boxMin[2]};
    for (int dim = 0; dim < 3; dim++) {
      delta[dim] = pos2[dim] - pos1[dim];
      if (delta[dim] > 0.5 * boxLen[dim]) delta[dim] -= boxLen[dim];
      if (delta[dim] <= -0.5 * boxLen[dim]) delta[dim] += boxLen[dim];
    }
    return delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
  };

  auto getPeriodicOverlapBoxes = [&](const std::array<double, 3> &lower, const std::array<double, 3> &upper) {
    std::vector<std::pair<std::array<double, 3>, std::array<double, 3>>> overlapBoxes;

    const std::array<double, 3> boxLength = boxMax - boxMin;

    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          std::array<double, 3> shiftVec = {dx * boxLength[0], dy * boxLength[1], dz * boxLength[2]};

          std::array<double, 3> shiftedLower = lower + shiftVec;
          std::array<double, 3> shiftedUpper = upper + shiftVec;

          bool overlaps = true;
          for (int d = 0; d < 3; ++d) {
            if (shiftedUpper[d] <= boxMin[d] || shiftedLower[d] >= boxMax[d]) {
              overlaps = false;
              break;
            }
          }

          if (overlaps) {
            std::array<double, 3> clampedLower = autopas::utils::ArrayMath::max(shiftedLower, boxMin);
            std::array<double, 3> clampedUpper = autopas::utils::ArrayMath::min(shiftedUpper, boxMax);
            overlapBoxes.emplace_back(clampedLower, clampedUpper);
          }
        }
      }
    }

    return overlapBoxes;
  };

  int numThreads = autopas::autopas_get_max_threads();
  std::vector<std::vector<double>> odfBinsThread(numThreads, std::vector<double>(_numBins, 0));
  std::vector<std::vector<size_t>> odfBinsCountThread(numThreads, std::vector<size_t>(_numBins, 0));

  const double radiusSquared = _radiusMax * _radiusMax;

  // get lower and upper bounds
  const std::array<double, 3UL> offset = {_radiusMax, _radiusMax, _radiusMax};

  AUTOPAS_OPENMP(parallel reduction(+ : totalNumParticlesInInnerArea)) {
    int threadID = autopas::autopas_get_thread_num();
    auto &localBins = odfBinsThread[threadID];
    auto &localBinsCount = odfBinsCountThread[threadID];

    for (auto particle1 = _autoPasContainer->begin(autopas::IteratorBehavior::owned); particle1.isValid();
         ++particle1) {
      const auto pos1 = particle1->getR();
      const auto q1 = particle1->getQuaternion();

      // Skip particles in the guard area
      if (isInGuardArea(pos1)) continue;
      totalNumParticlesInInnerArea++;

      const auto lower = pos1 - offset;
      const auto upper = pos1 + offset;

      const auto overlapBoxes = getPeriodicOverlapBoxes(lower, upper);

      for (const auto &[lo, up] : overlapBoxes) {
        for (auto particle2 = _autoPasContainer->getRegionIterator(
                 lo, up, autopas::IteratorBehavior::owned | autopas::IteratorBehavior::forceSequential);
             particle2.isValid(); ++particle2) {
          if (*particle1 == *particle2) {
            continue;
          }

          const auto pos2 = particle2->getR();
          const auto q2 = particle2->getQuaternion();

          double r2 = _periodicBoundaries ? calcMinImageDistSquared(pos1, pos2)
                                          : ((std::pow(pos2[0] - pos1[0], 2) + std::pow(pos2[1] - pos1[1], 2) +
                                              std::pow(pos2[2] - pos1[2], 2)));

          if (r2 < radiusSquared) {
            const auto shellSize = _radiusMax / _numBins;
            const auto bin = static_cast<int>(std::sqrt(r2) / shellSize);

            // calculate |u1 (dot) u2|
            const auto u1 = autopas::utils::ArrayMath::normalize(
                autopas::utils::quaternion::rotatePosition(q1, std::array{0.0, 0.0, 1.0}));
            const auto u2 = autopas::utils::ArrayMath::normalize(
                autopas::utils::quaternion::rotatePosition(q2, std::array{0.0, 0.0, 1.0}));
            const auto odfVal = std::abs(autopas::utils::ArrayMath::dot(u1, u2));

            localBins[bin] += odfVal;
            localBinsCount[bin] += 1;
          }
        }
      }
    }
  }

  // Reduce local histograms into the final odfBins
  for (const auto &localBins : odfBinsThread) {
    for (size_t i = 0; i < _numBins; ++i) {
      odfBins[i] += localBins[i];
    }
  }
  for (const auto &localBinsCount : odfBinsCountThread) {
    for (size_t i = 0; i < _numBins; ++i) {
      odfBinsCount[i] += localBinsCount[i];
    }
  }

  // compute ODF
  for (size_t i = 0; i < _numBins; ++i) {
    if (odfBins[i] > 0) {
      _finalOdf[i].second += odfBins[i] / static_cast<double>(odfBinsCount[i]);
    }
  }

  _numODFs++;
}

void ODF::computeFinalODF() {
  for (size_t i = 0; i < _numBins; ++i) {
    _finalOdf[i].second /= static_cast<double>(_numODFs);
  }

  _odfFinished = true;
}

std::vector<std::pair<double, double>> &ODF::getFinalODF() {
  if (_odfFinished) {
    return _finalOdf;
  } else {
    throw std::runtime_error("Final ODF is not computed yet! Call computeFinalODF() before getFinalODF()");
  }
}

void ODF::writeToCSV(std::string outputFolder, std::string filename) {
  // Ensure the output folder exists
  std::filesystem::create_directories(outputFolder);

  // Construct the full file path
  std::string filePath = outputFolder + "/" + filename + ".csv";

  std::ofstream file(filePath);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filePath << " for writing." << std::endl;
    return;
  }

  // Write header
  file << "distance,value\n";

  // Write data
  for (const auto &[distance, value] : _finalOdf) {
    file << distance << "," << value << "\n";
  }

  file.close();
}