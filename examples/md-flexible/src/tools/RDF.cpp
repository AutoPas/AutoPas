/**
 * @file RDF.cpp
 * @author D. Martin
 * @date 07.03.2025
 */
#include "RDF.h"

RDF::RDF(const std::shared_ptr<autopas::AutoPas<ParticleType>> autoPasContainer, double radiusMin, double radiusMax,
         size_t numBins, double guardArea, bool periodicBoundaries)
    : _autoPasContainer{autoPasContainer},
      _radiusMin{radiusMin},
      _radiusMax{radiusMax},
      _numBins{numBins},
      _guardArea{guardArea},
      _periodicBoundaries{periodicBoundaries} {
  _finalRdf.resize(_numBins);
  for (size_t i = 0; i < _numBins; ++i) {
    _finalRdf[i] =
        std::pair<double, double>{static_cast<double>(i + 1) / static_cast<double>(_numBins) * _radiusMax, 0};
  }
}

RDF::RDF(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << " for reading." << std::endl;
    return;
  }

  _finalRdf.clear();
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
      _finalRdf.emplace_back(r, g);
    } catch (const std::exception &e) {
      std::cerr << "Error: Could not convert values in line: " << line << "\n";
    }
  }
  file.close();

  _rdfFinished = true;
}

void RDF::reset() {
  _rdfFinished = false;
  _finalRdf.clear();
  _finalRdf.resize(_numBins);
  _numRDFs = 0;
  for (size_t i = 0; i < _numBins; ++i) {
    _finalRdf[i] =
        std::pair<double, double>{static_cast<double>(i + 1) / static_cast<double>(_numBins) * _radiusMax, 0};
  }
}

void RDF::captureRDF() {
  using namespace autopas::utils::ArrayMath::literals;

  if (_rdfFinished) {
    throw std::runtime_error("Final RDF is already computed!");
  }

  // Parameters for RDF calculation
  const double binWidth = (_radiusMax - _radiusMin) / _numBins;

  const auto &boxMin = _autoPasContainer->getBoxMin();
  const auto &boxMax = _autoPasContainer->getBoxMax();

  std::vector<int> rdfBins(_numBins, 0);
  std::vector<double> rdf(_numBins, 0);
  std::vector<double> shellVolumes(_numBins, 0);
  size_t totalNumParticlesInInnerArea = 0;

  auto calcVolumeInnerArea = [&]() {
    const std::array<double, 3> boxMinNew = {boxMin[0] + _guardArea, boxMin[1] + _guardArea, boxMin[2] + _guardArea};
    const std::array<double, 3> boxMaxNew = {boxMax[0] - _guardArea, boxMax[1] - _guardArea, boxMax[2] - _guardArea};
    const auto lenX = std::abs(boxMaxNew[0] - boxMinNew[0]);
    const auto lenY = std::abs(boxMaxNew[1] - boxMinNew[1]);
    const auto lenZ = std::abs(boxMaxNew[2] - boxMinNew[2]);
    return lenX * lenY * lenZ;
  };

  auto calcShellVolume = [](const double radius, const double thickness) {
    const auto outerSphereVolume =
        (4.0 / 3.0) * M_PI * (radius + thickness) * (radius + thickness) * (radius + thickness);
    const auto innerSphereVolume = (4.0 / 3.0) * M_PI * radius * radius * radius;
    return outerSphereVolume - innerSphereVolume;
  };

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

  const auto volumeInnerArea = calcVolumeInnerArea();

  // calculate the shell volumes
  for (size_t i = 0; i < _numBins; ++i) {
    const double t = _radiusMax / _numBins;
    const double r = i * t;
    shellVolumes[i] = calcShellVolume(r, t);
  }

  int numThreads = autopas::autopas_get_max_threads();
  std::vector<std::vector<int>> rdfBinsThread(numThreads, std::vector<int>(_numBins, 0));

  const double radiusSquared = _radiusMax * _radiusMax;

  // get lower and upper bounds
  const std::array<double, 3UL> offset = {_radiusMax, _radiusMax, _radiusMax};

  AUTOPAS_OPENMP(parallel reduction(+ : totalNumParticlesInInnerArea)) {
    int threadID = autopas::autopas_get_thread_num();
    auto &localBins = rdfBinsThread[threadID];

    for (auto particle1 = _autoPasContainer->begin(autopas::IteratorBehavior::owned); particle1.isValid();
         ++particle1) {
      const auto pos1 = particle1->getR();

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

          double r2 = _periodicBoundaries ? calcMinImageDistSquared(pos1, pos2)
                                          : ((std::pow(pos2[0] - pos1[0], 2) + std::pow(pos2[1] - pos1[1], 2) +
                                              std::pow(pos2[2] - pos1[2], 2)));

          if (r2 < radiusSquared) {
            const auto shellSize = _radiusMax / _numBins;
            const auto bin = static_cast<int>(std::sqrt(r2) / shellSize);
            localBins[bin] += 1;
          }
        }
      }
    }
  }

  // Reduce local histograms into the final rdfBins
  for (const auto &localBins : rdfBinsThread) {
    for (size_t i = 0; i < _numBins; ++i) {
      rdfBins[i] += localBins[i];
    }
  }

  // compute RDF
  for (size_t i = 0; i < _numBins; ++i) {
    _finalRdf[i].second +=
        (static_cast<double>(rdfBins[i]) / static_cast<double>(totalNumParticlesInInnerArea)) /
        (static_cast<double>(totalNumParticlesInInnerArea - 1) * (shellVolumes[i] / volumeInnerArea));
  }

  _numRDFs++;
}

void RDF::computeFinalRDF() {
  for (size_t i = 0; i < _numBins; ++i) {
    _finalRdf[i].second /= static_cast<double>(_numRDFs);
  }

  _rdfFinished = true;
}

std::vector<std::pair<double, double>> &RDF::getFinalRDF() {
  if (_rdfFinished) {
    return _finalRdf;
  } else {
    throw std::runtime_error("Final RDF is not computed yet! Call computeFinalRDF() before getFinalRDF()");
  }
}

void RDF::writeToCSV(std::string outputFolder, std::string filename) {
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
  for (const auto &[distance, value] : _finalRdf) {
    file << distance << "," << value << "\n";
  }

  file.close();
}