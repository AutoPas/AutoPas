

#include "src/zonalMethods/Midpoint.h"

#include "autopas/utils/ArrayMath.h"
#include "src/ParticleCommunicator.h"

Midpoint::Midpoint(double cutoff, double verletSkinWidth, int ownRank, RectRegion homeBoxRegion,
                   RectRegion globalBoxRegion, bool useNewton3, bool pairwiseInteraction,
                   autopas::AutoPas_MPI_Comm comm, std::array<int, 26> allNeighbourIndices,
                   std::array<options::BoundaryTypeOption, 3> boundaryType)
    : ZonalMethod(26, ownRank, homeBoxRegion, globalBoxRegion, comm, allNeighbourIndices, boundaryType),
      _cutoff(cutoff),
      _useNewton3(useNewton3),
      _pairwiseInteraction(pairwiseInteraction) {
  _exportRegions.reserve(_regionCount);
  _importRegions.reserve(_regionCount);

  auto hsCondition = [](const int d[3]) {
    /**
     * Stencil:
     *  from every side
     */
    return true;
  };

  auto identifyZone = [this](const int d[3]) {
    return std::to_string(convRelNeighboursToIndex(std::array<int, 3>{d[0], d[1], d[2]}));
  };

  // calculate exportRegions
  getRectRegionsConditional(_homeBoxRegion, cutoff / 2, verletSkinWidth, _exportRegions, hsCondition, identifyZone,
                            false);

  // calculate importRegions
  getRectRegionsConditional(_homeBoxRegion, cutoff / 2, verletSkinWidth, _importRegions, hsCondition, identifyZone,
                            true);

  std::reverse(_importRegions.begin(), _importRegions.end());

  // calculate interaction schedules
  calculateInteractionSchedule(identifyZone);
}

Midpoint::~Midpoint() = default;

void Midpoint::collectParticles(AutoPasType &autoPasContainer) {
  size_t index = 0;
  for (auto &region : _exportRegions) {
    // NOTE Optimization: Could reserve buffer in advance
    _regionBuffers.at(index).clear();
    if (needToCollectParticles(region.getNeighbour())) {
      region.collectParticles(autoPasContainer, _regionBuffers.at(index));
      wrapAroundPeriodicBoundary(region.getNeighbour(), _regionBuffers.at(index));
    }
    ++index;
  }
}

void Midpoint::SendAndReceiveExports(AutoPasType &autoPasContainer) {
  ParticleCommunicator particleCommunicator(_comm);
  size_t bufferIndex = 0;
  for (auto &exRegion : _exportRegions) {
    auto index = convRelNeighboursToIndex(exRegion.getNeighbour());
    auto neighbourRank = _allNeighbourIndices.at(index);
    if (neighbourRank != _ownRank) {
      particleCommunicator.sendParticles(_regionBuffers.at(bufferIndex), neighbourRank);
    }
    ++bufferIndex;
  }
  // receive
  // NOTE Optimization: Could reserve buffer in advance
  bufferIndex = 0;
  for (auto &imRegion : _importRegions) {
    _importBuffers.at(bufferIndex).clear();
    auto index = convRelNeighboursToIndex(imRegion.getNeighbour());
    auto neighbourRank = _allNeighbourIndices.at(index);
    if (neighbourRank != _ownRank) {
      particleCommunicator.receiveParticles(_importBuffers.at(bufferIndex), neighbourRank);
    } else {
      _importBuffers.at(bufferIndex)
          .insert(_importBuffers.at(bufferIndex).end(), _regionBuffers.at(bufferIndex).begin(),
                  _regionBuffers.at(bufferIndex).end());
    }
    autoPasContainer.addHaloParticles(_importBuffers.at(bufferIndex));
    ++bufferIndex;
  }
  particleCommunicator.waitForSendRequests();
}

void Midpoint::SendAndReceiveResults(AutoPasType &autoPasContainer) {
  ParticleCommunicator particleCommunicator(_comm);
  // send results
  size_t bufferIndex = 0;
  for (auto &imRegion : _importRegions) {
    auto index = convRelNeighboursToIndex(imRegion.getNeighbour());
    auto neighbourRank = _allNeighbourIndices.at(index);
    if (neighbourRank != _ownRank) {
      particleCommunicator.sendParticles(_importBuffers.at(bufferIndex), neighbourRank);
    } else {
      _regionBuffers.at(bufferIndex).clear();
      // NOTE: We can only add the results inside the container if
      // we do not have sent exported in to the home box from both directions
      // <- which is guaranteed not the case for Midpoint
      _regionBuffers.at(bufferIndex)
          .insert(_regionBuffers.at(bufferIndex).end(), _importBuffers.at(bufferIndex).begin(),
                  _importBuffers.at(bufferIndex).end());
    }
    ++bufferIndex;
  }

  // receive reseults
  bufferIndex = 0;
  for (auto &exRegion : _exportRegions) {
    auto index = convRelNeighboursToIndex(exRegion.getNeighbour());
    auto neighbourRank = _allNeighbourIndices.at(index);
    if (neighbourRank != _ownRank) {
      auto size = _regionBuffers.at(bufferIndex).size();
      _regionBuffers.at(bufferIndex).clear();
      _regionBuffers.at(bufferIndex).reserve(size);
      particleCommunicator.receiveParticles(_regionBuffers.at(bufferIndex), neighbourRank);
    }
    ++bufferIndex;
  }
  particleCommunicator.waitForSendRequests();

  // save results to container
  using namespace autopas::utils::ArrayMath::literals;
  bufferIndex = 0;
  // for all exported regions
  for (auto &exRegion : _exportRegions) {
    // go over all exported particles in the container
    for (auto particleIter = autoPasContainer.getRegionIterator(exRegion._origin, exRegion._origin + exRegion._size,
                                                                autopas::IteratorBehavior::owned);
         particleIter.isValid(); ++particleIter) {
      // find the corresponding result in the buffer
      size_t result_index = 0;
      for (auto &result : _regionBuffers.at(bufferIndex)) {
        if (particleIter->getID() == result.getID()) {
          // if found, add the result and delete from buffer
          particleIter->addF(result.getF());
          _regionBuffers.at(bufferIndex).erase(_regionBuffers.at(bufferIndex).begin() + result_index);
          break;
        }
        ++result_index;
      }
    }
    // sanity check
    if (_regionBuffers.at(bufferIndex).size()) {
      throw std::runtime_error("Midpoint: Not all results were found in the container - Something went wrong!");
    }
    ++bufferIndex;
  }
}

void Midpoint::recollectResultsFromContainer(AutoPasType &autoPasContainer) {
  // clear and reserve space
  for (auto &buffer : _importBuffers) {
    auto size = buffer.size();
    buffer.clear();
    buffer.reserve(size);
  }
  // iterate over halo particles and insert into respecitve buffer
  for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::halo); iter.isValid(); ++iter) {
    size_t bufferIndex = 0;
    for (auto &imRegion : _importRegions) {
      if (imRegion.contains(iter->getR())) {
        _importBuffers.at(bufferIndex).push_back(*iter);
        break;
      }
      ++bufferIndex;
    }
  }
}

void Midpoint::calculateZonalInteractionPairwise(std::string zone1, std::string zone2,
                                                 std::function<void(ParticleType &, ParticleType &)> aosFunctor) {
  auto index1 = (_regionCount - std::stoi(zone1)) - 1;
  auto index2 = (_regionCount - std::stoi(zone2)) - 1;
  // calculate forces between the collected results
  for (auto &p1 : _importBuffers.at(index1)) {
    for (auto &p2 : _importBuffers.at(index2)) {
      // check midpoint
      using namespace autopas::utils::ArrayMath::literals;
      auto midpoint = (p1.getR() + p2.getR()) * 0.5;
      if (!_homeBoxRegion.contains(midpoint)) {
        continue;
      }
      aosFunctor(p1, p2);
    }
  }
}

void Midpoint::calculateInteractionSchedule(std::function<std::string(const int[3])> identifyZone) {
  // calculate interaction schedule for triwise
  if (!_pairwiseInteraction) {
    int d[3];
    for (d[0] = -1; d[0] <= 1; d[0]++) {
      for (d[1] = -1; d[1] <= 1; d[1]++) {
        for (d[2] = -1; d[2] <= 1; d[2]++) {
          if (d[0] == 0 && d[1] == 0 && d[2] == 0) {
            continue;
          }
          auto zone = identifyZone(d);
          _interactionZones.push_back(zone);
          for (int i = 0; i < 26; i++) {
            if (std::stoi(zone) == i) continue;
            auto d2 = convZoneStringIntoRelNeighbourIndex(std::to_string(i));
            bool isNeighbour =
                (std::abs(d[0] - d2[0]) <= 1) && (std::abs(d[1] - d2[1]) <= 1) && (std::abs(d[2] - d2[2]) <= 1);
            if (!isNeighbour) continue;
            _interactionSchedule[zone].push_back(std::to_string(i));
          }
        }
      }
    }
    return;
  }

  /**
   * Calculate interaction schedule for pairwise
   * See:
   * Evaluation of Zonal Methods for Small
   * Molecular Systems
   * Julian Spahl
   * p.10
   * NOTE: Actually, there are some interactions missing in this paper, which are
   * the interactions between two center neighbours
   * NOTE: Almost all of these interactions can be skipped, if the dimensions of the
   * interaction boxes are greater equal than the cutoff
   */
  bool optimize = true;
  for (size_t i = 0; i < 3; i++) {
    if (_homeBoxRegion._size[i] < _cutoff) {
      optimize = false;
      break;
    }
  }
  std::cout << std::string("Midpoint interaction schedule optimization: ") +
                   std::string(optimize ? "enabled" : "disabled")
            << std::endl;
  int d[3];
  for (d[0] = -1; d[0] <= 1; d[0]++) {
    for (d[1] = -1; d[1] <= 1; d[1]++) {
      for (d[2] = -1; d[2] <= 1; d[2]++) {
        if (d[0] == 0 && d[1] == 0 && d[2] == 0) {
          continue;
        }
        auto zone = identifyZone(d);
        // add the zone to the zones vector
        _interactionZones.push_back(zone);
        _interactionSchedule.insert_or_assign(zone, std::vector<std::string>{});
        /* if neighbourIndex is smaller than 13, interact with opposite box
         * Fig 3.3
         */
        if (!optimize && convRelNeighboursToIndex({d[0], d[1], d[2]}) < 13) {
          _interactionSchedule[zone].push_back(identifyZone(new int[3]{-d[0], -d[1], -d[2]}));
        }
        // distinguish neighbour types
        std::vector<size_t> zeroIndices;
        std::vector<size_t> nonZeroIndices;
        for (size_t i = 0; i < 3; i++) {
          if (d[i] == 0) {
            zeroIndices.push_back(i);
          } else {
            nonZeroIndices.push_back(i);
          }
        }
        /* neighbour type center:
         * interact with opposite ring (Fig 3.5) +
         * with the "next" center neighbour for two directions +
         */
        if (zeroIndices.size() == 2) {
          std::vector<std::string> oppositeRing;
          oppositeRing.reserve(8);
          int opp[3];
          opp[nonZeroIndices[0]] = -d[nonZeroIndices[0]];
          for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
              if (i == 0 && j == 0) {
                continue;
              }
              opp[zeroIndices[0]] = i;
              opp[zeroIndices[1]] = j;
              oppositeRing.push_back(identifyZone(opp));
            }
          }
          if (!optimize) {
            _interactionSchedule.at(zone).insert(_interactionSchedule.at(zone).end(), oppositeRing.begin(),
                                                 oppositeRing.end());
          }
          // for two directions, interact with the next center piece
          for (size_t i = 0; i < 2; i++) {
            int opp[3] = {d[0], d[1], d[2]};
            if (nonZeroIndices[0] > zeroIndices[i]) {
              if (opp[nonZeroIndices[0]] == 1) {
                opp[zeroIndices[i]] = 1;
              } else {
                opp[zeroIndices[i]] = -1;
              }
            } else {
              if (opp[nonZeroIndices[0]] == 1) {
                opp[zeroIndices[i]] = -1;
              } else {
                opp[zeroIndices[i]] = 1;
              }
            }
            opp[nonZeroIndices[0]] = 0;
            _interactionSchedule.at(zone).push_back(identifyZone(opp));
          }
        }
        /* neighbour type edge
         * interact with opposite wing
         * Fig 3.6
         */
        else if (!optimize && zeroIndices.size() == 1) {
          std::vector<std::string> oppositeWing;
          oppositeWing.reserve(2);
          int opp[3];
          opp[nonZeroIndices[0]] = -d[nonZeroIndices[0]];
          opp[nonZeroIndices[1]] = -d[nonZeroIndices[1]];
          opp[zeroIndices[0]] = -1;
          _interactionSchedule.at(zone).push_back(identifyZone(opp));
          opp[zeroIndices[0]] = 1;
          _interactionSchedule.at(zone).push_back(identifyZone(opp));
        }
        // neighbour type corner
        else {
          // do nothing
        }
      }
    }
  }
}

void Midpoint::calculateZonalInteractionTriwise(
    std::string zone, std::function<void(ParticleType &, ParticleType &, ParticleType &)> aosFunctor) {
  CombinedBuffer combinedBuffer;
  auto zoneIndex = (_regionCount - std::stoi(zone)) - 1;
  auto &zoneBuffer = _importBuffers.at(zoneIndex);

  // initiate combined buffer
  for (auto otherZone : _interactionSchedule.at(zone)) {
    auto otherIndex = (_regionCount - std::stoi(otherZone)) - 1;
    combinedBuffer.add_buffer(_importBuffers.at(otherIndex));
  }

  // calculate interaction of triples
  // for triples with 1 particle in original zone
  for (auto &p1 : zoneBuffer) {
    for (size_t i = 0; i < combinedBuffer.size(); i++) {
      for (size_t j = i + 1; j < combinedBuffer.size(); j++) {
        ParticleType &p2 = combinedBuffer.at(i);
        ParticleType &p3 = combinedBuffer.at(j);
        using namespace autopas::utils::ArrayMath::literals;
        if (!isMidpointInsideDomain(p1.getR(), p2.getR(), p3.getR(), _homeBoxRegion._origin,
                                    _homeBoxRegion._origin + _homeBoxRegion._size)) {
          continue;
        }
        aosFunctor(p1, p2, p3);
      }
    }
  }
  // for triples with 2 particles in original zone
  for (size_t i = 0; i < zoneBuffer.size(); i++) {
    for (size_t j = i + 1; j < zoneBuffer.size(); j++) {
      ParticleType &p1 = zoneBuffer.at(i);
      ParticleType &p2 = zoneBuffer.at(j);
      for (size_t k = 0; k < combinedBuffer.size(); k++) {
        ParticleType &p3 = combinedBuffer.at(k);
        using namespace autopas::utils::ArrayMath::literals;
        if (!isMidpointInsideDomain(p1.getR(), p2.getR(), p3.getR(), _homeBoxRegion._origin,
                                    _homeBoxRegion._origin + _homeBoxRegion._size)) {
          continue;
        }
        aosFunctor(p1, p2, p3);
        aosFunctor(p2, p1, p3);
      }
    }
  }
}

std::array<int, 3> Midpoint::convZoneStringIntoRelNeighbourIndex(std::string s) {
  std::array<int, 3> relNeighbour;
  auto index = std::stoi(s);
  if (index >= 13) index++;
  relNeighbour[2] = index % 3;
  index -= relNeighbour[2];
  relNeighbour[1] = (index % 9) / 3;
  index -= relNeighbour[1] * 3;
  relNeighbour[0] = index / 9;

  using namespace autopas::utils::ArrayMath::literals;
  relNeighbour -= 1;

  return relNeighbour;
}

const std::vector<RectRegion> Midpoint::getExportRegions() { return _exportRegions; }
const std::vector<RectRegion> Midpoint::getImportRegions() { return _importRegions; }
