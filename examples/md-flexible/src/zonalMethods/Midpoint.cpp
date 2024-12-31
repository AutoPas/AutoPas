

#include "src/zonalMethods/Midpoint.h"

#include "autopas/utils/ArrayMath.h"
#include "src/ParticleCommunicator.h"

Midpoint::Midpoint(double cutoff, double verletSkinWidth, int ownRank, RectRegion homeBoxRegion,
                   RectRegion globalBoxRegion, autopas::AutoPas_MPI_Comm comm, std::array<int, 26> allNeighbourIndices,
                   std::array<options::BoundaryTypeOption, 3> boundaryType)
    : ZonalMethod(26, ownRank, homeBoxRegion, globalBoxRegion, comm, allNeighbourIndices, boundaryType) {
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
  /**
   * See:
   * Evaluation of Zonal Methods for Small
   * Molecular Systems
   * Julian Spahl
   * p.10
   */
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
        if (convRelNeighboursToIndex({d[0], d[1], d[2]}) < 13) {
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
         * interact with opposite ring
         * Fig 3.5
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
          _interactionSchedule.at(zone).insert(_interactionSchedule.at(zone).end(), oppositeRing.begin(),
                                               oppositeRing.end());
        }
        /* neighbour type edge
         * interact with opposite wing
         * Fig 3.6
         */
        else if (zeroIndices.size() == 1) {
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


const std::vector<RectRegion> Midpoint::getExportRegions() { return _exportRegions; }
const std::vector<RectRegion> Midpoint::getImportRegions() { return _importRegions; }
