

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

  // calculate interaction schedules
  calculateInteractionSchedule(identifyZone);
}

Midpoint::~Midpoint() = default;

void Midpoint::collectParticles(AutoPasType &autoPasContainer) {
  size_t index = 0;
  for (auto &region : _exportRegions) {
    // NOTE Optimization: Could reserve buffer in advance
    _regionBuffers[index].clear();
    if (needToCollectParticles(region.getNeighbour())) {
      region.collectParticles(autoPasContainer, _regionBuffers[index]);
      wrapAroundPeriodicBoundary(region.getNeighbour(), _regionBuffers[index]);
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
      particleCommunicator.sendParticles(_regionBuffers[bufferIndex], neighbourRank);
    }
    ++bufferIndex;
  }
  // receive
  // NOTE Optimization: Could reserve buffer in advance
  bufferIndex = 0;
  for (auto &imRegion : _importRegions) {
    _importBuffers[bufferIndex].clear();
    auto index = convRelNeighboursToIndex(imRegion.getNeighbour());
    auto neighbourRank = _allNeighbourIndices.at(index);
    if (neighbourRank != _ownRank) {
      particleCommunicator.receiveParticles(_importBuffers[bufferIndex], neighbourRank);
    } else {
      _importBuffers[bufferIndex].insert(_importBuffers[bufferIndex].end(), _regionBuffers[bufferIndex].begin(),
                                         _regionBuffers[bufferIndex].end());
    }
    autoPasContainer.addHaloParticles(_importBuffers[bufferIndex]);
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
      particleCommunicator.sendParticles(_importBuffers[bufferIndex], neighbourRank);
    } else {
      // NOTE: We can only add the results inside the container if
      // we do not have sent exported in to the home box from both directions
      // <- which is guaranteed no the case for Midpoint
      _regionBuffers[bufferIndex].insert(_regionBuffers[bufferIndex].end(), _importBuffers[bufferIndex].begin(),
                                         _importBuffers[bufferIndex].end());
    }
    ++bufferIndex;
  }

  // receive reseults
  bufferIndex = 0;
  for (auto &exRegion : _exportRegions) {
    auto index = convRelNeighboursToIndex(exRegion.getNeighbour());
    auto neighbourRank = _allNeighbourIndices.at(index);
    if (neighbourRank != _ownRank) {
      auto size = _regionBuffers[bufferIndex].size();
      _regionBuffers[bufferIndex].clear();
      _regionBuffers[bufferIndex].reserve(size);
      particleCommunicator.receiveParticles(_regionBuffers[bufferIndex], neighbourRank);
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
      for (auto &result : _regionBuffers[bufferIndex]) {
        if (particleIter->getID() == result.getID()) {
          // if found, add the result and delete from buffer
          particleIter->addF(result.getF());
          _regionBuffers[bufferIndex].erase(_regionBuffers[bufferIndex].begin() + result_index);
          break;
        }
        ++result_index;
      }
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
        _importBuffers[bufferIndex].push_back(*iter);
        break;
      }
      ++bufferIndex;
    }
  }
}

void Midpoint::calculateZonalInteractionPairwise(std::string zone1, std::string zone2,
                                                 std::function<void(ParticleType &, ParticleType &)> aosFunctor) {}

void Midpoint::calculateInteractionSchedule(std::function<std::string(const int[3])> identifyZone) {
  /**
   * See:
   * The midpoint method for parallelization of particle
   * simulations
   * Kevin J. Bowers; Ron O. Dror; David E. Shaw
   * Figure 4
   */
  std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> interactionSchedule2D;
  interactionSchedule2D.insert_or_assign({1, 0},
                                         std::vector<std::pair<int, int>>{{0, -1}, {-1, 0}, {0, 1}, {-1, -1}, {-1, 1}});
  interactionSchedule2D.insert_or_assign({0, -1}, std::vector<std::pair<int, int>>{{-1, 0}, {0, 1}, {1, 1}, {-1, 1}});
  interactionSchedule2D.insert_or_assign({-1, 0}, std::vector<std::pair<int, int>>{{0, 1}, {1, 1}, {-1, 1}});
  interactionSchedule2D.insert_or_assign({0, 1}, std::vector<std::pair<int, int>>{{-1, 1}, {-1, -1}});
  interactionSchedule2D.insert_or_assign({1, 1}, std::vector<std::pair<int, int>>{{-1, -1}});
  interactionSchedule2D.insert_or_assign({-1, 1}, std::vector<std::pair<int, int>>{{-1, 1}});

  auto get2DCoord = [](int d[3], size_t dim) -> std::pair<int, int> {
    if (dim == 0) {
      return {d[1], d[2]};
    } else if (dim == 1) {
      return {d[0], d[2]};
    }
    return {d[0], d[1]};
  };

  auto restore3DCoord = [](std::pair<int, int> pair, size_t dim, int coord) -> std::array<int, 3> {
    if (dim == 0) {
      return {coord, pair.first, pair.second};
    } else if (dim == 1) {
      return {pair.first, coord, pair.second};
    }
    return {pair.first, pair.second, coord};
  };

  // reserve in advance
  _interactionZones.reserve(_zoneCount);

  // for each relative neighbour
  int d[3];
  for (d[0] = -1; d[0] <= 1; d[0]++) {
    for (d[1] = -1; d[1] <= 1; d[1]++) {
      for (d[2] = -1; d[2] <= 1; d[2]++) {
        auto zone = identifyZone(d);
        // add the zone to the zones vecotr
        _interactionZones.push_back(zone);
        // apply the 2D interaction into each direction
        for (size_t dim = 0; dim < 3; dim++) {
          // if 
          if (d[dim] == 0) {
            auto coord = get2DCoord(d, dim);
            auto interactions = interactionSchedule2D.at(coord);
            // convert 2D interactions to respective 3D neighbours
            std::vector<std::string> interactions3D;
            for (auto inter : interactions) {
              interactions3D.push_back(identifyZone(restore3DCoord(inter, dim, d[dim]).data()));
            }
            _interactionSchedule.at(zone).insert(_interactionSchedule.at(zone).end(), interactions3D.begin(),
                                                 interactions3D.end());
          }
        }
      }
    }
  }
}
