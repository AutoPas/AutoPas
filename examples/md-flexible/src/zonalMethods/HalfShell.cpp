
#include "src/zonalMethods/HalfShell.h"

#include "autopas/utils/ArrayMath.h"
#include "src/ParticleCommunicator.h"

HalfShell::HalfShell(double cutoff, double verletSkinWidth, int ownRank, RectRegion homeBoxRegion,
                     RectRegion globalBoxRegion, bool useNewton3, bool pairwiseInteraction,
                     autopas::AutoPas_MPI_Comm comm, std::array<int, 26> allNeighbourIndices,
                     std::array<options::BoundaryTypeOption, 3> boundaryType)
    : ZonalMethod(1, ownRank, homeBoxRegion, globalBoxRegion, comm, allNeighbourIndices, boundaryType),
      _useNewton3(useNewton3),
      _pairwiseInteraction(pairwiseInteraction),
      _cutoff(cutoff),
      _verletSkinWidth(verletSkinWidth) {
  _exportRegions.reserve(_regionCount);
  _importRegions.reserve(_regionCount);

  auto hsCondition = [](const int d[3]) {
    /**
     * Stencil:
     *  z > 0 +
     *  z == 0 and y > 0 +
     *  z == 0 and y == 0 and x > 0
     */
    return d[2] > 0 or (d[2] == 0 and (d[1] > 0 or (d[1] == 0 and d[0] > 0)));
  };

  auto importCondition = [hsCondition](const int d[3]) { return !hsCondition(d); };

  auto identifyZone = [](const int d[3]) { return "A"; };

  // calculate exportRegions
  getRectRegionsConditional(_homeBoxRegion, cutoff, verletSkinWidth, _exportRegions, hsCondition, identifyZone, false);

  // calculate importRegions
  getRectRegionsConditional(_homeBoxRegion, cutoff, verletSkinWidth, _importRegions, importCondition, identifyZone,
                            true);

  std::reverse(_importRegions.begin(), _importRegions.end());

  _interactionZones.push_back("A");
  _interactionSchedule.insert_or_assign("A", std::vector<std::string>{});
}

HalfShell::~HalfShell() = default;

void HalfShell::collectParticles(AutoPasType &autoPasContainer) {
  _exportParticleMap.clear();
  size_t index = 0;
  for (auto &region : _exportRegions) {
    // NOTE Optimization: Could reserve buffer in advance
    _regionBuffers.at(index).clear();
    if (needToCollectParticles(region.getNeighbour())) {
      region.collectParticles(autoPasContainer, _regionBuffers.at(index), std::ref(_exportParticleMap));
      wrapAroundPeriodicBoundary(region.getNeighbour(), _regionBuffers.at(index));
    }
    ++index;
  }
}

void HalfShell::SendAndReceiveImports(AutoPasType &autoPasContainer) {
  ParticleCommunicator particleCommunicator(_comm);
  size_t bufferIndex = 0;
  for (auto &exRegion : _exportRegions) {
    auto index = convRelNeighboursToIndex(exRegion.getNeighbour());
    auto neighbourRank = _allNeighbourIndices.at(index);
    if (neighbourRank != _ownRank) {
      particleCommunicator.sendParticlePositions(_regionBuffers.at(bufferIndex), neighbourRank);
    }
    ++bufferIndex;
  }
  // receive
  // NOTE Optimization: Could reserve buffer in advance
  bufferIndex = 0;
  _importParticleMap.clear();
  for (auto &imRegion : _importRegions) {
    _importBuffers.at(bufferIndex).clear();
    auto index = convRelNeighboursToIndex(imRegion.getNeighbour());
    auto neighbourRank = _allNeighbourIndices.at(index);
    if (neighbourRank != _ownRank) {
      particleCommunicator.receiveParticlePositions(_importBuffers.at(bufferIndex), neighbourRank);
    } else {
      _importBuffers.at(bufferIndex)
          .insert(_importBuffers.at(bufferIndex).end(), _regionBuffers.at(bufferIndex).begin(),
                  _regionBuffers.at(bufferIndex).end());
    }
    for (size_t i = 0; i < _importBuffers.at(bufferIndex).size(); ++i) {
      ParticleType &particle = _importBuffers.at(bufferIndex).at(i);
      _importParticleMap.insert_or_assign({particle.getID(), particle.getR()}, std::ref(particle));
    }
    autoPasContainer.addHaloParticles(_importBuffers.at(bufferIndex));
    ++bufferIndex;
  }
  particleCommunicator.waitForSendRequests();
}

void HalfShell::SendAndReceiveResults(AutoPasType &autoPasContainer) {
  ParticleCommunicator particleCommunicator(_comm);
  // send results
  size_t bufferIndex = 0;
  for (auto &imRegion : _importRegions) {
    auto index = convRelNeighboursToIndex(imRegion.getNeighbour());
    auto neighbourRank = _allNeighbourIndices.at(index);
    if (neighbourRank != _ownRank) {
      particleCommunicator.sendParticleForces(_importBuffers.at(bufferIndex), neighbourRank);
    } else {
      _regionBuffers.at(bufferIndex).clear();
      // NOTE: We can only add the results inside the container if
      // we do not have sent results in to the home box from both directions
      // <- which is guaranteed no the case for HalfShell
      _regionBuffers.at(bufferIndex)
          .insert(_regionBuffers.at(bufferIndex).end(), _importBuffers.at(bufferIndex).begin(),
                  _importBuffers.at(bufferIndex).end());
    }
    ++bufferIndex;
  }

  // receive results
  bufferIndex = 0;
  for (auto &exRegion : _exportRegions) {
    auto index = convRelNeighboursToIndex(exRegion.getNeighbour());
    auto neighbourRank = _allNeighbourIndices.at(index);
    if (neighbourRank != _ownRank) {
      auto size = _regionBuffers.at(bufferIndex).size();
      _regionBuffers.at(bufferIndex).clear();
      _regionBuffers.at(bufferIndex).reserve(size);
      particleCommunicator.receiveParticleForces(_regionBuffers.at(bufferIndex), neighbourRank);
    }
    ++bufferIndex;
  }
  particleCommunicator.waitForSendRequests();

  // save results to container
  using namespace autopas::utils::ArrayMath::literals;
  bufferIndex = 0;
  // for all exported regions
  for (auto &exRegion : _exportRegions) {
    for (auto &result : _regionBuffers.at(bufferIndex)) {
      auto &particle = _exportParticleMap.at(result.getID());
      particle.get().addF(result.getF());
    }
    ++bufferIndex;
  }
}

void HalfShell::recollectResultsFromContainer(AutoPasType &autoPasContainer) {
  // iterate over halo particles and overwrite them in the respective buffer
  for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::halo); iter.isValid(); ++iter) {
    auto &b = _importParticleMap.at({iter->getID(), iter->getR()});
    b.get() = *iter;
  }
}

void HalfShell::calculateZonalInteractionPairwise(std::string zone1, std::string zone2,
                                                  std::function<void(ParticleType &, ParticleType &)> aosFunctor) {}

void HalfShell::calculateZonalInteractionTriwise(
    std::string zone, std::function<void(ParticleType &, ParticleType &, ParticleType &, bool)> aosFunctor) {}

void HalfShell::resizeHomeBoxRegion(RectRegion homeBoxRegion) {
  _homeBoxRegion = homeBoxRegion;
  _exportRegions.clear();
  _importRegions.clear();
  _exportRegions.reserve(_regionCount);
  _importRegions.reserve(_regionCount);

  auto hsCondition = [](const int d[3]) {
    /**
     * Stencil:
     *  z > 0 +
     *  z == 0 and y > 0 +
     *  z == 0 and y == 0 and x > 0
     */
    return d[2] > 0 or (d[2] == 0 and (d[1] > 0 or (d[1] == 0 and d[0] > 0)));
  };

  auto importCondition = [hsCondition](const int d[3]) { return !hsCondition(d); };

  auto identifyZone = [](const int d[3]) { return "A"; };

  // calculate exportRegions
  getRectRegionsConditional(_homeBoxRegion, _cutoff, _verletSkinWidth, _exportRegions, hsCondition, identifyZone,
                            false);

  // calculate importRegions
  getRectRegionsConditional(_homeBoxRegion, _cutoff, _verletSkinWidth, _importRegions, importCondition, identifyZone,
                            true);

  std::reverse(_importRegions.begin(), _importRegions.end());
}

const std::vector<RectRegion> HalfShell::getExportRegions() { return _exportRegions; }

const std::vector<RectRegion> HalfShell::getImportRegions() { return _importRegions; }
