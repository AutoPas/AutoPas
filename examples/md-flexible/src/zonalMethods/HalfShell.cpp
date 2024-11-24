
#include "src/zonalMethods/HalfShell.h"

#include "src/ParticleCommunicator.h"

HalfShell::HalfShell(double cutoff, double verletSkinWidth, int ownRank, RectRegion homeBoxRegion,
                     RectRegion globalBoxRegion, autopas::AutoPas_MPI_Comm comm,
                     std::array<int, 26> allNeighbourIndices, std::array<options::BoundaryTypeOption, 3> boundaryType)
    : ZonalMethod(1, ownRank, homeBoxRegion, globalBoxRegion, comm, allNeighbourIndices, boundaryType) {
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

  auto identifyZone = [](const int d[3]) { return 'A'; };

  // calculate exportRegions
  getRectRegionsConditional(_homeBoxRegion, cutoff, verletSkinWidth, _exportRegions, hsCondition, identifyZone, false);

  // calculate importRegions
  getRectRegionsConditional(_homeBoxRegion, cutoff, verletSkinWidth, _importRegions, hsCondition, identifyZone, true);

  _interactionZones.push_back('A');
}

HalfShell::~HalfShell() = default;

void HalfShell::collectParticles(AutoPasType &autoPasContainer) {
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

void HalfShell::SendAndReceiveExports(AutoPasType &autoPasContainer) {
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
  _importParticles.clear();
  for (auto &imRegion : _importRegions) {
    auto index = convRelNeighboursToIndex(imRegion.getNeighbour());
    auto neighbourRank = _allNeighbourIndices.at(index);
    if (neighbourRank != _ownRank) {
      particleCommunicator.receiveParticles(_importParticles, neighbourRank);
    } else {
      _importParticles.insert(_importParticles.end(), _regionBuffers[bufferIndex].begin(),
                              _regionBuffers[bufferIndex].end());
    }
    ++bufferIndex;
  }
  particleCommunicator.waitForSendRequests();

  autoPasContainer.addHaloParticles(_importParticles);
}

void HalfShell::SendAndReceiveResults(AutoPasType &autoPasContainer) {
  // get cell
}

void HalfShell::calculateZonalInteractionPairwise(char zone1, char zone2,
                                                  std::function<void(ParticleType &, ParticleType &)> aosFunctor) {}
