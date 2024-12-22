
#include "src/zonalMethods/FullShell.h"

#include "autopas/AutoPas.h"
#include "src/ParticleCommunicator.h"

FullShell::FullShell(double cutoff, double verletSkinWidth, int ownRank, RectRegion homeBoxRegion,
                     RectRegion globalBoxRegion, autopas::AutoPas_MPI_Comm comm,
                     std::array<int, 26> allNeighbourIndices, std::array<options::BoundaryTypeOption, 3> boundaryType)
    : ZonalMethod(1, ownRank, homeBoxRegion, globalBoxRegion, comm, allNeighbourIndices, boundaryType) {
  _exportRegions.reserve(_regionCount);
  _importRegions.reserve(_regionCount);

  auto fsCondition = [](const int d[3]) {
    /**
     * Stencil:
     *   - export to all
     *   - no need to import
     */
    return true;
  };

  auto identifyZone = [](const int d[3]) { return "A"; };

  // calculate exportRegions
  getRectRegionsConditional(_homeBoxRegion, cutoff, verletSkinWidth, _exportRegions, fsCondition, identifyZone, false);

  // calculate importRegions
  getRectRegionsConditional(_homeBoxRegion, cutoff, verletSkinWidth, _importRegions, fsCondition, identifyZone, true);

  std::reverse(_importRegions.begin(), _importRegions.end());

  _interactionZones.push_back("A");
  _interactionSchedule.insert_or_assign("A", std::vector<std::string>{});
}

FullShell::~FullShell() = default;

void FullShell::collectParticles(AutoPasType &autoPasContainer) {
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

void FullShell::SendAndReceiveExports(AutoPasType &autoPasContainer) {
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
  _importParticles.clear();
  for (auto &imRegion : _importRegions) {
    auto index = convRelNeighboursToIndex(imRegion.getNeighbour());
    auto neighbourRank = _allNeighbourIndices.at(index);
    if (neighbourRank != _ownRank) {
      particleCommunicator.receiveParticles(_importParticles, neighbourRank);
    } else {
      _importParticles.insert(_importParticles.end(), _regionBuffers.at(bufferIndex).begin(),
                              _regionBuffers.at(bufferIndex).end());
    }

    ++bufferIndex;
  }
  particleCommunicator.waitForSendRequests();

  autoPasContainer.addHaloParticles(_importParticles);
}

void FullShell::SendAndReceiveResults(AutoPasType &autoPasContainer) {}

void FullShell::recollectResultsFromContainer(AutoPasType &autoPasContainer) {}

void FullShell::calculateZonalInteractionPairwise(std::string zone1, std::string zone2,
                                                  std::function<void(ParticleType &, ParticleType &)> aosFunctor) {}

const std::vector<RectRegion> FullShell::getExportRegions() { return _exportRegions; }

const std::vector<RectRegion> FullShell::getImportRegions() { return _importRegions; }
