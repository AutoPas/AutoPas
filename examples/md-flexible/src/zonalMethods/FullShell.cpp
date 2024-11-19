
#include "src/zonalMethods/FullShell.h"

#include "autopas/AutoPas.h"
#include "src/ParticleCommunicator.h"

FullShell::FullShell(RectRegion homeBoxRegion, double cutoff, double verletSkinWidth) : ZonalMethod(1) {
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

  auto identifyZone = [](const int d[3]) { return 'A'; };

  // calculate exportRegions
  getRectRegionsConditional(homeBoxRegion, cutoff, verletSkinWidth, _exportRegions, fsCondition, identifyZone, false);

  // skip calculation import regions since not needed

  _interactionSchedule.push_back('A');
}

FullShell::~FullShell() = default;

void FullShell::collectParticles(AutoPasType &autoPasContainer) {
  size_t index = 0;
  for (auto &region : _exportRegions) {
    // NOTE Optimization: Could reserve buffer in advance
    _regionBuffers[index].clear();
    region.collectParticles(autoPasContainer, _regionBuffers[index++]);
  }
}

void FullShell::SendAndReceiveExports(AutoPasType &autoPasContainer, autopas::AutoPas_MPI_Comm comm,
                                      std::array<int, 26> allNeighbourIndices, int ownRank) {
  ParticleCommunicator particleCommunicator(comm);
  size_t index = 0;
  for (auto &exRegion : _exportRegions) {
    auto index = convRelNeighboursToIndex(exRegion.getNeighbour());
    auto neighbourRank = allNeighbourIndices.at(index);
    if (neighbourRank == ownRank) continue;
    particleCommunicator.sendParticles(_regionBuffers[index++], neighbourRank);
  }
  // receive
  // NOTE Optimization: Could reserve buffer in advance
  _importParticles.clear();
  for (auto &imRegion : _importRegions) {
    auto index = convRelNeighboursToIndex(imRegion.getNeighbour());
    auto neighbourRank = allNeighbourIndices.at(index);
    if (neighbourRank == ownRank) continue;
    particleCommunicator.receiveParticles(_importParticles, neighbourRank);
  }
  particleCommunicator.waitForSendRequests();

  autoPasContainer.addHaloParticles(_importParticles);
}

void FullShell::SendAndReceiveResults(AutoPasType &autoPasContainer, autopas::AutoPas_MPI_Comm comm,
                                      std::array<int, 26> allNeighbourIndices, int ownRank) {}
