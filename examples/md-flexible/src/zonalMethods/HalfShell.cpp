
#include "src/zonalMethods/HalfShell.h"

#include "src/ParticleCommunicator.h"

HalfShell::HalfShell(RectRegion homeBoxRegion, double cutoff, double verletSkinWidth) : ZonalMethod(1) {
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

  // calculate exportRegions
  getRectRegionsConditional(homeBoxRegion, cutoff, verletSkinWidth, _exportRegions, hsCondition, false);

  // calculate importRegions
  getRectRegionsConditional(homeBoxRegion, cutoff, verletSkinWidth, _importRegions, hsCondition);

  // for the RegionIterator to operate properly, normalize the regions
  for (auto &ex : _exportRegions) ex.normalize();

  for (auto &im : _importRegions) im.normalize();
}

HalfShell::~HalfShell() = default;

void HalfShell::collectParticles(AutoPasType &autoPasContainer) {
  size_t index = 0;
  for (auto &region : _exportRegions) {
    // NOTE Optimization: Could reserve buffer in advance
    _regionBuffers[index].clear();
    region.collectParticles(autoPasContainer, _regionBuffers[index++]);
  }
}

void HalfShell::SendAndReceiveExports(AutoPasType &autoPasContainer, autopas::AutoPas_MPI_Comm comm,
                                      std::array<int, 26> allNeighbourIndices) {
  ParticleCommunicator particleCommunicator(comm);
  size_t index = 0;
  for (auto &exRegion : _exportRegions) {
    auto index = convRelNeighboursToIndex(exRegion.getNeighbour());
    auto neighbourRank = allNeighbourIndices.at(index);
    particleCommunicator.sendParticles(_regionBuffers[index++], neighbourRank);
  }
  // receive
  // NOTE Optimization: Could reserve buffer in advance
  _importParticles.clear();
  for (auto &imRegion : _importRegions) {
    auto index = convRelNeighboursToIndex(imRegion.getNeighbour());
    auto neighbourRank = allNeighbourIndices.at(index);
    particleCommunicator.receiveParticles(_importParticles, neighbourRank);
  }
  particleCommunicator.waitForSendRequests();

  autoPasContainer.addHaloParticles(_importParticles);
}

void HalfShell::SendAndReceiveResults(AutoPasType &autoPasContainer, autopas::AutoPas_MPI_Comm comm,
                                      std::array<int, 26> allNeighbourIndices) {}
