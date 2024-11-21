
#include "src/zonalMethods/FullShell.h"

#include "autopas/AutoPas.h"
#include "src/ParticleCommunicator.h"

FullShell::FullShell(double cutoff, double verletSkinWidth, int ownRank,
                     RectRegion homeBoxRegion, RectRegion globalBoxRegion, autopas::AutoPas_MPI_Comm comm,
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

  auto identifyZone = [](const int d[3]) { return 'A'; };

  // calculate exportRegions
  getRectRegionsConditional(_homeBoxRegion, cutoff, verletSkinWidth, _exportRegions, fsCondition, identifyZone, false);

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

void FullShell::SendAndReceiveExports(AutoPasType &autoPasContainer) {
  ParticleCommunicator particleCommunicator(_comm);
  size_t index = 0;
  for (auto &exRegion : _exportRegions) {
    auto index = convRelNeighboursToIndex(exRegion.getNeighbour());
    auto neighbourRank = _allNeighbourIndices.at(index);
    if (neighbourRank == _ownRank) {
      // if periodic boundary, add halo particles manually
      bool periodic = false;
      // check if we are at a periodic border
      for (int i = 0; i < 3; i++) {
        bool atBorder = false;

        if (i != 0 && _boundaryType[i] == options::BoundaryTypeOption::periodic) {
          periodic = true;
        }
      }
      if (periodic) {
      }
      continue;
    }
    particleCommunicator.sendParticles(_regionBuffers[index++], neighbourRank);
  }
  // receive
  // NOTE Optimization: Could reserve buffer in advance
  _importParticles.clear();
  for (auto &imRegion : _importRegions) {
    auto index = convRelNeighboursToIndex(imRegion.getNeighbour());
    auto neighbourRank = _allNeighbourIndices.at(index);
    if (neighbourRank == _ownRank) {
      // NOTE: if periodic boundary, add halo particles
      continue;
    }
    particleCommunicator.receiveParticles(_importParticles, neighbourRank);
  }
  particleCommunicator.waitForSendRequests();

  autoPasContainer.addHaloParticles(_importParticles);
}

void FullShell::SendAndReceiveResults(AutoPasType &autoPasContainer) {}
