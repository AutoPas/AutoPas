
#include "src/zonalMethods/FullShell.h"

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

  // calculate exportRegions
  getRectRegionsConditional(homeBoxRegion, cutoff, verletSkinWidth, _exportRegions, fsCondition, false);

  // calculate importRegions - NOTE: no real need to calculate
  getRectRegionsConditional(homeBoxRegion, cutoff, verletSkinWidth, _importRegions, fsCondition);

  // for the RegionIterator to operate properly, normalize the regions
  for (auto &ex : _exportRegions) ex.normalize();

  for (auto &im : _importRegions) im.normalize();
}

FullShell::~FullShell() = default;

void FullShell::collectParticles(AutoPasType &autoPasContainer) {
  size_t index = 0;
  for (auto &region : _exportRegions) {
    region.collectParticles(autoPasContainer, _regionBuffers[index++]);
  }
}
