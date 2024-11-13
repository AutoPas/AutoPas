
#include "src/zonalMethods/HalfShell.h"

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
    region.collectParticles(autoPasContainer, _regionBuffers[index++]);
  }
}
