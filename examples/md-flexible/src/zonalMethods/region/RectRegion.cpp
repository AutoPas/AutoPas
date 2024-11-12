
#include "src/zonalMethods/region/RectRegion.h"

#include "autopas/utils/ArrayMath.h"

void RectRegion::collectParticles(AutoPasType &autoPasContainer, std::vector<ParticleType> &buffer) {
  using namespace autopas::utils::ArrayMath::literals;
  auto boxMax = _origin + _size;
  for (auto particleIter = autoPasContainer.getRegionIterator(_origin, boxMax, autopas::IteratorBehavior::owned);
       particleIter.isValid(); ++particleIter) {
    buffer.push_back(*particleIter);
  }
}
