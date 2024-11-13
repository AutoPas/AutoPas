
#include "src/zonalMethods/region/RectRegion.h"

#include "autopas/utils/ArrayMath.h"

RectRegion::RectRegion(std::array<double, 3> origin, std::array<double, 3> size, char zoneID)
    : _origin(origin), _size(size), Region(zoneID) {
  normalize();
}

void RectRegion::collectParticles(AutoPasType &autoPasContainer, std::vector<ParticleType> &buffer) {
  using namespace autopas::utils::ArrayMath::literals;
  auto boxMax = _origin + _size;
  for (auto particleIter = autoPasContainer.getRegionIterator(_origin, boxMax, autopas::IteratorBehavior::owned);
       particleIter.isValid(); ++particleIter) {
    buffer.push_back(*particleIter);
  }
}

bool RectRegion::operator==(const RectRegion &other) const { return _origin == other._origin and _size == other._size; }

void RectRegion::normalize() {
  for (int i = 0; i < 3; ++i) {
    if (_size[i] < 0) {
      _origin[i] += _size[i];
      _size[i] *= -1;
    }
  }
}
