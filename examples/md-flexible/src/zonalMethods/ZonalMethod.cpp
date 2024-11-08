
#include "src/zonalMethods/ZonalMethod.h"

ZonalMethod::ZonalMethod(unsigned int zoneCount) : _zoneCount(zoneCount) {
  _zones = std::make_unique<std::unique_ptr<Zone>[]>(zoneCount);
}
