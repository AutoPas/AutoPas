
#include "src/zonalMethods/ZonalMethod.h"

ZonalMethod::ZonalMethod(unsigned int zoneCount) : _zoneCount(zoneCount) {}

ZonalMethod::~ZonalMethod() = default;

void ZonalMethod::getRectRegionsConditional(RectRegion &homeBoxRegion, double cutoffRadius, double verletSkinWidth,
                                            std::vector<RectRegion> &regions,
                                            const std::function<bool(const int[3])> &condition, bool calcImports) {
  // factor for calculating import or export regions
  double factor = 1.0;
  if (!calcImports) {
    factor = -1.0;
  }

  // iterate over all neighbours
  int d[3];
  for (d[0] = -1; d[0] <= 1; d[0]++) {
    for (d[1] = -1; d[1] <= 1; d[1]++) {
      for (d[2] = -1; d[2] <= 1; d[2]++) {
        if ((d[0] || d[1] || d[2]) == 0)  // if all are 0 (false), then continue
          continue;

        // we don't include anything, that does not conform with the condition
        if (!condition(d)) {
          continue;
        }

        RectRegion tmp = homeBoxRegion;
        for (unsigned int dimension = 0; dimension < 3; dimension++) {
          if (d[dimension] == 0) {
            tmp._origin[dimension] = homeBoxRegion._origin[dimension];
            tmp._size[dimension] = homeBoxRegion._size[dimension];
          } else if (d[dimension] == -1) {  // LOWER
            tmp._origin[dimension] = homeBoxRegion._origin[dimension];
            tmp._size[dimension] = factor * (-cutoffRadius - verletSkinWidth);
          } else {  //  UPPER
            tmp._origin[dimension] = homeBoxRegion._origin[dimension] + homeBoxRegion._size[dimension];
            tmp._size[dimension] = factor * (cutoffRadius + verletSkinWidth);
          }
        }

        tmp.setNeighbour({d[0], d[1], d[2]});

        regions.push_back(tmp);
      }
    }
  }
}
