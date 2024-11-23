
#include "src/zonalMethods/ZonalMethod.h"

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Math.h"

ZonalMethod::ZonalMethod(unsigned int zoneCount, int ownRank, RectRegion homeBoxRegion, RectRegion globalBoxRegion,
                         autopas::AutoPas_MPI_Comm comm, std::array<int, 26> allNeighbourIndices,
                         std::array<options::BoundaryTypeOption, 3> boundaryType)
    : _zoneCount(zoneCount),
      _ownRank(ownRank),
      _homeBoxRegion(homeBoxRegion),
      _globalBoxRegion(globalBoxRegion),
      _comm(comm),
      _allNeighbourIndices(allNeighbourIndices),
      _boundaryType(boundaryType),
      _interactionSchedule({}),
      _interactionZones({}) {
  using namespace autopas::utils::Math;
  // intialize _isAtGlobalBoundaryMin and _isAtGlobalBoundaryMax
  for (size_t i = 0; i < 3; i++) {
    if (isNearRel(_homeBoxRegion._origin.at(i), _globalBoxRegion._origin.at(i))) {
      _isAtGlobalBoundaryMin.at(i) = 1;
    } else {
      _isAtGlobalBoundaryMin.at(i) = 0;
    }
    if (isNearRel(_homeBoxRegion._origin.at(i) + _homeBoxRegion._size.at(i),
                  _globalBoxRegion._origin.at(i) + _globalBoxRegion._size.at(i))) {
      _isAtGlobalBoundaryMax.at(i) = 1;
    } else {
      _isAtGlobalBoundaryMax.at(i) = 0;
    }
  }
}

ZonalMethod::~ZonalMethod() = default;

void ZonalMethod::getRectRegionsConditional(RectRegion &homeBoxRegion, double cutoffRadius, double verletSkinWidth,
                                            std::vector<RectRegion> &regions,
                                            const std::function<bool(const int[3])> &condition,
                                            const std::function<char(const int[3])> &identifyZone, bool calcImports) {
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

        tmp.setZoneID(identifyZone(d));

        tmp.normalize();

        regions.push_back(tmp);
      }
    }
  }
}

size_t ZonalMethod::convRelNeighboursToIndex(std::array<int, 3> relNeighbour) {
  auto index = (relNeighbour[0] + 1) * 9 + (relNeighbour[1] + 1) * 3 + (relNeighbour[2] + 1);
  if (index > 13) {
    index--;
  }
  return index;
}

void ZonalMethod::wrapAroundPeriodicBoundary(std::array<int, 3> relNeighbour, std::vector<ParticleType> &particles) {
  using namespace autopas::utils::ArrayMath::literals;
  // check if neighbour is over a global boundary (all directions)
  std::array<double, 3> overGlobalBoundary = {0, 0, 0};
  for (size_t i = 0; i < 3; i++) {
    // if the boundary in this direction is periodic
    if (_boundaryType.at(i) == options::BoundaryTypeOption::periodic) {
      // check if the neighbour is over the boundary
      // and set the respective values
      if (_isAtGlobalBoundaryMin.at(i) && relNeighbour.at(i) < 0) {
        overGlobalBoundary.at(i) = 1;
      }
      if (_isAtGlobalBoundaryMax.at(i) && relNeighbour.at(i) > 0) {
        overGlobalBoundary.at(i) = -1;
      }
    }
  }

  // if there is no wrapping todo, return
  if (overGlobalBoundary == std::array<double, 3>{0, 0, 0}) {
    return;
  }

  // wrap around global box
  for (ParticleType &p : particles) {
    auto pos = p.getR();
    pos += overGlobalBoundary * _globalBoxRegion._size;
    p.setR(pos);
  }
}

bool ZonalMethod::needToCollectParticles(std::array<int, 3> relNeighbour) {
  bool allOverBorderArePeriodic = true;
  for (size_t i = 0; i < 3; i++) {
    // check if we are going over a global border
    bool overBoder = (relNeighbour.at(i) < 0 and _isAtGlobalBoundaryMin.at(i)) or
                     (relNeighbour.at(i) > 0 and _isAtGlobalBoundaryMax.at(i));
    // if overBorder and not periodic
    if (overBoder and _boundaryType.at(i) != options::BoundaryTypeOption::periodic) {
      allOverBorderArePeriodic = false;
    }
  }
  return allOverBorderArePeriodic;
}
