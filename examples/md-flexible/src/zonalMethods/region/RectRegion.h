
#pragma once

#include <array>

#include "src/zonalMethods/region/Region.h"

/*
 * This class represents the geographic regions of a import / export zone of a zonal method.
 * The described zone is a rectangular parallelepiped.
 * */
class RectRegion : public Region {
 public:
  /**
   * Stores the origin of the zone.
   * */
  std::array<double, 3> _origin;

  /**
   * Stores the size of the zone.
   * */
  std::array<double, 3> _size;

  /**
   * Constructor
   * @param origin
   * @param size
   * @param zoneID
   */
  RectRegion(std::array<double, 3> origin, std::array<double, 3> size, char zoneID = 0);

  /**
   * Constructor
   */
  RectRegion() = default;

  /**
   * Collect particles from the AutoPas container and store them in the given buffer.
   * @param autoPasContainer
   * @param buffer
   */
  void collectParticles(AutoPasType &autoPasContainer, std::vector<ParticleType> &buffer) override;

  /**
   * Equality operator
   * @param other The other RectRegion to compare with.
   * @return true if the origins and sizes of both RectRegions are equal.
   */
  bool operator==(const RectRegion &other) const;

  /**
   * Normalizes the origin and size of the region such that all components of
   * the size vector are positive,
   * without changing the specified region.
   * */
  void normalize();
};

#include "autopas/utils/ArrayUtils.h"
using namespace autopas::utils::ArrayUtils;

inline std::ostream &operator<<(std::ostream &strm, const RectRegion &r) {
  return strm << "Region( origin: " << to_string(r._origin) << " | size: " << to_string(r._size) << " )";
}
