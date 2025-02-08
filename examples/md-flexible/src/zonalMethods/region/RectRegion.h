
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
   * Constructor of RectRegion
   * Normalizes origin and size when called.
   * @param origin
   * @param size
   * @param zoneID
   */
  RectRegion(std::array<double, 3> origin, std::array<double, 3> size, std::string zoneID = "I");

  /**
   * Constructor
   */
  RectRegion() = default;

  /**
   * Collect particles from the AutoPas container and store them in the given buffer.
   * @param autoPasContainer
   * @param buffer
   */
  void collectParticles(AutoPasType &autoPasContainer, std::vector<ParticleType> &buffer,
                        std::optional<ParticleMap> = std::nullopt) override;

  /**
   * Check if the given position is inside the region
   * @param position
   * @return
   */
  bool contains(std::array<double, 3> position) override;

  /**
   * Equality operator
   * @param other The other RectRegion to compare with.
   * @return true if the origins and sizes of both RectRegions are equal.
   */
  bool operator==(const RectRegion &other) const;

  /**
   * Normalizes the origin and size of the region such that all components of
   * the _size vector are positive,
   * without changing the specified region.
   * */
  void normalize();

  /**
   * Returns a string representation of the region.
   * @return
   */
  std::string toString() const;
};

#include "autopas/utils/ArrayUtils.h"
using namespace autopas::utils::ArrayUtils;

inline std::ostream &operator<<(std::ostream &strm, const RectRegion &r) {
  return strm << "Region( origin: " << to_string(r._origin) << " | size: " << to_string(r._size) << " )";
}
