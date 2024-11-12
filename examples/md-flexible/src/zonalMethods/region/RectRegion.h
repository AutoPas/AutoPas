
#pragma once

#include <array>

#include "src/zonalMethods/region/Region.h"

/*
 * This class represents the geographic regions of a import / export zone of a zonal method.
 * The described zone is a rectangular parallelepiped.
 * */
class RectRegion : Region {

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
  RectRegion(std::array<double, 3> origin, std::array<double, 3> size, char zoneID = 0)
      : _origin(origin), _size(size), Region(zoneID) {}

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
};
