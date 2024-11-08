

#include <array>

#include "src/zonalMethods/zone/Zone.h"

/*
 * This class represents the geographic regions of a import / export zone of a zonal method.
 * The described zone is a rectangular parallelepiped.
 * */
class RectZone : Zone {
 private:
  /**
   * Stores the origin of the zone.
   * */
  std::array<double, 3> _origin;

  /**
   * Stores the size of the zone.
   * */
  std::array<double, 3> _size;

 public:
  /**
   * Constructor
   * @param origin
   * @param size
   */
  RectZone(std::array<double, 3> origin, std::array<double, 3> size) : _origin(origin), _size(size) {}

  /**
   * Getter for the origin
   * @return
   */
  inline std::array<double, 3> getOrigin() const { return _origin; }

  /**
   * Getter for the size
   * @return
   */
  inline std::array<double, 3> getSize() const { return _size; }

  /**
   * Collect particles from the AutoPas container and store them in the given buffer.
   * @param autoPasContainer
   * @param buffer
   */
  void collectParticles(AutoPasType &autoPasContainer, std::vector<ParticleType> &buffer) override;
};
