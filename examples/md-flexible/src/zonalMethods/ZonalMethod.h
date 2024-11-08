
#include "src/zonalMethods/zone/Zone.h"

/**
 * This abstract class represents a zonal method.
 * All specific zonal method classes should extend this class.
 * */
class ZonalMethod {
 private:
  /**
   * Stores the number of zones.
   */
  unsigned int _zoneCount;

  /**
   * Stores a pointer to a dynamically allocated array of Zones.
   * The array is allocated in the constructor with the required number of zones
   * for the specific method.
   */
  std::unique_ptr<std::unique_ptr<Zone>[]> _zones;

 public:
  /**
   * Type for the AutoPas container
   */
  using AutoPasType = autopas::AutoPas<ParticleType>;

  /**
   * Constructor
   * @param zoneCount
   */
  ZonalMethod(unsigned int zoneCount);

  /**
   * Collect particles from the AutoPas container and store them internally.
   * @param autoPasContainer
   */
  virtual void collectParticles(AutoPasType &autoPasContainer) = 0;
};
