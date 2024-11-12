
#pragma once
#include "autopas/AutoPasDecl.h"
#include "src/TypeDefinitions.h"

/*
 * This abstract class represents a geographic region of a import / export zone of a zonal method.
 * All specific region classes should extend this class.
 * */
class Region {
 public:
  /**
   * Type for the AutoPas container
   */
  using AutoPasType = autopas::AutoPas<ParticleType>;

  /**
   * Collect particles from the AutoPas container and store them in the given buffer.
   * @param autoPasContainer
   * @param buffer
   */
  virtual void collectParticles(AutoPasType &autoPasContainer, std::vector<ParticleType> &buffer) = 0;

  /**
   * Constructor of Region
   * @param zoneID
   */
  Region(size_t zoneID = 0) : _zoneID(zoneID) {}

  /**
   * Destructor
   */
  virtual ~Region() = default;

 protected:
  /**
   * Stores the id of which zone it belongs to
   * */
  char _zoneID;

  /**
   * Stores the relative postitioning of the communication partner of this zone
   * */
  std::array<int, 3> _neighbour;
};
