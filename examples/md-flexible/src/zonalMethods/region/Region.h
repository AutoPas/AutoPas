
#pragma once
#include <spdlog/spdlog.h>

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
   * Check if the given position is inside the region
   * @param position
   * @return
   */
  virtual bool contains(std::array<double, 3> position) = 0;

  /**
   * Constructor of Region
   * @param zoneID
   */
  Region(std::string zoneID = "I") : _zoneID(zoneID) {}

  /**
   * Destructor
   */
  virtual ~Region() = default;

 protected:
  /**
   * Stores the id of which zone it belongs to
   * */
  std::string _zoneID;

  /**
   * Stores the relative postitioning of the communication partner of this zone
   * */
  std::array<int, 3> _neighbour;

 public:
  /**
   * Returns the id of the zone
   */
  inline std::string getZoneID() const { return _zoneID; }

  /**
   * Sets the id of the zone
   */
  inline void setZoneID(std::string zoneID) { _zoneID = zoneID; }

  /**
   * Returns the relative postitioning of the communication partner of this zone
   */
  inline std::array<int, 3> getNeighbour() const { return _neighbour; }

  /**
   * Sets the relative postitioning of the communication partner of this zone
   */
  inline void setNeighbour(std::array<int, 3> neighbour) { _neighbour = neighbour; }
};
