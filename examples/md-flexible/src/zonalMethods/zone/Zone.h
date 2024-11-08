
#include "autopas/AutoPasDecl.h"
#include "src/TypeDefinitions.h"

/*
 * This abstract class represents the geographic regions of a import / export zone of a zonal method.
 * All specific zonal method classes should extend this class.
 * */
class Zone {
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
};
