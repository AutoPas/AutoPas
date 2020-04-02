/**
 * @file ParticleDeletedObserver.h
 * @author seckler
 * @date 27.03.20
 */

#pragma once

namespace autopas::internal {

/**
 * Class that is notified when a particle is deleted.
 * For VerletClusterLists
 */
class ParticleDeletedObserver {
 public:
  /**
   * This function is called when a particle is deleted.
   */
  virtual void notifyParticleDeleted() = 0;
};

}  // namespace autopas::internal