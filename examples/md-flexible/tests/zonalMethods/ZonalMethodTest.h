
#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "src/zonalMethods/ZonalMethod.h"

/**
 * Test class for the Zone classes
 */
class ZonalMethodTest : public ZonalMethod, public AutoPasTestBase {
 public:
  /**
   * Constructor.
   */
  ZonalMethodTest()
      : ZonalMethod(1, 0, RectRegion({0, 0, 0}, {1, 1, 1}), RectRegion({0, 0, 0}, {1, 1, 1})), AutoPasTestBase() {}

  void collectParticles(AutoPasType &autoPasContainer) override;

  void SendAndReceiveExports(AutoPasType &autoPasContainer) override;

  void SendAndReceiveResults(AutoPasType &autoPasContainer) override;

  void calculateZonalInteractionPairwise(std::string zone1, std::string zone2,
                                         std::function<void(ParticleType &, ParticleType &)> aosFunctor) override;
  void calculateZonalInteractionTriwise(
      std::string zone, std::function<void(ParticleType &, ParticleType &, ParticleType &, bool)> aosFunctor) override;

  void recollectResultsFromContainer(AutoPasType &autoPasContainer) override;

  void resizeHomeBoxRegion(RectRegion homeBoxRegion) override;

};
