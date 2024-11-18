
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
  ZonalMethodTest() : ZonalMethod(1), AutoPasTestBase() {}

  void collectParticles(AutoPasType &autoPasContainer) override;

  void SendAndReceiveExports(AutoPasType &autoPasContainer, autopas::AutoPas_MPI_Comm comm,
                             std::array<int, 26> allNeighbourIndices, int ownRank) override;

  void SendAndReceiveResults(AutoPasType &autoPasContainer, autopas::AutoPas_MPI_Comm comm,
                             std::array<int, 26> allNeighbourIndices, int ownRank) override;
};
