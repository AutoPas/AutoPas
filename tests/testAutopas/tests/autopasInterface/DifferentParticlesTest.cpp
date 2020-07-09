/**
 * @file DifferentParticlesTest.cpp
 * @author seckler
 * @date 20.02.2020
 */

#include "DifferentParticlesTest.h"

#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * A particle class without an actual constructor (only copy, etc.).
 */
class NonConstructibleParticle : public Particle {
 public:
  /**
   * Default constructor.
   */
  NonConstructibleParticle() = default;
};

/**
 * Tests if AutoPas still compiles with a Particle that implements the normal interface, BUT no constructor.
 */
TEST_F(DifferentParticlesTest, testNonConstructibleParticle) {
  autopas::AutoPas<NonConstructibleParticle> autoPas;
  autoPas.setBoxMax({10., 10., 10.});
  autoPas.init();

  // We also check if iteratePairwise can be instantiated.
  MockFunctor<NonConstructibleParticle, FPCell> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(false));
  autoPas.iteratePairwise(&functor);
}