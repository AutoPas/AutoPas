/**
 * @file DifferentParticlesTest.cpp
 * @author seckler
 * @date 20.02.2020
 */

#include "DifferentParticlesTest.h"

#include "autopas/AutoPasDecl.h"
#include "testingHelpers/NonConstructibleParticle.h"
#include "testingHelpers/commonTypedefs.h"

extern template class autopas::AutoPas<NonConstructibleParticle>;
extern template bool autopas::AutoPas<NonConstructibleParticle>::computeInteractions(
    MockPairwiseFunctor<NonConstructibleParticle> *);

/**
 * Tests if AutoPas still compiles with a Particle that implements the normal interface, BUT no constructor.
 */
TEST_F(DifferentParticlesTest, testNonConstructibleParticle) {
  autopas::AutoPas<NonConstructibleParticle> autoPas;
  autoPas.setBoxMax({10., 10., 10.});
  autoPas.init();

  // We also check if computeInteractions can be instantiated.
  MockPairwiseFunctor<NonConstructibleParticle> functor;
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, isVecPatternAllowed(::testing::_)).WillRepeatedly(::testing::Return(true));
  autoPas.computeInteractions(&functor);
}