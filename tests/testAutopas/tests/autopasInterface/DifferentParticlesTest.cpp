/**
 * @file DifferentParticlesTest.cpp
 * @author seckler
 * @date 20.02.2020
 */

#include "DifferentParticlesTest.h"

#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * A particle class with only a default constructor.
 */
class OnlyDefaultConstructibleParticle : public Particle {
 public:
  /**
   * Default constructor.
   */
  OnlyDefaultConstructibleParticle() = default;
};

/**
 * Tests if AutoPas still compiles with a Particle that implements the normal interface, BUT implements a different
 * compiler.
 */
TEST_F(DifferentParticlesTest, testOnlyDefaultConstructibleParticle) {
  autopas::AutoPas<OnlyDefaultConstructibleParticle, autopas::FullParticleCell<OnlyDefaultConstructibleParticle>>
      autoPas;
  autoPas.setBoxMin({0., 0., 0.});
  autoPas.setBoxMax({10., 10., 10.});
  autoPas.setCutoff(1.);
  autoPas.init();
}

/**
 * A particle class with only a non-default constructor.
 */
class NonDefaultConstructibleParticle : public Particle {
 public:
  /**
   * Non-default constructor.
   */
  NonDefaultConstructibleParticle(int, int){};
};

/**
 * Tests if AutoPas still compiles with a Particle that implements the normal interface, BUT implements a different
 * compiler.
 */
TEST_F(DifferentParticlesTest, testNonDefaultConstructibleParticle) {
  autopas::AutoPas<NonDefaultConstructibleParticle, autopas::FullParticleCell<NonDefaultConstructibleParticle>> autoPas;
  autoPas.setBoxMin({0., 0., 0.});
  autoPas.setBoxMax({10., 10., 10.});
  autoPas.setCutoff(1.);
  autoPas.init();
}