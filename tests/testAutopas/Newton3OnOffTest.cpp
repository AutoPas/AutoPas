/**
 * @file Newton3OnOffTest.cpp
 * @author seckler
 * @date 18.04.18
 */

#include "Newton3OnOffTest.h"
#include "mocks/MockFunctor.h"

TEST_F(Newton3OnOffTest, test) {
  MockFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>>
      f;
}