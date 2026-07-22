/**
 * @file NeighborIdentificationFunctorAutoPasCompileTest.cpp
 * @author S. Newcome
 * @date 22.07.2026
 */

#include <gtest/gtest.h>

#include "autopas/AutoPas.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "neighborIdentification/NeighborIdentificationFunctor.h"

namespace {

/**
 * Checks that AutoPas compiles with the public-facing NeighborIdentificationFunctor: calling computeInteractions() with
 * it forces AutoPas to instantiate its full container/traversal/tuning machinery for the functor type. The container is
 * intentionally left empty, so the value of this test is the compile-time instantiation rather than any runtime result.
 */
TEST(NeighborIdentificationFunctor_Test, CompilesWithAutoPas) {
  autopas::AutoPas<autopas::ParticleBaseFP64> autoPas;
  autoPas.setBoxMin({0., 0., 0.});
  autoPas.setBoxMax({10., 10., 10.});
  autoPas.setCutoff(1.0);
  autoPas.setVerletSkin(0.2);
  autoPas.init();

  autopas::NeighborIdentificationFunctor<autopas::ParticleBaseFP64>::NeighborListAoSType neighborLists;
  autopas::NeighborIdentificationFunctor<autopas::ParticleBaseFP64> functor(neighborLists, /*interactionLength*/ 1.0,
                                                                            /*gatherNewton3Lists*/ false);

  autoPas.computeInteractions(&functor);
}

}  // namespace