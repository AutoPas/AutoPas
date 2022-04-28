/**
 * @file ReflectiveBoundaryConditionTest.cpp
 * @author F. Gratl
 * @date 26.04.22
 */

#include "ReflectiveBoundaryConditionTest.h"

#include "autopas/AutoPasDecl.h"
#include "autopas/utils/ArrayMath.h"
#include "src/domainDecomposition/RegularGridDecomposition.h"
#include "src/options/BoundaryTypeOption.h"

extern template class autopas::AutoPas<ParticleType>;

/**
 * Very simple test of reflective boundaries in all 3 dimension. Places 2 particles on every face, with one travelling
 * towards the boundary and the other away.
 */
TEST_P(ReflectiveBoundaryConditionTest, simpleReflectionTest) {
  // initialise AutoPas container & domainDecomposition
  const std::array<double, 3> boxMin = {0., 0., 0.};
  const std::array<double, 3> boxMax = {5., 5., 5.};
  const std::array<double, 3> boxLength = autopas::utils::ArrayMath::sub(boxMax, boxMin);
  const std::array<bool, 3> subdivideDimension = {true, true, true};
  const double cutoffWidth = 0.3;
  const double skinWidth = 0.2;
  const std::array<options::BoundaryTypeOption, 3> boundaryConditions = {options::BoundaryTypeOption::reflective,
                                                                         options::BoundaryTypeOption::reflective,
                                                                         options::BoundaryTypeOption::reflective};

  RegularGridDecomposition domainDecomposition(boxMin, boxMax, subdivideDimension, cutoffWidth, skinWidth,
                                               boundaryConditions);

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);

  autoPasContainer->setBoxMin(domainDecomposition.getLocalBoxMin());
  autoPasContainer->setBoxMax(domainDecomposition.getLocalBoxMax());
  autoPasContainer->setCutoff(cutoffWidth);
  autoPasContainer->setVerletSkin(skinWidth);
  autoPasContainer->init();

  // get particle properties
  const std::array<double, 3> particlePosition = std::get<0>(GetParam());
  const std::array<double, 3> particleVelocity = std::get<1>(GetParam());

  // derive expected position
  auto calcReflectedVelocity = [](const std::array<double, 3> velocities, const std::array<bool, 3> isReflected) {
    auto reflVel = velocities;
    for (int i = 0; i < 3; ++i) {
      if (isReflected[i]) {
        reflVel[i] *= -1;
      }
    }
    return reflVel;
  };

  const std::array<double, 3> expectedPosition = particlePosition;
  const std::array<double, 3> expectedVelocity = calcReflectedVelocity(particleVelocity, std::get<2>(GetParam()));

  // create particle and add to container
  // in the MPI case this is expected to be wrong for all but one rank
  if (domainDecomposition.isInsideLocalDomain(particlePosition)) {
    ParticleType particle;
    particle.setID(0);
    particle.setR(particlePosition);
    particle.setV(particleVelocity);
    autoPasContainer->addParticle(particle);
#if not defined(AUTOPAS_INCLUDE_MPI)
  } else {
    GTEST_FAIL() << "Test Particle is not in the box -> setup is wrong!";
#endif
  }

  // apply reflective BCs (as periodic + domain exchange)
  domainDecomposition.exchangeMigratingParticles(autoPasContainer);
  domainDecomposition.reflectParticlesAtBoundaries(autoPasContainer);
  domainDecomposition.exchangeHaloParticles(autoPasContainer);

  if (domainDecomposition.isInsideLocalDomain(expectedPosition)) {
    EXPECT_EQ(1, autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned));
    // check particles have been successfully reflected (and not translated)
    const auto &reflectedParticle = autoPasContainer->begin(autopas::IteratorBehavior::owned);
    const auto &reflectedPosition = reflectedParticle->getR();
    const auto &reflectedVelocity = reflectedParticle->getV();
    for (size_t i = 0; i < 3; ++i) {
      EXPECT_NEAR(reflectedPosition[i], expectedPosition[i], 1e-13) << "Unexpected position[" << i << "]";
      EXPECT_NEAR(reflectedVelocity[i], expectedVelocity[i], 1e-13) << "Unexpected velocity[" << i << "]";
    }
#if not defined(AUTOPAS_INCLUDE_MPI)
  } else {
    GTEST_FAIL() << "Expected position is not in the box -> setup is wrong!";
  }
  // check that there are no halo particles. We can only guarantee this in the non-MPI case
  EXPECT_EQ(0, autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::halo));
#else
  }
#endif
}

INSTANTIATE_TEST_SUITE_P(
    TestSimpleReflections, ReflectiveBoundaryConditionTest,
    testing::Values(/*position*/ /*velocity*/ /*is reflected*/
                    std::make_tuple(std::array<double, 3>{0.05, 2.50, 2.50}, std::array<double, 3>{1, 1, -1},
                                    std::array<bool, 3>{false, false, false}),
                    std::make_tuple(std::array<double, 3>{0.05, 2.50, 2.50}, std::array<double, 3>{-1, 1, -1},
                                    std::array<bool, 3>{true, false, false}),
                    std::make_tuple(std::array<double, 3>{4.95, 2.50, 2.50}, std::array<double, 3>{1, 1, -1},
                                    std::array<bool, 3>{true, false, false}),
                    std::make_tuple(std::array<double, 3>{4.95, 2.50, 2.50}, std::array<double, 3>{-1, 1, -1},
                                    std::array<bool, 3>{false, false, false}),

                    std::make_tuple(std::array<double, 3>{2.50, 0.05, 2.50}, std::array<double, 3>{1, 1, -1},
                                    std::array<bool, 3>{false, false, false}),
                    std::make_tuple(std::array<double, 3>{2.50, 0.05, 2.50}, std::array<double, 3>{1, -1, -1},
                                    std::array<bool, 3>{false, true, false}),
                    std::make_tuple(std::array<double, 3>{2.50, 4.95, 2.50}, std::array<double, 3>{1, 1, -1},
                                    std::array<bool, 3>{false, true, false}),
                    std::make_tuple(std::array<double, 3>{2.50, 4.95, 2.50}, std::array<double, 3>{1, -1, -1},
                                    std::array<bool, 3>{false, false, false}),

                    std::make_tuple(std::array<double, 3>{2.50, 2.50, 0.05}, std::array<double, 3>{1, -1, 1},
                                    std::array<bool, 3>{false, false, false}),
                    std::make_tuple(std::array<double, 3>{2.50, 2.50, 0.05}, std::array<double, 3>{1, -1, -1},
                                    std::array<bool, 3>{false, false, true}),
                    std::make_tuple(std::array<double, 3>{2.50, 2.50, 4.95}, std::array<double, 3>{1, -1, 1},
                                    std::array<bool, 3>{false, false, true}),
                    std::make_tuple(std::array<double, 3>{2.50, 2.50, 4.95}, std::array<double, 3>{1, -1, -1},
                                    std::array<bool, 3>{false, false, false}))
    //    ,ReflectiveBoundaryConditionTest::PrintToStringParamName());
);