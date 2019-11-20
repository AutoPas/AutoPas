/**
 * @file TimeDiscretizationTest.cpp
 * @author N. Fottner
 * @date 05/22/19.
 */

#include "TimeDiscretizationTest.h"

void TimeDiscretizationTest::fillWithParticlesAndInit(autopas::AutoPas<PrintableMolecule,
                                                                       autopas::FullParticleCell<PrintableMolecule>> &autopas) {
  autopas.setBoxMin({0., 0., 0.});
  autopas.setBoxMax({5., 5., 5.});
  autopas.init();
  PrintableMolecule dummy;
  dummy.setF({0., 0., 1.});
  dummy.setV({0., 0., 1.});
  GridGenerator::fillWithParticles(autopas, {2, 2, 2}, dummy, {1, 1, 1}, {0., 0., 0.});
}

TEST_F(TimeDiscretizationTest, calcVelocities) {
  auto autoPas = autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>();
  fillWithParticlesAndInit(autoPas);

  TimeDiscretization<decltype(autoPas), decltype(_particlePropertiesLibrary)> timeDiscretization(0.1, _particlePropertiesLibrary);

  timeDiscretization.calculateVelocities(autoPas);
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
      // only velocity in one direction is expected
      EXPECT_EQ(iter->getV()[0], 0);
      EXPECT_EQ(iter->getV()[1], 0);
      // Strömer-Verlet: (0+1)/2 * 0.1 + 0= 0.05
      EXPECT_NEAR(iter->getV()[2], 0.05, 1e-13);

      // set force for next iteration
      iter->setOldF(iter->getF());
      iter->setF({0, 0, 2});
  }

  timeDiscretization.calculateVelocities(autoPas);
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    // only velocity in one direction is expected
    EXPECT_EQ(iter->getV()[0], 0);
    EXPECT_EQ(iter->getV()[1], 0);
    // Strömer-Verlet: (1+2)/2 * 0.1 + 0.05 = 0.2
    EXPECT_NEAR(iter->getV()[2], 0.2, 1e-13);
  }
}

TEST_F(TimeDiscretizationTest, calcPositions) {
  auto autoPas = autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>();
  auto autoPasRef = autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>();
  fillWithParticlesAndInit(autoPas);
  fillWithParticlesAndInit(autoPasRef);

  TimeDiscretization<decltype(autoPas), decltype(_particlePropertiesLibrary)> timeDiscretization(0.1, _particlePropertiesLibrary);

  timeDiscretization.calculatePositions(autoPas);
  for (auto iter = autoPas.begin(), iterRef = autoPasRef.begin(); iter.isValid(); ++iter, ++iterRef) {
    // only change in one direction is expected
    EXPECT_EQ(iter->getR()[0], iterRef->getR()[0]);
    EXPECT_EQ(iter->getR()[1], iterRef->getR()[1]);
    // Strömer-Verlet: 0.1 * 1 + 0.1^2 * (1 / 2) = 0.105
    EXPECT_NEAR(iter->getR()[2], iterRef->getR()[2] + 0.105, 1e-13);

    // expect force to be reset
    EXPECT_THAT(iter->getF(), ::testing::ElementsAreArray({0,0,0}));

    // set force and velocity for next iteration
    iter->setF({0, 0, 2});
    iter->setV({0, 0, .5});

    // update reference position
    iterRef->setR(iter->getR());
  }

  timeDiscretization.calculatePositions(autoPas);
  for (auto iter = autoPas.begin(), iterRef = autoPasRef.begin(); iter.isValid(); ++iter, ++iterRef) {
    // only velocity in one direction is expected
    EXPECT_EQ(iter->getR()[0], iterRef->getR()[0]);
    EXPECT_EQ(iter->getR()[1], iterRef->getR()[1]);
    // Strömer-Verlet: 0.1 * .5 + 0.1^2 * (2 / 2) = 0.06
    EXPECT_NEAR(iter->getR()[2], iterRef->getR()[2] + 0.06, 1e-13);
  }
}