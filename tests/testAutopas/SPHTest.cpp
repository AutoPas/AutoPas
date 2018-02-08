//
// Created by seckler on 22.01.18.
//

#include "SPHTest.h"
#include "autopas.h"
#include "sph/autopassph.h"

TEST_F(SPHTest, testW) {
  double value = autopas::sph::SPHKernels::W({1., 1., 1.}, 1.);
  double should_be_value = 0.00944773;
  EXPECT_NEAR(value, should_be_value, 1e-8);

  value = autopas::sph::SPHKernels::W({1., .5, .25}, .5);
  should_be_value = 0.00151727;
  EXPECT_NEAR(value, should_be_value, 1e-8);
}

TEST_F(SPHTest, testGradW) {
  std::array<double, 3> value =
      autopas::sph::SPHKernels::gradW({1., 1., 1.}, 1.);
  std::array<double, 3> should_be_value = {-0.0213086, -0.0213086, -0.0213086};
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(value[i], should_be_value[i], 1e-7);
  }

  value = autopas::sph::SPHKernels::gradW({1., .5, .25}, .5);
  should_be_value = {-0.038073, -0.0190365, -0.00951825};
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(value[i], should_be_value[i], 1e-7);
  }
}

TEST_F(SPHTest, testSPHCalcDensityFunctor) {
  autopas::sph::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5,
                                         0.7, 0.6);
  autopas::sph::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5,
                                         1.3, 0.8);

  autopas::sph::SPHCalcDensityFunctor densityFunctor;
  densityFunctor.AoSFunctor(sphParticle1, sphParticle2);
  densityFunctor.AoSFunctor(sphParticle2, sphParticle1);

  EXPECT_NEAR(sphParticle1.getDensity(), 0.559026, 1e-6);
  EXPECT_NEAR(sphParticle2.getDensity(), 0.172401, 1e-6);
}

TEST_F(SPHTest, testSPHCalcPressure) {
  autopas::sph::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5,
                                         0.7, 0.6);
  autopas::sph::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5,
                                         1.3, 0.8);

  // simulate density functor call:
  sphParticle1.addDensity(0.559026);
  sphParticle2.addDensity(0.172401);

  // set pressure:
  sphParticle1.setEnergy(2.5);
  sphParticle2.setEnergy(2.5);
  sphParticle1.calcPressure();
  sphParticle2.calcPressure();

  EXPECT_NEAR(sphParticle1.getPressure(), 0.559026, 1e-6);
  EXPECT_NEAR(sphParticle2.getPressure(), 0.172401, 1e-6);
  EXPECT_NEAR(sphParticle1.getSoundSpeed(), 1.18322, 1e-5);
  EXPECT_NEAR(sphParticle2.getSoundSpeed(), 1.18322, 1e-5);
}

TEST_F(SPHTest, testSPHCalcHydroForceFunctor) {
  autopas::sph::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5,
                                         0.7, 0.6);
  autopas::sph::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5,
                                         1.3, 0.8);

  // simulate density functor call:
  sphParticle1.addDensity(0.559026);
  sphParticle2.addDensity(0.172401);

  // set pressure:
  sphParticle1.setEnergy(2.5);
  sphParticle2.setEnergy(2.5);
  sphParticle1.setPressure(0.559026);
  sphParticle2.setPressure(0.172401);
  sphParticle1.setSoundSpeed(1.18322);
  sphParticle2.setSoundSpeed(1.18322);

  autopas::sph::SPHCalcHydroForceFunctor hydroForceFunctor;
  hydroForceFunctor.AoSFunctor(sphParticle1, sphParticle2);
  hydroForceFunctor.AoSFunctor(sphParticle2, sphParticle1);

  EXPECT_NEAR(sphParticle1.getAcceleration().at(0), -2.26921, 1e-5);
  EXPECT_NEAR(sphParticle1.getAcceleration().at(1), -4.53843, 1e-5);
  EXPECT_NEAR(sphParticle1.getAcceleration().at(2), -6.80764, 1e-5);
  EXPECT_NEAR(sphParticle2.getAcceleration().at(0), 3.78202, 1e-5);
  EXPECT_NEAR(sphParticle2.getAcceleration().at(1), 7.56405, 1e-5);
  EXPECT_NEAR(sphParticle2.getAcceleration().at(2), 11.3461, 1e-4);
  EXPECT_NEAR(sphParticle1.getEngDot(), 5.46311, 1e-5);
  EXPECT_NEAR(sphParticle2.getEngDot(), 13.0197, 1e-4);

  sphParticle1.calcDt();
  sphParticle2.calcDt();
  EXPECT_NEAR(sphParticle1.getDt(), 0.0595165, 1e-7);
  EXPECT_NEAR(sphParticle2.getDt(), 0.110531, 1e-6);
}