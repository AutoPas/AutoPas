/**
 * @file ParticlePropertiesLibraryTest.cpp
 * @author N. Fottner
 * @date 02/08/19
 */
#include "ParticlePropertiesLibraryTest.h"

#include "autopas/molecularDynamics/LJFunctor.h"
#include "src/configuration/MDFlexConfig.h"
#include "testingHelpers/commonTypedefs.h"

//double ParticlePropertiesLibraryTest::mixingE(double e1, double e2) { return std::sqrt(e1 * e2); }
//double ParticlePropertiesLibraryTest::mixingS(double s1, double s2) { return ((s1 + s2) / 2); }

TEST_F(ParticlePropertiesLibraryTest, massTest) {
//  //ASSERT_EQ(PPL.getTypes.size(), masses.size()); // todo this
//  EXPECT_EQ(PPL.getSiteMass(0), masses[0]);
//  EXPECT_EQ(PPL.getSiteMass(1), masses[1]);
}

TEST_F(ParticlePropertiesLibraryTest, epsilonTest) {
//  //ASSERT_EQ(PPL.getTypes().size(), epsilons.size()); // todo this
//
//  EXPECT_EQ(PPL.get24Epsilon(0), epsilons[0] * 24);
//  EXPECT_EQ(PPL.get24Epsilon(1), epsilons[1] * 24);
//
//  for (unsigned int i = 0; i < epsilons.size(); ++i) {
//    for (unsigned int j = 0; j < epsilons.size(); ++j) {
//      auto expectedVal = mixingE(epsilons[i], epsilons[j]) * 24;
//      EXPECT_EQ(PPL.mixing24Epsilon(i, j), expectedVal) << "For i=" << i << " j=" << j;
//    }
//  }
}

TEST_F(ParticlePropertiesLibraryTest, sigmaTest) {
//  //ASSERT_EQ(PPL.getTypes().size(), sigmas.size()); // todo this
//
//  EXPECT_EQ(PPL.getSigmaSquare(0), sigmas[0] * sigmas[0]);
//  EXPECT_EQ(PPL.getSigmaSquare(1), sigmas[1] * sigmas[1]);
//
//  for (unsigned int i = 0; i < sigmas.size(); ++i) {
//    for (unsigned int j = 0; j < sigmas.size(); ++j) {
//      auto expectedVal = mixingS(sigmas[i], sigmas[j]);
//      expectedVal *= expectedVal;
//      EXPECT_EQ(PPL.mixingSigmaSquare(i, j), expectedVal) << "For i=" << i << " j=" << j;
//    }
//  }
}

TEST_F(ParticlePropertiesLibraryTest, shiftTest) {
//  //ASSERT_EQ(PPL.getTypes().size(), shifts.size()); // todo this
//
//  EXPECT_EQ(PPL.mixingShift6(0, 0), shifts[0]);
//  EXPECT_EQ(PPL.mixingShift6(1, 1), shifts[1]);
}

/**
 * Idea: Two particles with distance of (almost) cutoff should produce (almost) zero shifted potential.
 */
TEST_F(ParticlePropertiesLibraryTest, mixedShiftTestUpot) {
//  Molecule m1({0, 0, 0}, {0, 0, 0}, 0, 0);
//  Molecule m2({cutoff - 1e-14, 0, 0}, {0, 0, 0}, 1, 1);
//
//  autopas::LJFunctor<Molecule, /* shifting */ true, /*mixing*/ true, autopas::FunctorN3Modes::Both,
//                     /*globals*/ true>
//      functor(cutoff, PPL);
//
//  functor.initTraversal();
//  functor.AoSFunctor(m1, m2, true);
//  functor.endTraversal(true);
//  EXPECT_NE(functor.getUpot(), 0);
//  EXPECT_NEAR(functor.getUpot(), 0, 1e-10);
}

// todo: move this to yaml test
TEST_F(ParticlePropertiesLibraryTest, ParticlePropertiesInitialization) {
//  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
//                                        std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml"};
//
//  char *argv[3] = {arguments[0].data(), arguments[1].data(), arguments[2].data()};
//
//  MDFlexConfig configuration(3, argv);
//
//  EXPECT_EQ(configuration.getParticlePropertiesLibrary()->getSiteMass(0), 1.0);
//  EXPECT_EQ(configuration.getParticlePropertiesLibrary()->get24Epsilon(0), 24.0);
//  EXPECT_EQ(configuration.getParticlePropertiesLibrary()->getSigmaSquare(0), 1.0);
//  EXPECT_EQ(configuration.getParticlePropertiesLibrary()->getSiteMass(1), 2.);
//  EXPECT_EQ(configuration.getParticlePropertiesLibrary()->get24Epsilon(1), 48.0);
//  EXPECT_EQ(configuration.getParticlePropertiesLibrary()->getSigmaSquare(1), 4.0);
//  EXPECT_EQ(configuration.getParticlePropertiesLibrary()->getSiteMass(2), 3.0);
//  EXPECT_EQ(configuration.getParticlePropertiesLibrary()->get24Epsilon(2), 72.0);
//  EXPECT_EQ(configuration.getParticlePropertiesLibrary()->getSigmaSquare(2), 9.0);
//  EXPECT_EQ(configuration.getParticlePropertiesLibrary()->getSiteMass(3), 4.0);
//  EXPECT_EQ(configuration.getParticlePropertiesLibrary()->get24Epsilon(3), 96.0);
//  EXPECT_EQ(configuration.getParticlePropertiesLibrary()->getSigmaSquare(3), 16.0);
}

/**
 * Initializes a ParticleProperties Library, adds two sites, and tests that the getters for site values return correct
 * site information. Then tests that an error is thrown if a site with a non-consecutive site Id is added.
 */
TEST_F(ParticlePropertiesLibraryTest, SitePropertiesSettingAndGettingTest) {
  const double cutoff = 0.1;
  std::shared_ptr<ParticlePropertiesLibrary<double, unsigned int>> PPL =
      std::make_shared<ParticlePropertiesLibrary<double, unsigned int>>(cutoff);

  // Check that PPL is empty
  EXPECT_EQ(PPL->getNumberRegisteredSiteTypes(), 0);

  // Add site 0
  const double epsilon0 = 1.;
  const double sigma0 = 1.;
  const double mass0 = 1.;
  PPL->addSiteType(0, epsilon0, sigma0, mass0);

  // Check successfully getting of information
  EXPECT_EQ(PPL->getNumberRegisteredMolTypes(), 1);
  EXPECT_EQ(PPL->get24Epsilon(0), 24.*epsilon0);
  EXPECT_EQ(PPL->getSigmaSquared(0), sigma0*sigma0);
  EXPECT_EQ(PPL->getSiteMass(0), mass0);

  // Add site 1
  const double epsilon1 = 0.2;
  const double sigma1 = 0.7;
  const double mass1 = 1.2;
  PPL->addSiteType(1, epsilon1, sigma1, mass1);

  // Check successfully getting of information
  EXPECT_EQ(PPL->getNumberRegisteredMolTypes(), 2);
  EXPECT_EQ(PPL->get24Epsilon(0), 24.*epsilon0);
  EXPECT_EQ(PPL->getSigmaSquared(0), sigma0*sigma0);
  EXPECT_EQ(PPL->getSiteMass(0), mass0);
  EXPECT_EQ(PPL->get24Epsilon(1), 24.*epsilon1);
  EXPECT_EQ(PPL->getSigmaSquared(1), sigma1*sigma1);
  EXPECT_EQ(PPL->getSiteMass(1), mass1);

  // Check addSiteType with an inappropriate siteId throws an error.
  EXPECT_ANY_THROW(PPL->addSiteType(1, 1., 1., 1.));
  EXPECT_ANY_THROW(PPL->addSiteType(5, 1., 1., 1.));
}