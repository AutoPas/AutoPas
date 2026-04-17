/**
 * @file ParticlePropertiesLibraryTest.cpp
 * @author N. Fottner & S. Newcome
 * @date 02/08/19
 */
#include "ParticlePropertiesLibraryTest.h"

#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * Initializes a ParticleProperties Library, adds two sites, and tests that the getters for site values return correct
 * site information. Then tests that an error is thrown if a site with a non-consecutive site Id is added.
 */
TEST_F(ParticlePropertiesLibraryTest, SitePropertiesAddingAndGettingTest) {
  const double cutoff = 0.1;
  std::shared_ptr<ParticlePropertiesLibrary<double, unsigned int>> PPL =
      std::make_shared<ParticlePropertiesLibrary<double, unsigned int>>(cutoff);

  // Check that PPL is empty
  EXPECT_EQ(PPL->getNumberRegisteredSiteTypes(), 0);

  // Add site 0
  const double epsilon0 = 1.;
  const double sigma0 = 1.;
  const double nu0 = 0.1;
  const double mass0 = 1.;
  PPL->addSiteType(0, mass0);
  PPL->addLJParametersToSite(0, epsilon0, sigma0);
  PPL->addATMParametersToSite(0, nu0);

  // Check successfully getting of information
  EXPECT_EQ(PPL->getNumberRegisteredSiteTypes(), 1);
  EXPECT_EQ(PPL->getEpsilon(0), epsilon0);
  EXPECT_EQ(PPL->getSigma(0), sigma0);
  EXPECT_EQ(PPL->getNu(0), nu0);
  EXPECT_EQ(PPL->getSiteMass(0), mass0);

  // Add site 1
  const double epsilon1 = 0.2;
  const double sigma1 = 0.7;
  const double nu1 = 0.2;
  const double mass1 = 1.2;
  PPL->addSiteType(1, mass1);
  PPL->addLJParametersToSite(1, epsilon1, sigma1);
  PPL->addATMParametersToSite(1, nu1);

  // Check successfully getting of information
  EXPECT_EQ(PPL->getNumberRegisteredSiteTypes(), 2);
  EXPECT_EQ(PPL->getEpsilon(0), epsilon0);
  EXPECT_EQ(PPL->getSigma(0), sigma0);
  EXPECT_EQ(PPL->getNu(0), nu0);
  EXPECT_EQ(PPL->getSiteMass(0), mass0);
  EXPECT_EQ(PPL->getEpsilon(1), epsilon1);
  EXPECT_EQ(PPL->getSigma(1), sigma1);
  EXPECT_EQ(PPL->getNu(1), nu1);
  EXPECT_EQ(PPL->getSiteMass(1), mass1);

  // Check addSiteType with an inappropriate siteId throws an error.
  EXPECT_ANY_THROW(PPL->addSiteType(1, 1.));
  EXPECT_ANY_THROW(PPL->addSiteType(5, 1.));
}

/**
 * Adds a LJ site and a AT site to a PPL, and tests that the getters for site values return correct
 * site information being 0 for not explicitly set parameters. Then tests that an error is thrown if LJ or AT
 * parameters are added for a not yet initialized site.
 */
TEST_F(ParticlePropertiesLibraryTest, AddingDifferentSitesTest) {
  const double cutoff = 0.1;
  std::shared_ptr<ParticlePropertiesLibrary<double, unsigned int>> PPL =
      std::make_shared<ParticlePropertiesLibrary<double, unsigned int>>(cutoff);

  // Add a site 0 as a pure LJ site
  const double epsilon0 = 1.;
  const double sigma0 = 1.;
  const double mass0 = 1.;
  PPL->addSiteType(0, mass0);
  PPL->addLJParametersToSite(0, epsilon0, sigma0);

  // Add site 1 as a pure AT site
  const double nu1 = 0.1;
  const double mass1 = 1.2;
  PPL->addSiteType(1, mass1);
  PPL->addATMParametersToSite(1, nu1);

  // Check successfully getting of information
  // Uninitialized LJ or AT parameters should be 0
  EXPECT_EQ(PPL->getNumberRegisteredSiteTypes(), 2);
  EXPECT_EQ(PPL->getEpsilon(0), epsilon0);
  EXPECT_EQ(PPL->getSigma(0), sigma0);
  EXPECT_EQ(PPL->getNu(0), 0.0);
  EXPECT_EQ(PPL->getSiteMass(0), mass0);
  EXPECT_EQ(PPL->getEpsilon(1), 0.0);
  EXPECT_EQ(PPL->getSigma(1), 0.0);
  EXPECT_EQ(PPL->getNu(1), nu1);
  EXPECT_EQ(PPL->getSiteMass(1), mass1);

  // Calculate Mixing Data
  PPL->calculateMixingCoefficients();

  // Check mixing data with zero-initialized sites
  EXPECT_EQ(PPL->getMixing24Epsilon(0, 1), 0.0);
  EXPECT_EQ(PPL->getMixingSigmaSquared(1, 0), 0.25);
  EXPECT_EQ(PPL->getMixingNu(0, 0, 1), 0.0);
  EXPECT_EQ(PPL->getMixingNu(1, 1, 0), 0.0);

  // Check that add<LJ/AT>Site for a non-initialized siteId throws an error.
  EXPECT_ANY_THROW(PPL->addLJParametersToSite(2, 1., 1.));
  EXPECT_ANY_THROW(PPL->addATMParametersToSite(3, 1.));
}

/**
 * Tests adding multi-site molecule types and getting this information.
 *
 * Initializes a ParticleProperties Library and adds two sites. Then adds two molecule types, testing that the getters
 * return correct information.
 *
 * When md-flexible is compiled without support for multi-site molecules, this test checks that attempting to add
 * multi-site molecule types results in an exception being thrown.
 *
 * @note This molecular information is not intended to be mathematically sound.
 */
TEST_F(ParticlePropertiesLibraryTest, MolPropertiesAddingAndGettingTest) {
  const double cutoff = 0.1;
  std::shared_ptr<ParticlePropertiesLibrary<double, unsigned int>> PPL =
      std::make_shared<ParticlePropertiesLibrary<double, unsigned int>>(cutoff);

#if MD_FLEXIBLE_MODE == MULTISITE
  // add two site types
  PPL->addSiteType(0, 1.);
  PPL->addSiteType(1, 1.2);
  PPL->addLJParametersToSite(0, 1., 1.);
  PPL->addLJParametersToSite(1, 0.2, 0.7);

  // Check that PPL is empty of molecule types
  EXPECT_EQ(PPL->getNumberRegisteredMolTypes(), 0);

  // Add Molecule Type 0
  const std::vector<unsigned int> siteIds0 = {0};
  const std::vector<std::array<double, 3>> sitePositions0 = {{0., 0., 0.}};
  const std::array<double, 3> MoI0 = {1., 1., 1.};
  PPL->addMolType(0, siteIds0, sitePositions0, MoI0);

  // Check getters
  EXPECT_EQ(PPL->getNumberRegisteredMolTypes(), 1);
  EXPECT_THAT(PPL->getSiteTypes(0), ::testing::ElementsAreArray(siteIds0));
  EXPECT_THAT(PPL->getSitePositions(0), ::testing::ElementsAreArray(sitePositions0));
  EXPECT_THAT(PPL->getMomentOfInertia(0), ::testing::ElementsAreArray(MoI0));
  EXPECT_EQ(PPL->getNumSites(0), 1);

  // Add Molecule Type 1
  const std::vector<unsigned int> siteIds1 = {0, 1, 1};
  const std::vector<std::array<double, 3>> sitePositions1 = {{1., 0., 0.}, {-0.5, 0., 0.}, {0.5, 0., 0.}};
  const std::array<double, 3> MoI1 = {1., -1., 0.5};
  PPL->addMolType(1, siteIds1, sitePositions1, MoI1);

  // Check getters
  EXPECT_EQ(PPL->getNumberRegisteredMolTypes(), 2);
  EXPECT_THAT(PPL->getSiteTypes(0), ::testing::ElementsAreArray(siteIds0));
  EXPECT_THAT(PPL->getSitePositions(0), ::testing::ElementsAreArray(sitePositions0));
  EXPECT_THAT(PPL->getMomentOfInertia(0), ::testing::ElementsAreArray(MoI0));
  EXPECT_EQ(PPL->getNumSites(0), 1);
  EXPECT_THAT(PPL->getSiteTypes(1), ::testing::ElementsAreArray(siteIds1));
  EXPECT_THAT(PPL->getSitePositions(1), ::testing::ElementsAreArray(sitePositions1));
  EXPECT_THAT(PPL->getMomentOfInertia(1), ::testing::ElementsAreArray(MoI1));
  EXPECT_EQ(PPL->getNumSites(1), 3);

  // Try adding molecules with inappropriate IDs.
  EXPECT_ANY_THROW(PPL->addMolType(1, {0}, {{0., 0., 0.}}, {1., 1., 1.}););
  EXPECT_ANY_THROW(PPL->addMolType(5, {0}, {{0., 0., 0.}}, {1., 1., 1.}););

  // Try adding molecules with non-matching sizes of site type Ids and site position vectors
  EXPECT_ANY_THROW(PPL->addMolType(2, {0}, {{0., 0., 0.}, {0., 0., 0.}}, {1., 1., 1.}););
  EXPECT_ANY_THROW(PPL->addMolType(2, {0, 0}, {{0., 0., 0.}}, {1., 1., 1.}););

  // Try adding molecules with non-existant site Ids
  EXPECT_ANY_THROW(PPL->addMolType(2, {2}, {{0., 0., 0.}}, {1., 1., 1.}););
#else
  // Add Molecule Type
  EXPECT_ANY_THROW(PPL->addMolType(0, {0}, {{0., 0., 0.}}, {1., 1., 1.}));
#endif
}

/**
 * Simple unit test for PPL's calculation of shift6. Compares that the calculation is the exact same as the calculation
 * in this unit test.
 */
TEST_F(ParticlePropertiesLibraryTest, LennardJonesCalculateShift6Simple) {
  const double cutoff = 1.1;
  const double sigma = 0.5;
  const double epsilon = 1.3;
  std::shared_ptr<ParticlePropertiesLibrary<double, unsigned int>> PPL =
      std::make_shared<ParticlePropertiesLibrary<double, unsigned int>>(cutoff);

  const auto cutoffSquared = cutoff * cutoff;
  const auto sigmaSquared = sigma * sigma;
  const auto epsilon24 = 24. * epsilon;

  // Calculate expected shift6
  const auto sigmaDivCutoffPow2 = sigmaSquared / cutoffSquared;
  const auto sigmaDivCutoffPow6 = sigmaDivCutoffPow2 * sigmaDivCutoffPow2 * sigmaDivCutoffPow2;
  const auto expectedShift6 = epsilon24 * (sigmaDivCutoffPow6 - sigmaDivCutoffPow6 * sigmaDivCutoffPow6);

  // Calculate shift6
  const double shift6 = PPL->calcShift6(epsilon24, sigmaSquared, cutoffSquared);

  // Compare
  EXPECT_DOUBLE_EQ(expectedShift6, shift6);
}

/**
 * Tests correctness of calculated shift6 by testing that the potential energy at the cutoff, plus the shift, is zero.
 */
TEST_F(ParticlePropertiesLibraryTest, LennardJonesTestShiftGivesCorrectEnergyAtCutoff) {
  const double cutoff = 1.1;
  const double sigma = 0.5;
  const double epsilon = 1.3;
  std::shared_ptr<ParticlePropertiesLibrary<double, size_t>> PPL =
      std::make_shared<ParticlePropertiesLibrary<double, size_t>>(cutoff);

  PPL->addSiteType(0, 1.);
  PPL->addLJParametersToSite(0, epsilon, sigma);
  PPL->calculateMixingCoefficients();

  const auto cutoffSquared = cutoff * cutoff;
  const auto sigmaSquared = sigma * sigma;
  const auto epsilon24 = 24. * epsilon;

  // Calculate shift6
  const double shift6 = PPL->calcShift6(epsilon24, sigmaSquared, cutoffSquared);
  const auto shift = shift6 / 6.;

  // Create two LJ Molecules that are cutoff apart
  mdLib::MoleculeLJ molA({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  mdLib::MoleculeLJ molB({cutoff, 0., 0.}, {0., 0., 0.}, 1, 0);

  // Create LJ Functor class and use it to calculate the potential energy between the two
  mdLib::LJFunctor<mdLib::MoleculeLJ, /* shifting */ false, /*mixing*/ true, autopas::FunctorN3Modes::Both,
                   /*globals*/ true>
      ljFunctor(cutoff, *PPL);

  ljFunctor.initTraversal();
  ljFunctor.AoSFunctor(molA, molB, true);
  ljFunctor.endTraversal(true);

  EXPECT_DOUBLE_EQ(ljFunctor.getPotentialEnergy() + shift, 0.);
}

/**
 * Tests that getting the Lennard-Jones mixing data works correctly and that the mixing rules applied in PPL are
 * correct. In addition, tests that all of the mixing data getters are consistent with each other.
 */
TEST_F(ParticlePropertiesLibraryTest, LennardJonesMixingTest) {
  const double cutoff = 1.1;
  const double cutoffSquared = cutoff * cutoff;
  std::shared_ptr<ParticlePropertiesLibrary<double, unsigned int>> PPL =
      std::make_shared<ParticlePropertiesLibrary<double, unsigned int>>(cutoff);

  // Add three sites
  const double epsilon0 = 0.6;
  const double sigma0 = 1.2;
  const double epsilon1 = 0.7;
  const double sigma1 = 1.4;
  const double epsilon2 = 1.;
  const double sigma2 = 1.;
  PPL->addSiteType(0, 1.);
  PPL->addLJParametersToSite(0, epsilon0, sigma0);
  PPL->addSiteType(1, 1.);
  PPL->addLJParametersToSite(1, epsilon1, sigma1);
  PPL->addSiteType(2, 1.);
  PPL->addLJParametersToSite(2, epsilon2, sigma2);

  // Calculate mixing coefficients
  PPL->calculateMixingCoefficients();

  // Calculate expected values. Mixing rules are commutative so only half are calculated.
  const auto epsilon01 = std::sqrt(epsilon0 * epsilon1);
  const auto epsilon02 = std::sqrt(epsilon0 * epsilon2);
  const auto epsilon12 = std::sqrt(epsilon1 * epsilon2);

  const auto sigma01 = (sigma0 + sigma1) / 2.;
  const auto sigma02 = (sigma0 + sigma2) / 2.;
  const auto sigma12 = (sigma1 + sigma2) / 2.;

  const auto shift6_00 = PPL->calcShift6(24. * epsilon0, sigma0 * sigma0, cutoffSquared);
  const auto shift6_01 = PPL->calcShift6(24. * epsilon01, sigma01 * sigma01, cutoffSquared);
  const auto shift6_02 = PPL->calcShift6(24. * epsilon02, sigma02 * sigma02, cutoffSquared);
  const auto shift6_11 = PPL->calcShift6(24. * epsilon1, sigma1 * sigma1, cutoffSquared);
  const auto shift6_12 = PPL->calcShift6(24. * epsilon12, sigma12 * sigma12, cutoffSquared);
  const auto shift6_22 = PPL->calcShift6(24. * epsilon2, sigma2 * sigma2, cutoffSquared);

  // Compare PPL's calculated mixing coefficients against the expected values
  EXPECT_DOUBLE_EQ(PPL->getMixing24Epsilon(0, 0), 24. * epsilon0);
  EXPECT_DOUBLE_EQ(PPL->getMixing24Epsilon(0, 1), 24. * epsilon01);
  EXPECT_DOUBLE_EQ(PPL->getMixing24Epsilon(0, 2), 24. * epsilon02);
  EXPECT_DOUBLE_EQ(PPL->getMixing24Epsilon(1, 0), 24. * epsilon01);
  EXPECT_DOUBLE_EQ(PPL->getMixing24Epsilon(1, 1), 24. * epsilon1);
  EXPECT_DOUBLE_EQ(PPL->getMixing24Epsilon(1, 2), 24. * epsilon12);
  EXPECT_DOUBLE_EQ(PPL->getMixing24Epsilon(2, 0), 24. * epsilon02);
  EXPECT_DOUBLE_EQ(PPL->getMixing24Epsilon(2, 1), 24. * epsilon12);
  EXPECT_DOUBLE_EQ(PPL->getMixing24Epsilon(2, 2), 24. * epsilon2);

  EXPECT_DOUBLE_EQ(PPL->getMixingSigmaSquared(0, 0), sigma0 * sigma0);
  EXPECT_DOUBLE_EQ(PPL->getMixingSigmaSquared(0, 1), sigma01 * sigma01);
  EXPECT_DOUBLE_EQ(PPL->getMixingSigmaSquared(0, 2), sigma02 * sigma02);
  EXPECT_DOUBLE_EQ(PPL->getMixingSigmaSquared(1, 0), sigma01 * sigma01);
  EXPECT_DOUBLE_EQ(PPL->getMixingSigmaSquared(1, 1), sigma1 * sigma1);
  EXPECT_DOUBLE_EQ(PPL->getMixingSigmaSquared(1, 2), sigma12 * sigma12);
  EXPECT_DOUBLE_EQ(PPL->getMixingSigmaSquared(2, 0), sigma02 * sigma02);
  EXPECT_DOUBLE_EQ(PPL->getMixingSigmaSquared(2, 1), sigma12 * sigma12);
  EXPECT_DOUBLE_EQ(PPL->getMixingSigmaSquared(2, 2), sigma2 * sigma2);

  EXPECT_DOUBLE_EQ(PPL->getMixingShift6(0, 0), shift6_00);
  EXPECT_DOUBLE_EQ(PPL->getMixingShift6(0, 1), shift6_01);
  EXPECT_DOUBLE_EQ(PPL->getMixingShift6(0, 2), shift6_02);
  EXPECT_DOUBLE_EQ(PPL->getMixingShift6(1, 0), shift6_01);
  EXPECT_DOUBLE_EQ(PPL->getMixingShift6(1, 1), shift6_11);
  EXPECT_DOUBLE_EQ(PPL->getMixingShift6(1, 2), shift6_12);
  EXPECT_DOUBLE_EQ(PPL->getMixingShift6(2, 0), shift6_02);
  EXPECT_DOUBLE_EQ(PPL->getMixingShift6(2, 1), shift6_12);
  EXPECT_DOUBLE_EQ(PPL->getMixingShift6(2, 2), shift6_22);

  // Confirm that PPL's individual get sime mixing data functions match PPL's get all mixing data for a site-type pair.
  for (unsigned int i = 0; i < PPL->getNumberRegisteredSiteTypes(); i++) {
    for (unsigned int j = 0; j < PPL->getNumberRegisteredSiteTypes(); j++) {
      const auto allMixingData = PPL->getLJMixingData(i, j);
      const auto epsilon24 = PPL->getMixing24Epsilon(i, j);
      const auto sigmaSquared = PPL->getMixingSigmaSquared(i, j);
      const auto shift6 = PPL->getMixingShift6(i, j);

      EXPECT_DOUBLE_EQ(allMixingData.epsilon24, epsilon24);
      EXPECT_DOUBLE_EQ(allMixingData.sigmaSquared, sigmaSquared);
      EXPECT_DOUBLE_EQ(allMixingData.shift6, shift6);
    }
  }
}

/**
 * Tests that getting the Axilrod-Teller-Muto mixing data works correctly and that the mixing rules applied in PPL are
 * correct. In addition, tests that the mixing data getters are consistent with each other.
 */
TEST_F(ParticlePropertiesLibraryTest, AxilrodTellerMixingTest) {
  const double cutoff = 1.1;
  const double cutoffSquared = cutoff * cutoff;
  std::shared_ptr<ParticlePropertiesLibrary<double, unsigned int>> PPL =
      std::make_shared<ParticlePropertiesLibrary<double, unsigned int>>(cutoff);

  // Add three sites
  const double nu0 = 0.1;
  const double nu1 = 0.7;
  const double nu2 = 0.3;

  PPL->addSiteType(0, 1.);
  PPL->addATMParametersToSite(0, nu0);
  PPL->addSiteType(1, 1.);
  PPL->addATMParametersToSite(1, nu1);
  PPL->addSiteType(2, 1.);
  PPL->addATMParametersToSite(2, nu2);

  // Calculate mixing coefficients
  PPL->calculateMixingCoefficients();

  // Calculate expected values. Mixing rules are commutative so only unique values are calculated.
  const auto nu001 = std::cbrt(nu0 * nu0 * nu1);
  const auto nu011 = std::cbrt(nu0 * nu1 * nu1);
  const auto nu112 = std::cbrt(nu1 * nu1 * nu2);
  const auto nu122 = std::cbrt(nu1 * nu2 * nu2);
  const auto nu012 = std::cbrt(nu0 * nu1 * nu2);

  // Compare PPL's calculated mixing coefficients against the expected values
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(0, 0, 0), nu0);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(0, 0, 1), nu001);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(0, 1, 0), nu001);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(1, 0, 0), nu001);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(0, 1, 1), nu011);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(1, 0, 1), nu011);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(1, 1, 0), nu011);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(1, 1, 1), nu1);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(1, 1, 2), nu112);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(1, 2, 1), nu112);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(2, 1, 1), nu112);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(1, 2, 2), nu122);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(2, 1, 2), nu122);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(2, 2, 1), nu122);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(2, 2, 2), nu2);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(0, 1, 2), nu012);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(0, 2, 1), nu012);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(1, 0, 2), nu012);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(1, 2, 0), nu012);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(2, 0, 1), nu012);
  EXPECT_DOUBLE_EQ(PPL->getMixingNu(2, 1, 0), nu012);

  // Confirm that PPL's individual get same mixing data functions match PPL's get all mixing data for a site-type pair.
  for (unsigned int i = 0; i < PPL->getNumberRegisteredSiteTypes(); i++) {
    for (unsigned int j = 0; j < PPL->getNumberRegisteredSiteTypes(); j++) {
      for (unsigned int k = 0; k < PPL->getNumberRegisteredSiteTypes(); k++) {
        const auto allMixingData = PPL->getATMMixingData(i, j, k);
        const auto nu = PPL->getMixingNu(i, j, k);

        EXPECT_DOUBLE_EQ(allMixingData.nu, nu);
      }
    }
  }
}