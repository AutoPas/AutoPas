/**
 * @file ATMFunctorTestNoGlobals.cpp
 * @author muehlhaeusser
 * @date 26.09.23
 */

#include "ATMFunctorTestNoGlobals.h"

TYPED_TEST_SUITE_P(ATMFunctorTestNoGlobals);

TYPED_TEST_P(ATMFunctorTestNoGlobals, testAoSNoGlobalsATM) {
  using FuncType = typename TypeParam::FuncType;
  constexpr bool mixing = FuncType::getMixing();
  constexpr bool newton3 = TypeParam::newton3;

  ParticlePropertiesLibrary<double, size_t> particlePropertiesLibrary(this->cutoff);
  std::unique_ptr<FuncType> functor;

  particlePropertiesLibrary.addSiteType(0, 1.0);
  particlePropertiesLibrary.addATMParametersToSite(0, this->nu);
  if constexpr (mixing) {
    functor = std::make_unique<FuncType>(this->cutoff, particlePropertiesLibrary);
    particlePropertiesLibrary.addSiteType(1, 1.0);
    particlePropertiesLibrary.addATMParametersToSite(1, this->nu2);
    particlePropertiesLibrary.addSiteType(2, 1.0);
    particlePropertiesLibrary.addATMParametersToSite(2, this->nu3);

  } else {
    functor = std::make_unique<FuncType>(this->cutoff);
    functor->setParticleProperties(this->nu);
  }
  particlePropertiesLibrary.calculateMixingCoefficients();

  Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, (mixing) ? 1 : 0);
  Molecule p3({0.3, 0.2, 0.1}, {0., 0., 0.}, 2, (mixing) ? 2 : 0);

  if (auto msg = this->shouldSkipIfNotImplemented([&]() { functor->AoSFunctor(p1, p3, p2, newton3); }); msg != "") {
    GTEST_SKIP() << msg;
  }

  const auto f1one = p1.getF();
  const auto f2one = p2.getF();
  const auto f3one = p3.getF();

  if (mixing) {
    EXPECT_NEAR(f1one[0], this->expectedForceMixingP1[0], this->absDelta);
    EXPECT_NEAR(f1one[1], this->expectedForceMixingP1[1], this->absDelta);
    EXPECT_NEAR(f1one[2], this->expectedForceMixingP1[2], this->absDelta);
  } else {
    EXPECT_NEAR(f1one[0], this->expectedForceNonMixingP1[0], this->absDelta);
    EXPECT_NEAR(f1one[1], this->expectedForceNonMixingP1[1], this->absDelta);
    EXPECT_NEAR(f1one[2], this->expectedForceNonMixingP1[2], this->absDelta);
  }
  if (newton3) {
    if (mixing) {
      EXPECT_NEAR(f2one[0], this->expectedForceMixingP2[0], this->absDelta);
      EXPECT_NEAR(f2one[1], this->expectedForceMixingP2[1], this->absDelta);
      EXPECT_NEAR(f2one[2], this->expectedForceMixingP2[2], this->absDelta);
      EXPECT_NEAR(f3one[0], this->expectedForceMixingP3[0], this->absDelta);
      EXPECT_NEAR(f3one[1], this->expectedForceMixingP3[1], this->absDelta);
      EXPECT_NEAR(f3one[2], this->expectedForceMixingP3[2], this->absDelta);
    } else {
      EXPECT_NEAR(f2one[0], this->expectedForceNonMixingP2[0], this->absDelta);
      EXPECT_NEAR(f2one[1], this->expectedForceNonMixingP2[1], this->absDelta);
      EXPECT_NEAR(f2one[2], this->expectedForceNonMixingP2[2], this->absDelta);
      EXPECT_NEAR(f3one[0], this->expectedForceNonMixingP3[0], this->absDelta);
      EXPECT_NEAR(f3one[1], this->expectedForceNonMixingP3[1], this->absDelta);
      EXPECT_NEAR(f3one[2], this->expectedForceNonMixingP3[2], this->absDelta);
    }
  } else {
    EXPECT_DOUBLE_EQ(f2one[0], 0);
    EXPECT_DOUBLE_EQ(f2one[1], 0);
    EXPECT_DOUBLE_EQ(f2one[2], 0);
    EXPECT_DOUBLE_EQ(f3one[0], 0);
    EXPECT_DOUBLE_EQ(f3one[1], 0);
    EXPECT_DOUBLE_EQ(f3one[2], 0);
  }

  functor->AoSFunctor(p2, p1, p3, newton3);

  const auto f1two = p1.getF();
  const auto f2two = p2.getF();
  const auto f3two = p3.getF();

  double factor = newton3 ? 2. : 1.;
  if (mixing) {
    EXPECT_NEAR(f1two[0], factor * this->expectedForceMixingP1[0], this->absDelta);
    EXPECT_NEAR(f1two[1], factor * this->expectedForceMixingP1[1], this->absDelta);
    EXPECT_NEAR(f1two[2], factor * this->expectedForceMixingP1[2], this->absDelta);

    EXPECT_NEAR(f2two[0], factor * this->expectedForceMixingP2[0], this->absDelta);
    EXPECT_NEAR(f2two[1], factor * this->expectedForceMixingP2[1], this->absDelta);
    EXPECT_NEAR(f2two[2], factor * this->expectedForceMixingP2[2], this->absDelta);
  } else {
    EXPECT_NEAR(f1two[0], factor * this->expectedForceNonMixingP1[0], this->absDelta);
    EXPECT_NEAR(f1two[1], factor * this->expectedForceNonMixingP1[1], this->absDelta);
    EXPECT_NEAR(f1two[2], factor * this->expectedForceNonMixingP1[2], this->absDelta);

    EXPECT_NEAR(f2two[0], factor * this->expectedForceNonMixingP2[0], this->absDelta);
    EXPECT_NEAR(f2two[1], factor * this->expectedForceNonMixingP2[1], this->absDelta);
    EXPECT_NEAR(f2two[2], factor * this->expectedForceNonMixingP2[2], this->absDelta);
  }
  if (newton3) {
    if (mixing) {
      EXPECT_NEAR(f3two[0], factor * this->expectedForceMixingP3[0], this->absDelta);
      EXPECT_NEAR(f3two[1], factor * this->expectedForceMixingP3[1], this->absDelta);
      EXPECT_NEAR(f3two[2], factor * this->expectedForceMixingP3[2], this->absDelta);
    } else {
      EXPECT_NEAR(f3two[0], factor * this->expectedForceNonMixingP3[0], this->absDelta);
      EXPECT_NEAR(f3two[1], factor * this->expectedForceNonMixingP3[1], this->absDelta);
      EXPECT_NEAR(f3two[2], factor * this->expectedForceNonMixingP3[2], this->absDelta);
    }
  } else {
    EXPECT_DOUBLE_EQ(f3two[0], 0);
    EXPECT_DOUBLE_EQ(f3two[1], 0);
    EXPECT_DOUBLE_EQ(f3two[2], 0);
  }

  functor->AoSFunctor(p3, p1, p2, newton3);

  const auto f1three = p1.getF();
  const auto f2three = p2.getF();
  const auto f3three = p3.getF();

  factor = newton3 ? 3. : 1.;
  if (mixing) {
    EXPECT_NEAR(f1three[0], factor * this->expectedForceMixingP1[0], this->absDelta);
    EXPECT_NEAR(f1three[1], factor * this->expectedForceMixingP1[1], this->absDelta);
    EXPECT_NEAR(f1three[2], factor * this->expectedForceMixingP1[2], this->absDelta);

    EXPECT_NEAR(f2three[0], factor * this->expectedForceMixingP2[0], this->absDelta);
    EXPECT_NEAR(f2three[1], factor * this->expectedForceMixingP2[1], this->absDelta);
    EXPECT_NEAR(f2three[2], factor * this->expectedForceMixingP2[2], this->absDelta);

    EXPECT_NEAR(f3three[0], factor * this->expectedForceMixingP3[0], this->absDelta);
    EXPECT_NEAR(f3three[1], factor * this->expectedForceMixingP3[1], this->absDelta);
    EXPECT_NEAR(f3three[2], factor * this->expectedForceMixingP3[2], this->absDelta);
  } else {
    EXPECT_NEAR(f1three[0], factor * this->expectedForceNonMixingP1[0], this->absDelta);
    EXPECT_NEAR(f1three[1], factor * this->expectedForceNonMixingP1[1], this->absDelta);
    EXPECT_NEAR(f1three[2], factor * this->expectedForceNonMixingP1[2], this->absDelta);

    EXPECT_NEAR(f2three[0], factor * this->expectedForceNonMixingP2[0], this->absDelta);
    EXPECT_NEAR(f2three[1], factor * this->expectedForceNonMixingP2[1], this->absDelta);
    EXPECT_NEAR(f2three[2], factor * this->expectedForceNonMixingP2[2], this->absDelta);

    EXPECT_NEAR(f3three[0], factor * this->expectedForceNonMixingP3[0], this->absDelta);
    EXPECT_NEAR(f3three[1], factor * this->expectedForceNonMixingP3[1], this->absDelta);
    EXPECT_NEAR(f3three[2], factor * this->expectedForceNonMixingP3[2], this->absDelta);
  }

  // order of second and third particle should not matter
  p1.setF({0., 0., 0.});
  p2.setF({0., 0., 0.});
  p3.setF({0., 0., 0.});
  functor->AoSFunctor(p1, p3, p2, newton3);

  const auto f1four = p1.getF();
  const auto f2four = p2.getF();
  const auto f3four = p3.getF();

  if (mixing) {
    EXPECT_NEAR(f1one[0], this->expectedForceMixingP1[0], this->absDelta);
    EXPECT_NEAR(f1one[1], this->expectedForceMixingP1[1], this->absDelta);
    EXPECT_NEAR(f1one[2], this->expectedForceMixingP1[2], this->absDelta);
  } else {
    EXPECT_NEAR(f1one[0], this->expectedForceNonMixingP1[0], this->absDelta);
    EXPECT_NEAR(f1one[1], this->expectedForceNonMixingP1[1], this->absDelta);
    EXPECT_NEAR(f1one[2], this->expectedForceNonMixingP1[2], this->absDelta);
  }
  if (newton3) {
    if (mixing) {
      EXPECT_NEAR(f2one[0], this->expectedForceMixingP2[0], this->absDelta);
      EXPECT_NEAR(f2one[1], this->expectedForceMixingP2[1], this->absDelta);
      EXPECT_NEAR(f2one[2], this->expectedForceMixingP2[2], this->absDelta);
      EXPECT_NEAR(f3one[0], this->expectedForceMixingP3[0], this->absDelta);
      EXPECT_NEAR(f3one[1], this->expectedForceMixingP3[1], this->absDelta);
      EXPECT_NEAR(f3one[2], this->expectedForceMixingP3[2], this->absDelta);
    } else {
      EXPECT_NEAR(f2one[0], this->expectedForceNonMixingP2[0], this->absDelta);
      EXPECT_NEAR(f2one[1], this->expectedForceNonMixingP2[1], this->absDelta);
      EXPECT_NEAR(f2one[2], this->expectedForceNonMixingP2[2], this->absDelta);
      EXPECT_NEAR(f3one[0], this->expectedForceNonMixingP3[0], this->absDelta);
      EXPECT_NEAR(f3one[1], this->expectedForceNonMixingP3[1], this->absDelta);
      EXPECT_NEAR(f3one[2], this->expectedForceNonMixingP3[2], this->absDelta);
    }
  } else {
    EXPECT_DOUBLE_EQ(f2one[0], 0);
    EXPECT_DOUBLE_EQ(f2one[1], 0);
    EXPECT_DOUBLE_EQ(f2one[2], 0);
    EXPECT_DOUBLE_EQ(f3one[0], 0);
    EXPECT_DOUBLE_EQ(f3one[1], 0);
    EXPECT_DOUBLE_EQ(f3one[2], 0);
  }
}

TYPED_TEST_P(ATMFunctorTestNoGlobals, testSoANoGlobalsATM) {
  using namespace autopas::utils::ArrayMath::literals;
  using FuncType = typename TypeParam::FuncType;
  using TestType = ATMFunctorTestNoGlobals<FuncType>;
  constexpr bool mixing = FuncType::getMixing();
  constexpr bool newton3 = TypeParam::newton3;

  // Helper lambdas to check the forces
  auto testNonZeroForce = [&](const std::array<double, 3> &particleForce, const std::array<double, 3> &expectedForce,
                              const std::string &type) {
    for (size_t i = 0; i < 3; i++) {
      EXPECT_NEAR(particleForce[i], expectedForce[i], this->absDelta)
          << "Failed for SoAFunctor: " << type << ", Newton3: " << newton3;
    }
  };
  auto testZeroForce = [&](const std::array<double, 3> &particleForce, const std::string &type) {
    for (size_t i = 0; i < 3; i++) {
      EXPECT_DOUBLE_EQ(particleForce[i], 0) << "Failed for SoAFunctor: " << type << ", Newton3: " << newton3;
    }
  };

  for (typename TestType::SoAFunctorType soaFunctorType :
       {TestType::SoAFunctorType::triple, TestType::SoAFunctorType::pair12, TestType::SoAFunctorType::pair21,
        TestType::SoAFunctorType::single, TestType::SoAFunctorType::verlet}) {
    ParticlePropertiesLibrary<double, size_t> particlePropertiesLibrary(this->cutoff);
    std::unique_ptr<FuncType> functor;

    if constexpr (mixing) {
      particlePropertiesLibrary.addSiteType(0, 1.0);
      particlePropertiesLibrary.addATMParametersToSite(0, this->nu);
      particlePropertiesLibrary.addSiteType(1, 1.0);
      particlePropertiesLibrary.addATMParametersToSite(1, this->nu2);
      particlePropertiesLibrary.addSiteType(2, 1.0);
      particlePropertiesLibrary.addATMParametersToSite(2, this->nu3);
      particlePropertiesLibrary.calculateMixingCoefficients();
      functor = std::make_unique<FuncType>(this->cutoff, particlePropertiesLibrary);
    } else {
      functor = std::make_unique<FuncType>(this->cutoff);
      functor->setParticleProperties(this->nu);
    }

    const auto expectedForceP1 = mixing ? this->expectedForceMixingP1 : this->expectedForceNonMixingP1;
    const auto expectedForceP2 = mixing ? this->expectedForceMixingP2 : this->expectedForceNonMixingP2;
    const auto expectedForceP3 = mixing ? this->expectedForceMixingP3 : this->expectedForceNonMixingP3;

    // Set up the particles
    // This array stores pointers to the 3 particles that might reside in the same or different cells
    std::array<Molecule *, 3> moleculePtrs{};

    FMCell cell1, cell2, cell3;
    // Add particle 1 to cell1
    Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0, 0);
    cell1.addParticle(p1);

    // Distribute particles 2 and 3 depending on the SoA functor to test
    // Gather pointers to the 3 molecules for easier testing
    Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, (mixing) ? 1 : 0);
    Molecule p3({0.3, 0.2, 0.1}, {0., 0., 0.}, 2, (mixing) ? 2 : 0);
    switch (soaFunctorType) {
      case TestType::SoAFunctorType::single:
        cell1.addParticle(p2);
        cell1.addParticle(p3);
        moleculePtrs[1] = &cell1._particles[1];
        moleculePtrs[2] = &cell1._particles[2];
        break;
      case TestType::SoAFunctorType::pair21:
        cell1.addParticle(p2);
        cell2.addParticle(p3);
        moleculePtrs[1] = &cell1._particles[1];
        moleculePtrs[2] = &cell2._particles[0];
        break;
      case TestType::SoAFunctorType::pair12:
        cell2.addParticle(p2);
        cell2.addParticle(p3);
        moleculePtrs[1] = &cell2._particles[0];
        moleculePtrs[2] = &cell2._particles[1];
        break;
      case TestType::SoAFunctorType::triple:
        cell2.addParticle(p2);
        cell3.addParticle(p3);
        moleculePtrs[1] = &cell2._particles[0];
        moleculePtrs[2] = &cell3._particles[0];
        break;
      case TestType::SoAFunctorType::verlet:
        cell1.addParticle(p2);
        cell1.addParticle(p3);
        moleculePtrs[1] = &cell1._particles[1];
        moleculePtrs[2] = &cell1._particles[2];
        break;
      default:
        FAIL();
    }
    moleculePtrs[0] = &cell1._particles[0];

    // Create helper lambdas to load/extract particles to/from the SoAs
    auto cells = std::array{&cell1, &cell2, &cell3};
    auto loadParticlesToSoA = [&]() {
      for (auto &cell : cells) {
        functor->SoALoader(*cell, cell->_particleSoABuffer, 0, /*skipSoAResize*/ false);
      }
    };
    auto extractParticlesFromSoA = [&]() {
      for (auto &cell : cells) {
        functor->SoAExtractor(*cell, cell->_particleSoABuffer, 0);
      }
    };

    // Lambda helper to call the corresponding SoA functor and perform the force interactions
    auto computeFunctorInteractions = [&](auto &c1, auto &c2, auto &c3) {
      switch (soaFunctorType) {
        case TestType::SoAFunctorType::single:
          functor->SoAFunctorSingle(c1._particleSoABuffer, newton3);
          break;
        case TestType::SoAFunctorType::pair12:
        case TestType::SoAFunctorType::pair21:
          functor->SoAFunctorPair(c1._particleSoABuffer, c2._particleSoABuffer, newton3);
          break;
        case TestType::SoAFunctorType::triple:
          functor->SoAFunctorTriple(c1._particleSoABuffer, c2._particleSoABuffer, c3._particleSoABuffer, newton3);
          break;
        case TestType::SoAFunctorType::verlet:
          // Build verlet list
          std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> neighborList(3);
          neighborList[0] = {1, 2};
          if (not newton3) {
            neighborList[1] = {0, 2};
            neighborList[2] = {0, 1};
          }
          functor->SoAFunctorVerlet(c1._particleSoABuffer, 0, neighborList[0], newton3);
          functor->SoAFunctorVerlet(c1._particleSoABuffer, 1, neighborList[1], newton3);
          functor->SoAFunctorVerlet(c1._particleSoABuffer, 2, neighborList[2], newton3);
      }
    };

    ///// First Iteration /////

    loadParticlesToSoA();
    // Perform the interaction based on the soa functor type - Skip if not yet implemented
    if (auto msg = this->shouldSkipIfNotImplemented([&]() { computeFunctorInteractions(cell1, cell2, cell3); });
        msg != "") {
      std::cout << msg << std::endl;
      continue;
    }
    extractParticlesFromSoA();

    // Test the forces
    // Particle 1
    testNonZeroForce(moleculePtrs[0]->getF(), expectedForceP1, ATMFunctorTest::to_string(soaFunctorType));

    // Particle 2
    if (newton3 or soaFunctorType == TestType::SoAFunctorType::single or
        soaFunctorType == TestType::SoAFunctorType::pair21 or soaFunctorType == TestType::SoAFunctorType::verlet) {
      // If particle 2 is in the same cell as particle 1, its forces should already be calculated
      testNonZeroForce(moleculePtrs[1]->getF(), expectedForceP2, ATMFunctorTest::to_string(soaFunctorType));
    } else {  // Particle 2 is in a different cell, so the forces are zero when newton3 == false
      testZeroForce(moleculePtrs[1]->getF(), ATMFunctorTest::to_string(soaFunctorType));
    }

    // Particle 3
    if (newton3 or soaFunctorType == TestType::SoAFunctorType::single or
        soaFunctorType == TestType::SoAFunctorType::verlet) {
      testNonZeroForce(moleculePtrs[2]->getF(), expectedForceP3, ATMFunctorTest::to_string(soaFunctorType));
    } else {
      testZeroForce(moleculePtrs[2]->getF(), ATMFunctorTest::to_string(soaFunctorType));
    }

    ///// Second Iteration /////

    // Only done for types with multiple cells
    if (soaFunctorType != TestType::SoAFunctorType::single and soaFunctorType != TestType::SoAFunctorType::verlet) {
      loadParticlesToSoA();
      computeFunctorInteractions(cell2, cell1, cell3);
      extractParticlesFromSoA();

      // Test the forces after the second interaction
      double factor = newton3 ? 2. : 1.;
      testNonZeroForce(moleculePtrs[0]->getF(), expectedForceP1 * factor, ATMFunctorTest::to_string(soaFunctorType));
      testNonZeroForce(moleculePtrs[1]->getF(), expectedForceP2 * factor, ATMFunctorTest::to_string(soaFunctorType));

      if (newton3 or soaFunctorType != TestType::SoAFunctorType::triple) {
        testNonZeroForce(moleculePtrs[2]->getF(), expectedForceP3 * factor, ATMFunctorTest::to_string(soaFunctorType));
      } else {
        testZeroForce(moleculePtrs[2]->getF(), ATMFunctorTest::to_string(soaFunctorType));
      }

      ///// Third Iteration /////

      // Only needed for SoAFunctorTriple
      if (soaFunctorType == TestType::SoAFunctorType::triple) {
        loadParticlesToSoA();
        computeFunctorInteractions(cell3, cell1, cell2);
        extractParticlesFromSoA();

        factor = newton3 ? 3. : 1.;
        testNonZeroForce(moleculePtrs[0]->getF(), expectedForceP1 * factor, ATMFunctorTest::to_string(soaFunctorType));
        testNonZeroForce(moleculePtrs[1]->getF(), expectedForceP2 * factor, ATMFunctorTest::to_string(soaFunctorType));
        testNonZeroForce(moleculePtrs[2]->getF(), expectedForceP3 * factor, ATMFunctorTest::to_string(soaFunctorType));
      }
    }

    if (::testing::Test::HasFailure()) {
      std::cerr << "Failures for options: " << std::endl
                << "\tSoAFunctorType: " << soaFunctorType << std::endl
                << "\tnewton3: " << newton3 << std::endl;
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P(ATMFunctorTestNoGlobals, testAoSNoGlobalsATM, testSoANoGlobalsATM);

template <class Func, bool n3>
struct TypeWrapper {
  using FuncType = Func;
  constexpr static bool newton3 = n3;
};

// struct aliasing for readable names
template <class FuncType>
struct Newton3True : public TypeWrapper<FuncType, true> {};
template <class FuncType>
struct Newton3False : public TypeWrapper<FuncType, false> {};

using MyTypes =
    ::testing::Types<Newton3True<ATMFunMixNoGlob>, Newton3False<ATMFunMixNoGlob>, Newton3True<ATMFunNoMixNoGlob>,
                     Newton3False<ATMFunNoMixNoGlob>, Newton3True<ATMFunHWYMixNoGlob>, Newton3False<ATMFunHWYMixNoGlob>,
                     Newton3True<ATMFunHWYNoMixNoGlob>, Newton3False<ATMFunHWYNoMixNoGlob>>;

template <typename TypeList>
struct ATMFunctorTypeNames {
  template <typename T>
  static std::string GetName(int) {
    using F = typename T::FuncType;

    const std::string n3 = T::newton3 ? "Newton3" : "NoNewton3";

    if constexpr (std::is_same_v<F, ATMFunHWYMixNoGlob>) {
      return "HWY_Mix_" + n3;
    } else if constexpr (std::is_same_v<F, ATMFunHWYNoMixNoGlob>) {
      return "HWY_NoMix_" + n3;
    } else if constexpr (std::is_same_v<F, ATMFunMixNoGlob>) {
      return "Mix_" + n3;
    } else if constexpr (std::is_same_v<F, ATMFunNoMixNoGlob>) {
      return "NoMix_" + n3;
    }
  }
};

INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, ATMFunctorTestNoGlobals, MyTypes, ATMFunctorTypeNames<MyTypes>);