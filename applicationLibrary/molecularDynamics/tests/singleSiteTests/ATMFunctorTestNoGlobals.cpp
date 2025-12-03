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
    EXPECT_NEAR(f1one[0], this->expectedForceP1[0], this->absDelta);
    EXPECT_NEAR(f1one[1], this->expectedForceP1[1], this->absDelta);
    EXPECT_NEAR(f1one[2], this->expectedForceP1[2], this->absDelta);
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
      EXPECT_NEAR(f2one[0], this->expectedForceP2[0], this->absDelta);
      EXPECT_NEAR(f2one[1], this->expectedForceP2[1], this->absDelta);
      EXPECT_NEAR(f2one[2], this->expectedForceP2[2], this->absDelta);
      EXPECT_NEAR(f3one[0], this->expectedForceP3[0], this->absDelta);
      EXPECT_NEAR(f3one[1], this->expectedForceP3[1], this->absDelta);
      EXPECT_NEAR(f3one[2], this->expectedForceP3[2], this->absDelta);
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
    EXPECT_NEAR(f1two[0], factor * this->expectedForceP1[0], this->absDelta);
    EXPECT_NEAR(f1two[1], factor * this->expectedForceP1[1], this->absDelta);
    EXPECT_NEAR(f1two[2], factor * this->expectedForceP1[2], this->absDelta);

    EXPECT_NEAR(f2two[0], factor * this->expectedForceP2[0], this->absDelta);
    EXPECT_NEAR(f2two[1], factor * this->expectedForceP2[1], this->absDelta);
    EXPECT_NEAR(f2two[2], factor * this->expectedForceP2[2], this->absDelta);
  }
  if (newton3) {
    if (mixing) {
      EXPECT_NEAR(f3two[0], factor * this->expectedForceMixingP3[0], this->absDelta);
      EXPECT_NEAR(f3two[1], factor * this->expectedForceMixingP3[1], this->absDelta);
      EXPECT_NEAR(f3two[2], factor * this->expectedForceMixingP3[2], this->absDelta);
    } else {
      EXPECT_NEAR(f3two[0], factor * this->expectedForceP3[0], this->absDelta);
      EXPECT_NEAR(f3two[1], factor * this->expectedForceP3[1], this->absDelta);
      EXPECT_NEAR(f3two[2], factor * this->expectedForceP3[2], this->absDelta);
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
    EXPECT_NEAR(f1three[0], factor * this->expectedForceP1[0], this->absDelta);
    EXPECT_NEAR(f1three[1], factor * this->expectedForceP1[1], this->absDelta);
    EXPECT_NEAR(f1three[2], factor * this->expectedForceP1[2], this->absDelta);

    EXPECT_NEAR(f2three[0], factor * this->expectedForceP2[0], this->absDelta);
    EXPECT_NEAR(f2three[1], factor * this->expectedForceP2[1], this->absDelta);
    EXPECT_NEAR(f2three[2], factor * this->expectedForceP2[2], this->absDelta);

    EXPECT_NEAR(f3three[0], factor * this->expectedForceP3[0], this->absDelta);
    EXPECT_NEAR(f3three[1], factor * this->expectedForceP3[1], this->absDelta);
    EXPECT_NEAR(f3three[2], factor * this->expectedForceP3[2], this->absDelta);
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
    EXPECT_NEAR(f1one[0], this->expectedForceP1[0], this->absDelta);
    EXPECT_NEAR(f1one[1], this->expectedForceP1[1], this->absDelta);
    EXPECT_NEAR(f1one[2], this->expectedForceP1[2], this->absDelta);
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
      EXPECT_NEAR(f2one[0], this->expectedForceP2[0], this->absDelta);
      EXPECT_NEAR(f2one[1], this->expectedForceP2[1], this->absDelta);
      EXPECT_NEAR(f2one[2], this->expectedForceP2[2], this->absDelta);
      EXPECT_NEAR(f3one[0], this->expectedForceP3[0], this->absDelta);
      EXPECT_NEAR(f3one[1], this->expectedForceP3[1], this->absDelta);
      EXPECT_NEAR(f3one[2], this->expectedForceP3[2], this->absDelta);
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

  for (typename TestType::InteractionType interactionType :
       {TestType::InteractionType::triple, TestType::InteractionType::pair12, TestType::InteractionType::pair21,
        TestType::InteractionType::verlet, TestType::InteractionType::own}) {
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

    FMCell cell1, cell2, cell3;
    {
      // particle 1 is always in cell1
      Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0, 0);
      cell1.addParticle(p1);

      // The cell of particle 2 depends on the InteractionType.
      Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, (mixing) ? 1 : 0);
      switch (interactionType) {
        case TestType::InteractionType::verlet:
          // same as for own
        case TestType::InteractionType::pair21:
          // same as for own
        case TestType::InteractionType::own:
          // If we interact one cell with itself, it should be in cell1 as well.
          cell1.addParticle(p2);
          break;
        case TestType::InteractionType::pair12:
          // same as for cell triple
        case TestType::InteractionType::triple:
          // Second particle should be in cell2
          cell2.addParticle(p2);
          break;
        default:
          FAIL();
      }

      // The cell of particle 3 depends on the InteractionType.
      Molecule p3({0.3, 0.2, 0.1}, {0., 0., 0.}, 2, (mixing) ? 2 : 0);
      switch (interactionType) {
        case TestType::InteractionType::verlet:
          // same as for own
        case TestType::InteractionType::own:
          // If we interact one cell with itself, it should be in cell1 as well.
          cell1.addParticle(p3);
          break;
        case TestType::InteractionType::pair21:
          // same as for pair12
        case TestType::InteractionType::pair12:
          // third particle is always in cell2 for pair tests
          cell2.addParticle(p3);
          break;
        case TestType::InteractionType::triple:
          // third particle should be in cell3
          cell3.addParticle(p3);
          break;
        default:
          FAIL();
      }
    }
    // Load the particles into the soa.
    functor->SoALoader(cell1, cell1._particleSoABuffer, 0, /*skipSoAResize*/ false);
    functor->SoALoader(cell2, cell2._particleSoABuffer, 0, /*skipSoAResize*/ false);
    functor->SoALoader(cell3, cell3._particleSoABuffer, 0, /*skipSoAResize*/ false);

    if (auto msg = this->shouldSkipIfNotImplemented([&]() {
          switch (interactionType) {
            case TestType::InteractionType::own:
              // Interaction of one cell with itself
              functor->SoAFunctorSingle(cell1._particleSoABuffer, newton3);
              break;
            case TestType::InteractionType::pair12:
              // Interaction of a cell pair
              functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
              break;
            case TestType::InteractionType::pair21:
              // Interaction of a cell pair
              functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
              break;
            case TestType::InteractionType::triple:
              // Interaction of a cell pair
              functor->SoAFunctorTriple(cell1._particleSoABuffer, cell2._particleSoABuffer, cell3._particleSoABuffer,
                                        newton3);
              break;
            case TestType::InteractionType::verlet:
              // Build verlet list
              std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> neighborList(3);
              neighborList[0].push_back(1);
              neighborList[0].push_back(2);
              if (not newton3) {
                neighborList[1].push_back(0);
                neighborList[1].push_back(2);
                neighborList[2].push_back(0);
                neighborList[2].push_back(1);
              }
              functor->SoAFunctorVerlet(cell1._particleSoABuffer, 0, neighborList[0], newton3);
              functor->SoAFunctorVerlet(cell1._particleSoABuffer, 1, neighborList[1], newton3);
              functor->SoAFunctorVerlet(cell1._particleSoABuffer, 2, neighborList[2], newton3);
          }
        });
        msg != "") {
      GTEST_SKIP() << msg;
    }

    // Extract the particles from the soa
    functor->SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor->SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    functor->SoAExtractor(cell3, cell3._particleSoABuffer, 0);

    // force of particle 1
    auto f1 = cell1.begin()->getF();

    if (mixing) {
      EXPECT_NEAR(f1[0], this->expectedForceMixingP1[0], this->absDelta);
      EXPECT_NEAR(f1[1], this->expectedForceMixingP1[1], this->absDelta);
      EXPECT_NEAR(f1[2], this->expectedForceMixingP1[2], this->absDelta);
    } else {
      EXPECT_NEAR(f1[0], this->expectedForceP1[0], this->absDelta);
      EXPECT_NEAR(f1[1], this->expectedForceP1[1], this->absDelta);
      EXPECT_NEAR(f1[2], this->expectedForceP1[2], this->absDelta);
    }

    // force of particle 2
    std::array<double, 3> f2 = {0., 0., 0.};
    switch (interactionType) {
      case TestType::InteractionType::verlet:
      case TestType::InteractionType::pair21:
      case TestType::InteractionType::own:
        f2 = (++cell1.begin())->getF();
        break;
      case TestType::InteractionType::pair12:
      case TestType::InteractionType::triple:
        f2 = cell2.begin()->getF();
        break;
    }

    // if the interactionType is own, then the forces of the second particle should always be calculated!
    if (newton3 or interactionType == TestType::InteractionType::own or
        interactionType == TestType::InteractionType::pair21 or interactionType == TestType::InteractionType::verlet) {
      if (mixing) {
        EXPECT_NEAR(f2[0], this->expectedForceMixingP2[0], this->absDelta);
        EXPECT_NEAR(f2[1], this->expectedForceMixingP2[1], this->absDelta);
        EXPECT_NEAR(f2[2], this->expectedForceMixingP2[2], this->absDelta);
      } else {
        EXPECT_NEAR(f2[0], this->expectedForceP2[0], this->absDelta);
        EXPECT_NEAR(f2[1], this->expectedForceP2[1], this->absDelta);
        EXPECT_NEAR(f2[2], this->expectedForceP2[2], this->absDelta);
      }
    } else {
      EXPECT_DOUBLE_EQ(f2[0], 0);
      EXPECT_DOUBLE_EQ(f2[1], 0);
      EXPECT_DOUBLE_EQ(f2[2], 0);
    }

    // force of particle 3
    std::array<double, 3> f3 = {0., 0., 0.};
    switch (interactionType) {
      case TestType::InteractionType::verlet:
      case TestType::InteractionType::own: {
        auto iter = ++cell1.begin();
        f3 = (++iter)->getF();
        break;
      }
      case TestType::InteractionType::pair21:
        f3 = cell2.begin()->getF();
        break;
      case TestType::InteractionType::pair12:
        f3 = (++cell2.begin())->getF();
        break;
      case TestType::InteractionType::triple:
        f3 = cell3.begin()->getF();
        break;
    }

    // if the interactionType is own, then the forces of the third particle should always be calculated!
    if (newton3 or interactionType == TestType::InteractionType::own or
        interactionType == TestType::InteractionType::verlet) {
      if (mixing) {
        EXPECT_NEAR(f3[0], this->expectedForceMixingP3[0], this->absDelta);
        EXPECT_NEAR(f3[1], this->expectedForceMixingP3[1], this->absDelta);
        EXPECT_NEAR(f3[2], this->expectedForceMixingP3[2], this->absDelta);
      } else {
        EXPECT_NEAR(f3[0], this->expectedForceP3[0], this->absDelta);
        EXPECT_NEAR(f3[1], this->expectedForceP3[1], this->absDelta);
        EXPECT_NEAR(f3[2], this->expectedForceP3[2], this->absDelta);
      }
    } else {
      EXPECT_DOUBLE_EQ(f3[0], 0);
      EXPECT_DOUBLE_EQ(f3[1], 0);
      EXPECT_DOUBLE_EQ(f3[2], 0);
    }

    // Second cell interaction for types with multiple cells
    if (interactionType != TestType::InteractionType::own and interactionType != TestType::InteractionType::verlet) {
      functor->SoALoader(cell1, cell1._particleSoABuffer, 0, /*skipSoAResize*/ false);
      functor->SoALoader(cell2, cell2._particleSoABuffer, 0, /*skipSoAResize*/ false);
      functor->SoALoader(cell3, cell3._particleSoABuffer, 0, /*skipSoAResize*/ false);
      switch (interactionType) {
        case TestType::InteractionType::pair12:
          // Interaction of a cell pair
          functor->SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, newton3);
          break;
        case TestType::InteractionType::pair21:
          // Interaction of a cell pair
          functor->SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, newton3);
          break;
        case TestType::InteractionType::triple:
          // Interaction of a cell triplet
          functor->SoAFunctorTriple(cell2._particleSoABuffer, cell1._particleSoABuffer, cell3._particleSoABuffer,
                                    newton3);
          break;
        default:
          break;
      }
      functor->SoAExtractor(cell1, cell1._particleSoABuffer, 0);
      functor->SoAExtractor(cell2, cell2._particleSoABuffer, 0);
      functor->SoAExtractor(cell3, cell3._particleSoABuffer, 0);

      // Force on particles 1 and 2
      f1 = cell1.begin()->getF();
      std::array<double, 3> f2 = {0., 0., 0.};
      switch (interactionType) {
        case TestType::InteractionType::pair21:
          f2 = (++cell1.begin())->getF();
          break;
        case TestType::InteractionType::pair12:
        case TestType::InteractionType::triple:
          f2 = cell2.begin()->getF();
          break;
        default:
          break;
      }

      // Test the forces after the second interaction
      double factor = newton3 ? 2. : 1.;
      if (mixing) {
        EXPECT_NEAR(f1[0], factor * this->expectedForceMixingP1[0], this->absDelta);
        EXPECT_NEAR(f1[1], factor * this->expectedForceMixingP1[1], this->absDelta);
        EXPECT_NEAR(f1[2], factor * this->expectedForceMixingP1[2], this->absDelta);

        EXPECT_NEAR(f2[0], factor * this->expectedForceMixingP2[0], this->absDelta);
        EXPECT_NEAR(f2[1], factor * this->expectedForceMixingP2[1], this->absDelta);
        EXPECT_NEAR(f2[2], factor * this->expectedForceMixingP2[2], this->absDelta);
      } else {
        EXPECT_NEAR(f1[0], factor * this->expectedForceP1[0], this->absDelta);
        EXPECT_NEAR(f1[1], factor * this->expectedForceP1[1], this->absDelta);
        EXPECT_NEAR(f1[2], factor * this->expectedForceP1[2], this->absDelta);

        EXPECT_NEAR(f2[0], factor * this->expectedForceP2[0], this->absDelta);
        EXPECT_NEAR(f2[1], factor * this->expectedForceP2[1], this->absDelta);
        EXPECT_NEAR(f2[2], factor * this->expectedForceP2[2], this->absDelta);
      }

      // Force on particle 3
      std::array<double, 3> f3 = {0., 0., 0.};
      switch (interactionType) {
        case TestType::InteractionType::pair21:
          f3 = cell2.begin()->getF();
          break;
        case TestType::InteractionType::pair12:
          f3 = (++cell2.begin())->getF();
          break;
        case TestType::InteractionType::triple:
          f3 = cell3.begin()->getF();
          break;
        default:
          break;
      }

      if (newton3 or interactionType != TestType::InteractionType::triple) {
        if (mixing) {
          EXPECT_NEAR(f3[0], factor * this->expectedForceMixingP3[0], this->absDelta);
          EXPECT_NEAR(f3[1], factor * this->expectedForceMixingP3[1], this->absDelta);
          EXPECT_NEAR(f3[2], factor * this->expectedForceMixingP3[2], this->absDelta);
        } else {
          EXPECT_NEAR(f3[0], factor * this->expectedForceP3[0], this->absDelta);
          EXPECT_NEAR(f3[1], factor * this->expectedForceP3[1], this->absDelta);
          EXPECT_NEAR(f3[2], factor * this->expectedForceP3[2], this->absDelta);
        }
      } else {
        EXPECT_DOUBLE_EQ(f3[0], 0);
        EXPECT_DOUBLE_EQ(f3[1], 0);
        EXPECT_DOUBLE_EQ(f3[2], 0);
      }

      // Third interaction for 3 cells
      if (interactionType == TestType::InteractionType::triple) {
        functor->SoALoader(cell1, cell1._particleSoABuffer, 0, /*skipSoAResize*/ false);
        functor->SoALoader(cell2, cell2._particleSoABuffer, 0, /*skipSoAResize*/ false);
        functor->SoALoader(cell3, cell3._particleSoABuffer, 0, /*skipSoAResize*/ false);

        functor->SoAFunctorTriple(cell3._particleSoABuffer, cell1._particleSoABuffer, cell2._particleSoABuffer,
                                  newton3);

        functor->SoAExtractor(cell1, cell1._particleSoABuffer, 0);
        functor->SoAExtractor(cell2, cell2._particleSoABuffer, 0);
        functor->SoAExtractor(cell3, cell3._particleSoABuffer, 0);

        double factor = newton3 ? 3. : 1.;
        if (mixing) {
          EXPECT_NEAR(f1[0], factor * this->expectedForceMixingP1[0], this->absDelta);
          EXPECT_NEAR(f1[1], factor * this->expectedForceMixingP1[1], this->absDelta);
          EXPECT_NEAR(f1[2], factor * this->expectedForceMixingP1[2], this->absDelta);

          EXPECT_NEAR(f2[0], factor * this->expectedForceMixingP2[0], this->absDelta);
          EXPECT_NEAR(f2[1], factor * this->expectedForceMixingP2[1], this->absDelta);
          EXPECT_NEAR(f2[2], factor * this->expectedForceMixingP2[2], this->absDelta);

          EXPECT_NEAR(f3[0], factor * this->expectedForceMixingP3[0], this->absDelta);
          EXPECT_NEAR(f3[1], factor * this->expectedForceMixingP3[1], this->absDelta);
          EXPECT_NEAR(f3[2], factor * this->expectedForceMixingP3[2], this->absDelta);
        } else {
          EXPECT_NEAR(f1[0], factor * this->expectedForceP1[0], this->absDelta);
          EXPECT_NEAR(f1[1], factor * this->expectedForceP1[1], this->absDelta);
          EXPECT_NEAR(f1[2], factor * this->expectedForceP1[2], this->absDelta);

          EXPECT_NEAR(f2[0], factor * this->expectedForceP2[0], this->absDelta);
          EXPECT_NEAR(f2[1], factor * this->expectedForceP2[1], this->absDelta);
          EXPECT_NEAR(f2[2], factor * this->expectedForceP2[2], this->absDelta);

          EXPECT_NEAR(f3[0], factor * this->expectedForceP3[0], this->absDelta);
          EXPECT_NEAR(f3[1], factor * this->expectedForceP3[1], this->absDelta);
          EXPECT_NEAR(f3[2], factor * this->expectedForceP3[2], this->absDelta);
        }
      }
    }

    if (::testing::Test::HasFailure()) {
      std::cerr << "Failures for options: " << std::endl
                << "\tInteractionType: " << interactionType << std::endl
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

using MyTypes = ::testing::Types<Newton3True<ATMFunMixNoGlob>, Newton3False<ATMFunMixNoGlob>,
                                 Newton3True<ATMFunNoMixNoGlob>, Newton3False<ATMFunNoMixNoGlob>>;

INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, ATMFunctorTestNoGlobals, MyTypes);