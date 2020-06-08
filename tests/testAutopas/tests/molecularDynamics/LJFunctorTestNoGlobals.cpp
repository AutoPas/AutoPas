/**
 * @file LJFunctorTestNoGlobals.cpp
 * @author F. Gratl
 * @date 20.03.20
 */

#include "LJFunctorTestNoGlobals.h"

TYPED_TEST_SUITE_P(LJFunctorTestNoGlobals);

TYPED_TEST_P(LJFunctorTestNoGlobals, testAoSNoGlobals) {
  using FuncType = typename TypeParam::FuncType;
  constexpr bool mixing = FuncType::getMixing();
  constexpr bool newton3 = TypeParam::newton3;

  ParticlePropertiesLibrary<double, size_t> particlePropertiesLibrary(this->cutoff);
  std::unique_ptr<FuncType> functor;

  particlePropertiesLibrary.addType(0, this->epsilon, this->sigma, 1.0);
  if constexpr (mixing) {
    functor = std::make_unique<FuncType>(this->cutoff, particlePropertiesLibrary);
    particlePropertiesLibrary.addType(1, this->epsilon2, this->sigma2, 1.0);
  } else {
    functor = std::make_unique<FuncType>(this->cutoff);
    functor->setParticleProperties(this->epsilon * 24, 1);
  }

  Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, (mixing) ? 1 : 0);

  if (auto msg = this->shouldSkipIfNotImplemented([&]() { functor->AoSFunctor(p1, p2, newton3); }); msg != "") {
    GTEST_SKIP() << msg;
  }

  auto f1one = p1.getF();
  auto f2one = p2.getF();

  if (mixing) {
    EXPECT_NEAR(f1one[0], this->expectedForceMixing[0], this->absDelta);
    EXPECT_NEAR(f1one[1], this->expectedForceMixing[1], this->absDelta);
    EXPECT_NEAR(f1one[2], this->expectedForceMixing[2], this->absDelta);
  } else {
    EXPECT_NEAR(f1one[0], this->expectedForce[0], this->absDelta);
    EXPECT_NEAR(f1one[1], this->expectedForce[1], this->absDelta);
    EXPECT_NEAR(f1one[2], this->expectedForce[2], this->absDelta);
  }
  if (newton3) {
    if (mixing) {
      EXPECT_NEAR(f2one[0], -this->expectedForceMixing[0], this->absDelta);
      EXPECT_NEAR(f2one[1], -this->expectedForceMixing[1], this->absDelta);
      EXPECT_NEAR(f2one[2], -this->expectedForceMixing[2], this->absDelta);
    } else {
      EXPECT_NEAR(f2one[0], -this->expectedForce[0], this->absDelta);
      EXPECT_NEAR(f2one[1], -this->expectedForce[1], this->absDelta);
      EXPECT_NEAR(f2one[2], -this->expectedForce[2], this->absDelta);
    }
  } else {
    EXPECT_DOUBLE_EQ(f2one[0], 0);
    EXPECT_DOUBLE_EQ(f2one[1], 0);
    EXPECT_DOUBLE_EQ(f2one[2], 0);
  }

  functor->AoSFunctor(p2, p1, newton3);

  auto f1two = p1.getF();
  auto f2two = p2.getF();

  double factor = newton3 ? 2. : 1.;
  if (mixing) {
    EXPECT_NEAR(f1two[0], factor * this->expectedForceMixing[0], this->absDelta);
    EXPECT_NEAR(f1two[1], factor * this->expectedForceMixing[1], this->absDelta);
    EXPECT_NEAR(f1two[2], factor * this->expectedForceMixing[2], this->absDelta);

    EXPECT_NEAR(f2two[0], -factor * this->expectedForceMixing[0], this->absDelta);
    EXPECT_NEAR(f2two[1], -factor * this->expectedForceMixing[1], this->absDelta);
    EXPECT_NEAR(f2two[2], -factor * this->expectedForceMixing[2], this->absDelta);
  } else {
    EXPECT_NEAR(f1two[0], factor * this->expectedForce[0], this->absDelta);
    EXPECT_NEAR(f1two[1], factor * this->expectedForce[1], this->absDelta);
    EXPECT_NEAR(f1two[2], factor * this->expectedForce[2], this->absDelta);

    EXPECT_NEAR(f2two[0], -factor * this->expectedForce[0], this->absDelta);
    EXPECT_NEAR(f2two[1], -factor * this->expectedForce[1], this->absDelta);
    EXPECT_NEAR(f2two[2], -factor * this->expectedForce[2], this->absDelta);
  }
}

TYPED_TEST_P(LJFunctorTestNoGlobals, testSoANoGlobals) {
  using FuncType = typename TypeParam::FuncType;
  using TestType = LJFunctorTestNoGlobals<FuncType>;
  constexpr bool mixing = FuncType::getMixing();
  constexpr bool newton3 = TypeParam::newton3;

  for (typename TestType::InteractionType interactionType :
       {TestType::InteractionType::pair, TestType::InteractionType::verlet, TestType::InteractionType::own}) {
    ParticlePropertiesLibrary<double, size_t> particlePropertiesLibrary(this->cutoff);
    std::unique_ptr<FuncType> functor;

    if constexpr (mixing) {
      particlePropertiesLibrary.addType(0, this->epsilon, this->sigma, 1.0);
      particlePropertiesLibrary.addType(1, this->epsilon2, this->sigma2, 1.0);
      functor = std::make_unique<FuncType>(this->cutoff, particlePropertiesLibrary);
    } else {
      functor = std::make_unique<FuncType>(this->cutoff);
      functor->setParticleProperties(this->epsilon * 24, 1);
    }

    FMCell cell1, cell2;
    {
      // particle 1 is always in cell1
      Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0, 0);
      cell1.addParticle(p1);

      // The cell of particle 2 depends on the InteractionType.
      Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, (mixing) ? 1 : 0);
      switch (interactionType) {
        case TestType::InteractionType::verlet:
          // same as for own
        case TestType::InteractionType::own:
          // If we interact one cell with itself, it should be in cell1 as well.
          cell1.addParticle(p2);
          break;
        case TestType::InteractionType::pair:
          // If we interact a cell pair, it should be in cell2.
          cell2.addParticle(p2);
          break;
        default:
          FAIL();
      }
    }
    // Load the particles into the soa.
    functor->SoALoader(cell1, cell1._particleSoABuffer, 0);
    functor->SoALoader(cell2, cell2._particleSoABuffer, 0);

    if (auto msg = this->shouldSkipIfNotImplemented([&]() {
          switch (interactionType) {
            case TestType::InteractionType::own:
              // Interation of one cell with itself
              functor->SoAFunctorSingle(cell1._particleSoABuffer, newton3);
              break;
            case TestType::InteractionType::pair:
              // Interation of a cell pair
              functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
              break;
            case TestType::InteractionType::verlet:
              // Build verlet list
              std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> neighborList(2);
              neighborList[0].push_back(1);
              if (not newton3) {
                neighborList[1].push_back(0);
              }
              functor->SoAFunctorVerlet(cell1._particleSoABuffer, 0, neighborList[0], newton3);
              functor->SoAFunctorVerlet(cell1._particleSoABuffer, 1, neighborList[1], newton3);
          }
        });
        msg != "") {
      GTEST_SKIP() << msg;
    }

    // Extract the particles from the soa
    functor->SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor->SoAExtractor(cell2, cell2._particleSoABuffer, 0);

    // force of particle 1
    auto f1 = cell1.begin()->getF();

    if (mixing) {
      EXPECT_NEAR(f1[0], this->expectedForceMixing[0], this->absDelta);
      EXPECT_NEAR(f1[1], this->expectedForceMixing[1], this->absDelta);
      EXPECT_NEAR(f1[2], this->expectedForceMixing[2], this->absDelta);
    } else {
      EXPECT_NEAR(f1[0], this->expectedForce[0], this->absDelta);
      EXPECT_NEAR(f1[1], this->expectedForce[1], this->absDelta);
      EXPECT_NEAR(f1[2], this->expectedForce[2], this->absDelta);
    }

    // force of particle 2
    std::array<double, 3> f2 = {0., 0., 0.};
    switch (interactionType) {
      case TestType::InteractionType::verlet:
      case TestType::InteractionType::own:
        f2 = (++cell1.begin())->getF();
        break;
      case TestType::InteractionType::pair:
        f2 = cell2.begin()->getF();
        break;
    }
    // if the interactiontype is own, then the forces of the second particle should always be calculated!
    if (newton3 or interactionType != TestType::InteractionType::pair) {
      if (mixing) {
        EXPECT_NEAR(f2[0], -this->expectedForceMixing[0], this->absDelta);
        EXPECT_NEAR(f2[1], -this->expectedForceMixing[1], this->absDelta);
        EXPECT_NEAR(f2[2], -this->expectedForceMixing[2], this->absDelta);
      } else {
        EXPECT_NEAR(f2[0], -this->expectedForce[0], this->absDelta);
        EXPECT_NEAR(f2[1], -this->expectedForce[1], this->absDelta);
        EXPECT_NEAR(f2[2], -this->expectedForce[2], this->absDelta);
      }
    } else {
      EXPECT_DOUBLE_EQ(f2[0], 0);
      EXPECT_DOUBLE_EQ(f2[1], 0);
      EXPECT_DOUBLE_EQ(f2[2], 0);
    }

    if (interactionType == TestType::InteractionType::pair) {
      functor->SoALoader(cell1, cell1._particleSoABuffer, 0);
      functor->SoALoader(cell2, cell2._particleSoABuffer, 0);
      functor->SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, newton3);
      functor->SoAExtractor(cell1, cell1._particleSoABuffer, 0);
      functor->SoAExtractor(cell2, cell2._particleSoABuffer, 0);

      f1 = cell1.begin()->getF();
      f2 = cell2.begin()->getF();

      double factor = newton3 ? 2. : 1.;
      if (mixing) {
        EXPECT_NEAR(f1[0], factor * this->expectedForceMixing[0], this->absDelta);
        EXPECT_NEAR(f1[1], factor * this->expectedForceMixing[1], this->absDelta);
        EXPECT_NEAR(f1[2], factor * this->expectedForceMixing[2], this->absDelta);

        EXPECT_NEAR(f2[0], -factor * this->expectedForceMixing[0], this->absDelta);
        EXPECT_NEAR(f2[1], -factor * this->expectedForceMixing[1], this->absDelta);
        EXPECT_NEAR(f2[2], -factor * this->expectedForceMixing[2], this->absDelta);
      } else {
        EXPECT_NEAR(f1[0], factor * this->expectedForce[0], this->absDelta);
        EXPECT_NEAR(f1[1], factor * this->expectedForce[1], this->absDelta);
        EXPECT_NEAR(f1[2], factor * this->expectedForce[2], this->absDelta);

        EXPECT_NEAR(f2[0], -factor * this->expectedForce[0], this->absDelta);
        EXPECT_NEAR(f2[1], -factor * this->expectedForce[1], this->absDelta);
        EXPECT_NEAR(f2[2], -factor * this->expectedForce[2], this->absDelta);
      }
    }
    if (::testing::Test::HasFailure()) {
      std::cerr << "Failures for options: " << std::endl
                << "\tInteractionType: " << interactionType << std::endl
                << "\tnewton3: " << newton3 << std::endl;
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P(LJFunctorTestNoGlobals, testAoSNoGlobals, testSoANoGlobals);

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

using MyTypes = ::testing::Types<Newton3True<LJFunShiftMixNoGlob>, Newton3False<LJFunShiftMixNoGlob>,
                                 Newton3True<LJFunShiftNoMixNoGlob>, Newton3False<LJFunShiftNoMixNoGlob>
#ifdef __AVX__
                                 ,
                                 Newton3True<LJFunAVXShiftMixNoGlob>, Newton3False<LJFunAVXShiftMixNoGlob>,
                                 Newton3True<LJFunAVXShiftNoMixNoGlob>, Newton3False<LJFunAVXShiftNoMixNoGlob>
#endif
                                 >;
INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, LJFunctorTestNoGlobals, MyTypes);