/**
 * @file ATFunctorTestGlobals.cpp
 * @author muehlhaeusser
 * @date 29.08.23
 */

#include "ATFunctorTestGlobals.h"

TYPED_TEST_SUITE_P(ATFunctorTestGlobals);

template <class FuncType>
void ATFunctorTestGlobals<FuncType>::ATFunctorTestGlobalsNoMixing(ATFunctorTestGlobals<FuncType>::where_type where,
                                                                  bool newton3) {
  FuncType functor(cutoff);
  functor.setParticleProperties(nu);

  double whereFactor;
  std::string where_str;
  bool owned1, owned2, owned3;
  switch (where) {
    case inside:
      whereFactor = 1.;
      where_str = "inside";
      owned1 = owned2 = owned3 = true;
      break;
    case ininout:
      whereFactor = 2. / 3.;
      where_str = "ininout";
      owned1 = owned2 = true;
      owned3 = false;
      break;
    case inoutout:
      whereFactor = 1. / 3.;
      where_str = "inoutout";
      owned2 = true;
      owned1 = owned3 = false;
      break;
    case outside:
      whereFactor = 0.;
      where_str = "outside";
      owned1 = owned2 = owned3 = false;
      break;
    default:
      FAIL() << "not in enum where_type";
  }

  const std::array<double, 3> p1Pos{0., 0., 0.};
  const std::array<double, 3> p2Pos{0.1, 0., 0.};
  const std::array<double, 3> p3Pos{0., 0.2, 0.3};

  Molecule p1(p1Pos, {0., 0., 0.}, 0, 0);
  p1.setOwnershipState(owned1 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
  Molecule p2(p2Pos, {0., 0., 0.}, 1, 0);
  p2.setOwnershipState(owned2 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
  Molecule p3(p3Pos, {0., 0., 0.}, 2, 0);
  p3.setOwnershipState(owned3 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);

  functor.initTraversal();
  functor.AoSFunctor(p1, p2, p3, newton3);
  if (not newton3) {
    functor.AoSFunctor(p2, p1, p3, newton3);
    functor.AoSFunctor(p3, p1, p2, newton3);
  }
  functor.endTraversal(newton3);

  const double potentialEnergy = functor.getPotentialEnergy();
  const double virial = functor.getVirial();

  const double expectedEnergy = calculateATPotential(p1Pos, p2Pos, p3Pos, cutoff, nu);
  const auto [virial1, virial2, virial3] = calculateATVirialTotalPerParticle(p1Pos, p2Pos, p3Pos, cutoff, nu);
  const double expectedVirial = virial1 * owned1 + virial2 * owned2 + virial3 * owned3;

  EXPECT_NEAR(potentialEnergy, whereFactor * expectedEnergy, absDelta)
      << "where: " << where_str << ", newton3: " << newton3;
  EXPECT_NEAR(virial, expectedVirial, absDelta) << "where: " << where_str << ", newton3: " << newton3;
}

TYPED_TEST_P(ATFunctorTestGlobals, testAoSATFunctorGlobalsOpenMPParallel) {
  using FuncType = TypeParam;
  using TestType = ATFunctorTestGlobals<FuncType>;

  const bool newton3 = true;

  constexpr std::array<double, 3> p1Pos{0., 0., 0.};
  constexpr std::array<double, 3> p2Pos{0.1, 0., 0.};
  constexpr std::array<double, 3> p3Pos{0., 0.2, 0.3};
  constexpr std::array<double, 3> p4Pos{0., 2., 0.};
  constexpr std::array<double, 3> p5Pos{0.1, 2., 0.};
  constexpr std::array<double, 3> p6Pos{0., 2.2, 0.3};

  Molecule p1(p1Pos, {0., 0., 0.}, 0, 0);
  Molecule p2(p2Pos, {0., 0., 0.}, 1, 0);
  Molecule p3(p3Pos, {0., 0., 0.}, 2, 0);

  Molecule p4(p4Pos, {0., 0., 0.}, 0, 0);
  Molecule p5(p5Pos, {0., 0., 0.}, 1, 0);
  Molecule p6(p6Pos, {0., 0., 0.}, 2, 0);

  FuncType functor(this->cutoff);
  functor.setParticleProperties(this->nu);

  functor.initTraversal();

  std::string msg = "";
  // This is a basic check for the global calculations, by checking the handling of two particle interactions in
  // parallel. If interactions are dangerous, archer will complain.
#if defined(AUTOPAS_OPENMP)
// reduction for appending strings: "abc" + "def" -> "abcdef"
#pragma omp declare reduction(stringAppend : std::string : omp_out.append(omp_in))

#pragma omp parallel reduction(stringAppend : msg)
#endif
  {
#if defined(AUTOPAS_OPENMP)
#pragma omp sections
#endif
    {
#if defined(AUTOPAS_OPENMP)
#pragma omp section
#endif
      {
        msg += this->shouldSkipIfNotImplemented([&]() { functor.AoSFunctor(p1, p2, p3, newton3); });
      }  // pragma omp section
#if defined(AUTOPAS_OPENMP)
#pragma omp section
#endif
      {
        msg += this->shouldSkipIfNotImplemented([&]() { functor.AoSFunctor(p4, p5, p6, newton3); });
      }  // pragma omp section
    }    // pragma omp sections
  }      // pragma omp parallel

  if (not msg.empty()) {
    GTEST_SKIP() << msg;
  }

  functor.endTraversal(newton3);

  const double potentialEnergy = functor.getPotentialEnergy();
  const double virial = functor.getVirial();

  const double expectedEnergyTriplet1 = calculateATPotential(p1Pos, p2Pos, p3Pos, this->cutoff, this->nu);
  const double expectedVirialTriplet1 = calculateATVirialTotal(p1Pos, p2Pos, p3Pos, this->cutoff, this->nu);
  const double expectedEnergyTriplet2 = calculateATPotential(p4Pos, p5Pos, p6Pos, this->cutoff, this->nu);
  const double expectedVirialTriplet2 = calculateATVirialTotal(p4Pos, p5Pos, p6Pos, this->cutoff, this->nu);
  const double expectedEnergy = expectedEnergyTriplet1 + expectedEnergyTriplet2;
  const double expectedVirial = expectedVirialTriplet1 + expectedVirialTriplet2;

  EXPECT_NEAR(potentialEnergy, expectedEnergy, this->absDelta) << "newton3: " << newton3;
  EXPECT_NEAR(virial, expectedVirial, this->absDelta) << "newton3: " << newton3;
}

TYPED_TEST_P(ATFunctorTestGlobals, testATFunctorGlobalsThrowBad) {
  using exception_type = autopas::utils::ExceptionHandler::AutoPasException;

  using FuncType = TypeParam;
  FuncType functor(this->cutoff);

  // getPotentialEnergy without postprocessing is not allowed
  EXPECT_THROW(functor.getPotentialEnergy(), exception_type);
  EXPECT_THROW(functor.getVirial(), exception_type);

  EXPECT_NO_THROW(functor.initTraversal());

  EXPECT_NO_THROW(functor.endTraversal(true));
  EXPECT_NO_THROW(functor.initTraversal());
  EXPECT_NO_THROW(functor.endTraversal(true));
  // repeated postprocessing is not allowed
  EXPECT_THROW(functor.endTraversal(true), exception_type);

  EXPECT_NO_THROW(functor.initTraversal());
  EXPECT_NO_THROW(functor.endTraversal(true));
}

TYPED_TEST_P(ATFunctorTestGlobals, testAoSATFunctorGlobals) {
  using FuncType = TypeParam;
  using TestType = ATFunctorTestGlobals<FuncType>;

  for (typename TestType::where_type where : {TestType::where_type::inside, TestType::where_type::ininout,
                                              TestType::where_type::inoutout, TestType::where_type::outside}) {
    for (bool newton3 : {false, true}) {
      if (auto msg = this->shouldSkipIfNotImplemented([&]() { this->ATFunctorTestGlobalsNoMixing(where, newton3); });
          msg != "") {
        GTEST_SKIP() << msg;
      }
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P(ATFunctorTestGlobals, testAoSATFunctorGlobals, testATFunctorGlobalsThrowBad,
                            testAoSATFunctorGlobalsOpenMPParallel);

using MyTypes = ::testing::Types<ATFunNoMixGlob
#ifdef __AVX__

// TODO: Add AVX Functor
#endif
                                 >;

INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, ATFunctorTestGlobals, MyTypes);
