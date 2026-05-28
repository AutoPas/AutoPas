/**
 * @file ATMFunctorTestGlobals.cpp
 * @author muehlhaeusser
 * @date 29.08.23
 */

#include "ATMFunctorTestGlobals.h"

#include "testingHelpers/ATMPotential.h"

TYPED_TEST_SUITE_P(ATMFunctorTestGlobals);

/*
 * Tests potential energy and virial calculations the AoSFunctor
 * @tparam FuncType The type of the functor to test
 * @param where The triplet configuration which is either all halo, 1 halo and 2 owned, 2 halo and 1 owned or 3 owned
 * @param newton3 Whether to use newton 3
 */
template <class FuncType>
void ATMFunctorTestGlobals<FuncType>::ATMFunctorTestGlobalsNoMixingAoS(where_type where, bool newton3) {
  FuncType functor(cutoff);
  functor.setParticleProperties(nu);

  // Map `where` to configuration values. The tuple consists of a factor for the expected pot. energy and booleans for
  // the owned status of the 3 particles.
  const std::map<where_type, std::tuple<double, bool, bool, bool>> whereConfig = {
      {allInside, {1., true, true, true}},
      {ininout, {2. / 3., true, true, false}},
      {inoutout, {1. / 3., false, true, false}},
      {allOutside, {0., false, false, false}}};
  const auto &[whereFactor, owned1, owned2, owned3] = whereConfig.at(where);

  constexpr std::array<double, 3> p1Pos{0., 0., 0.};
  constexpr std::array<double, 3> p2Pos{0.1, 0., 0.};
  constexpr std::array<double, 3> p3Pos{0., 0.2, 0.3};

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

  constexpr double expectedEnergy = calculateATMPotential(p1Pos, p2Pos, p3Pos, cutoff, nu);
  const auto [virial1, virial2, virial3] = calculateATMVirialTotalPerParticle(p1Pos, p2Pos, p3Pos, cutoff, nu);
  const double expectedVirial = virial1 * owned1 + virial2 * owned2 + virial3 * owned3;

  EXPECT_NEAR(potentialEnergy, whereFactor * expectedEnergy, absDelta)
      << "where: " << to_string(where) << ", newton3: " << newton3;
  EXPECT_NEAR(virial, expectedVirial, absDelta) << "where: " << to_string(where) << ", newton3: " << newton3;
}

/*
 * Tests potential energy and virial calculations the AoSFunctor
 * @tparam FuncType The type of the functor to test
 * @param where The triplet configuration which is either all halo, 1 halo and 2 owned, 2 halo and 1 owned or 3 owned
 * @param soaFunctorType Test for SoAFunctorSingle, SoAFunctorPair12, SoAFunctorPair21 or SoAFunctorTriple
 * @param newton3 Whether to use newton 3
 */
template <class FuncType>
void ATMFunctorTestGlobals<FuncType>::ATMFunctorTestSoAGlobals(where_type where, SoAFunctorType soaFunctorType,
                                                               bool newton3) {
  FuncType functor(cutoff);
  functor.setParticleProperties(nu);

  // Map `where` to configuration values
  const std::map<where_type, std::tuple<double, bool, bool, bool>> whereConfig = {
      {allInside, {1., true, true, true}},
      {ininout, {2. / 3., true, true, false}},
      {inoutout, {1. / 3., false, true, false}},
      {allOutside, {0., false, false, false}}};
  const auto &[whereFactor, owned1, owned2, owned3] = whereConfig.at(where);

  constexpr std::array<double, 3> p1Pos{0., 0., 0.};
  constexpr std::array<double, 3> p2Pos{0.1, 0., 0.};
  constexpr std::array<double, 3> p3Pos{0., 0.2, 0.3};
  Molecule p1(p1Pos, {0., 0., 0.}, 0, 0);
  p1.setOwnershipState(owned1 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
  Molecule p2(p2Pos, {0., 0., 0.}, 1, 0);
  p2.setOwnershipState(owned2 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
  Molecule p3(p3Pos, {0., 0., 0.}, 2, 0);
  p3.setOwnershipState(owned3 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);

  // Construct the cells
  FMCell cell1, cell2, cell3;
  std::vector<FMCell *> cells;

  // Distribute the particles according to the SoAFunctorType to test
  switch (soaFunctorType) {
    case single:
      cell1.addParticle(p1);
      cell1.addParticle(p2);
      cell1.addParticle(p3);
      cells = {&cell1};
      break;
    case pair12:
      cell1.addParticle(p1);
      cell2.addParticle(p2);
      cell2.addParticle(p3);
      cells = {&cell1, &cell2};
      break;
    case pair21:
      cell1.addParticle(p1);
      cell1.addParticle(p2);
      cell2.addParticle(p3);
      cells = {&cell1, &cell2};
      break;
    case triple:
      cell1.addParticle(p1);
      cell2.addParticle(p2);
      cell3.addParticle(p3);
      cells = {&cell1, &cell2, &cell3};
      break;
    default:
      FAIL() << "Not a valid SoAFunctorType";
  }

  // Fill SoA for all involved cells
  for (auto *c : cells) {
    functor.SoALoader(*c, c->_particleSoABuffer, 0, false);
  }
  functor.initTraversal();

  // Invoke the appropriate SoA functor based on the SoAFunctorType and newton3
  switch (soaFunctorType) {
    case single:
      functor.SoAFunctorSingle(cells[0]->_particleSoABuffer, newton3);
      break;
    case pair12:
    case pair21:
      functor.SoAFunctorPair(cells[0]->_particleSoABuffer, cells[1]->_particleSoABuffer, newton3);
      if (not newton3) {
        functor.SoAFunctorPair(cells[1]->_particleSoABuffer, cells[0]->_particleSoABuffer, newton3);
      }
      break;
    case triple:
      functor.SoAFunctorTriple(cells[0]->_particleSoABuffer, cells[1]->_particleSoABuffer, cells[2]->_particleSoABuffer,
                               newton3);
      if (not newton3) {
        functor.SoAFunctorTriple(cells[1]->_particleSoABuffer, cells[0]->_particleSoABuffer,
                                 cells[2]->_particleSoABuffer, newton3);
        functor.SoAFunctorTriple(cells[2]->_particleSoABuffer, cells[0]->_particleSoABuffer,
                                 cells[1]->_particleSoABuffer, newton3);
      }
      break;
    default:
      FAIL() << "Not a valid SoAFunctorType";
  }

  functor.endTraversal(newton3);

  const double potentialEnergy = functor.getPotentialEnergy();
  const double virial = functor.getVirial();

  constexpr double expectedEnergy = calculateATMPotential(p1Pos, p2Pos, p3Pos, cutoff, nu);
  const auto [virial1, virial2, virial3] = calculateATMVirialTotalPerParticle(p1Pos, p2Pos, p3Pos, cutoff, nu);
  const double expectedVirial = virial1 * owned1 + virial2 * owned2 + virial3 * owned3;

  EXPECT_NEAR(potentialEnergy, whereFactor * expectedEnergy, absDelta)
      << "SoAFunctorType: " << to_string(soaFunctorType) << ", particles are: " << to_string(where)
      << ", newton3: " << newton3;
  EXPECT_NEAR(virial, expectedVirial, absDelta) << "SoAFunctorType: " << to_string(soaFunctorType)
                                                << ", particles are: " << to_string(where) << ", newton3: " << newton3;
}

/*
 * Tests potential energy and virial calculations for a triplet split over periodic boundary conditions.
 * This especially interesting for the virial
 * @tparam FuncType The type of the functor to test
 * @param functorToTest Test either AoS, SoAPair or SoATriple
 * @param newton3 Whether to use newton 3
 */
template <class FuncType>
void ATMFunctorTestGlobals<FuncType>::ATMFunctorTestGlobalsPeriodicBCs(const std::string &functorToTest, bool newton3) {
  using namespace autopas::utils::ArrayMath::literals;
  FuncType functor(5.0);
  functor.setParticleProperties(nu);

  // Imaginary box is from {0, 0, 0} to {10, 10, 10}
  // p1 is close to the left border, p2 and p3 close to the right border
  // They will interact as a triplet with their corresponding halo particles
  constexpr std::array<double, 3> p1OwnedPos{1., 4., 5.};
  constexpr std::array<double, 3> p2OwnedPos{9., 5., 4.};
  constexpr std::array<double, 3> p3OwnedPos{9., 4., 4.};

  constexpr std::array<double, 3> p1HaloPos{11., 4., 5.};
  constexpr std::array<double, 3> p2HaloPos{-1., 5., 4.};
  constexpr std::array<double, 3> p3HaloPos{-1., 4., 4.};

  Molecule p1(p1OwnedPos, {0., 0., 0.}, 0, 0);
  p1.setOwnershipState(autopas::OwnershipState::owned);
  Molecule p2(p2OwnedPos, {0., 0., 0.}, 1, 0);
  p2.setOwnershipState(autopas::OwnershipState::owned);
  Molecule p3(p3OwnedPos, {0., 0., 0.}, 2, 0);
  p3.setOwnershipState(autopas::OwnershipState::owned);

  Molecule p1Halo(p1HaloPos, {0., 0., 0.}, 3, 0);
  p1Halo.setOwnershipState(autopas::OwnershipState::halo);
  Molecule p2Halo(p2HaloPos, {0., 0., 0.}, 4, 0);
  p2Halo.setOwnershipState(autopas::OwnershipState::halo);
  Molecule p3Halo(p3HaloPos, {0., 0., 0.}, 5, 0);
  p3Halo.setOwnershipState(autopas::OwnershipState::halo);

  if (functorToTest == "AoS") {
    functor.initTraversal();
    functor.AoSFunctor(p1, p2Halo, p3Halo, newton3);
    functor.AoSFunctor(p1Halo, p2, p3, newton3);
    if (not newton3) {
      functor.AoSFunctor(p2, p3, p1Halo, newton3);
      functor.AoSFunctor(p3, p2, p1Halo, newton3);
    }
    functor.endTraversal(newton3);
  } else if (functorToTest == "SoAPair") {
    FMCell ownedCell{}, haloCellLeft{}, haloCellRight{};
    haloCellLeft.addParticle(p2Halo);
    haloCellLeft.addParticle(p3Halo);
    ownedCell.addParticle(p1);
    ownedCell.addParticle(p2);
    ownedCell.addParticle(p3);
    haloCellRight.addParticle(p1Halo);
    functor.SoALoader(ownedCell, ownedCell._particleSoABuffer, 0, false);
    functor.SoALoader(haloCellLeft, haloCellLeft._particleSoABuffer, 0, false);
    functor.SoALoader(haloCellRight, haloCellRight._particleSoABuffer, 0, false);

    functor.initTraversal();
    functor.SoAFunctorPair(ownedCell._particleSoABuffer, haloCellLeft._particleSoABuffer, newton3);
    functor.SoAFunctorPair(ownedCell._particleSoABuffer, haloCellRight._particleSoABuffer, newton3);
    functor.endTraversal(newton3);
  } else {
    FMCell ownedCell1{}, ownedCell2{}, ownedCell3{}, haloCellLeft1{}, haloCellLeft2{}, haloCellRight{};
    haloCellLeft1.addParticle(p2Halo);
    haloCellLeft2.addParticle(p3Halo);
    ownedCell1.addParticle(p1);
    ownedCell2.addParticle(p2);
    ownedCell3.addParticle(p3);
    haloCellRight.addParticle(p1Halo);
    functor.SoALoader(ownedCell1, ownedCell1._particleSoABuffer, 0, false);
    functor.SoALoader(ownedCell2, ownedCell2._particleSoABuffer, 0, false);
    functor.SoALoader(ownedCell3, ownedCell3._particleSoABuffer, 0, false);
    functor.SoALoader(haloCellLeft1, haloCellLeft1._particleSoABuffer, 0, false);
    functor.SoALoader(haloCellLeft2, haloCellLeft2._particleSoABuffer, 0, false);
    functor.SoALoader(haloCellRight, haloCellRight._particleSoABuffer, 0, false);

    functor.initTraversal();
    functor.SoAFunctorTriple(ownedCell1._particleSoABuffer, haloCellLeft1._particleSoABuffer,
                             haloCellLeft2._particleSoABuffer, newton3);
    functor.SoAFunctorTriple(ownedCell2._particleSoABuffer, ownedCell3._particleSoABuffer,
                             haloCellRight._particleSoABuffer, newton3);
    if (not newton3) {
      functor.SoAFunctorTriple(ownedCell3._particleSoABuffer, ownedCell2._particleSoABuffer,
                               haloCellRight._particleSoABuffer, newton3);
    }
    functor.endTraversal(newton3);
  }

  // Get reference globals for both mixed owned/halo triplets
  // Note on virial: W = f_12 * r_12 + f_13 * r_13 + f_23 * r_23
  // Only for non-periodic BCs, this is equal to W = f_1 * r_1 + f_2 * r_2 + f_3 + r_3
  const auto expectedVirial1 = calculateATMVirialTotal(p1OwnedPos, p2HaloPos, p3HaloPos, cutoff, nu);
  const auto expectedVirial2 = calculateATMVirialTotal(p1HaloPos, p2OwnedPos, p3OwnedPos, cutoff, nu);

  constexpr double expectedEnergy1 = calculateATMPotential(p1OwnedPos, p2HaloPos, p3HaloPos, cutoff, nu);
  constexpr double expectedEnergy2 = calculateATMPotential(p1HaloPos, p2OwnedPos, p3OwnedPos, cutoff, nu);

  // Sanity check; Globals should only depend on relative positions
  EXPECT_NEAR(expectedVirial1, expectedVirial2, absDelta);
  EXPECT_NEAR(expectedEnergy1, expectedEnergy2, absDelta);

  // Get calculated globals from functor
  const double potentialEnergy = functor.getPotentialEnergy();
  const double globalVirial = functor.getVirial();

  EXPECT_NEAR(potentialEnergy, expectedEnergy1, absDelta)
      << "FunctorType: " << functorToTest << ", newton3: " << newton3;
  EXPECT_NEAR(globalVirial, expectedVirial1, absDelta) << "FunctorType: " << functorToTest << ", newton3: " << newton3;
}

// Tests global values for the AoS functor with various owned/halo configurations
TYPED_TEST_P(ATMFunctorTestGlobals, testAoSATMFunctorGlobals) {
  using FuncType = TypeParam;
  using TestType = ATMFunctorTestGlobals<FuncType>;

  for (typename TestType::where_type where : {TestType::where_type::allInside, TestType::where_type::ininout,
                                              TestType::where_type::inoutout, TestType::where_type::allOutside}) {
    for (bool newton3 : {false, true}) {
      if (auto msg =
              this->shouldSkipIfNotImplemented([&]() { this->ATMFunctorTestGlobalsNoMixingAoS(where, newton3); });
          msg != "") {
        GTEST_SKIP() << msg;
      }
    }
  }
}

// Tests global values for all the SoA functor variants with various owned/halo configurations
TYPED_TEST_P(ATMFunctorTestGlobals, testSoAATMFunctorGlobals) {
  using FuncType = TypeParam;
  using TestType = ATMFunctorTestGlobals<FuncType>;
  for (typename TestType::SoAFunctorType soaFunctorType :
       {TestType::single, TestType::pair12, TestType::pair21, TestType::triple, /*TestType::verlet*/}) {
    for (bool newton3 : {false, true}) {
      for (typename TestType::where_type where : {TestType::where_type::allInside, TestType::where_type::ininout,
                                                  TestType::where_type::inoutout, TestType::where_type::allOutside}) {
        this->ATMFunctorTestSoAGlobals(where, soaFunctorType, newton3);
      }
    }
  }
}

// Tests global values for a triplet that interacts over periodic boundaries.
TYPED_TEST_P(ATMFunctorTestGlobals, testGlobalsWithPeriodicBCs) {
  for (const auto &functorToTest : {"AoS", "SoAPair", "SoATriple"}) {
    for (bool newton3 : {false, true}) {
      this->ATMFunctorTestGlobalsPeriodicBCs(functorToTest, newton3);
    }
  }
}

TYPED_TEST_P(ATMFunctorTestGlobals, testAoSATMFunctorGlobalsOpenMPParallel) {
  using FuncType = TypeParam;

  constexpr bool newton3 = true;

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

  std::string msg;
  // This is a basic check for the global calculations, by checking the handling of two particle interactions in
  // parallel. If interactions are dangerous, archer will complain.
  // reduction for appending strings: "abc" + "def" -> "abcdef"
  AUTOPAS_OPENMP(declare reduction(stringAppend : std::string : omp_out.append(omp_in)))
  AUTOPAS_OPENMP(parallel reduction(stringAppend : msg)) {
    AUTOPAS_OPENMP(sections) {
      AUTOPAS_OPENMP(section) {
        msg += this->shouldSkipIfNotImplemented([&]() { functor.AoSFunctor(p1, p2, p3, newton3); });
      }  // pragma omp section
      AUTOPAS_OPENMP(section) {
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

  const double expectedEnergyTriplet1 = calculateATMPotential(p1Pos, p2Pos, p3Pos, this->cutoff, this->nu);
  const double expectedVirialTriplet1 = calculateATMVirialTotal(p1Pos, p2Pos, p3Pos, this->cutoff, this->nu);
  const double expectedEnergyTriplet2 = calculateATMPotential(p4Pos, p5Pos, p6Pos, this->cutoff, this->nu);
  const double expectedVirialTriplet2 = calculateATMVirialTotal(p4Pos, p5Pos, p6Pos, this->cutoff, this->nu);
  const double expectedEnergy = expectedEnergyTriplet1 + expectedEnergyTriplet2;
  const double expectedVirial = expectedVirialTriplet1 + expectedVirialTriplet2;

  EXPECT_NEAR(potentialEnergy, expectedEnergy, this->absDelta) << "newton3: " << newton3;
  EXPECT_NEAR(virial, expectedVirial, this->absDelta) << "newton3: " << newton3;
}

TYPED_TEST_P(ATMFunctorTestGlobals, testATMFunctorGlobalsThrowBad) {
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

REGISTER_TYPED_TEST_SUITE_P(ATMFunctorTestGlobals, testAoSATMFunctorGlobals, testATMFunctorGlobalsThrowBad,
                            testAoSATMFunctorGlobalsOpenMPParallel, testSoAATMFunctorGlobals,
                            testGlobalsWithPeriodicBCs);

using MyTypes = ::testing::Types<ATMFunNoMixGlob
#ifdef __AVX__

// TODO: Add AVX Functor
#endif
                                 >;

INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, ATMFunctorTestGlobals, MyTypes);
