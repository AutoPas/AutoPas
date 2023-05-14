/**
 * @file LJFunctorTestGlobals.cpp
 * @author F. Gratl
 * @date 20.03.20
 */

#include "LJFunctorTestGlobals.h"

TYPED_TEST_SUITE_P(LJFunctorTestGlobals);

template <class FuncType>
void LJFunctorTestGlobals<FuncType>::testAoSGlobals(LJFunctorTestGlobals<FuncType>::where_type where, bool newton3) {
  double p1X = 0.;
  double p1Y = 0.;
  double p1Z = 0.;
  double p2X = 0.1;
  double p2Y = 0.2;
  double p2Z = 0.3;

  FuncType functor(cutoff);
  functor.setParticleProperties(epsilon * 24, sigma);
  double xOffset;
  double whereFactor;
  std::string where_str;
  bool owned1, owned2;
  switch (where) {
    case inside:
      xOffset = 0.;
      whereFactor = 1.;
      where_str = "inside";
      owned1 = owned2 = true;
      break;
    case boundary:
      xOffset = 4.9;
      // if there are no duplicated calculations all calculations count, therefore factor = 1
      // if there are duplicated calculations there shouldn't be only a partial (factor 0.5) contribution to the energy
      // if one particle is inside and one outside
      whereFactor = 0.5;
      where_str = "boundary";
      owned1 = true;
      owned2 = false;
      break;
    case outside:
      xOffset = 5.0;
      // if there are no duplicated calculations all calculations count, therefore factor = 1
      // if there are duplicated calculations there shouldn't be any contribution to the energy if both particles are
      // outside
      whereFactor = 0.;
      where_str = "outside";
      owned1 = owned2 = false;
      break;
    default:
      FAIL() << "not in enum where_type";
  }

  std::array<double, 3> p1Pos{p1X + xOffset, p1Y, p1Z};
  std::array<double, 3> p2Pos{p2X + xOffset, p2Y, p2Z};

  Molecule p1(p1Pos, {0., 0., 0.}, 0, 0);
  p1.setOwnershipState(owned1 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
  Molecule p2(p2Pos, {0., 0., 0.}, 1, 0);
  p2.setOwnershipState(owned2 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
  functor.initTraversal();

  functor.AoSFunctor(p1, p2, newton3);
  if (not newton3) {
    functor.AoSFunctor(p2, p1, newton3);
  }
  functor.endTraversal(newton3);

  double upot = functor.getUpot();
  double virial = functor.getVirial();

  double expectedEnergy = calculateLJPotential(p1Pos, p2Pos, cutoff, sigma, epsilon);
  double expectedVirial = calculateLJVirialTotal(p1Pos, p2Pos, cutoff, sigma, epsilon);

  EXPECT_NEAR(upot, whereFactor * expectedEnergy, absDelta) << "where: " << where_str << ", newton3: " << newton3;
  EXPECT_NEAR(virial, whereFactor * expectedVirial, absDelta) << "where: " << where_str << ", newton3: " << newton3;
}

template <class FuncType>
void LJFunctorTestGlobals<FuncType>::testAoSGlobalsMixedN3(LJFunctorTestGlobals<FuncType>::where_type where) {
  double xOffset;
  std::string where_str;
  bool owned1, owned2, owned3;
  switch (where) {
    case inside:
      xOffset = 0.;
      where_str = "inside";
      owned1 = owned2 = owned3 = true;
      break;
    case boundary:
      xOffset = 4.1;
      where_str = "boundary";
      owned1 = true;
      owned2 = false;
      owned3 = true;
      break;
    case outside:
      xOffset = 5.;
      where_str = "outside";
      owned1 = owned2 = owned3 = false;
      break;
    default:
      FAIL() << "not in enum where_type";
  }

  double refVirial = 0;
  double refUpot = 0;

  bool firstIterationDone = false;

  for (int i = 0; i < 8; i++) {
    for (bool newton3 : {false, true}) {
      FuncType functor(cutoff);
      functor.setParticleProperties(epsilon * 24, sigma);

      Molecule p1({0. + xOffset, 0., 0.}, {0., 0., 0.}, 0, 0);
      p1.setOwnershipState(owned1 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
      Molecule p2({0.9 + xOffset, 0, 0}, {0., 0., 0.}, 1, 0);
      p2.setOwnershipState(owned2 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
      Molecule p3({0 + xOffset, 0.4, 0}, {0., 0., 0.}, 2, 0);
      p3.setOwnershipState(owned3 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);

      functor.initTraversal();

      const bool p1p2 = (bool)((i >> 0) & 1);
      const bool p1p3 = (bool)((i >> 1) & 1);
      const bool p2p3 = (bool)((i >> 2) & 1);

      functor.AoSFunctor(p1, p2, p1p2);
      if (not p1p2) {
        functor.AoSFunctor(p2, p1, p1p2);
      }
      functor.AoSFunctor(p1, p3, p1p3);
      if (not p1p3) {
        functor.AoSFunctor(p3, p1, p1p3);
      }
      functor.AoSFunctor(p2, p3, p2p3);
      if (not p2p3) {
        functor.AoSFunctor(p3, p2, p2p3);
      }

      functor.endTraversal(newton3);

      auto finalUpot = functor.getUpot();
      auto finalVirial = functor.getVirial();

      // We assume the globals in the first iteration with p1p2 = p1p3 = p2p3 = false and global newton3 = false
      // are calculated correctly, as this should be ensured by the testAoSGlobals tests. We compare the globals
      // calculated with different newton3 configurations per pairwise interaction agaist the first result.
      if (firstIterationDone) {
        EXPECT_NEAR(finalUpot, refUpot, absDelta)
            << "newton3-configuration: (p1<->p2: " << p1p2 << ", p1<->p3: " << p1p3 << ", p2<->p3: " << p2p3
            << "), where: " << where_str << ", global newton3: " << newton3;
        EXPECT_NEAR(finalVirial, refVirial, absDelta)
            << "newton3-configuration: (p1<->p2: " << p1p2 << ", p1<->p3: " << p1p3 << ", p2<->p3: " << p2p3
            << "), where: " << where_str << ", global newton3: " << newton3;
      } else {
        refUpot = finalUpot;
        refVirial = finalVirial;
        firstIterationDone = true;
      }
    }
  }
}

template <class FuncType>
void LJFunctorTestGlobals<FuncType>::testSoAGlobals(LJFunctorTestGlobals<FuncType>::where_type where, bool newton3,
                                                    InteractionType interactionType,
                                                    size_t additionalParticlesToVerletNumber,
                                                    uint64_t numParticleReplicas, bool mixedNewton3FunctorCalls) {
  constexpr bool shifting = true;
  constexpr bool mixing = false;

  // static coords for the particles
  constexpr double p1X = 0.;
  constexpr double p1Y = 0.;
  constexpr double p1Z = 0.;
  constexpr double p2X = 0.1;
  constexpr double p2Y = 0.2;
  constexpr double p2Z = 0.3;
  constexpr double pAddX = 1.2;
  constexpr double pAddY = 0.;
  constexpr double pAddZ = 0.;

  // calculate the reference values for the globals
  std::array<double, 3> p1Pos{p1X, p1Y, p1Z};
  std::array<double, 3> p2Pos{p2X, p2Y, p2Z};
  double expectedEnergy = calculateLJPotential(p1Pos, p2Pos, cutoff, sigma, epsilon);
  double expectedVirial = calculateLJVirialTotal(p1Pos, p2Pos, cutoff, sigma, epsilon);

  autopas::LJFunctor<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> functor(cutoff);
  functor.setParticleProperties(epsilon * 24, sigma);
  double xOffset;
  double whereFactor = 0.;
  std::string where_str;
  bool owned1, owned2;
  switch (where) {
    case inside:
      xOffset = 0.;
      where_str = "inside";
      owned1 = owned2 = true;
      break;
    case boundary:
      xOffset = 4.9;
      where_str = "boundary";
      owned1 = true;
      owned2 = false;
      break;
    case outside:
      xOffset = 5.0;
      where_str = "outside";
      owned1 = owned2 = false;
      break;
    default:
      FAIL() << "not in enum where_type";
  }
  FMCell cell1, cell2, cell3, cell4;
  for (uint64_t replicaID = 0; replicaID < numParticleReplicas; ++replicaID) {
    if (replicaID > 0) {
      if (owned1 && owned2) {
        owned1 = false;
        owned2 = false;
      } else {
        owned1 = true;
        owned2 = true;
      }
    }

    Molecule p1({p1X + xOffset, p1Y + 2. * replicaID, p1Z}, {0., 0., 0.}, 2 * replicaID, 0);
    p1.setOwnershipState(owned1 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
    Molecule p2({p2X + xOffset, p2Y + 2. * replicaID, p2Z}, {0., 0., 0.}, 2 * replicaID + 1, 0);
    p2.setOwnershipState(owned2 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);

    // calculate whereFactor:
    double currentWhereFactor = 0.;
    currentWhereFactor += p1.isOwned() ? .5 : 0.;
    currentWhereFactor += p2.isOwned() ? .5 : 0.;
    whereFactor += currentWhereFactor;

    switch (interactionType) {
      case InteractionType::verlet:
        // same as own: add to cell1.
      case InteractionType::own:
        cell1.addParticle(p1);
        cell1.addParticle(p2);

        cell2.addParticle(p1);
        cell2.addParticle(p2);

        break;
      case InteractionType::pair:
        cell1.addParticle(p1);
        cell2.addParticle(p2);

        cell3.addParticle(p1);
        cell4.addParticle(p2);

        break;
      default:
        FAIL();
    }
  }

  if (interactionType == InteractionType::verlet) {
    Molecule pAdditional({pAddX + xOffset, pAddY, pAddZ}, {0., 0., 0.}, std::numeric_limits<uint64_t>::max(), 0);
    pAdditional.setOwnershipState(owned2 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
    // add dummy particles outside of the cutoff. this will only change the number of particles in the verlet lists,
    // but will leave the desired result unchanged. the higher number of particles is useful to test the soa
    // functor version of verlet lists.
    for (size_t i = 0; i < additionalParticlesToVerletNumber; ++i) {
      cell1.addParticle(pAdditional);

      cell2.addParticle(pAdditional);
    }
  }

  functor.initTraversal();

  functor.SoALoader(cell1, cell1._particleSoABuffer, 0);
  functor.SoALoader(cell2, cell2._particleSoABuffer, 0);
  functor.SoALoader(cell3, cell3._particleSoABuffer, 0);
  functor.SoALoader(cell4, cell4._particleSoABuffer, 0);

  switch (interactionType) {
    case InteractionType::verlet: {
      // Build verlet list
      std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> neighborListNonN3(2 * numParticleReplicas);
      for (uint64_t replicaID = 0; replicaID < numParticleReplicas; ++replicaID) {
        neighborListNonN3[2 * replicaID].push_back(2 * replicaID + 1);
        for (size_t i = 0; i < additionalParticlesToVerletNumber; ++i) {
          neighborListNonN3[2 * replicaID].push_back(2 * numParticleReplicas + i);
        }
        neighborListNonN3[2 * replicaID + 1].push_back(2 * replicaID);
        for (size_t i = 0; i < additionalParticlesToVerletNumber; ++i) {
          neighborListNonN3[2 * replicaID + 1].push_back(2 * numParticleReplicas + i);
        }
      }
      std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> neighborListN3(2 * numParticleReplicas);
      for (uint64_t replicaID = 0; replicaID < numParticleReplicas; ++replicaID) {
        neighborListN3[2 * replicaID].push_back(2 * replicaID + 1);
        for (size_t i = 0; i < additionalParticlesToVerletNumber; ++i) {
          neighborListN3[2 * replicaID].push_back(2 * numParticleReplicas + i);
        }
      }
      for (uint64_t i = 0; i < 2 * numParticleReplicas; ++i) {
        if (not mixedNewton3FunctorCalls) {
          functor.SoAFunctorVerlet(cell1._particleSoABuffer, i, newton3 ? neighborListN3[i] : neighborListNonN3[i],
                                   newton3);
          functor.SoAFunctorVerlet(cell2._particleSoABuffer, i, newton3 ? neighborListN3[i] : neighborListNonN3[i],
                                   newton3);
        } else {
          functor.SoAFunctorVerlet(cell1._particleSoABuffer, i, newton3 ? neighborListN3[i] : neighborListNonN3[i],
                                   newton3);
          functor.SoAFunctorVerlet(cell2._particleSoABuffer, i,
                                   (not newton3) ? neighborListN3[i] : neighborListNonN3[i], not newton3);
        }
      }
    } break;
    case InteractionType::own:
      if (not mixedNewton3FunctorCalls) {
        functor.SoAFunctorSingle(cell1._particleSoABuffer, newton3);
        functor.SoAFunctorSingle(cell2._particleSoABuffer, newton3);
      } else {
        functor.SoAFunctorSingle(cell1._particleSoABuffer, newton3);
        functor.SoAFunctorSingle(cell2._particleSoABuffer, not newton3);
      }
      break;
    case InteractionType::pair:
      if (not mixedNewton3FunctorCalls) {
        functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
        functor.SoAFunctorPair(cell3._particleSoABuffer, cell4._particleSoABuffer, newton3);
        if (not newton3) {
          functor.SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, newton3);
          functor.SoAFunctorPair(cell4._particleSoABuffer, cell3._particleSoABuffer, newton3);
        }
      } else {
        if (newton3) {
          functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
          functor.SoAFunctorPair(cell3._particleSoABuffer, cell4._particleSoABuffer, not newton3);
          functor.SoAFunctorPair(cell4._particleSoABuffer, cell3._particleSoABuffer, not newton3);
        } else {
          functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
          functor.SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, newton3);
          functor.SoAFunctorPair(cell3._particleSoABuffer, cell4._particleSoABuffer, not newton3);
        }
      }
      break;
  }
  functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
  functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
  functor.SoAExtractor(cell3, cell3._particleSoABuffer, 0);
  functor.SoAExtractor(cell4, cell4._particleSoABuffer, 0);

  functor.endTraversal(newton3);

  double upot = functor.getUpot();
  double virial = functor.getVirial();

  // we multiply the expected values by two as we have always p1 and p2 both in cell1 and cell2
  EXPECT_NEAR(upot, whereFactor * expectedEnergy * 2., absDelta)
      << "where: " << where_str << ", newton3: " << newton3
      << ", interactionType: " << (interactionType == pair ? "pair" : (interactionType == own ? "own" : "verlet"))
      << ", additionalVerletDummyParticles: " << additionalParticlesToVerletNumber
      << ", numParticleReplicas: " << numParticleReplicas;
  EXPECT_NEAR(virial, whereFactor * expectedVirial * 2., absDelta)
      << "where: " << where_str << ", newton3: " << newton3
      << ", interactionType: " << (interactionType == pair ? "pair" : (interactionType == own ? "own" : "verlet"))
      << ", additionalVerletDummyParticles: " << additionalParticlesToVerletNumber
      << ", numParticleReplicas: " << numParticleReplicas;
}

TYPED_TEST_P(LJFunctorTestGlobals, testAoSFunctorGlobals) {
  using FuncType = TypeParam;
  using TestType = LJFunctorTestGlobals<FuncType>;

  for (typename TestType::where_type where :
       {TestType::where_type::inside, TestType::where_type::boundary, TestType::where_type::outside}) {
    for (bool newton3 : {false, true}) {
      if (auto msg = this->shouldSkipIfNotImplemented([&]() { this->testAoSGlobals(where, newton3); }); msg != "") {
        GTEST_SKIP() << msg;
      }
    }
  }
}

TYPED_TEST_P(LJFunctorTestGlobals, testAoSFunctorGlobalsMixedN3) {
  using FuncType = TypeParam;
  using TestType = LJFunctorTestGlobals<FuncType>;

  for (typename TestType::where_type where :
       {TestType::where_type::inside, TestType::where_type::boundary, TestType::where_type::outside}) {
    if (auto msg = this->shouldSkipIfNotImplemented([&]() { this->testAoSGlobalsMixedN3(where); }); msg != "") {
      GTEST_SKIP() << msg;
    }
  }
}

TYPED_TEST_P(LJFunctorTestGlobals, testSoAFunctorGlobalsOwn) {
  using FuncType = TypeParam;
  using TestType = LJFunctorTestGlobals<FuncType>;

  for (bool mixedNewton3FunctorCalls : {false, true}) {
    // the own functor can only be called for inner or outside pairs! (if two particles lie in one cell they can be
    // either both inside the process or neither of them is)
    for (typename TestType::where_type where : {TestType::inside, TestType::outside}) {
      for (bool newton3 : {false, true}) {
        for (uint64_t numParticleReplicas : {1, 2}) {
          this->testSoAGlobals(where, newton3, TestType::own, 0, numParticleReplicas, mixedNewton3FunctorCalls);
        }
      }
    }
  }
}

TYPED_TEST_P(LJFunctorTestGlobals, testSoAFunctorGlobalsVerlet) {
  using FuncType = TypeParam;
  using TestType = LJFunctorTestGlobals<FuncType>;

  for (bool mixedNewton3FunctorCalls : {false, true}) {
    for (size_t additionalDummyParticles = 0; additionalDummyParticles < 30; additionalDummyParticles += 5) {
      for (typename TestType::where_type where : {TestType::inside, TestType::boundary, TestType::outside}) {
        for (bool newton3 : {false, true}) {
          for (uint64_t numParticleReplicas : {1, 2}) {
            this->testSoAGlobals(where, newton3, TestType::verlet, additionalDummyParticles, numParticleReplicas,
                                 mixedNewton3FunctorCalls);
          }
        }
      }
    }
  }
}

TYPED_TEST_P(LJFunctorTestGlobals, testSoAFunctorGlobalsPair) {
  using FuncType = TypeParam;
  using TestType = LJFunctorTestGlobals<FuncType>;

  for (bool mixedNewton3FunctorCalls : {false, true}) {
    for (typename TestType::where_type where : {TestType::inside, TestType::boundary, TestType::outside}) {
      for (bool newton3 : {false, true}) {
        for (uint64_t numParticleReplicas : {1, 2}) {
          this->testSoAGlobals(where, newton3, TestType::pair, 0, numParticleReplicas, mixedNewton3FunctorCalls);
        }
      }
    }
  }
}

TYPED_TEST_P(LJFunctorTestGlobals, testAoSFunctorGlobalsOpenMPParallel) {
  using FuncType = TypeParam;
  using TestType = LJFunctorTestGlobals<FuncType>;

  constexpr bool newton3 = true;
  constexpr double whereFactor = 1.;  // all inside, so factor 1
  std::string where_str = "inside";

  std::array<double, 3> p1Pos{0., 0., 0.};
  std::array<double, 3> p2Pos{0.1, 0.2, 0.3};
  std::array<double, 3> p3Pos{0., 2., 0.};
  std::array<double, 3> p4Pos{0.1, 2.2, 0.3};

  Molecule p1(p1Pos, {0., 0., 0.}, 0, 0);
  Molecule p2(p2Pos, {0., 0., 0.}, 1, 0);

  Molecule p3(p3Pos, {0., 0., 0.}, 0, 0);
  Molecule p4(p4Pos, {0., 0., 0.}, 1, 0);
  FuncType functor(this->cutoff);
  functor.setParticleProperties(this->epsilon * 24, 1);

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
        msg += this->shouldSkipIfNotImplemented([&]() { functor.AoSFunctor(p1, p2, newton3); });
      }  // pragma omp section
#if defined(AUTOPAS_OPENMP)
#pragma omp section
#endif
      {
        msg += this->shouldSkipIfNotImplemented([&]() { functor.AoSFunctor(p3, p4, newton3); });
      }  // pragma omp section
    }    // pragma omp sections
  }      // pragma omp parallel

  if (not msg.empty()) {
    GTEST_SKIP() << msg;
  }

  functor.endTraversal(newton3);

  double upot = functor.getUpot();
  double virial = functor.getVirial();

  double expectedEnergy = calculateLJPotential(p1Pos, p2Pos, this->cutoff, this->sigma, this->epsilon);
  double expectedVirial = calculateLJVirialTotal(p1Pos, p2Pos, this->cutoff, this->sigma, this->epsilon);
  expectedEnergy += calculateLJPotential(p3Pos, p4Pos, this->cutoff, this->sigma, this->epsilon);
  expectedVirial += calculateLJVirialTotal(p3Pos, p4Pos, this->cutoff, this->sigma, this->epsilon);

  EXPECT_NEAR(upot, whereFactor * expectedEnergy, this->absDelta) << "where: " << where_str << ", newton3: " << newton3;
  EXPECT_NEAR(virial, whereFactor * expectedVirial, this->absDelta)
      << "where: " << where_str << ", newton3: " << newton3;
}

TYPED_TEST_P(LJFunctorTestGlobals, testFunctorGlobalsThrowBad) {
  using exception_type = autopas::utils::ExceptionHandler::AutoPasException;

  using FuncType = TypeParam;
  FuncType functor(this->cutoff);

  // getupot without postprocessing is not allowed
  EXPECT_THROW(functor.getUpot(), exception_type);
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

REGISTER_TYPED_TEST_SUITE_P(LJFunctorTestGlobals, testAoSFunctorGlobals, testAoSFunctorGlobalsOpenMPParallel,
                            testSoAFunctorGlobalsOwn, testSoAFunctorGlobalsPair, testSoAFunctorGlobalsVerlet,
                            testFunctorGlobalsThrowBad, testAoSFunctorGlobalsMixedN3);

using MyTypes = ::testing::Types<LJFunShiftNoMixGlob
#ifdef __AVX__
                                 ,
                                 LJFunAVXShiftNoMixGlob
#endif
                                 >;
INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, LJFunctorTestGlobals, MyTypes);
