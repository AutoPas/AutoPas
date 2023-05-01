/**
 * @file LJFunctorTestGlobals.cpp
 * @author F. Gratl
 * @date 20.03.20
 */

#include "LJFunctorTestGlobals.h"

TYPED_TEST_SUITE_P(LJFunctorTestGlobals);

template <class FuncType>
void LJFunctorTestGlobals<FuncType>::testAoSGlobals(LJFunctorTestGlobals<FuncType>::where_type where, bool newton3) {
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
  Molecule p1({0. + xOffset, 0., 0.}, {0., 0., 0.}, 0, 0);
  p1.setOwnershipState(owned1 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
  Molecule p2({0.1 + xOffset, 0.2, 0.3}, {0., 0., 0.}, 1, 0);
  p2.setOwnershipState(owned2 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
  functor.initTraversal();

  functor.AoSFunctor(p1, p2, newton3);
  if (not newton3) {
    functor.AoSFunctor(p2, p1, newton3);
  }
  functor.endTraversal(newton3);

  double upot = functor.getUpot();
  double virial = functor.getVirial();

  EXPECT_NEAR(upot, whereFactor * expectedEnergy, absDelta) << "where: " << where_str << ", newton3: " << newton3;
  EXPECT_NEAR(virial, whereFactor * expectedVirial, absDelta) << "where: " << where_str << ", newton3: " << newton3;
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

template <class FuncType>
void LJFunctorTestGlobals<FuncType>::testSoAGlobals(LJFunctorTestGlobals<FuncType>::where_type where, bool newton3,
                                                    InteractionType interactionType,
                                                    size_t additionalParticlesToVerletNumber,
                                                    uint64_t numParticleReplicas) {
  constexpr bool shifting = true;
  constexpr bool mixing = false;
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
  FMCell cell1, cell2;
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

    Molecule p1({0. + xOffset, 0. + 2. * replicaID, 0.}, {0., 0., 0.}, 2 * replicaID, 0);
    p1.setOwnershipState(owned1 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
    cell1.addParticle(p1);
    Molecule p2({0.1 + xOffset, 0.2 + 2. * replicaID, 0.3}, {0., 0., 0.}, 2 * replicaID + 1, 0);
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
        cell1.addParticle(p2);
        break;
      case InteractionType::pair:
        cell2.addParticle(p2);
        break;
      default:
        FAIL();
    }
  }

  if (interactionType == InteractionType::verlet) {
    Molecule pAdditional({1.2 + xOffset, 0., 0.}, {0., 0., 0.}, std::numeric_limits<uint64_t>::max(), 0);
    pAdditional.setOwnershipState(owned2 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
    // add dummy particles outside of the cutoff. this will only change the number of particles in the verlet lists,
    // but will leave the desired result unchanged. the higher number of particles is useful to test the soa
    // functor version of verlet lists.
    for (size_t i = 0; i < additionalParticlesToVerletNumber; ++i) {
      cell1.addParticle(pAdditional);
    }
  }

  functor.initTraversal();

  functor.SoALoader(cell1, cell1._particleSoABuffer, 0);
  functor.SoALoader(cell2, cell2._particleSoABuffer, 0);

  switch (interactionType) {
    case InteractionType::verlet: {
      // Build verlet list
      std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> neighborList(2 * numParticleReplicas);
      for (uint64_t replicaID = 0; replicaID < numParticleReplicas; ++replicaID) {
        neighborList[2 * replicaID].push_back(2 * replicaID + 1);
        for (size_t i = 0; i < additionalParticlesToVerletNumber; ++i) {
          neighborList[2 * replicaID].push_back(2 * numParticleReplicas + i);
        }
        if (not newton3) {
          neighborList[2 * replicaID + 1].push_back(2 * replicaID);
          for (size_t i = 0; i < additionalParticlesToVerletNumber; ++i) {
            neighborList[2 * replicaID + 1].push_back(2 * numParticleReplicas + i);
          }
        }
      }
      for (uint64_t i = 0; i < 2 * numParticleReplicas; ++i) {
        functor.SoAFunctorVerlet(cell1._particleSoABuffer, i, neighborList[i], newton3);
      }
    } break;
    case InteractionType::own:
      functor.SoAFunctorSingle(cell1._particleSoABuffer, newton3);
      break;
    case InteractionType::pair:
      functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
      if (not newton3) {
        functor.SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, newton3);
      }
      break;
  }
  functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
  functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);

  functor.endTraversal(newton3);

  double upot = functor.getUpot();
  double virial = functor.getVirial();

  EXPECT_NEAR(upot, whereFactor * expectedEnergy, absDelta)
      << "where: " << where_str << ", newton3: " << newton3
      << ", interactionType: " << (interactionType == pair ? "pair" : (interactionType == own ? "own" : "verlet"))
      << ", additionalVerletDummyParticles: " << additionalParticlesToVerletNumber
      << ", numParticleReplicas: " << numParticleReplicas;
  EXPECT_NEAR(virial, whereFactor * expectedVirial, absDelta)
      << "where: " << where_str << ", newton3: " << newton3
      << ", interactionType: " << (interactionType == pair ? "pair" : (interactionType == own ? "own" : "verlet"))
      << ", additionalVerletDummyParticles: " << additionalParticlesToVerletNumber
      << ", numParticleReplicas: " << numParticleReplicas;
}

TYPED_TEST_P(LJFunctorTestGlobals, testSoAFunctorGlobalsOwn) {
  using FuncType = TypeParam;
  using TestType = LJFunctorTestGlobals<FuncType>;

  // the own functor can only be called for inner or outside pairs! (if two particles lie in one cell they can be
  // either both inside the process or neither of them is)
  for (typename TestType::where_type where : {TestType::inside, TestType::outside}) {
    for (bool newton3 : {false, true}) {
      for (uint64_t numParticleReplicas : {1, 2}) {
        this->testSoAGlobals(where, newton3, TestType::own, 0, numParticleReplicas);
      }
    }
  }
}

TYPED_TEST_P(LJFunctorTestGlobals, testSoAFunctorGlobalsVerlet) {
  using FuncType = TypeParam;
  using TestType = LJFunctorTestGlobals<FuncType>;

  for (size_t additionalDummyParticles = 0; additionalDummyParticles < 30; additionalDummyParticles += 5) {
    // the own functor can only be called for inner or outside pairs! (if two particles lie in one cell they can be
    // either both inside the process or neither of them is)
    for (typename TestType::where_type where : {TestType::inside, TestType::boundary, TestType::outside}) {
      for (bool newton3 : {false, true}) {
        for (uint64_t numParticleReplicas : {1, 2}) {
          this->testSoAGlobals(where, newton3, TestType::verlet, additionalDummyParticles, numParticleReplicas);
        }
      }
    }
  }
}

TYPED_TEST_P(LJFunctorTestGlobals, testSoAFunctorGlobalsPair) {
  using FuncType = TypeParam;
  using TestType = LJFunctorTestGlobals<FuncType>;

  for (typename TestType::where_type where : {TestType::inside, TestType::boundary, TestType::outside}) {
    for (bool newton3 : {false, true}) {
      for (uint64_t numParticleReplicas : {1, 2}) {
        this->testSoAGlobals(where, newton3, TestType::pair, 0, numParticleReplicas);
      }
    }
  }
}

TYPED_TEST_P(LJFunctorTestGlobals, testAoSFunctorGlobalsOpenMPParallel) {
  using FuncType = TypeParam;
  using TestType = LJFunctorTestGlobals<FuncType>;

  constexpr bool newton3 = true;
  constexpr double multiParticleFactor = 2.;  // two particles, so factor 2
  constexpr double whereFactor = 1.;          // all inside, so factor 1
  std::string where_str = "inside";
  Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p3({0., 2., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p4({0.1, 2.2, 0.3}, {0., 0., 0.}, 1, 0);
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

  EXPECT_NEAR(upot, whereFactor * multiParticleFactor * this->expectedEnergy, this->absDelta)
      << "where: " << where_str << ", newton3: " << newton3;
  EXPECT_NEAR(virial, whereFactor * multiParticleFactor * this->expectedVirial, this->absDelta)
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
      xOffset = 4.9;
      // if there are no duplicated calculations all calculations count, therefore factor = 1
      // if there are duplicated calculations there shouldn't be only a partial (factor 0.5) contribution to the energy
      // if one particle is inside and one outside
      where_str = "boundary";
      owned1 = true;
      owned2 = false;
      owned3 = false;
      break;
    case outside:
      xOffset = 5.0;
      // if there are no duplicated calculations all calculations count, therefore factor = 1
      // if there are duplicated calculations there shouldn't be any contribution to the energy if both particles are
      // outside
      where_str = "outside";
      owned1 = owned2 = owned3 = false;
      break;
    default:
      FAIL() << "not in enum where_type";
  }

  std::vector<double> allVirials;
  std::vector<double> allUpots;

  for (int i = 0; i < 8; i++) {
    for (bool newton3 : {false, true}) {
      FuncType functor(cutoff);
      functor.setParticleProperties(epsilon * 24, sigma);

      Molecule p1({0. + xOffset, 0., 0.}, {0., 0., 0.}, 0, 0);
      p1.setOwnershipState(owned1 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
      Molecule p2({0.1 + xOffset, 0.2, 0.3}, {0., 0., 0.}, 1, 0);
      p2.setOwnershipState(owned2 ? autopas::OwnershipState::owned : autopas::OwnershipState::halo);
      Molecule p3({0.1 + xOffset, 0.3, 0.4}, {0., 0., 0.}, 2, 0);
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

      allUpots.push_back(functor.getUpot());
      allVirials.push_back(functor.getVirial());
    }
  }

  const auto checkEquality = [](const std::vector<double> v) {
    if (v.size() == 0) return false;
    for (int i = 1; i < v.size(); i++) {
      if (std::fabs(v[0] - v[i]) > 1e-3) {
        return false;
      }
    }
    return true;
  };

  const bool virialsOK = checkEquality(allVirials);
  const bool upotsOK = checkEquality(allUpots);

  EXPECT_TRUE(virialsOK);
  EXPECT_TRUE(upotsOK);
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
