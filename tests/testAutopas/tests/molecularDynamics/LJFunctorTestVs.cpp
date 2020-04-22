/**
 * @file LJFunctorTestVs.cpp
 * @author F. Gratl
 * @date 20.03.20
 */

#include "LJFunctorTestVs.h"

TYPED_TEST_SUITE_P(LJFunctorTestVs);

TYPED_TEST_P(LJFunctorTestVs, testSetPropertiesVSPPLSoA) {
  using FunPPL = typename std::tuple_element_t<0, TypeParam>;
  using FunNoPPL = typename std::tuple_element_t<1, TypeParam>;

  FunNoPPL funNoPPL(this->cutoff);
  funNoPPL.setParticleProperties(24 * this->epsilon, this->sigma);

  ParticlePropertiesLibrary<double, size_t> particlePropertiesLibrary(this->cutoff);
  particlePropertiesLibrary.addType(0, this->epsilon, this->sigma, 1);
  FunPPL funPPL(this->cutoff, particlePropertiesLibrary);

  size_t numParticlesPerCell = 9;

  Molecule defaultParticle;
  FMCell cell1NoPPL;
  FMCell cell2NoPPL;
  autopasTools::generators::RandomGenerator::fillWithParticles(cell1NoPPL, defaultParticle, {0, 0, 0}, {5, 5, 5},
                                                               numParticlesPerCell, 42);
  autopasTools::generators::RandomGenerator::fillWithParticles(cell2NoPPL, defaultParticle, {0, 0, 0}, {5, 5, 5},
                                                               numParticlesPerCell, 43);

  funNoPPL.SoALoader(cell1NoPPL, cell1NoPPL._particleSoABuffer);
  funNoPPL.SoALoader(cell2NoPPL, cell2NoPPL._particleSoABuffer);

  FMCell cell1PPL(cell1NoPPL);
  FMCell cell2PPL(cell2NoPPL);

  constexpr bool newton3 = true;
  constexpr bool cellWiseOwnedState = true;
  funNoPPL.SoAFunctorPair(cell1NoPPL._particleSoABuffer, cell2NoPPL._particleSoABuffer, newton3, cellWiseOwnedState);
  funPPL.SoAFunctorPair(cell1PPL._particleSoABuffer, cell2PPL._particleSoABuffer, newton3, cellWiseOwnedState);

  funPPL.SoAExtractor(cell1PPL, cell1PPL._particleSoABuffer);
  funPPL.SoAExtractor(cell2PPL, cell2PPL._particleSoABuffer);
  funNoPPL.SoAExtractor(cell1NoPPL, cell1NoPPL._particleSoABuffer);
  funNoPPL.SoAExtractor(cell2NoPPL, cell2NoPPL._particleSoABuffer);

  for (size_t i = 0; i < numParticlesPerCell; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      EXPECT_EQ(cell1NoPPL[i], cell1PPL[i]) << "cell1NoPPL[i] = " << cell1NoPPL[i].toString() << std::endl
                                            << "cell1PPL[i] = " << cell1PPL[i].toString();
      EXPECT_EQ(cell2NoPPL[i], cell2PPL[i]) << "cell2NoPPL[i] = " << cell2NoPPL[i].toString() << std::endl
                                            << "cell2PPL[i] = " << cell2PPL[i].toString();
    }
  }
}

TYPED_TEST_P(LJFunctorTestVs, testSetPropertiesVSPPLAoS) {
  using FunPPL = typename std::tuple_element_t<0, TypeParam>;
  using FunNoPPL = typename std::tuple_element_t<1, TypeParam>;

  FunNoPPL funNoPPL(this->cutoff);
  funNoPPL.setParticleProperties(24 * this->epsilon, this->sigma);

  ParticlePropertiesLibrary<double, size_t> particlePropertiesLibrary(this->cutoff);
  particlePropertiesLibrary.addType(0, this->epsilon, this->sigma, 1);
  FunPPL funPPL(this->cutoff, particlePropertiesLibrary);

  std::vector<Molecule> moleculesNoPPL = {Molecule({0, 0, 0}, {0, 0, 0}, 0, 0), Molecule({0, 0, 1}, {0, 0, 0}, 1, 0)};
  std::vector<Molecule> moleculesPPL(moleculesNoPPL);

  constexpr bool newton3 = false;
  if (auto msg = this->shouldSkipIfNotImplemented([&]() {
        funPPL.AoSFunctor(moleculesPPL[0], moleculesPPL[1], newton3);
        funNoPPL.AoSFunctor(moleculesNoPPL[0], moleculesNoPPL[1], newton3);
      });
      msg != "") {
    GTEST_SKIP() << msg;
  }

  // sanity check
  ASSERT_GT(moleculesPPL.size(), 0);
  // Molecules should be exactly the same
  EXPECT_THAT(moleculesNoPPL, testing::ElementsAreArray(moleculesPPL));
}

REGISTER_TYPED_TEST_SUITE_P(LJFunctorTestVs, testSetPropertiesVSPPLSoA, testSetPropertiesVSPPLAoS);

/**
 * Compare:
 * - Mixing vs not mixing (-> using ppl vs not using ppl when only one type of particle exists)
 * - AVX vs not AVX
 * - combinations of the above
 */
using MyTypes = ::testing::Types<
    // LJFunctor<mixing> VS LJFunctor<not mixing>
    std::tuple<LJFunShiftMixGlob, LJFunShiftNoMixGlob>
#ifdef __AVX__
    ,
    // LJFunctorAVX<mixing> VS LJFunctorAVX<not mixing>
    std::tuple<LJFunAVXShiftMixGlob, LJFunAVXShiftNoMixGlob>,
    // LJFunctor<mixing> VS LJFunctorAVX<not mixing>
    std::tuple<LJFunShiftMixGlob, LJFunAVXShiftNoMixGlob>,
    // LJFunctorAVX<mixing> VS LJFunctor<not mixing>
    std::tuple<LJFunAVXShiftMixGlob, LJFunShiftNoMixGlob>
#endif
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, LJFunctorTestVs, MyTypes);