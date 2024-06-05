/**
 * @file LJFunctorTestVs.cpp
 * @author F. Gratl
 * @date 20.03.20
 */

#include "LJFunctorTestVs.h"

TYPED_TEST_SUITE_P(LJFunctorTestVs);

/**
 * Test LJ Functor with mixing parameters stored with molecule and non-mixing parameters stored in functor
 */
TYPED_TEST_P(LJFunctorTestVs, testMixingVsNoMixingSoA) {
  using FunMixing = typename std::tuple_element_t<0, TypeParam>;
  using FunNoMixing = typename std::tuple_element_t<1, TypeParam>;

  FunNoMixing funNoMixing(this->cutoff);
  funNoMixing.setParticleProperties(24 * this->epsilon, this->sigma);

  FunMixing funMixing(this->cutoff);

  size_t numParticlesPerCell = 9;

  Molecule defaultParticle({0.,0.,0.}, {0.,0.,0.}, 0, std::sqrt(this->epsilon), this->sigma/2);
  FMCell cell1NoPPL;
  FMCell cell2NoPPL;
  autopasTools::generators::RandomGenerator::fillWithParticles(cell1NoPPL, defaultParticle, {0, 0, 0}, {5, 5, 5},
                                                               numParticlesPerCell, 42);
  autopasTools::generators::RandomGenerator::fillWithParticles(cell2NoPPL, defaultParticle, {0, 0, 0}, {5, 5, 5},
                                                               numParticlesPerCell, 43);

  funNoMixing.SoALoader(cell1NoPPL, cell1NoPPL._particleSoABuffer, 0, /*skipSoAResize*/ false);
  funNoMixing.SoALoader(cell2NoPPL, cell2NoPPL._particleSoABuffer, 0, /*skipSoAResize*/ false);

  FMCell cell1PPL(cell1NoPPL);
  FMCell cell2PPL(cell2NoPPL);

  constexpr bool newton3 = true;
  funNoMixing.SoAFunctorPair(cell1NoPPL._particleSoABuffer, cell2NoPPL._particleSoABuffer, newton3);
  funMixing.SoAFunctorPair(cell1PPL._particleSoABuffer, cell2PPL._particleSoABuffer, newton3);

  funMixing.SoAExtractor(cell1PPL, cell1PPL._particleSoABuffer, 0);
  funMixing.SoAExtractor(cell2PPL, cell2PPL._particleSoABuffer, 0);
  funNoMixing.SoAExtractor(cell1NoPPL, cell1NoPPL._particleSoABuffer, 0);
  funNoMixing.SoAExtractor(cell2NoPPL, cell2NoPPL._particleSoABuffer, 0);

  for (size_t i = 0; i < numParticlesPerCell; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      EXPECT_EQ(cell1NoPPL[i], cell1PPL[i]) << "cell1NoPPL[i] = " << cell1NoPPL[i].toString() << std::endl
                                            << "cell1PPL[i] = " << cell1PPL[i].toString();
      EXPECT_EQ(cell2NoPPL[i], cell2PPL[i]) << "cell2NoPPL[i] = " << cell2NoPPL[i].toString() << std::endl
                                            << "cell2PPL[i] = " << cell2PPL[i].toString();
    }
  }
}

TYPED_TEST_P(LJFunctorTestVs, testMixingVsNoMixingAoS) {
  using FunMixing = typename std::tuple_element_t<0, TypeParam>;
  using FunNoMixing = typename std::tuple_element_t<1, TypeParam>;

  FunNoMixing funNoMixing(this->cutoff);
  funNoMixing.setParticleProperties(24 * this->epsilon, this->sigma);

  FunMixing funMixing(this->cutoff);

  std::vector<Molecule> moleculesNoPPL = {Molecule({0, 0, 0}, {0, 0, 0}, 0, std::sqrt(this->epsilon), this->sigma/2), Molecule({0, 0, 1}, {0, 0, 0}, 1, std::sqrt(this->epsilon), this->sigma/2)};
  std::vector<Molecule> moleculesPPL(moleculesNoPPL);

  constexpr bool newton3 = false;
  if (auto msg = this->shouldSkipIfNotImplemented([&]() {
        funMixing.AoSFunctor(moleculesPPL[0], moleculesPPL[1], newton3);
        funNoMixing.AoSFunctor(moleculesNoPPL[0], moleculesNoPPL[1], newton3);
      });
      msg != "") {
    GTEST_SKIP() << msg;
  }

  // sanity check
  ASSERT_GT(moleculesPPL.size(), 0);
  // Molecules should be exactly the same
  EXPECT_THAT(moleculesNoPPL, testing::ElementsAreArray(moleculesPPL));
}

REGISTER_TYPED_TEST_SUITE_P(LJFunctorTestVs, testMixingVsNoMixingSoA, testMixingVsNoMixingAoS);

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
#ifdef __ARM_FEATURE_SVE
    ,
    // LJFunctorSVE<mixing> VS LJFunctorSVE<not mixing>
    std::tuple<LJFunSVEShiftMixGlob, LJFunSVEShiftNoMixGlob>,
    // LJFunctor<mixing> VS LJFunctorSVE<not mixing>
    std::tuple<LJFunShiftMixGlob, LJFunSVEShiftNoMixGlob>,
    // LJFunctorSVE<mixing> VS LJFunctor<not mixing>
    std::tuple<LJFunSVEShiftMixGlob, LJFunShiftNoMixGlob>
#endif
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, LJFunctorTestVs, MyTypes);