/**
 * @file DEMFunctorTestNoGlobals.cpp
 * @author Joon Kim
 * @date 15.11.2024
 */

#include "DEMFunctorTestNoGlobals.h"

TYPED_TEST_SUITE_P(DEMFunctorTestNoGlobals);

TYPED_TEST_P(DEMFunctorTestNoGlobals, testAoSNoGlobals) {
  using namespace autopas::utils::ArrayMath::literals;
  using FuncType = typename TypeParam::FuncType;
  constexpr bool mixing = FuncType::getMixing();
  constexpr bool newton3 = TypeParam::newton3;

  // Given
  ParticlePropertiesLibrary<double, size_t> particlePropertiesLibrary(this->cutoff);
  std::unique_ptr<FuncType> functor;

  particlePropertiesLibrary.addSiteType(0, this->mass);
  particlePropertiesLibrary.addDEMParametersToSite(0, this->radius1, this->specificHeat);
  demLib::DEMParameters demParameters{
      this->elasticStiffness,     this->normalViscosity,     this->frictionViscosity,    this->rollingViscosity,
      this->torsionViscosity,     this->staticFrictionCoeff, this->dynamicFrictionCoeff, this->rollingFrictionCoeff,
      this->torsionFrictionCoeff, this->heatConductivity,    this->heatGenerationFactor};

  if constexpr (mixing) {
    functor = std::make_unique<FuncType>(this->cutoff, demParameters, particlePropertiesLibrary);
    particlePropertiesLibrary.addSiteType(1, this->mass);
    particlePropertiesLibrary.addDEMParametersToSite(1, this->radius2, this->specificHeat);
  } else {
    functor = std::make_unique<FuncType>(this->cutoff, demParameters);
    functor->setParticleProperties(this->radius1);
  }
  particlePropertiesLibrary.calculateMixingCoefficients();

  Granular p1(this->startingPosOverlap[0], this->startingVel[0], this->startingAngVel[0], 0, 0);
  Granular p2(this->startingPosOverlap[1], this->startingVel[1], this->startingAngVel[1], 1, (mixing) ? 1 : 0);

  // When
  if (auto msg = this->shouldSkipIfNotImplemented([&]() { functor->AoSFunctor(p1, p2, newton3); }); msg != "") {
    GTEST_SKIP() << msg;
  }

  const auto f1one = p1.getF();
  const auto f2one = p2.getF();
  const auto q1one = p1.getTorque();
  const auto q2one = p2.getTorque();
  const auto hf1one = p1.getHeatFlux();
  const auto hf2one = p2.getHeatFlux();

  // Then
  if (mixing) {
    EXPECT_NEAR(f1one[0], this->expectedForceMixing[0], this->absDelta);
    EXPECT_NEAR(f1one[1], this->expectedForceMixing[1], this->absDelta);
    EXPECT_NEAR(f1one[2], this->expectedForceMixing[2], this->absDelta);

    EXPECT_NEAR(q1one[0], this->expectedTorqueMixing[0], this->absDelta);
    EXPECT_NEAR(q1one[1], this->expectedTorqueMixing[1], this->absDelta);
    EXPECT_NEAR(q1one[2], this->expectedTorqueMixing[2], this->absDelta);

    EXPECT_NEAR(hf1one, this->expectedHeatFluxMixing, this->absDelta);
  } else {
    EXPECT_NEAR(f1one[0], this->expectedForce[0], this->absDelta);
    EXPECT_NEAR(f1one[1], this->expectedForce[1], this->absDelta);
    EXPECT_NEAR(f1one[2], this->expectedForce[2], this->absDelta);

    EXPECT_NEAR(q1one[0], this->expectedTorque[0], this->absDelta);
    EXPECT_NEAR(q1one[1], this->expectedTorque[1], this->absDelta);
    EXPECT_NEAR(q1one[2], this->expectedTorque[2], this->absDelta);

    EXPECT_NEAR(hf1one, this->expectedHeatFlux, this->absDelta);
  }
  if (newton3) {
    if (mixing) {
      EXPECT_NEAR(f2one[0], -this->expectedForceMixing[0], this->absDelta);
      EXPECT_NEAR(f2one[1], -this->expectedForceMixing[1], this->absDelta);
      EXPECT_NEAR(f2one[2], -this->expectedForceMixing[2], this->absDelta);

      EXPECT_NEAR(q2one[0], this->expectedTorqueMixingNewton3[0], this->absDelta);
      EXPECT_NEAR(q2one[1], this->expectedTorqueMixingNewton3[1], this->absDelta);
      EXPECT_NEAR(q2one[2], this->expectedTorqueMixingNewton3[2], this->absDelta);

      EXPECT_NEAR(hf2one, this->expectedHeatFluxMixing, this->absDelta);
    } else {
      EXPECT_NEAR(f2one[0], -this->expectedForce[0], this->absDelta);
      EXPECT_NEAR(f2one[1], -this->expectedForce[1], this->absDelta);
      EXPECT_NEAR(f2one[2], -this->expectedForce[2], this->absDelta);

      EXPECT_NEAR(q2one[0], this->expectedTorqueNewton3[0], this->absDelta);
      EXPECT_NEAR(q2one[1], this->expectedTorqueNewton3[1], this->absDelta);
      EXPECT_NEAR(q2one[2], this->expectedTorqueNewton3[2], this->absDelta);

      EXPECT_NEAR(hf2one, this->expectedHeatFlux, this->absDelta);
    }
  } else {
    EXPECT_DOUBLE_EQ(f2one[0], 0);
    EXPECT_DOUBLE_EQ(f2one[1], 0);
    EXPECT_DOUBLE_EQ(f2one[2], 0);

    EXPECT_DOUBLE_EQ(q2one[0], 0);
    EXPECT_DOUBLE_EQ(q2one[1], 0);
    EXPECT_DOUBLE_EQ(q2one[2], 0);

    EXPECT_DOUBLE_EQ(hf2one, 0);
  }
}

REGISTER_TYPED_TEST_SUITE_P(DEMFunctorTestNoGlobals, testAoSNoGlobals);

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

using MyTypes = ::testing::Types<Newton3True<DEMFunMixGlob>, Newton3False<DEMFunMixGlob>,
                                 Newton3True<DEMFunNoMixGlob>, Newton3False<DEMFunNoMixGlob>>;

INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, DEMFunctorTestNoGlobals, MyTypes);