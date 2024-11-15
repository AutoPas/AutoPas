/**
 * @file DEMFunctorTestNoGlobals.cpp
 * @author Joon Kim
 * @date 15.11.2024
 */

#include "DEMFunctorTestNoGlobals.h"

TYPED_TEST_SUITE_P(DEMFunctorTestNoGlobals);

TYPED_TEST_P(DEMFunctorTestNoGlobals, testAoSNoGlobalsOverlap) {
  using FuncType = typename TypeParam::FuncType;
  constexpr bool mixing = FuncType::getMixing();
  constexpr bool newton3 = TypeParam::newton3;

  // Given
  ParticlePropertiesLibrary<double, size_t> particlePropertiesLibrary(this->cutoff);
  std::unique_ptr<FuncType> functor;

  particlePropertiesLibrary.addSiteType(0, this->epsilon, this->sigma, 1.0, this->radius);
  if constexpr (mixing) {
    functor = std::make_unique<FuncType>(
        this->cutoff, this->elasticStiffness, this->adhesiveStiffness, this->frictionStiffness, this->normalViscosity,
        this->frictionViscosity, this->staticFrictionCoeff, this->dynamicFrictionCoeff, particlePropertiesLibrary);
    particlePropertiesLibrary.addSiteType(1, this->epsilon2, this->sigma2, 1.0, this->radius);
  } else {
    functor = std::make_unique<FuncType>(this->cutoff, this->elasticStiffness, this->adhesiveStiffness,
                                         this->frictionStiffness, this->normalViscosity, this->frictionViscosity,
                                         this->staticFrictionCoeff, this->dynamicFrictionCoeff);
    functor->setParticleProperties(this->epsilon * 24, 1, this->radius);
  }
  particlePropertiesLibrary.calculateMixingCoefficients();

  GranularParticle p1(this->startingPosOverlap[0], this->startingVel[0], this->startingAngVel[0], 0, 0);
  GranularParticle p2(this->startingPosOverlap[1], this->startingVel[1], this->startingAngVel[1], 1, (mixing) ? 1 : 0);

  // When
  if (auto msg = this->shouldSkipIfNotImplemented([&]() { functor->AoSFunctor(p1, p2, newton3); }); msg != "") {
    GTEST_SKIP() << msg;
  }

  const auto f1one = p1.getF();
  const auto f2one = p2.getF();

  // Then
  const std::array<double, 3> expectedForce = {
      this->expectedNormalContactForceOverlap[0] + this->expectedNormalVdWForceOverlap[0] +
          this->expectedFrictionForceOverlap[0],
      this->expectedNormalContactForceOverlap[1] + this->expectedNormalVdWForceOverlap[1] +
          this->expectedFrictionForceOverlap[1],
      this->expectedNormalContactForceOverlap[2] + this->expectedNormalVdWForceOverlap[2] +
          this->expectedFrictionForceOverlap[2]};

  const std::array<double, 3> expectedForceMixing = {
      this->expectedNormalContactForceMixingOverlap[0] + this->expectedNormalVdWForceMixingOverlap[0] +
          this->expectedFrictionForceMixingOverlap[0],
      this->expectedNormalContactForceMixingOverlap[1] + this->expectedNormalVdWForceMixingOverlap[1] +
          this->expectedFrictionForceMixingOverlap[1],
      this->expectedNormalContactForceMixingOverlap[2] + this->expectedNormalVdWForceMixingOverlap[2] +
          this->expectedFrictionForceMixingOverlap[2]};

  if (mixing) {
    EXPECT_NEAR(f1one[0], expectedForceMixing[0], this->absDelta);
    EXPECT_NEAR(f1one[1], expectedForceMixing[1], this->absDelta);
    EXPECT_NEAR(f1one[2], expectedForceMixing[2], this->absDelta);
  } else {
    EXPECT_NEAR(f1one[0], expectedForce[0], this->absDelta);
    EXPECT_NEAR(f1one[1], expectedForce[1], this->absDelta);
    EXPECT_NEAR(f1one[2], expectedForce[2], this->absDelta);
  }
  if (newton3) {
    if (mixing) {
      EXPECT_NEAR(f2one[0], -expectedForceMixing[0], this->absDelta);
      EXPECT_NEAR(f2one[1], -expectedForceMixing[1], this->absDelta);
      EXPECT_NEAR(f2one[2], -expectedForceMixing[2], this->absDelta);
    } else {
      EXPECT_NEAR(f2one[0], -expectedForce[0], this->absDelta);
      EXPECT_NEAR(f2one[1], -expectedForce[1], this->absDelta);
      EXPECT_NEAR(f2one[2], -expectedForce[2], this->absDelta);
    }
  } else {
    EXPECT_DOUBLE_EQ(f2one[0], 0);
    EXPECT_DOUBLE_EQ(f2one[1], 0);
    EXPECT_DOUBLE_EQ(f2one[2], 0);
  }

  // When
  functor->AoSFunctor(p2, p1, newton3);

  const auto f1two = p1.getF();
  const auto f2two = p2.getF();

  // Then
  const double factor = newton3 ? 2. : 1.;
  if (mixing) {
    EXPECT_NEAR(f1two[0], factor * expectedForceMixing[0], this->absDelta);
    EXPECT_NEAR(f1two[1], factor * expectedForceMixing[1], this->absDelta);
    EXPECT_NEAR(f1two[2], factor * expectedForceMixing[2], this->absDelta);

    EXPECT_NEAR(f2two[0], -factor * expectedForceMixing[0], this->absDelta);
    EXPECT_NEAR(f2two[1], -factor * expectedForceMixing[1], this->absDelta);
    EXPECT_NEAR(f2two[2], -factor * expectedForceMixing[2], this->absDelta);
  } else {
    EXPECT_NEAR(f1two[0], factor * expectedForce[0], this->absDelta);
    EXPECT_NEAR(f1two[1], factor * expectedForce[1], this->absDelta);
    EXPECT_NEAR(f1two[2], factor * expectedForce[2], this->absDelta);

    EXPECT_NEAR(f2two[0], -factor * expectedForce[0], this->absDelta);
    EXPECT_NEAR(f2two[1], -factor * expectedForce[1], this->absDelta);
    EXPECT_NEAR(f2two[2], -factor * expectedForce[2], this->absDelta);
  }
}

TYPED_TEST_P(DEMFunctorTestNoGlobals, testAoSNoGlobalsNoOverlap) {
  using FuncType = typename TypeParam::FuncType;
  constexpr bool mixing = FuncType::getMixing();
  constexpr bool newton3 = TypeParam::newton3;

  // Given
  ParticlePropertiesLibrary<double, size_t> particlePropertiesLibrary(this->cutoff);
  std::unique_ptr<FuncType> functor;

  particlePropertiesLibrary.addSiteType(0, this->epsilon, this->sigma, 1.0, this->radius);
  if constexpr (mixing) {
    functor = std::make_unique<FuncType>(
        this->cutoff, this->elasticStiffness, this->adhesiveStiffness, this->frictionStiffness, this->normalViscosity,
        this->frictionViscosity, this->staticFrictionCoeff, this->dynamicFrictionCoeff, particlePropertiesLibrary);
    particlePropertiesLibrary.addSiteType(1, this->epsilon2, this->sigma2, 1.0, this->radius);
  } else {
    functor = std::make_unique<FuncType>(this->cutoff, this->elasticStiffness, this->adhesiveStiffness,
                                         this->frictionStiffness, this->normalViscosity, this->frictionViscosity,
                                         this->staticFrictionCoeff, this->dynamicFrictionCoeff);
    functor->setParticleProperties(this->epsilon * 24, 1, this->radius);
  }
  particlePropertiesLibrary.calculateMixingCoefficients();

  GranularParticle p1(this->startingPosNoOverlap[0], this->startingVel[0], this->startingAngVel[0], 0, 0);
  GranularParticle p2(this->startingPosNoOverlap[1], this->startingVel[1], this->startingAngVel[1], 1,
                      (mixing) ? 1 : 0);

  // When
  if (auto msg = this->shouldSkipIfNotImplemented([&]() { functor->AoSFunctor(p1, p2, newton3); }); msg != "") {
    GTEST_SKIP() << msg;
  }

  const auto f1one = p1.getF();
  const auto f2one = p2.getF();

  // Then
  const std::array<double, 3> expectedForce = {
      this->expectedNormalContactForceNoOverlap[0] + this->expectedNormalVdWForceNoOverlap[0] +
          this->expectedFrictionForceNoOverlap[0],
      this->expectedNormalContactForceNoOverlap[1] + this->expectedNormalVdWForceNoOverlap[1] +
          this->expectedFrictionForceNoOverlap[1],
      this->expectedNormalContactForceNoOverlap[2] + this->expectedNormalVdWForceNoOverlap[2] +
          this->expectedFrictionForceNoOverlap[2]};

  const std::array<double, 3> expectedForceMixing = {
      this->expectedNormalContactForceMixingNoOverlap[0] + this->expectedNormalVdWForceMixingNoOverlap[0] +
          this->expectedFrictionForceMixingNoOverlap[0],
      this->expectedNormalContactForceMixingNoOverlap[1] + this->expectedNormalVdWForceMixingNoOverlap[1] +
          this->expectedFrictionForceMixingNoOverlap[1],
      this->expectedNormalContactForceMixingNoOverlap[2] + this->expectedNormalVdWForceMixingNoOverlap[2] +
          this->expectedFrictionForceMixingNoOverlap[2]};

  if (mixing) {
    EXPECT_NEAR(f1one[0], expectedForceMixing[0], this->absDelta);
    EXPECT_NEAR(f1one[1], expectedForceMixing[1], this->absDelta);
    EXPECT_NEAR(f1one[2], expectedForceMixing[2], this->absDelta);
  } else {
    EXPECT_NEAR(f1one[0], expectedForce[0], this->absDelta);
    EXPECT_NEAR(f1one[1], expectedForce[1], this->absDelta);
    EXPECT_NEAR(f1one[2], expectedForce[2], this->absDelta);
  }
  if (newton3) {
    if (mixing) {
      EXPECT_NEAR(f2one[0], -expectedForceMixing[0], this->absDelta);
      EXPECT_NEAR(f2one[1], -expectedForceMixing[1], this->absDelta);
      EXPECT_NEAR(f2one[2], -expectedForceMixing[2], this->absDelta);
    } else {
      EXPECT_NEAR(f2one[0], -expectedForce[0], this->absDelta);
      EXPECT_NEAR(f2one[1], -expectedForce[1], this->absDelta);
      EXPECT_NEAR(f2one[2], -expectedForce[2], this->absDelta);
    }
  } else {
    EXPECT_DOUBLE_EQ(f2one[0], 0);
    EXPECT_DOUBLE_EQ(f2one[1], 0);
    EXPECT_DOUBLE_EQ(f2one[2], 0);
  }

  // When
  functor->AoSFunctor(p2, p1, newton3);

  const auto f1two = p1.getF();
  const auto f2two = p2.getF();

  // Then
  const double factor = newton3 ? 2. : 1.;
  if (mixing) {
    EXPECT_NEAR(f1two[0], factor * expectedForceMixing[0], this->absDelta);
    EXPECT_NEAR(f1two[1], factor * expectedForceMixing[1], this->absDelta);
    EXPECT_NEAR(f1two[2], factor * expectedForceMixing[2], this->absDelta);

    EXPECT_NEAR(f2two[0], -factor * expectedForceMixing[0], this->absDelta);
    EXPECT_NEAR(f2two[1], -factor * expectedForceMixing[1], this->absDelta);
    EXPECT_NEAR(f2two[2], -factor * expectedForceMixing[2], this->absDelta);
  } else {
    EXPECT_NEAR(f1two[0], factor * expectedForce[0], this->absDelta);
    EXPECT_NEAR(f1two[1], factor * expectedForce[1], this->absDelta);
    EXPECT_NEAR(f1two[2], factor * expectedForce[2], this->absDelta);

    EXPECT_NEAR(f2two[0], -factor * expectedForce[0], this->absDelta);
    EXPECT_NEAR(f2two[1], -factor * expectedForce[1], this->absDelta);
    EXPECT_NEAR(f2two[2], -factor * expectedForce[2], this->absDelta);
  }
}

REGISTER_TYPED_TEST_SUITE_P(DEMFunctorTestNoGlobals, testAoSNoGlobalsOverlap, testAoSNoGlobalsNoOverlap);

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

using MyTypes = ::testing::Types<Newton3True<DEMFunMixNoGlob>, Newton3False<DEMFunMixNoGlob>,
                                 Newton3True<DEMFunNoMixNoGlob>, Newton3False<DEMFunNoMixNoGlob>>;

INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, DEMFunctorTestNoGlobals, MyTypes);