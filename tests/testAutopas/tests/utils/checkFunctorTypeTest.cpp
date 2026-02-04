/**
 * @file checkFunctorTypeTest.cpp
 * @author muehlhaeusser
 * @date 19.09.24
 */

#include <gtest/gtest.h>

#include "autopas/utils/checkFunctorType.h"
#include "testingHelpers/commonTypedefs.h"

using namespace autopas;

// Define different kind of Functors for testing
class PairwiseTestFunctor : public PairwiseFunctor<ParticleFP64, PairwiseTestFunctor> {
 public:
  PairwiseTestFunctor() : PairwiseFunctor(1.0){};
  std::string getName() override { return "PairwiseTestFunctor"; };
  bool allowsNewton3() override { return true; };
  bool allowsNonNewton3() override { return true; };
  bool isRelevantForTuning() override { return true; };
};
class TriwiseTestFunctor : public TriwiseFunctor<ParticleFP64, TriwiseTestFunctor> {
 public:
  TriwiseTestFunctor() : TriwiseFunctor(1.0){};
  std::string getName() override { return "TriwiseTestFunctor"; };
  bool allowsNewton3() override { return true; };
  bool allowsNonNewton3() override { return true; };
  bool isRelevantForTuning() override { return true; };
};

class ChildPairwiseTestFunctor : public PairwiseTestFunctor {};
class ChildTriwiseTestFunctor : public PairwiseTestFunctor {};
class BaseTestFunctor : public Functor<ParticleFP64, BaseTestFunctor> {};
class InvalidFunctor {};

// Test that only a pairwise functor is recognized as PairwiseFunctor
TEST(checkFunctorTypeTest, testIsPairwiseFunctor) {
  // This should be true because PairwiseTestFunctor inherits from PairwiseFunctor
  EXPECT_TRUE((std::is_same_v<utils::isPairwiseFunctor<PairwiseTestFunctor>, std::true_type>));

  // This should be true because ChildPairwiseTestFunctor inherits from PairwiseTestFunctor which inherits from
  // PairwiseFunctor
  EXPECT_TRUE((std::is_same_v<utils::isPairwiseFunctor<ChildPairwiseTestFunctor>, std::true_type>));

  // This should be false because TriwiseTestFunctor does not inherit from PairwiseFunctor
  EXPECT_TRUE((std::is_same_v<utils::isPairwiseFunctor<TriwiseTestFunctor>, std::false_type>));

  // This should be false because BaseTestFunctor does not inherit from PairwiseFunctor
  EXPECT_TRUE((std::is_same_v<utils::isPairwiseFunctor<BaseTestFunctor>, std::false_type>));

  // This should be false because InvalidFunctor does not inherit from PairwiseFunctor
  EXPECT_TRUE((std::is_same_v<utils::isPairwiseFunctor<InvalidFunctor>, std::false_type>));

  // This should be false because int is not related to PairwiseFunctor
  EXPECT_TRUE((std::is_same_v<utils::isPairwiseFunctor<int>, std::false_type>));

  // Test by infering the type from a constructed object
  auto pairTestFunctor = PairwiseTestFunctor();
  EXPECT_TRUE((std::is_same_v<utils::isPairwiseFunctor<decltype(pairTestFunctor)>, std::true_type>));
}

// Test that only a triwise functor is recognized as TriwiseFunctor
TEST(checkFunctorTypeTest, testIsTriwiseFunctor) {
  // This should be true because TriwiseTestFunctor inherits from TriwiseFunctor
  EXPECT_TRUE((std::is_same_v<utils::isTriwiseFunctor<TriwiseTestFunctor>, std::true_type>));

  // This should be true because ChildTriwiseTestFunctor inherits from TriwiseTestFunctor which inherits from
  // TriwiseFunctor
  EXPECT_TRUE((std::is_same_v<utils::isPairwiseFunctor<ChildTriwiseTestFunctor>, std::true_type>));

  // This should be false because PairwiseTestFunctor does not inherit from TriwiseFunctor
  EXPECT_TRUE((std::is_same_v<utils::isTriwiseFunctor<PairwiseTestFunctor>, std::false_type>));

  // This should be false because BaseTestFunctor does not inherit from TriwiseFunctor
  EXPECT_TRUE((std::is_same_v<utils::isTriwiseFunctor<BaseTestFunctor>, std::false_type>));

  // This should be false because InvalidFunctor does not inherit from TriwiseFunctor
  EXPECT_TRUE((std::is_same_v<utils::isTriwiseFunctor<InvalidFunctor>, std::false_type>));

  // This should be false because int is not related to TriwiseFunctor
  EXPECT_TRUE((std::is_same_v<utils::isTriwiseFunctor<int>, std::false_type>));

  // Test by infering the type from a constructed object
  auto triTestFunctor = TriwiseTestFunctor();
  EXPECT_TRUE((std::is_same_v<utils::isTriwiseFunctor<decltype(triTestFunctor)>, std::true_type>));
}
