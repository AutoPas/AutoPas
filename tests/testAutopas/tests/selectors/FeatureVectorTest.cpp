/**
 * @file FeatureVectorTest.cpp
 * @author Jan Nguyen
 * @date 11.08.19
 */

#include "FeatureVectorTest.h"
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

TEST_F(FeatureVectorTest, lhsSample) {
  autopas::Random rand(42);
  size_t n = 100;

  auto vecList = autopas::FeatureVector::lhsSampleFeatures(n, rand, autopas::NumberInterval<double>(1., 2.),
                                                           autopas::allTraversalOptions, autopas::allDataLayoutOptions,
                                                           autopas::allNewton3Options);

  EXPECT_EQ(vecList.size(), n);
}

TEST_F(FeatureVectorTest, onehot) {
  autopas::Random rand(42);
  auto vecList = autopas::FeatureVector::lhsSampleFeatures(100, rand, autopas::NumberInterval<double>(1., 2.),
                                                           autopas::allTraversalOptions, autopas::allDataLayoutOptions,
                                                           autopas::allNewton3Options);

  for (auto fv : vecList) {
    auto vec = fv.oneHotEncode();
    ASSERT_EQ(vec.size(), autopas::FeatureVector::oneHotDims);

    auto decode = autopas::FeatureVector::oneHotDecode(vec);
    EXPECT_EQ(decode, fv);
  }
}
