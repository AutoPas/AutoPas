/**
 * @file FeatureVectorTest.cpp
 * @author Jan Nguyen
 * @date 11.08.19
 */

#include "FeatureVectorTest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

using namespace autopas;

TEST_F(FeatureVectorTest, lhsSample) {
  autopas::Random rand(42);
  size_t n = 100;

  auto vecList = autopas::FeatureVector::lhsSampleFeatures(
      n, rand, autopas::NumberInterval<double>(1., 2.), autopas::TraversalOption::getAllOptions(),
      autopas::LoadEstimatorOption::getAllOptions(), autopas::DataLayoutOption::getAllOptions(),
      autopas::Newton3Option::getAllOptions());

  EXPECT_EQ(vecList.size(), n);
}

TEST_F(FeatureVectorTest, distanceTest) {
  autopas::FeatureVector f1(ContainerOption::linkedCells, 1., TraversalOption::c01, LoadEstimatorOption::none,
                            DataLayoutOption::aos, Newton3Option::enabled);
  autopas::FeatureVector f2(ContainerOption::linkedCells, 1., TraversalOption::c08, LoadEstimatorOption::none,
                            DataLayoutOption::aos, Newton3Option::enabled);
  autopas::FeatureVector f3(ContainerOption::linkedCells, 1., TraversalOption::c08, LoadEstimatorOption::none,
                            DataLayoutOption::soa, Newton3Option::enabled);
  autopas::FeatureVector f4(ContainerOption::linkedCells, 1., TraversalOption::c08, LoadEstimatorOption::none,
                            DataLayoutOption::soa, Newton3Option::disabled);

  EXPECT_EQ(static_cast<Eigen::VectorXd>(f1 - f1).squaredNorm(), 0);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f2 - f2).squaredNorm(), 0);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f3 - f3).squaredNorm(), 0);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f4 - f4).squaredNorm(), 0);

  EXPECT_EQ(static_cast<Eigen::VectorXd>(f1 - f2).squaredNorm(), 1);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f2 - f3).squaredNorm(), 1);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f3 - f4).squaredNorm(), 1);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f4 - f3).squaredNorm(), 1);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f3 - f2).squaredNorm(), 1);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f2 - f1).squaredNorm(), 1);

  EXPECT_EQ(static_cast<Eigen::VectorXd>(f1 - f3).squaredNorm(), 2);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f2 - f4).squaredNorm(), 2);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f4 - f2).squaredNorm(), 2);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f3 - f1).squaredNorm(), 2);

  EXPECT_EQ(static_cast<Eigen::VectorXd>(f1 - f4).squaredNorm(), 3);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f4 - f1).squaredNorm(), 3);
}

TEST_F(FeatureVectorTest, onehot) {
  autopas::Random rand(42);
  auto vecList = autopas::FeatureVector::lhsSampleFeatures(
      100, rand, autopas::NumberInterval<double>(1., 2.), autopas::TraversalOption::getAllOptions(),
      autopas::LoadEstimatorOption::getAllOptions(), autopas::DataLayoutOption::getAllOptions(),
      autopas::Newton3Option::getAllOptions());

  for (auto fv : vecList) {
    auto vec = fv.oneHotEncode();
    ASSERT_EQ(vec.size(), autopas::FeatureVector::oneHotDims);

    auto decode = autopas::FeatureVector::oneHotDecode(vec);
    EXPECT_EQ(decode, fv);
  }
}
