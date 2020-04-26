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
      autopas::DataLayoutOption::getAllOptions(), autopas::Newton3Option::getAllOptions());

  EXPECT_EQ(vecList.size(), n);
}

TEST_F(FeatureVectorTest, distanceTest) {
  autopas::FeatureVector f1(ContainerOption::linkedCells, 1., TraversalOption::c01, DataLayoutOption::aos,
                            Newton3Option::enabled);
  autopas::FeatureVector f2(ContainerOption::linkedCells, 1., TraversalOption::c08, DataLayoutOption::aos,
                            Newton3Option::enabled);
  autopas::FeatureVector f3(ContainerOption::linkedCells, 1., TraversalOption::c08, DataLayoutOption::soa,
                            Newton3Option::enabled);
  autopas::FeatureVector f4(ContainerOption::linkedCells, 1., TraversalOption::c08, DataLayoutOption::soa,
                            Newton3Option::disabled);

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
  // check if applying one-hot encode and decode lead to the initial vector

  autopas::Random rand(42);
  auto vecList = autopas::FeatureVector::lhsSampleFeatures(
      100, rand, autopas::NumberInterval<double>(1., 2.), autopas::TraversalOption::getAllOptions(),
      autopas::DataLayoutOption::getAllOptions(), autopas::Newton3Option::getAllOptions());

  for (auto fv : vecList) {
    // encode
    auto vec = fv.oneHotEncode();

    // check if encoded vector has expected size
    ASSERT_EQ(vec.size(), autopas::FeatureVector::oneHotDims);

    // check if decoding leads to the inital vector
    auto decode = autopas::FeatureVector::oneHotDecode(vec);
    EXPECT_EQ(decode, fv);
  }
}

TEST_F(FeatureVectorTest, clusterEncode) {
  // check if cluster-encode and decode lead to the initial vector

  autopas::Random rand(91);
  auto cellSizes = autopas::NumberInterval<double>(1., 2.);

  auto traversals = autopas::TraversalOption::getAllOptions();
  auto dataLayouts = autopas::DataLayoutOption::getAllOptions();
  auto newtons = autopas::Newton3Option::getAllOptions();

  std::vector<TraversalOption> traversalsVec(traversals.begin(), traversals.end());
  std::vector<DataLayoutOption> dataLayoutsVec(dataLayouts.begin(), dataLayouts.end());
  std::vector<Newton3Option> newtonsVec(newtons.begin(), newtons.end());

  auto vecList =
      autopas::FeatureVector::lhsSampleFeatures(100, rand, cellSizes, traversalsVec, dataLayoutsVec, newtonsVec);

  for (auto fv : vecList) {
    // encode vector
    auto encoded = fv.clusterEncode(traversalsVec, dataLayoutsVec, newtonsVec);

    // check expected size of discrete and continuous tuples
    EXPECT_EQ(encoded.first.size(), 3);
    EXPECT_EQ(encoded.second.size(), 1);

    // check if decoding leads to intial vector
    auto decoded = autopas::FeatureVector::clusterDecode(encoded, traversalsVec, dataLayoutsVec, newtonsVec);
    EXPECT_EQ(decoded, fv);
  }
}

TEST_F(FeatureVectorTest, clusterNeighbours) {
  autopas::Random rand(91);
  auto cellSizes = autopas::NumberInterval<double>(1., 2.);

  auto traversals = autopas::TraversalOption::getAllOptions();
  auto dataLayouts = autopas::DataLayoutOption::getAllOptions();
  auto newtons = autopas::Newton3Option::getAllOptions();

  std::vector<TraversalOption> traversalsVec(traversals.begin(), traversals.end());
  std::vector<DataLayoutOption> dataLayoutsVec(dataLayouts.begin(), dataLayouts.end());
  std::vector<Newton3Option> newtonsVec(newtons.begin(), newtons.end());

  std::vector<int> dimRestriction = {static_cast<int>(traversals.size()), static_cast<int>(dataLayouts.size()),
                                     static_cast<int>(newtons.size())};

  auto vecList =
      autopas::FeatureVector::lhsSampleFeatures(100, rand, cellSizes, traversalsVec, dataLayoutsVec, newtonsVec);

  for (auto fv : vecList) {
    auto encoded = fv.clusterEncode(traversalsVec, dataLayoutsVec, newtonsVec);
    auto neighbours = autopas::FeatureVector::neighboursManhattan1(encoded.first, dimRestriction);

    // neighbours should contain all traversals + all datalayouts + all newtons - 3 (initial vector is counted trice)
    EXPECT_EQ(neighbours.size(), traversals.size() + dataLayouts.size() + newtons.size() - 3);

    // neighbours should be unique
    for (size_t i = 0; i < neighbours.size(); ++i) {
      for (size_t j = i + 1; j < neighbours.size(); ++j) {
        EXPECT_NE(neighbours[i], neighbours[j]);
      }
    }
  }
}
