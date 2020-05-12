/**
 * @file FeatureVectorTest.cpp
 * @author Jan Nguyen
 * @date 11.08.19
 */

#include "FeatureVectorTest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

using namespace autopas;

/**
 * Check if correct number of samples is generated.
 */
TEST_F(FeatureVectorTest, lhsSample) {
  autopas::Random rand;
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

/**
 * Check if applying one-hot encode and decode lead to the initial vector.
 */
TEST_F(FeatureVectorTest, onehot) {
  autopas::Random rand;
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

/**
 * Check if cluster-encode and decode lead to the initial vector.
 */
TEST_F(FeatureVectorTest, clusterEncode) {
  auto container = autopas::ContainerOption::linkedCells;
  auto cellSizeFactor = 1.0;
  auto traversals = autopas::TraversalOption::getAllOptions();
  auto dataLayouts = autopas::DataLayoutOption::getAllOptions();
  auto newtons = autopas::Newton3Option::getAllOptions();

  std::vector<std::pair<ContainerOption, TraversalOption>> containerTraversalsVec;
  containerTraversalsVec.reserve(traversals.size());
  for (auto traversal : traversals) {
    containerTraversalsVec.emplace_back(container, traversal);
  }
  std::vector<DataLayoutOption> dataLayoutsVec(dataLayouts.begin(), dataLayouts.end());
  std::vector<Newton3Option> newtonsVec(newtons.begin(), newtons.end());

  // generate all possible combinations
  std::vector<FeatureVector> vecList;
  for (auto [container, traversal] : containerTraversalsVec) {
    for (auto dataLayout : dataLayouts) {
      for (auto newton3 : newtons) {
        vecList.emplace_back(container, cellSizeFactor, traversal, dataLayout, newton3);
      }
    }
  }

  for (auto fv : vecList) {
    fv.container = container;

    // encode vector
    auto encoded = fv.clusterEncode(containerTraversalsVec, dataLayoutsVec, newtonsVec);

    // check expected size of discrete and continuous tuples
    EXPECT_EQ(encoded.first.size(), 3);
    EXPECT_EQ(encoded.second.size(), 1);

    // check if decoding leads to intial vector
    auto decoded = autopas::FeatureVector::clusterDecode(encoded, containerTraversalsVec, dataLayoutsVec, newtonsVec);
    EXPECT_EQ(decoded, fv);
  }
}

/**
 * Check neighboursManhattan1 generates unique and correct number of neighbours.
 */
TEST_F(FeatureVectorTest, clusterNeighboursManhattan1) {
  auto container = autopas::ContainerOption::linkedCells;
  auto cellSizeFactor = 1.0;
  auto traversals = autopas::TraversalOption::getAllOptions();
  auto dataLayouts = autopas::DataLayoutOption::getAllOptions();
  auto newtons = autopas::Newton3Option::getAllOptions();

  std::vector<std::pair<ContainerOption, TraversalOption>> containerTraversalsVec;
  containerTraversalsVec.reserve(traversals.size());
  for (auto traversal : traversals) {
    containerTraversalsVec.emplace_back(container, traversal);
  }
  std::vector<DataLayoutOption> dataLayoutsVec(dataLayouts.begin(), dataLayouts.end());
  std::vector<Newton3Option> newtonsVec(newtons.begin(), newtons.end());

  std::vector<int> dimRestriction = {static_cast<int>(traversals.size()), static_cast<int>(dataLayouts.size()),
                                     static_cast<int>(newtons.size())};

  // generate all possible combinations
  std::vector<FeatureVector> vecList;
  for (auto [container, traversal] : containerTraversalsVec) {
    for (auto dataLayout : dataLayouts) {
      for (auto newton3 : newtons) {
        vecList.emplace_back(container, cellSizeFactor, traversal, dataLayout, newton3);
      }
    }
  }

  for (auto fv : vecList) {
    // get neighbours of encoded vector
    auto [encodedDiscrete, encodedContinuous] = fv.clusterEncode(containerTraversalsVec, dataLayoutsVec, newtonsVec);
    auto neighbours = autopas::FeatureVector::neighboursManhattan1(encodedDiscrete, dimRestriction);

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
