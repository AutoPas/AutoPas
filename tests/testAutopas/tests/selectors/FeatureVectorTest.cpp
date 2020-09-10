/**
 * @file FeatureVectorTest.cpp
 * @author Jan Nguyen
 * @date 11.08.19
 */

#include "FeatureVectorTest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"

using namespace autopas;

FeatureVectorTest::FeatureVectorTest() {
  for (const auto &containerOption : ContainerOption::getAllOptions()) {
    for (const auto &traversalOption : compatibleTraversals::allCompatibleTraversals(containerOption)) {
      for (const auto &loadEstimatorOption : loadEstimators::getApplicableLoadEstimators(
               containerOption, traversalOption, LoadEstimatorOption::getAllOptions())) {
        allCompatibleContainerTraversalEstimators.emplace_back(containerOption, traversalOption, loadEstimatorOption);
      }
    }
  }
}

/**
 * Check if lhsSampleFeatures generates correct number of samples.
 */
TEST_F(FeatureVectorTest, lhsSample) {
  autopas::Random rand;
  size_t n = 100;

  auto vecList = autopas::FeatureVector::lhsSampleFeatures(
      n, rand, autopas::NumberInterval<double>(1., 2.), allCompatibleContainerTraversalEstimators,
      autopas::DataLayoutOption::getAllOptions(), autopas::Newton3Option::getAllOptions());

  EXPECT_EQ(vecList.size(), n);
}

/**
 * Check if lhsSampleFeatures generates correct number of samples.
 */
TEST_F(FeatureVectorTest, lhsSampleContinuous) {
  autopas::Random rand;
  size_t n = 100;

  auto vecList = autopas::FeatureVector::lhsSampleFeatureContinuous(n, rand, autopas::NumberInterval<double>(1., 2.));

  EXPECT_EQ(vecList.size(), n);
}

TEST_F(FeatureVectorTest, distanceTest) {
  autopas::FeatureVector f1(ContainerOption::linkedCells, 1., TraversalOption::lc_c01, LoadEstimatorOption::none,
                            DataLayoutOption::aos, Newton3Option::enabled);
  autopas::FeatureVector f2(ContainerOption::linkedCells, 1., TraversalOption::lc_c08, LoadEstimatorOption::none,
                            DataLayoutOption::aos, Newton3Option::enabled);
  autopas::FeatureVector f3(ContainerOption::linkedCells, 1., TraversalOption::lc_c08, LoadEstimatorOption::none,
                            DataLayoutOption::soa, Newton3Option::enabled);
  autopas::FeatureVector f4(ContainerOption::linkedCells, 1., TraversalOption::lc_c08, LoadEstimatorOption::none,
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

/**
 * Check if applying one-hot encode and decode lead to the initial vector.
 */
TEST_F(FeatureVectorTest, onehot) {
  autopas::Random rand;
  auto dataLayouts = autopas::DataLayoutOption::getAllOptions();
  auto newtons = autopas::Newton3Option::getAllOptions();

  std::vector<DataLayoutOption> dataLayoutsVec(dataLayouts.begin(), dataLayouts.end());
  std::vector<Newton3Option> newtonsVec(newtons.begin(), newtons.end());

  auto vecList =
      autopas::FeatureVector::lhsSampleFeatures(100, rand, autopas::NumberInterval<double>(1., 2.),
                                                allCompatibleContainerTraversalEstimators, dataLayoutsVec, newtonsVec);

  FeatureVectorEncoder encoder(allCompatibleContainerTraversalEstimators, dataLayoutsVec, newtonsVec);

  for (auto fv : vecList) {
    // encode
    auto vec = encoder.oneHotEncode(fv);

    // check if encoded vector has expected size
    ASSERT_EQ(vec.size(), encoder.getOneHotDims());

    // check if decoding leads to the inital vector
    auto decode = encoder.oneHotDecode(vec);
    EXPECT_EQ(decode, fv);
  }
}

/**
 * Check if cluster-encode and decode lead to the initial vector.
 */
TEST_F(FeatureVectorTest, clusterEncode) {
  auto cellSizeFactor = 1.0;
  auto dataLayouts = autopas::DataLayoutOption::getAllOptions();
  auto newtons = autopas::Newton3Option::getAllOptions();

  std::vector<DataLayoutOption> dataLayoutsVec(dataLayouts.begin(), dataLayouts.end());
  std::vector<Newton3Option> newtonsVec(newtons.begin(), newtons.end());

  FeatureVectorEncoder encoder(allCompatibleContainerTraversalEstimators, dataLayoutsVec, newtonsVec);

  // generate all possible combinations
  std::vector<FeatureVector> vecList;
  for (const auto &[container, traversal, estimator] : allCompatibleContainerTraversalEstimators) {
    for (const auto &dataLayout : dataLayouts) {
      for (const auto &newton3 : newtons) {
        vecList.emplace_back(container, cellSizeFactor, traversal, estimator, dataLayout, newton3);
      }
    }
  }

  for (auto fv : vecList) {
    // encode vector
    auto encoded = encoder.convertToCluster(fv);

    // check expected size of discrete and continuous tuples
    EXPECT_EQ(encoded.first.size(), 3);
    EXPECT_EQ(encoded.second.size(), 1);

    // check if decoding leads to intial vector
    auto decoded = encoder.convertFromCluster(encoded);
    EXPECT_EQ(decoded, fv);
  }
}

/**
 * Check clusterNeighboursManhattan1 generates unique and correct number of neighbours.
 */
TEST_F(FeatureVectorTest, clusterNeighboursManhattan1) {
  auto cellSizeFactor = 1.0;
  auto dataLayouts = autopas::DataLayoutOption::getAllOptions();
  auto newtons = autopas::Newton3Option::getAllOptions();

  std::vector<DataLayoutOption> dataLayoutsVec(dataLayouts.begin(), dataLayouts.end());
  std::vector<Newton3Option> newtonsVec(newtons.begin(), newtons.end());

  FeatureVectorEncoder encoder(allCompatibleContainerTraversalEstimators, dataLayoutsVec, newtonsVec);

  std::vector<int> dimRestriction = {static_cast<int>(allCompatibleContainerTraversalEstimators.size()),
                                     static_cast<int>(dataLayouts.size()), static_cast<int>(newtons.size())};

  // generate all possible combinations
  std::vector<FeatureVector> vecList;
  for (auto [container, traversal, estimator] : allCompatibleContainerTraversalEstimators) {
    for (auto dataLayout : dataLayouts) {
      for (auto newton3 : newtons) {
        vecList.emplace_back(container, cellSizeFactor, traversal, estimator, dataLayout, newton3);
      }
    }
  }

  for (auto fv : vecList) {
    // get neighbours of encoded vector
    auto [encodedDiscrete, encodedContinuous] = encoder.convertToCluster(fv);
    auto neighbours = encoder.clusterNeighboursManhattan1(encodedDiscrete);

    // neighbours should contain all container-traversals-estimator + all datalayouts + all newtons - 3 (initial vector
    // is counted trice)
    EXPECT_EQ(neighbours.size(),
              allCompatibleContainerTraversalEstimators.size() + dataLayouts.size() + newtons.size() - 3);

    // neighbours should be unique
    for (size_t i = 0; i < neighbours.size(); ++i) {
      for (size_t j = i + 1; j < neighbours.size(); ++j) {
        EXPECT_NE(neighbours[i].first, neighbours[j].first);
      }
    }
  }
}

/**
 * Check clusterNeighboursManhattan1Container generates unique and correct number of neighbours.
 */
TEST_F(FeatureVectorTest, clusterNeighboursManhattan1Container) {
  auto cellSizeFactor = 1.0;
  auto dataLayouts = autopas::DataLayoutOption::getAllOptions();
  auto newtons = autopas::Newton3Option::getAllOptions();

  std::vector<DataLayoutOption> dataLayoutsVec(dataLayouts.begin(), dataLayouts.end());
  std::vector<Newton3Option> newtonsVec(newtons.begin(), newtons.end());

  FeatureVectorEncoder encoder(allCompatibleContainerTraversalEstimators, dataLayoutsVec, newtonsVec);

  std::vector<int> dimRestriction = {static_cast<int>(allCompatibleContainerTraversalEstimators.size()),
                                     static_cast<int>(dataLayouts.size()), static_cast<int>(newtons.size())};

  // generate all possible combinations
  std::vector<FeatureVector> vecList;
  for (auto [container, traversal, estimator] : allCompatibleContainerTraversalEstimators) {
    for (auto dataLayout : dataLayouts) {
      for (auto newton3 : newtons) {
        vecList.emplace_back(container, cellSizeFactor, traversal, estimator, dataLayout, newton3);
      }
    }
  }

  for (auto fv : vecList) {
    // get neighbours of encoded vector
    auto [encodedDiscrete, encodedContinuous] = encoder.convertToCluster(fv);
    auto neighbours = encoder.clusterNeighboursManhattan1Container(encodedDiscrete);

    // neighbours should contain all container-traversals-estimator + all datalayouts + all newtons - 3 (initial vector
    // is counted trice)
    EXPECT_EQ(neighbours.size(),
              allCompatibleContainerTraversalEstimators.size() + dataLayouts.size() + newtons.size() - 3);

    // neighbours should be unique
    for (size_t i = 0; i < neighbours.size(); ++i) {
      for (size_t j = i + 1; j < neighbours.size(); ++j) {
        EXPECT_NE(neighbours[i].first, neighbours[j].first);
      }
    }
  }
}
