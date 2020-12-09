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

namespace FeatureVectorTest {

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

  for (const auto &dl : autopas::DataLayoutOption::getAllOptions()) {
    allDataLayouts.push_back(dl);
  }

  for (const auto &n3 : autopas::Newton3Option::getAllOptions()) {
    allNewton3.push_back(n3);
  }
}

/**
 * Check if lhsSampleFeatures generates correct number of samples.
 */
TEST_F(FeatureVectorTest, lhsSampleFeature) {
  autopas::Random rand;
  size_t n = 100;

  FeatureVectorEncoder encoder(allCompatibleContainerTraversalEstimators, allDataLayouts, allNewton3,
                               autopas::NumberInterval<double>(1., 2.));
  auto vecList = encoder.lhsSampleFeatures(n, rand);

  EXPECT_EQ(vecList.size(), n);
}

/**
 * Check if lhsSampleFeaturesCluster generates correct number of samples.
 */
TEST_F(FeatureVectorTest, lhsSampleFeatureCluster) {
  autopas::Random rand;
  size_t n = 100;
  double iteration = 0;

  FeatureVectorEncoder encoder(allCompatibleContainerTraversalEstimators, allDataLayouts, allNewton3,
                               autopas::NumberInterval<double>(1., 2.));
  auto vecList = encoder.lhsSampleFeatureCluster(n, rand, iteration);

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
  FeatureVectorEncoder encoder(allCompatibleContainerTraversalEstimators, allDataLayouts, allNewton3,
                               NumberInterval<double>(0., 1.));
  auto vecList = encoder.lhsSampleFeatures(100, rand);

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
  double iteration = 0.;
  auto cellSizeFactor = 1.0;
  auto dataLayouts = autopas::DataLayoutOption::getAllOptions();
  auto newtons = autopas::Newton3Option::getAllOptions();

  std::vector<DataLayoutOption> dataLayoutsVec(dataLayouts.begin(), dataLayouts.end());
  std::vector<Newton3Option> newtonsVec(newtons.begin(), newtons.end());

  FeatureVectorEncoder encoder(allCompatibleContainerTraversalEstimators, dataLayoutsVec, newtonsVec,
                               NumberSetFinite<double>({cellSizeFactor}));

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
    auto encoded = encoder.convertToCluster(fv, iteration);

    // check expected size of discrete and continuous tuples
    EXPECT_EQ(encoded.first.size(), FeatureVectorEncoder::tunableDiscreteDims);
    EXPECT_EQ(encoded.second.size(), FeatureVectorEncoder::tunableContinuousDims + 1);

    // check if decoding leads to intial vector
    auto decoded = encoder.convertFromCluster(encoded);
    EXPECT_EQ(decoded, fv);
  }
}

/**
 * Check clusterNeighboursManhattan1 generates unique and correct number of neighbours.
 */
TEST_F(FeatureVectorTest, clusterNeighboursManhattan1) {
  double iteration = 42.;
  auto cellSizeFactor = 1.0;
  auto dataLayouts = autopas::DataLayoutOption::getAllOptions();
  auto newtons = autopas::Newton3Option::getAllOptions();

  std::vector<DataLayoutOption> dataLayoutsVec(dataLayouts.begin(), dataLayouts.end());
  std::vector<Newton3Option> newtonsVec(newtons.begin(), newtons.end());

  FeatureVectorEncoder encoder(allCompatibleContainerTraversalEstimators, dataLayoutsVec, newtonsVec,
                               NumberSetFinite<double>({cellSizeFactor}));

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
    auto [encodedDiscrete, encodedContinuous] = encoder.convertToCluster(fv, iteration);
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
  double iteration = 123.;
  auto cellSizeFactor = 1.0;
  auto dataLayouts = autopas::DataLayoutOption::getAllOptions();
  auto newtons = autopas::Newton3Option::getAllOptions();

  std::vector<DataLayoutOption> dataLayoutsVec(dataLayouts.begin(), dataLayouts.end());
  std::vector<Newton3Option> newtonsVec(newtons.begin(), newtons.end());

  FeatureVectorEncoder encoder(allCompatibleContainerTraversalEstimators, dataLayoutsVec, newtonsVec,
                               NumberSetFinite<double>({cellSizeFactor}));

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
    auto [encodedDiscrete, encodedContinuous] = encoder.convertToCluster(fv, iteration);
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

}  // end namespace FeatureVectorTest
