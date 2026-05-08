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
    for (const auto &traversalOption :
         compatibleTraversals::allCompatibleTraversals(containerOption, autopas::InteractionTypeOption::pairwise)) {
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

  for (const auto &ompKind : OpenMPKindOption::getAllOptions()) {
    allOMPKinds.push_back(ompKind);
  }
}

/**
 * Check if lhsSampleFeatures generates correct number of samples.
 */
TEST_F(FeatureVectorTest, lhsSampleFeature) {
  autopas::Random rand;
  size_t n = 100;

  FeatureVectorEncoder encoder(allCompatibleContainerTraversalEstimators, allDataLayouts, allNewton3,
                               autopas::NumberInterval<double>(1., 2.), allOMPKinds, NumberSetFinite(std::set<size_t>({1, 2, 250})),
                               InteractionTypeOption::pairwise);
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
                               autopas::NumberInterval<double>(1., 2.), allOMPKinds, NumberSetFinite(std::set<size_t>({1, 2, 250})), InteractionTypeOption::pairwise);
  auto vecList = encoder.lhsSampleFeatureCluster(n, rand, iteration);

  EXPECT_EQ(vecList.size(), n);
}

TEST_F(FeatureVectorTest, distanceTest) {
  autopas::FeatureVector f1(ContainerOption::linkedCells, 1., TraversalOption::lc_c01, LoadEstimatorOption::none,
                            DataLayoutOption::aos, Newton3Option::enabled, OpenMPKindOption::omp_dynamic, 1, InteractionTypeOption::pairwise);
  autopas::FeatureVector f2(ContainerOption::linkedCells, 1., TraversalOption::lc_c08, LoadEstimatorOption::none,
                            DataLayoutOption::aos, Newton3Option::enabled, OpenMPKindOption::omp_dynamic, 1, InteractionTypeOption::pairwise);
  autopas::FeatureVector f3(ContainerOption::linkedCells, 1., TraversalOption::lc_c08, LoadEstimatorOption::none,
                            DataLayoutOption::soa, Newton3Option::enabled, OpenMPKindOption::omp_dynamic, 1,InteractionTypeOption::pairwise);
  autopas::FeatureVector f4(ContainerOption::linkedCells, 1., TraversalOption::lc_c08, LoadEstimatorOption::none,
                            DataLayoutOption::soa, Newton3Option::disabled, OpenMPKindOption::omp_dynamic, 1,InteractionTypeOption::pairwise);

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
                               NumberInterval<double>(0., 1.), allOMPKinds, NumberSetFinite(std::set<size_t>({1, 2, 250})), InteractionTypeOption::pairwise);
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
  const auto ompKinds = OpenMPKindOption::getAllOptions();
  const std::set<size_t> ompChunkSizes{1, 2, 250};

  FeatureVectorEncoder encoder(allCompatibleContainerTraversalEstimators, allDataLayouts, allNewton3,
                               NumberSetFinite<double>({cellSizeFactor}), allOMPKinds, NumberSetFinite(ompChunkSizes), InteractionTypeOption::pairwise);

  // generate all possible combinations
  std::vector<FeatureVector> vecList;
  for (const auto &[container, traversal, estimator] : allCompatibleContainerTraversalEstimators) {
    for (const auto &dataLayout : dataLayouts) {
      for (const auto &newton3 : newtons) {
        for (const auto &ompKind : ompKinds) {
          for (const auto &ompChunkSize : ompChunkSizes) {
            vecList.emplace_back(container, cellSizeFactor, traversal, estimator, dataLayout, newton3, ompKind, ompChunkSize,
                             InteractionTypeOption::pairwise);
          }
        }
      }
    }
  }

  for (auto fv : vecList) {
    // encode vector
    auto encoded = encoder.convertToCluster(fv, iteration);

    // check expected size of discrete and continuous tuples
    EXPECT_EQ(encoded.first.size(), FeatureVectorEncoder::numTunableOneHotIndices);
    EXPECT_EQ(encoded.second.size(), FeatureVectorEncoder::numTunableContinuousDims + 1);

    // check if decoding leads to intial vector
    auto decoded = encoder.convertFromCluster(encoded);
    EXPECT_EQ(decoded, fv);
  }
}

/**
 * Check clusterNeighborsManhattan1 generates unique and correct number of neighbors.
 */
TEST_F(FeatureVectorTest, clusterNeighborsManhattan1) {
  double iteration = 42.;
  auto cellSizeFactor = 1.0;
  auto dataLayouts = autopas::DataLayoutOption::getAllOptions();
  auto newtons = autopas::Newton3Option::getAllOptions();
  const auto ompKinds = OpenMPKindOption::getAllOptions();
  const std::set<size_t> ompChunkSizes{1, 2, 250};

  FeatureVectorEncoder encoder(allCompatibleContainerTraversalEstimators, allDataLayouts, allNewton3,
                               NumberSetFinite<double>({cellSizeFactor}), allOMPKinds, NumberSetFinite(ompChunkSizes), InteractionTypeOption::pairwise);

  std::vector<int> dimRestriction = {static_cast<int>(allCompatibleContainerTraversalEstimators.size()),
                                     static_cast<int>(dataLayouts.size()), static_cast<int>(newtons.size()),
                                     static_cast<int>(ompKinds.size())};

  // generate all possible combinations
  std::vector<FeatureVector> vecList;
  for (const auto &[container, traversal, estimator] : allCompatibleContainerTraversalEstimators) {
    for (const auto &dataLayout : dataLayouts) {
      for (const auto &newton3 : newtons) {
        for (const auto &ompKind : ompKinds) {
          for (const auto &ompChunkSize : ompChunkSizes) {
            vecList.emplace_back(container, cellSizeFactor, traversal, estimator, dataLayout, newton3, ompKind, ompChunkSize,
                             InteractionTypeOption::pairwise);
          }
        }
      }
    }
  }

  for (auto fv : vecList) {
    // get neighbors of encoded vector
    auto [encodedDiscrete, encodedContinuous] = encoder.convertToCluster(fv, iteration);
    auto neighbors = encoder.clusterNeighborsManhattan1(encodedDiscrete);

    // neighbors should contain all container-traversals-estimator + all datalayouts + all newtons + all OMP kinds - 3 (initial vector
    // is counted trice)
    EXPECT_EQ(neighbors.size(),
              allCompatibleContainerTraversalEstimators.size() + dataLayouts.size() + newtons.size() + ompKinds.size() - 3);

    // neighbors should be unique
    for (size_t i = 0; i < neighbors.size(); ++i) {
      for (size_t j = i + 1; j < neighbors.size(); ++j) {
        EXPECT_NE(neighbors[i].first, neighbors[j].first);
      }
    }
  }
}

/**
 * Check clusterNeighborsManhattan1Container generates unique and correct number of neighbors.
 */
TEST_F(FeatureVectorTest, clusterNeighborsManhattan1Container) {
  double iteration = 123.;
  auto cellSizeFactor = 1.0;
  auto dataLayouts = autopas::DataLayoutOption::getAllOptions();
  auto newtons = autopas::Newton3Option::getAllOptions();
  const auto ompKinds = OpenMPKindOption::getAllOptions();
  const std::set<size_t> ompChunkSizes{1, 2, 250};

  FeatureVectorEncoder encoder(allCompatibleContainerTraversalEstimators, allDataLayouts, allNewton3,
                               NumberSetFinite<double>({cellSizeFactor}), allOMPKinds, NumberSetFinite(ompChunkSizes), InteractionTypeOption::pairwise);

  std::vector<int> dimRestriction = {static_cast<int>(allCompatibleContainerTraversalEstimators.size()),
                                     static_cast<int>(dataLayouts.size()), static_cast<int>(newtons.size()),
                                     static_cast<int>(ompKinds.size())};

  // generate all possible combinations
  std::vector<FeatureVector> vecList;
  for (const auto &[container, traversal, estimator] : allCompatibleContainerTraversalEstimators) {
    for (const auto &dataLayout : dataLayouts) {
      for (const auto &newton3 : newtons) {
        for (const auto &ompKind : ompKinds) {
          for (const auto &ompChunkSize : ompChunkSizes) {
            vecList.emplace_back(container, cellSizeFactor, traversal, estimator, dataLayout, newton3, ompKind, ompChunkSize,
                             InteractionTypeOption::pairwise);
          }
        }
      }
    }
  }

  for (auto fv : vecList) {
    // get neighbors of encoded vector
    auto [encodedDiscrete, encodedContinuous] = encoder.convertToCluster(fv, iteration);
    auto neighbors = encoder.clusterNeighborsManhattan1Container(encodedDiscrete);

    // neighbors should contain all container-traversals-estimator + all datalayouts + all newtons - 3 (initial vector
    // is counted trice)
    EXPECT_EQ(neighbors.size(),
              allCompatibleContainerTraversalEstimators.size() + dataLayouts.size() + newtons.size() + ompKinds.size() - 3);

    // neighbours should be unique
    for (size_t i = 0; i < neighbors.size(); ++i) {
      for (size_t j = i + 1; j < neighbors.size(); ++j) {
        EXPECT_NE(neighbors[i].first, neighbors[j].first);
      }
    }
  }
}
