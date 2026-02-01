/**
 * @file TreeliteBasedDecisionTreeTuningTest.h
 * @author Elizaveta Polysaeva
 * @date 24.01.2026
 */

#include "TreeliteBasedDecisionTreeTuningTest.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <spdlog/sinks/null_sink.h>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <sstream>
#include <string>

#include "autopas/tuning/tuningStrategy/decisionTreeTuning/TreeliteBasedDecisionTreeTuning.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas {

namespace {
/**
 * Return the directory that contains decision tree tuning test artifacts.
 */
std::filesystem::path getTestDir() {
  return std::filesystem::path{AUTOPAS_SOURCE_DIR} / "tests/testAutopas/tests/tuning/tuningStrategy/decisionTreeTuning";
}

/**
 * Register spdlog logger "AutoPasLog".
 */
void ensureNullLogger() {
  if (!spdlog::get("AutoPasLog")) {
    auto null_logger = std::make_shared<spdlog::logger>("AutoPasLog", std::make_shared<spdlog::sinks::null_sink_mt>());
    spdlog::register_logger(null_logger);
  }
}
}  // namespace

/**
 * @test TestMissingModelFile
 *
 * Tests that when a non-existent model file is provided, the constructor throws.
 */
TEST(TreeliteBasedDecisionTreeTuningTest, TestMissingModelFile) {
  // #ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  std::set<Configuration> searchSpace;

  const auto baseDir = getTestDir();
  const std::string missingModelPath = (baseDir / "nonexistent_pairwise.tl").string();

  EXPECT_THROW(
      {
        TreeliteBasedDecisionTreeTuning tuningStrategy(searchSpace, missingModelPath, "" /*triwise*/, 0.0,
                                                       InteractionTypeOption::pairwise);
      },
      autopas::utils::ExceptionHandler::AutoPasException);
  // #else
  //   GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_TREELITE_BASED_TUNING=OFF";
  // #endif
}

/**
 * @test TestEmptyModelPath
 *
 * Tests that an empty model path is rejected.
 */
TEST(TreeliteBasedDecisionTreeTuningTest, TestEmptyModelPath) {
  // #ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  std::set<Configuration> searchSpace;

  EXPECT_THROW(
      {
        TreeliteBasedDecisionTreeTuning tuningStrategy(searchSpace, "" /*pairwise*/, "" /*triwise*/, 0.0,
                                                       InteractionTypeOption::pairwise);
      },
      autopas::utils::ExceptionHandler::AutoPasException);
  // #else
  //   GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_TREELITE_BASED_TUNING=OFF";
  // #endif
}

/**
 * @test TestModelWrongExtension
 *
 * Tests that a model path not ending in ".tl" is rejected by the constructor.
 */
TEST(TreeliteBasedDecisionTreeTuningTest, TestModelWrongExtension) {
  // #ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  std::set<Configuration> searchSpace;

  const auto baseDir = getTestDir();
  const std::string pklModelPath = (baseDir / "test_model.pkl").string();

  EXPECT_THROW(
      {
        TreeliteBasedDecisionTreeTuning tuningStrategy(searchSpace, pklModelPath, "" /*triwise*/, 0.0,
                                                       InteractionTypeOption::pairwise);
      },
      autopas::utils::ExceptionHandler::AutoPasException);
  // #else
  //   GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_TREELITE_BASED_TUNING=OFF";
  // #endif
}

/**
 * @test TestMissingClassesFile
 *
 * Tests that if the derived classes file (modelStem + "_classes.txt") is missing,
 * the constructor throws.
 */
TEST(TreeliteBasedDecisionTreeTuningTest, TestMissingClassesFile) {
  // #ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  std::set<Configuration> searchSpace;

  const auto baseDir = getTestDir() / "treelite_missing_classes_test";
  const auto modelPath = baseDir / "test_model_pairwise.tl";
  const auto featuresPath = baseDir / "features.json";
  const auto classesPath = baseDir / "test_model_pairwise_classes.txt";

  ASSERT_TRUE(std::filesystem::exists(modelPath)) << "Missing model file for test.";
  ASSERT_TRUE(std::filesystem::exists(featuresPath)) << "Missing features.json for test.";
  ASSERT_FALSE(std::filesystem::exists(classesPath)) << "Classes file should be missing for this test.";

  EXPECT_THROW(
      {
        TreeliteBasedDecisionTreeTuning tuningStrategy(searchSpace, modelPath.string(), "" /*triwise*/, 0.8,
                                                       InteractionTypeOption::pairwise);
      },
      autopas::utils::ExceptionHandler::AutoPasException);
  // #else
  //   GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_TREELITE_BASED_TUNING=OFF";
  // #endif
}

/**
 * @test TestMissingFeaturesFile
 *
 * Tests that if features.json is missing next to the model,
 * the constructor throws.
 */
TEST(TreeliteBasedDecisionTreeTuningTest, TestMissingFeaturesFile) {
  // #ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  std::set<Configuration> searchSpace;

  const auto baseDir = getTestDir() / "treelite_missing_features_test";
  const auto modelPath = baseDir / "test_model_pairwise.tl";
  const auto classesPath = baseDir / "test_model_pairwise_classes.txt";
  const auto featuresPath = baseDir / "features.json";

  ASSERT_TRUE(std::filesystem::exists(modelPath)) << "Missing model file for test.";
  ASSERT_TRUE(std::filesystem::exists(classesPath)) << "Missing classes file for test.";
  ASSERT_FALSE(std::filesystem::exists(featuresPath)) << "features.json should be missing for this test.";

  EXPECT_THROW(
      {
        TreeliteBasedDecisionTreeTuning tuningStrategy(searchSpace, modelPath.string(), "" /*triwise*/, 0.0,
                                                       InteractionTypeOption::pairwise);
      },
      autopas::utils::ExceptionHandler::AutoPasException);
  // #else
  //   GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_TREELITE_BASED_TUNING=OFF";
  // #endif
}

/**
 * @test TestValidPrediction
 *
 * Tests that a valid model produces a deterministic prediction and that
 * reset() updates the configuration queue accordingly.
 */
TEST(TreeliteBasedDecisionTreeTuningTest, TestValidPrediction) {
  // #ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  // Avoid failures when logging is enabled.
  ensureNullLogger();

  std::set<Configuration> searchSpace;

  const auto baseDir = getTestDir();
  const auto modelPath = baseDir / "test_model_pairwise.tl";
  const auto classesPath = baseDir / "test_model_pairwise_classes.txt";
  const auto featuresPath = baseDir / "features.json";

  ASSERT_TRUE(std::filesystem::exists(modelPath)) << "Test model file does not exist.";
  ASSERT_TRUE(std::filesystem::exists(classesPath)) << "Test classes file does not exist.";
  ASSERT_TRUE(std::filesystem::exists(featuresPath)) << "Test features file does not exist.";

  TreeliteBasedDecisionTreeTuning tuningStrategy(searchSpace, modelPath.string(), "" /*triwise*/, 0.0,
                                                 InteractionTypeOption::pairwise);

  double meanParticlesPerCell = 1.8148148;
  double medianParticlesPerCell = 0.0;
  double maxParticlesPerCell = 16.0;
  double relativeParticlesPerCellStdDev = 1.8209069;
  double threadCount = 1.0;
  double numCells = 972.0;
  double numEmptyCells = 602.0;
  double skin = 0.2;

  LiveInfo liveInfo;
  std::stringstream ss;
  ss << 8 << ' ' << "meanParticlesPerCell" << ' ' << 1 << ' ' << meanParticlesPerCell << ' ' << "medianParticlesPerCell"
     << ' ' << 1 << ' ' << medianParticlesPerCell << ' ' << "maxParticlesPerCell" << ' ' << 1 << ' '
     << maxParticlesPerCell << ' ' << "relativeParticlesPerCellStdDev" << ' ' << 1 << ' '
     << relativeParticlesPerCellStdDev << ' ' << "threadCount" << ' ' << 1 << ' ' << threadCount << ' ' << "numCells"
     << ' ' << 1 << ' ' << numCells << ' ' << "numEmptyCells" << ' ' << 1 << ' ' << numEmptyCells << ' ' << "skin"
     << ' ' << 1 << ' ' << skin << ' ';
  ss >> liveInfo;
  ASSERT_FALSE(liveInfo.get().empty()) << "LiveInfo deserialization produced an empty map.";

  tuningStrategy.receiveLiveInfo(liveInfo);

  std::vector<Configuration> configQueue;
  configQueue.push_back(Configuration(ContainerOption::verletClusterLists, 1.0, TraversalOption::lc_c18,
                                      LoadEstimatorOption::none, DataLayoutOption::aos, Newton3Option::disabled,
                                      InteractionTypeOption::pairwise));

  EvidenceCollection evidenceCollection;
  tuningStrategy.reset(0, 0, configQueue, evidenceCollection);

  ASSERT_FALSE(configQueue.empty());
  const Configuration &predictedConfig = configQueue.front();

  EXPECT_EQ(predictedConfig.container, ContainerOption::linkedCells);
  EXPECT_EQ(predictedConfig.traversal, TraversalOption::lc_c04_HCP);
  EXPECT_EQ(predictedConfig.dataLayout, DataLayoutOption::soa);
  EXPECT_EQ(predictedConfig.newton3, Newton3Option::enabled);
  EXPECT_EQ(predictedConfig.loadEstimator, LoadEstimatorOption::none);
  EXPECT_DOUBLE_EQ(predictedConfig.cellSizeFactor, 1.0);
  // #else
  //   GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_TREELITE_BASED_TUNING=OFF";
  // #endif
}

/**
 * @test TestConfidenceThresholdSkipsUpdate
 *
 * Tests that if the configured confidence threshold is not met, reset() does not replace
 * the configuration queue with a model prediction.
 */
TEST(TreeliteBasedDecisionTreeTuningTest, TestConfidenceThresholdSkipsUpdate) {
  // #ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  // Avoid failures when logging is enabled.
  ensureNullLogger();

  std::set<Configuration> searchSpace;

  const auto baseDir = getTestDir();
  const auto modelPath = baseDir / "test_model_pairwise.tl";
  const auto classesPath = baseDir / "test_model_pairwise_classes.txt";
  const auto featuresPath = baseDir / "features.json";

  ASSERT_TRUE(std::filesystem::exists(modelPath)) << "Test model file does not exist.";
  ASSERT_TRUE(std::filesystem::exists(classesPath)) << "Test classes file does not exist.";
  ASSERT_TRUE(std::filesystem::exists(featuresPath)) << "Test features file does not exist.";

  TreeliteBasedDecisionTreeTuning tuningStrategy(searchSpace, modelPath.string(), "" /*triwise*/, 1.1,
                                                 InteractionTypeOption::pairwise);

  double meanParticlesPerCell = 1.8148148;
  double medianParticlesPerCell = 0.0;
  double maxParticlesPerCell = 16.0;
  double relativeParticlesPerCellStdDev = 1.8209069;
  double threadCount = 1.0;
  double numCells = 972.0;
  double numEmptyCells = 602.0;
  double skin = 0.2;

  LiveInfo liveInfo;
  std::stringstream ss;
  ss << 8 << ' ' << "meanParticlesPerCell" << ' ' << 1 << ' ' << meanParticlesPerCell << ' ' << "medianParticlesPerCell"
     << ' ' << 1 << ' ' << medianParticlesPerCell << ' ' << "maxParticlesPerCell" << ' ' << 1 << ' '
     << maxParticlesPerCell << ' ' << "relativeParticlesPerCellStdDev" << ' ' << 1 << ' '
     << relativeParticlesPerCellStdDev << ' ' << "threadCount" << ' ' << 1 << ' ' << threadCount << ' ' << "numCells"
     << ' ' << 1 << ' ' << numCells << ' ' << "numEmptyCells" << ' ' << 1 << ' ' << numEmptyCells << ' ' << "skin"
     << ' ' << 1 << ' ' << skin << ' ';
  ss >> liveInfo;
  ASSERT_FALSE(liveInfo.get().empty()) << "LiveInfo deserialization produced an empty map.";

  tuningStrategy.receiveLiveInfo(liveInfo);

  std::vector<Configuration> configQueue;
  Configuration initialConfig(ContainerOption::verletClusterLists, 1.0, TraversalOption::lc_c18,
                              LoadEstimatorOption::none, DataLayoutOption::aos, Newton3Option::disabled,
                              InteractionTypeOption::pairwise);
  configQueue.push_back(initialConfig);

  EvidenceCollection evidenceCollection;
  tuningStrategy.reset(0, 0, configQueue, evidenceCollection);

  ASSERT_FALSE(configQueue.empty());
  const Configuration &predictedConfig = configQueue.front();

  EXPECT_EQ(predictedConfig.container, initialConfig.container);
  EXPECT_EQ(predictedConfig.traversal, initialConfig.traversal);
  EXPECT_EQ(predictedConfig.dataLayout, initialConfig.dataLayout);
  EXPECT_EQ(predictedConfig.newton3, initialConfig.newton3);
  // #else
  //   GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_TREELITE_BASED_TUNING=OFF";
  // #endif
}

/**
 * @test TestInvalidModelFile
 *
 * Tests that an invalid ".tl" file is rejected.
 */
TEST(TreeliteBasedDecisionTreeTuningTest, TestInvalidModelFile) {
  // #ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  std::set<Configuration> searchSpace;

  const auto baseDir = getTestDir() / "treelite_invalid_model_test";
  const auto modelPath = baseDir / "test_model_pairwise.tl";
  const auto classesPath = baseDir / "test_model_pairwise_classes.txt";
  const auto featuresPath = baseDir / "features.json";

  ASSERT_TRUE(std::filesystem::exists(modelPath)) << "Invalid test model file does not exist.";
  ASSERT_TRUE(std::filesystem::exists(classesPath)) << "Invalid test classes file does not exist.";
  ASSERT_TRUE(std::filesystem::exists(featuresPath)) << "Test features file does not exist.";

  EXPECT_THROW(
      {
        TreeliteBasedDecisionTreeTuning tuningStrategy(searchSpace, modelPath.string(), "" /*triwise*/, 0.8,
                                                       InteractionTypeOption::pairwise);
      },
      autopas::utils::ExceptionHandler::AutoPasException);
  // #else
  //   GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_TREELITE_BASED_TUNING=OFF";
  // #endif
}

/**
 * @test TestInvalidClassSize
 *
 * Tests that the constructor throws when the number of lables in a class is not as expected.
 */
TEST(TreeliteBasedDecisionTreeTuningTest, TestInvalidClassSize) {
  // #ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  std::set<Configuration> searchSpace;

  const auto baseDir = getTestDir() / "treelite_invalid_classes_test" / "treelite_invalid_class_size_test";
  const auto modelPath = baseDir / "test_model_pairwise.tl";
  const auto classesPath = baseDir / "test_model_pairwise_classes.txt";
  const auto featuresPath = baseDir / "features.json";

  ASSERT_TRUE(std::filesystem::exists(modelPath)) << "Test model file does not exist.";
  ASSERT_TRUE(std::filesystem::exists(classesPath)) << "Test classes file does not exist.";
  ASSERT_TRUE(std::filesystem::exists(featuresPath)) << "Test features file does not exist.";

  EXPECT_THROW(
      {
        TreeliteBasedDecisionTreeTuning tuningStrategy(searchSpace, modelPath.string(), "" /*triwise*/, 0.0,
                                                       InteractionTypeOption::pairwise);
      },
      autopas::utils::ExceptionHandler::AutoPasException);
  // #else
  //   GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_TREELITE_BASED_TUNING=OFF";
  // #endif
}

/**
 * @test TestInvalidLabelInClass
 *
 * Tests that the constructor throws when an invalid label is detected in a class.
 */
TEST(TreeliteBasedDecisionTreeTuningTest, TestInvalidLabelInClass) {
  // #ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  std::set<Configuration> searchSpace;

  const auto baseDir = getTestDir() / "treelite_invalid_classes_test" / "treelite_invalid_label_test";
  const auto modelPath = baseDir / "test_model_pairwise.tl";
  const auto classesPairwise = baseDir / "test_model_pairwise_classes.txt";
  const auto featuresPath = baseDir / "features.json";

  ASSERT_TRUE(std::filesystem::exists(modelPath)) << "Test model file does not exist.";
  ASSERT_TRUE(std::filesystem::exists(classesPairwise)) << "Test classes file does not exist.";
  ASSERT_TRUE(std::filesystem::exists(featuresPath)) << "Test features file does not exist.";

  EXPECT_THROW(
      {
        TreeliteBasedDecisionTreeTuning tuningStrategy(searchSpace, modelPath.string(), "" /*triwise*/, 0.0,
                                                       InteractionTypeOption::pairwise);
      },
      autopas::utils::ExceptionHandler::AutoPasException);
  // #else
  //   GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_TREELITE_BASED_TUNING=OFF";
  // #endif
}

/**
 * @test TestInvalidFeaturesSize
 *
 * Tests that the constructor throws when the number of features is not as expected.
 */
TEST(TreeliteBasedDecisionTreeTuningTest, TestInvalidFeaturesSize) {
  // #ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  std::set<Configuration> searchSpace;

  const auto baseDir = getTestDir() / "treelite_invalid_features_test" / "treelite_invalid_feature_list_size_test";
  const auto modelPath = baseDir / "test_model_pairwise.tl";
  const auto classesPairwise = baseDir / "test_model_pairwise_classes.txt";
  const auto featuresPath = baseDir / "features.json";

  ASSERT_TRUE(std::filesystem::exists(modelPath)) << "Test model file does not exist.";
  ASSERT_TRUE(std::filesystem::exists(classesPairwise)) << "Test classes file does not exist.";
  ASSERT_TRUE(std::filesystem::exists(featuresPath)) << "Test features file does not exist.";

  EXPECT_THROW(
      {
        TreeliteBasedDecisionTreeTuning tuningStrategy(searchSpace, modelPath.string(), "" /*triwise*/, 0.0,
                                                       InteractionTypeOption::pairwise);
      },
      autopas::utils::ExceptionHandler::AutoPasException);
  // #else
  //   GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_TREELITE_BASED_TUNING=OFF";
  // #endif
}

/**
 * @test TestInvalidFeatures
 *
 * Tests that the constructor throws when features.json doesn't match with known list of features.
 */
TEST(TreeliteBasedDecisionTreeTuningTest, TestInvalidFeatures) {
  // #ifdef AUTOPAS_ENABLE_TREELITE_BASED_TUNING
  std::set<Configuration> searchSpace;

  const auto baseDir = getTestDir() / "treelite_invalid_features_test" / "treelite_invalid_feature_list_test";
  const auto modelPath = baseDir / "test_model_pairwise.tl";
  const auto classesPairwise = baseDir / "test_model_pairwise_classes.txt";
  const auto featuresPath = baseDir / "features.json";

  ASSERT_TRUE(std::filesystem::exists(modelPath)) << "Test model file does not exist.";
  ASSERT_TRUE(std::filesystem::exists(classesPairwise)) << "Test classes file does not exist.";
  ASSERT_TRUE(std::filesystem::exists(featuresPath)) << "Test features file does not exist.";

  EXPECT_THROW(
      {
        TreeliteBasedDecisionTreeTuning tuningStrategy(searchSpace, modelPath.string(), "" /*triwise*/, 0.0,
                                                       InteractionTypeOption::pairwise);
      },
      autopas::utils::ExceptionHandler::AutoPasException);
  // #else
  //   GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_TREELITE_BASED_TUNING=OFF";
  // #endif
}

}  // namespace autopas
