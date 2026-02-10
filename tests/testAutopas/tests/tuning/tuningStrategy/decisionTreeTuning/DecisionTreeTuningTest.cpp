/**
 * @file DecisionTreeTuningTest.cpp
 * @author
 * Abdulkadir Pazar
 * @date 20.09.2024
 */

#include "DecisionTreeTuningTest.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
#include <pybind11/embed.h>
#endif

#include <filesystem>  // For checking the existence of test files

#include "autopas/tuning/tuningStrategy/decisionTreeTuning/DecisionTreeTuning.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/options/InteractionTypeOption.h"

namespace autopas {

#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
namespace py = pybind11;
#endif

/**
 * @class MockLiveInfo
 *
 * Mock class for LiveInfo to simulate system states for testing purposes.
 */
class MockLiveInfo : public LiveInfo {
 public:
  MOCK_METHOD((const std::map<std::string, InfoType> &), get, (), (const));
};

/**
 * @test TestScriptLoading
 *
 * Ensures that when a non-existent Python model file is provided, a runtime error is thrown.
 */
TEST(DecisionTreeTuningTest, TestNonExistantModelLoading) {
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  std::set<Configuration> searchSpace;

  EXPECT_THROW(
      {
        autopas::DecisionTreeTuning tuningStrategy(searchSpace, "/nonexistant_model.pkl", 0.8, InteractionTypeOption::pairwise);
        std::vector<autopas::Configuration> configQueue;
        autopas::EvidenceCollection evidenceCollection;
        tuningStrategy.reset(0, 0, configQueue, evidenceCollection);
      },
      utils::ExceptionHandler::AutoPasException);
#else
  GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_PYTHON_BASED_TUNING=OFF";
#endif
}

/**
 * @test TestValidPythonResponse
 *
 * Tests that a valid Python model response updates the configuration queue correctly.
 */
TEST(DecisionTreeTuningTest, TestValidPythonResponse) {
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  // Define the search space for configurations
  std::set<autopas::Configuration> searchSpace = {
      autopas::Configuration(autopas::ContainerOption::verletClusterLists, 1.0, autopas::TraversalOption::lc_c18,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                             autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise)};

  std::string modelPath = "/tests/testAutopas/tests/tuning/tuningStrategy/decisionTreeTuning/test_model.pkl";

  ASSERT_TRUE(std::filesystem::exists(std::string(AUTOPAS_SOURCE_DIR) + modelPath))
      << "Test model file does not exist.";

  autopas::DecisionTreeTuning tuningStrategy(searchSpace, modelPath, 0.8, InteractionTypeOption::pairwise);
  MockLiveInfo mockLiveInfo;
  std::map<std::string, LiveInfo::InfoType> liveInfoMap = {
      {"avgParticlesPerCell", 6.82}, {"maxParticlesPerCell", 33.0},    {"homogeneity", 0.42},
      {"maxDensity", 1.17},          {"particlesPerCellStdDev", 0.03}, {"threadCount", 1.0}};

  EXPECT_CALL(mockLiveInfo, get()).WillRepeatedly(::testing::ReturnRef(liveInfoMap));

  tuningStrategy.receiveLiveInfo(mockLiveInfo);
  std::vector<autopas::Configuration> configQueue;
  configQueue.push_back(autopas::Configuration(autopas::ContainerOption::verletClusterLists, 1.0,
                                               autopas::TraversalOption::lc_c18, autopas::LoadEstimatorOption::none,
                                               autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled,
                                               autopas::InteractionTypeOption::pairwise));

  autopas::EvidenceCollection evidenceCollection;

  tuningStrategy.reset(0, 0, configQueue, evidenceCollection);

  ASSERT_FALSE(configQueue.empty());

  const autopas::Configuration &predictedConfig = configQueue.front();

  EXPECT_EQ(predictedConfig.container, autopas::ContainerOption::linkedCells);
  EXPECT_EQ(predictedConfig.traversal, autopas::TraversalOption::lc_c08);
  EXPECT_EQ(predictedConfig.dataLayout, autopas::DataLayoutOption::soa);
  EXPECT_EQ(predictedConfig.newton3, autopas::Newton3Option::enabled);
#else
  GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_PYTHON_BASED_TUNING=OFF";
#endif
}

/**
 * @test TestInvalidModel
 *
 * Tests that a malformed model throws an exception during `reset()`.
 */
TEST(DecisionTreeTuningTest, TestInvalidModel) {
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  std::set<autopas::Configuration> searchSpace = {autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1.0, autopas::TraversalOption::lc_c08, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled, autopas::InteractionTypeOption::pairwise)};

  std::string modelPath = "/tests/testAutopas/tests/tuning/tuningStrategy/decisionTreeTuning/test_model_invalid.pkl";
  ASSERT_TRUE(std::filesystem::exists(std::string(AUTOPAS_SOURCE_DIR) + modelPath))
      << "Invalid response test model file does not exist.";
  autopas::DecisionTreeTuning tuningStrategy(searchSpace, modelPath, 0.8, InteractionTypeOption::pairwise);

  MockLiveInfo mockLiveInfo;
  std::map<std::string, LiveInfo::InfoType> liveInfoMap = {
      {"avgParticlesPerCell", 6.82}, {"maxParticlesPerCell", 33.0},    {"homogeneity", 0.42},
      {"maxDensity", 1.17},          {"particlesPerCellStdDev", 0.03}, {"threadCount", 1.0}};

  EXPECT_CALL(mockLiveInfo, get()).WillRepeatedly(::testing::ReturnRef(liveInfoMap));

  EXPECT_THROW(
      {
        std::vector<autopas::Configuration> configQueue;
        configQueue.push_back(
            autopas::Configuration(autopas::ContainerOption::verletClusterLists, 1.0, autopas::TraversalOption::lc_c18,
                                   autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                                   autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise));
        autopas::EvidenceCollection evidenceCollection;
        tuningStrategy.reset(0, 0, configQueue, evidenceCollection);
      },
      std::runtime_error);
#else
  GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_PYTHON_BASED_TUNING=OFF";
#endif
}

}  // namespace autopas
