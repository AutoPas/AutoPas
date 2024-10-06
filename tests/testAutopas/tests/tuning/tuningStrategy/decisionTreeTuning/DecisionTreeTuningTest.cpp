/**
 * @file DecisionTreeTuningTest.cpp
 * @author
 * Abdulkadir Pazar
 * @date 20.09.2024
 */

#include "DecisionTreeTuningTest.h"

#include <Python.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "autopas/tuning/searchSpace/Evidence.h"
#include "autopas/tuning/tuningStrategy/decisionTreeTuning/DecisionTreeTuning.h"

namespace autopas {

/**
 * @class MockLiveInfo
 *
 * The MockLiveInfo class inherits from LiveInfo and provides a mock implementation
 * of the get() method to simulate live info for DecisionTreeTuning tests.
 */
class MockLiveInfo : public LiveInfo {
public:
 MOCK_METHOD((const std::map<std::string, InfoType>&), get, (), (const));
};

/**
 * @test TestScriptLoading
 *
 * This test case ensures that when an invalid or non-existent Python model file is provided
 * to the DecisionTreeTuning constructor, a runtime error is thrown during the `reset()` call.
 */
TEST(DecisionTreeTuningTest, TestScriptLoading) {
 std::set<autopas::Configuration> searchSpace;

 EXPECT_THROW({
   autopas::DecisionTreeTuning tuningStrategy(searchSpace, "invalid_model.pkl");
   std::vector<autopas::Configuration> configQueue;
   autopas::EvidenceCollection evidenceCollection;
   tuningStrategy.reset(0, 0, configQueue, evidenceCollection);
 }, std::runtime_error);
}

/**
 * @test TestValidPythonResponse
 *
 * This test case simulates a valid configuration where the Python script successfully
 * returns the expected configuration, and the reset function updates the configuration queue.
 */
TEST(DecisionTreeTuningTest, TestValidPythonResponse) {
 std::set<autopas::Configuration> searchSpace = {
     autopas::Configuration(autopas::ContainerOption::linkedCells, 1.0, autopas::TraversalOption::lc_c08,
                            autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                            autopas::Newton3Option::enabled)
 };

 autopas::DecisionTreeTuning tuningStrategy(searchSpace, "test_model.pkl");


 MockLiveInfo mockLiveInfo;
 std::map<std::string, LiveInfo::InfoType> liveInfoMap = {
     {"avgParticlesPerCell", 6.82},
     {"maxParticlesPerCell", 33.0},
     {"homogeneity", 0.42},
     {"maxDensity", 1.17},
     {"particlesPerCellStdDev", 0.03},
     {"threadCount", 1.0}
 };

 EXPECT_CALL(mockLiveInfo, get()).WillRepeatedly(::testing::ReturnRef(liveInfoMap));

 tuningStrategy.receiveLiveInfo(mockLiveInfo);

 std::vector<autopas::Configuration> configQueue;
 autopas::EvidenceCollection evidenceCollection;

 tuningStrategy.reset(0, 0, configQueue, evidenceCollection);

 ASSERT_FALSE(configQueue.empty());
 const autopas::Configuration &predictedConfig = configQueue.front();

 EXPECT_EQ(predictedConfig.container, autopas::ContainerOption::linkedCells);
 EXPECT_EQ(predictedConfig.cellSizeFactor, 1.0);
 EXPECT_EQ(predictedConfig.traversal, autopas::TraversalOption::lc_c08);
 EXPECT_EQ(predictedConfig.loadEstimator, autopas::LoadEstimatorOption::none);
 EXPECT_EQ(predictedConfig.dataLayout, autopas::DataLayoutOption::soa);
 EXPECT_EQ(predictedConfig.newton3, autopas::Newton3Option::enabled);
}

/**
 * @test TestInvalidPythonResponse
 *
 * This test case verifies that when the Python script returns an invalid response (e.g., malformed JSON),
 * the `reset()` function throws a runtime error during the parsing of the configuration.
 */
TEST(DecisionTreeTuningTest, TestInvalidPythonResponse) {
 std::set<autopas::Configuration> searchSpace = {
     autopas::Configuration(autopas::ContainerOption::linkedCells, 1.0, autopas::TraversalOption::lc_c08,
                            autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                            autopas::Newton3Option::enabled)
 };

 autopas::DecisionTreeTuning tuningStrategy(searchSpace, "test_model_invalid_response.pkl");

 MockLiveInfo mockLiveInfo;
 std::map<std::string, LiveInfo::InfoType> liveInfoMap = {
     {"avgParticlesPerCell", 6.82},
     {"maxParticlesPerCell", 33.0},
     {"homogeneity", 0.42},
     {"maxDensity", 1.17},
     {"particlesPerCellStdDev", 0.03},
     {"threadCount", 1.0}
 };

 EXPECT_CALL(mockLiveInfo, get()).WillRepeatedly(::testing::ReturnRef(liveInfoMap));

 EXPECT_THROW({
   std::vector<autopas::Configuration> configQueue;
   autopas::EvidenceCollection evidenceCollection;
   tuningStrategy.reset(0, 0, configQueue, evidenceCollection);
 }, std::runtime_error);  // Expect an exception due to invalid JSON
}

/**
 * @test TestEmptyLiveInfo
 *
 * This test case verifies that even when the live info map is empty, the `reset()` function
 * can still successfully retrieve a configuration from the Python model.
 */
TEST(DecisionTreeTuningTest, TestEmptyLiveInfo) {
 std::set<autopas::Configuration> searchSpace = {
     autopas::Configuration(autopas::ContainerOption::linkedCells, 1.0, autopas::TraversalOption::lc_c08,
                            autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                            autopas::Newton3Option::enabled)
 };

 autopas::DecisionTreeTuning tuningStrategy(searchSpace, "test_model.pkl");

 MockLiveInfo mockLiveInfo;
 std::map<std::string, LiveInfo::InfoType> emptyLiveInfoMap = {};
 EXPECT_CALL(mockLiveInfo, get()).WillRepeatedly(::testing::ReturnRef(emptyLiveInfoMap));

 tuningStrategy.receiveLiveInfo(mockLiveInfo);

 std::vector<autopas::Configuration> configQueue;
 autopas::EvidenceCollection evidenceCollection;

 tuningStrategy.reset(0, 0, configQueue, evidenceCollection);

 ASSERT_FALSE(configQueue.empty());
 const autopas::Configuration &predictedConfig = configQueue.front();

 EXPECT_EQ(predictedConfig.container, autopas::ContainerOption::linkedCells);
 EXPECT_EQ(predictedConfig.cellSizeFactor, 1.0);
 EXPECT_EQ(predictedConfig.traversal, autopas::TraversalOption::lc_c08);
 EXPECT_EQ(predictedConfig.loadEstimator, autopas::LoadEstimatorOption::none);
 EXPECT_EQ(predictedConfig.dataLayout, autopas::DataLayoutOption::soa);
 EXPECT_EQ(predictedConfig.newton3, autopas::Newton3Option::enabled);
}

}  // namespace autopas
