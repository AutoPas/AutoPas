///**
//* @file DecisionTreeTuningTest.cpp
//* @author
//* Abdulkadir Pazar
//* @date 20.09.2024
//*/
//
//#include "DecisionTreeTuningTest.h"
//
//#include <gmock/gmock.h>
//#include <gtest/gtest.h>
//
//#include "autopas/tuning/tuningStrategy/decisionTreeTuning/DecisionTreeTuning.h"
//#include "autopas/utils/ExceptionHandler.h"
//#include <pybind11/embed.h>
//#include <filesystem>  // For checking the existence of test files
//
//namespace autopas {
//
//namespace py = pybind11;
//
//// Initialize the Python interpreter once for the entire test suite
//struct PythonInterpreter {
// PythonInterpreter() { py::initialize_interpreter(); }
// ~PythonInterpreter() { py::finalize_interpreter(); }
//};
//
//// Ensure the interpreter is initialized globally
//static PythonInterpreter globalPythonInterpreter;
//
///**
//* @class MockLiveInfo
//*
//* Mock class for LiveInfo to simulate system states for testing purposes.
//*/
//class MockLiveInfo : public LiveInfo {
//public:
// MOCK_METHOD((const std::map<std::string, InfoType> &), get, (), (const));
//};
//
///**
//* @test TestScriptLoading
//*
//* Ensures that when a non-existent Python model file is provided, a runtime error is thrown.
//*/
//TEST(DecisionTreeTuningTest, TestScriptLoading) {
// std::set<autopas::Configuration> searchSpace;
//
// EXPECT_THROW(
//     {
//       // Use a path that does not exist
//       autopas::DecisionTreeTuning tuningStrategy(searchSpace, "invalid_model.pkl", 0.8);
//       std::vector<autopas::Configuration> configQueue;
//       autopas::EvidenceCollection evidenceCollection;
//       tuningStrategy.reset(0, 0, configQueue, evidenceCollection);
//     },
//     std::runtime_error);  // Updated to match the exception type thrown for missing files
//}
//
///**
//* @test TestValidPythonResponse
//*
//* Tests that a valid Python model response updates the configuration queue correctly.
//*/
//TEST(DecisionTreeTuningTest, TestValidPythonResponse) {
// std::set<autopas::Configuration> searchSpace = {autopas::Configuration(
//     autopas::ContainerOption::linkedCells, 1.0, autopas::TraversalOption::lc_c08, autopas::LoadEstimatorOption::none,
//     autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled)};
//
// ASSERT_TRUE(std::filesystem::exists("test_model.pkl")) << "Test model file does not exist.";
//
// autopas::DecisionTreeTuning tuningStrategy(searchSpace, "test_model.pkl", 0.8);
//
// MockLiveInfo mockLiveInfo;
// std::map<std::string, LiveInfo::InfoType> liveInfoMap = {
//     {"avgParticlesPerCell", 6.82}, {"maxParticlesPerCell", 33.0},    {"homogeneity", 0.42},
//     {"maxDensity", 1.17},          {"particlesPerCellStdDev", 0.03}, {"threadCount", 1.0}};
//
// EXPECT_CALL(mockLiveInfo, get()).WillRepeatedly(::testing::ReturnRef(liveInfoMap));
//
// tuningStrategy.receiveLiveInfo(mockLiveInfo);
//
// std::vector<autopas::Configuration> configQueue;
// autopas::EvidenceCollection evidenceCollection;
//
// tuningStrategy.reset(0, 0, configQueue, evidenceCollection);
//
// ASSERT_FALSE(configQueue.empty());
// const autopas::Configuration &predictedConfig = configQueue.front();
//
// EXPECT_EQ(predictedConfig.container, autopas::ContainerOption::linkedCells);
// EXPECT_EQ(predictedConfig.traversal, autopas::TraversalOption::lc_c08);
// EXPECT_EQ(predictedConfig.dataLayout, autopas::DataLayoutOption::soa);
// EXPECT_EQ(predictedConfig.newton3, autopas::Newton3Option::enabled);
//}
//
///**
//* @test TestInvalidPythonResponse
//*
//* Tests that a malformed Python response throws an exception during `reset()`.
//*/
//TEST(DecisionTreeTuningTest, TestInvalidPythonResponse) {
// std::set<autopas::Configuration> searchSpace = {autopas::Configuration(
//     autopas::ContainerOption::linkedCells, 1.0, autopas::TraversalOption::lc_c08, autopas::LoadEstimatorOption::none,
//     autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled)};
//
// ASSERT_TRUE(std::filesystem::exists("test_model_invalid_response.pkl")) << "Invalid response test model file does not exist.";
//
// autopas::DecisionTreeTuning tuningStrategy(searchSpace, "test_model_invalid_response.pkl", 0.8);
//
// MockLiveInfo mockLiveInfo;
// std::map<std::string, LiveInfo::InfoType> liveInfoMap = {
//     {"avgParticlesPerCell", 6.82}, {"maxParticlesPerCell", 33.0},    {"homogeneity", 0.42},
//     {"maxDensity", 1.17},          {"particlesPerCellStdDev", 0.03}, {"threadCount", 1.0}};
//
// EXPECT_CALL(mockLiveInfo, get()).WillRepeatedly(::testing::ReturnRef(liveInfoMap));
//
// EXPECT_THROW(
//     {
//       std::vector<autopas::Configuration> configQueue;
//       autopas::EvidenceCollection evidenceCollection;
//       tuningStrategy.reset(0, 0, configQueue, evidenceCollection);
//     },
//     std::runtime_error);  // Updated to match exception type for invalid responses
//}
//
///**
//* @test TestEmptyLiveInfo
//*
//* Tests that an empty live info map does not prevent the model from providing a valid configuration.
//*/
//TEST(DecisionTreeTuningTest, TestEmptyLiveInfo) {
// std::set<autopas::Configuration> searchSpace = {autopas::Configuration(
//     autopas::ContainerOption::linkedCells, 1.0, autopas::TraversalOption::lc_c08, autopas::LoadEstimatorOption::none,
//     autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled)};
//
// ASSERT_TRUE(std::filesystem::exists("test_model.pkl")) << "Test model file does not exist.";
//
// autopas::DecisionTreeTuning tuningStrategy(searchSpace, "test_model.pkl", 0.8);
//
// MockLiveInfo mockLiveInfo;
// std::map<std::string, LiveInfo::InfoType> emptyLiveInfoMap = {};
// EXPECT_CALL(mockLiveInfo, get()).WillRepeatedly(::testing::ReturnRef(emptyLiveInfoMap));
//
// tuningStrategy.receiveLiveInfo(mockLiveInfo);
//
// std::vector<autopas::Configuration> configQueue;
// autopas::EvidenceCollection evidenceCollection;
//
// tuningStrategy.reset(0, 0, configQueue, evidenceCollection);
//
// ASSERT_FALSE(configQueue.empty());
// const autopas::Configuration &predictedConfig = configQueue.front();
//
// EXPECT_EQ(predictedConfig.container, autopas::ContainerOption::linkedCells);
// EXPECT_EQ(predictedConfig.traversal, autopas::TraversalOption::lc_c08);
// EXPECT_EQ(predictedConfig.dataLayout, autopas::DataLayoutOption::soa);
// EXPECT_EQ(predictedConfig.newton3, autopas::Newton3Option::enabled);
//}
//
//}  // namespace autopas
