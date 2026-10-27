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
#include <pybind11/stl.h>
#endif

#include <filesystem>
#include <memory>
#include <set>
#include <string>

#include "autopas/options/InteractionTypeOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/tuningStrategy/decisionTreeTuning/DecisionTreeTuning.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/StringUtils.h"
#include "testingHelpers/ArbitraryConfigurations.h"

namespace autopas {

#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
namespace py = pybind11;

/**
 * Imports generate_test_models.py, which writes the dummy models for these tests. Needs a running Python interpreter.
 * @return The imported module.
 */
py::module importTestModelGenerator() {
  py::module::import("sys").attr("path").attr("append")(
      std::string(AUTOPAS_SOURCE_DIR) + "/tests/testAutopas/tests/tuning/tuningStrategy/decisionTreeTuning");
  return py::module::import("generate_test_models");
}

/**
 * Get the directory for the dummy models of the current test, and create it if needed. Each test gets its own, so that
 * tests running in parallel do not overwrite each other's models.
 * @return Path of the directory.
 */
std::string getTestModelDir() {
  const auto dir =
      std::filesystem::temp_directory_path() /
      (std::string("AutoPasDecisionTreeTuningTest_") + ::testing::UnitTest::GetInstance()->current_test_info()->name());
  std::filesystem::create_directories(dir);
  return dir.string();
}
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
        autopas::DecisionTreeTuning tuningStrategy(searchSpace, "/nonexistant_model.pkl", 0.8,
                                                   InteractionTypeOption::pairwise);
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
  const auto pythonInterpreter = utils::PythonInterpreter::getSharedInterpreter();
  const auto modelDir = getTestModelDir();
  const auto predictedConfig = arbitraryConfigurations::_arbitrary_config_2B_2;
  importTestModelGenerator().attr("create_valid_model")(
      modelDir, utils::StringUtils::tokenize(autopas::Configuration().getCSVHeader(), ","),
      utils::StringUtils::tokenize(predictedConfig.getCSVLine(), ","));

  std::set<autopas::Configuration> searchSpace = {arbitraryConfigurations::_arbitrary_config_2B_0};
  autopas::DecisionTreeTuning tuningStrategy(searchSpace, modelDir + "/test_model.pkl", 0.8,
                                             InteractionTypeOption::pairwise);
  MockLiveInfo mockLiveInfo;
  std::map<std::string, LiveInfo::InfoType> liveInfoMap = {
      {"avgParticlesPerCell", 6.82}, {"maxParticlesPerCell", 33.0},    {"homogeneity", 0.42},
      {"maxDensity", 1.17},          {"particlesPerCellStdDev", 0.03}, {"threadCount", 1.0}};

  EXPECT_CALL(mockLiveInfo, get()).WillRepeatedly(::testing::ReturnRef(liveInfoMap));

  tuningStrategy.receiveLiveInfo(mockLiveInfo);
  std::vector<autopas::Configuration> configQueue{arbitraryConfigurations::_arbitrary_config_2B_0};

  autopas::EvidenceCollection evidenceCollection;

  tuningStrategy.reset(0, 0, configQueue, evidenceCollection);

  ASSERT_EQ(configQueue.size(), 1);
  EXPECT_EQ(configQueue.front(), predictedConfig);
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
  const auto pythonInterpreter = utils::PythonInterpreter::getSharedInterpreter();
  const auto modelDir = getTestModelDir();
  importTestModelGenerator().attr("create_invalid_model")(
      modelDir, utils::StringUtils::tokenize(autopas::Configuration().getCSVHeader(), ","));

  std::set<autopas::Configuration> searchSpace = {arbitraryConfigurations::_arbitrary_config_2B_2};
  autopas::DecisionTreeTuning tuningStrategy(searchSpace, modelDir + "/test_model_invalid.pkl", 0.8,
                                             InteractionTypeOption::pairwise);

  MockLiveInfo mockLiveInfo;
  std::map<std::string, LiveInfo::InfoType> liveInfoMap = {
      {"avgParticlesPerCell", 6.82}, {"maxParticlesPerCell", 33.0},    {"homogeneity", 0.42},
      {"maxDensity", 1.17},          {"particlesPerCellStdDev", 0.03}, {"threadCount", 1.0}};

  EXPECT_CALL(mockLiveInfo, get()).WillRepeatedly(::testing::ReturnRef(liveInfoMap));

  EXPECT_THROW(
      {
        std::vector<autopas::Configuration> configQueue{arbitraryConfigurations::_arbitrary_config_2B_0};
        autopas::EvidenceCollection evidenceCollection;
        tuningStrategy.reset(0, 0, configQueue, evidenceCollection);
      },
      utils::ExceptionHandler::AutoPasException);
#else
  GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_PYTHON_BASED_TUNING=OFF";
#endif
}

/**
 * @test ThrowsWithIncompleteModel
 *
 * Tests that a model which does not predict exactly the configuration components throws.
 * E.g. because it was trained before a component was added/removed to (or otherwise modified) the Configuration class.
 */
TEST(DecisionTreeTuningTest, ThrowsWithIncompleteModel) {
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  const auto pythonInterpreter = utils::PythonInterpreter::getSharedInterpreter();
  const auto modelDir = getTestModelDir();
  const auto configComponents = utils::StringUtils::tokenize(autopas::Configuration().getCSVHeader(), ",");
  const auto configValues =
      utils::StringUtils::tokenize(arbitraryConfigurations::_arbitrary_config_2B_2.getCSVLine(), ",");
  std::set<autopas::Configuration> searchSpace;

  // Writes a model that predicts the given components and values, and tries to construct a DecisionTreeTuning with it,
  // expecting that to throw.
  const auto expectThrowWithModel = [&](const std::vector<std::string> &modelComponents,
                                        const std::vector<std::string> &modelValues) {
    importTestModelGenerator().attr("create_valid_model")(modelDir, modelComponents, modelValues);
    EXPECT_THROW(
        {
          autopas::DecisionTreeTuning tuningStrategy(searchSpace, modelDir + "/test_model.pkl", 0.8,
                                                     InteractionTypeOption::pairwise);
        },
        utils::ExceptionHandler::AutoPasException);
  };

  {
    SCOPED_TRACE("Removed component");
    // Remove a component, e.g. like the model was trained before that component existed
    auto modelComponents = configComponents;
    auto modelValues = configValues;
    modelComponents.pop_back();
    modelValues.pop_back();
    expectThrowWithModel(modelComponents, modelValues);
  }
  {
    SCOPED_TRACE("Added component");
    // Add a fake component, e.g. like the model was trained before a component was removed
    auto modelComponents = configComponents;
    auto modelValues = configValues;
    modelComponents.push_back("FakeComponent");
    modelValues.push_back("FakeValue");
    expectThrowWithModel(modelComponents, modelValues);
  }
  {
    SCOPED_TRACE("Renamed component");
    // Replace a component with a fake one, e.g. like the model was trained before that component was renamed
    auto modelComponents = configComponents;
    modelComponents.back() = "FakeComponent";
    expectThrowWithModel(modelComponents, configValues);
  }
#else
  GTEST_SKIP() << "Skipping test as AUTOPAS_ENABLE_PYTHON_BASED_TUNING=OFF";
#endif
}

}  // namespace autopas