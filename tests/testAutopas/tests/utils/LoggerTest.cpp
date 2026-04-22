/**
 * @file LoggerTest.cpp
 * @author seckler
 * @date 16.04.18
 */

#include "LoggerTest.h"

#include <filesystem>
#include <fstream>

// Create a logger that logs to _stream
void LoggerTest::SetUp() { autopas::Logger::create(_stream); }

/**
 * Counts how many log levels print when the logger is set to the given level.
 * @param level
 * @param enabled
 * @return Number of levels >= the given level.
 */
int LoggerTest::testLevel(autopas::Logger::LogLevel level, bool enabled = true) {
  autopas::Logger::get()->set_level(level);
  if (not enabled) autopas::Logger::get()->set_level(autopas::Logger::LogLevel::off);

  _stream.str("");
  _stream.clear();

  AutoPasLog(TRACE, "trace");
  AutoPasLog(DEBUG, "debug");
  AutoPasLog(INFO, "info");
  AutoPasLog(WARN, "warn");
  AutoPasLog(ERROR, "error");
  AutoPasLog(CRITICAL, "critical");

  int lineCount = 0;
  std::string str;
  while (getline(_stream, str)) ++lineCount;

  return lineCount;
}

TEST_F(LoggerTest, LogLevelTest) {
  // Some log levels might be deactivated on compile time. Calculate how many we expect.
  constexpr int numLogLevelsTotal = 6;
  constexpr int numLogLevelsActive = numLogLevelsTotal - SPDLOG_ACTIVE_LEVEL;

  // The expected number of lines is limited by the active log levels
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::trace), std::clamp(6, 0, numLogLevelsActive));
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::debug), std::clamp(5, 0, numLogLevelsActive));
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::info), std::clamp(4, 0, numLogLevelsActive));
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::warn), std::clamp(3, 0, numLogLevelsActive));
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::err), std::clamp(2, 0, numLogLevelsActive));
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::critical), std::clamp(1, 0, numLogLevelsActive));
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::off), 0);
}

TEST_F(LoggerTest, LogLevelTestDisabled) {
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::trace, false), 0);
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::debug, false), 0);
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::info, false), 0);
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::warn, false), 0);
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::err, false), 0);
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::critical, false), 0);
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::off, false), 0);
}

/**
 * Tests that Logger::get() safely performs lazy initialization when called without a prior create().
 */
TEST_F(LoggerTest, GetWithoutCreate) {
  // SetUp() calls create(), so we unregister first to test lazy initialization
  autopas::Logger::unregister();
  // Verify that the logger is actually gone from the spdlog registry
  EXPECT_EQ(spdlog::get("AutoPasLog"), nullptr);

  EXPECT_NO_THROW({
    auto logger = autopas::Logger::get();
    EXPECT_NE(logger, nullptr);
    EXPECT_EQ(logger->name(), "AutoPasLog");
  });
}

/**
 * Tests that Logger::get() correctly returns the active logger when called after Logger::create().
 */
TEST_F(LoggerTest, GetAfterCreate) {
  // SetUp() already called create()
  EXPECT_NO_THROW({
    auto logger = autopas::Logger::get();
    EXPECT_NE(logger, nullptr);
  });
}

/**
 * Tests that calling Logger::create() multiple times safely drops and recreates the logger without errors.
 */
TEST_F(LoggerTest, CreateMultipleTimes) {
  EXPECT_NO_THROW({
    autopas::Logger::create();
    autopas::Logger::create();
  });
}

/**
 * Tests that the logger correctly routes its output messages to a provided std::ostream.
 */
TEST_F(LoggerTest, CreateWithStream) {
  EXPECT_NO_THROW({
    std::ostringstream oss;
    autopas::Logger::create(oss);

    auto logger = autopas::Logger::get();
    EXPECT_NE(logger, nullptr);

    // Verify that the stream actually receives the log messages
    logger->warn("Test warning message");
    logger->flush();
    EXPECT_TRUE(oss.str().find("Test warning message") != std::string::npos);
  });
}

/**
 * Tests that Logger::unregister() can be called safely regardless of whether a logger currently exists in the registry.
 */
TEST_F(LoggerTest, UnregisterMultipleTimes) {
  EXPECT_NO_THROW({
    autopas::Logger::unregister();  // Unregister existing (from SetUp)
    autopas::Logger::unregister();  // Unregister when it doesn't exist
    autopas::Logger::create();      // Create it
    autopas::Logger::unregister();  // Unregister when it does exist
  });
}

/**
 * Tests that the logger correctly routes its output messages to a file.
 */
TEST_F(LoggerTest, CreateWithFile) {
  const std::string filename = "test_autopas_logger.log";

  EXPECT_NO_THROW({
    autopas::Logger::create(filename);

    auto logger = autopas::Logger::get();
    EXPECT_NE(logger, nullptr);

    // Write a message and flush to ensure it is written to disk immediately
    logger->warn("Test file warning message");
    logger->flush();
  });

  // Verify file contents
  std::ifstream ifs(filename);
  ASSERT_TRUE(ifs.is_open());
  std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
  EXPECT_TRUE(content.find("Test file warning message") != std::string::npos);
  ifs.close();

  // Clean up: unregister to release the file lock, then delete the test file
  autopas::Logger::unregister();
  std::filesystem::remove(filename);
}