/**
 * @file LoggerTest.cpp
 * @author seckler
 * @date 16.04.18
 */

#include "LoggerTest.h"

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