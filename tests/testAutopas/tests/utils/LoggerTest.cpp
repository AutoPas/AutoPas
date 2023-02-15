/**
 * @file LoggerTest.cpp
 * @author seckler
 * @date 16.04.18
 */

#include "LoggerTest.h"

void LoggerTest::SetUp() { autopas::Logger::create(stream); }

void LoggerTest::TearDown() { autopas::Logger::unregister(); }

int LoggerTest::testLevel(autopas::Logger::LogLevel level, bool enabled = true) {
  autopas::Logger::get()->set_level(level);
  if (not enabled) autopas::Logger::get()->set_level(autopas::Logger::LogLevel::off);

  stream.flush();
  stream.clear();

  AutoPasLog(TRACE, "trace");
  AutoPasLog(DEBUG, "debug");
  AutoPasLog(INFO, "info");
  AutoPasLog(WARN, "warn");
  AutoPasLog(ERROR, "error");
  AutoPasLog(CRITICAL, "critical");

  int lineCount = 0;
  std::string str;
  while (getline(stream, str)) ++lineCount;

  return lineCount;
}

TEST_F(LoggerTest, LogLevelTest) {
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::trace), 6);
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::debug), 5);
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::info), 4);
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::warn), 3);
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::err), 2);
  EXPECT_EQ(testLevel(autopas::Logger::LogLevel::critical), 1);
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