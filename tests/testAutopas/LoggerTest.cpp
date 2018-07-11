/**
 * @file LoggerTest.cpp
 * @author seckler
 * @date 16.04.18
 */

#include "LoggerTest.h"

void LoggerTest::SetUp() { autopas::Logger::create(stream); }

void LoggerTest::TearDown() { autopas::Logger::unregister(); }

int LoggerTest::testLevel(spdlog::level::level_enum level, bool enabled = true) {
  AutoPasLogger->set_level(level);
  if (not enabled) AutoPasLogger->set_level(spdlog::level::off);

  stream.flush();
  stream.clear();

  AutoPasLogger->trace("trace");
  AutoPasLogger->debug("debug");
  AutoPasLogger->info("info");
  AutoPasLogger->warn("warn");
  AutoPasLogger->error("error");
  AutoPasLogger->critical("critical");

  int lineCount = 0;
  std::string str;
  while (getline(stream, str)) ++lineCount;

  return lineCount;
}

TEST_F(LoggerTest, LogLevelTest) {
  EXPECT_EQ(testLevel(spdlog::level::trace), 6);
  EXPECT_EQ(testLevel(spdlog::level::debug), 5);
  EXPECT_EQ(testLevel(spdlog::level::info), 4);
  EXPECT_EQ(testLevel(spdlog::level::warn), 3);
  EXPECT_EQ(testLevel(spdlog::level::err), 2);
  EXPECT_EQ(testLevel(spdlog::level::critical), 1);
  EXPECT_EQ(testLevel(spdlog::level::off), 0);
}

TEST_F(LoggerTest, LogLevelTestDisabled) {
  EXPECT_EQ(testLevel(spdlog::level::trace, false), 0);
  EXPECT_EQ(testLevel(spdlog::level::debug, false), 0);
  EXPECT_EQ(testLevel(spdlog::level::info, false), 0);
  EXPECT_EQ(testLevel(spdlog::level::warn, false), 0);
  EXPECT_EQ(testLevel(spdlog::level::err, false), 0);
  EXPECT_EQ(testLevel(spdlog::level::critical, false), 0);
  EXPECT_EQ(testLevel(spdlog::level::off, false), 0);
}