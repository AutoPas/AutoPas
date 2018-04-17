/**
 * @file LoggerTest.cpp
 * @author seckler
 * @date 16.04.18
 */

#include "LoggerTest.h"
#include "autopasIncludes.h"
#include "gtest/gtest.h"

using namespace autopas::log;

int testLevel(logLevel level, bool enabled = true) {
  std::stringstream stream;
  Logger log(level, &stream, &stream);
  log.setEnabled(enabled);

  log.debug() << "debug" << std::endl;
  log.info() << "info" << std::endl;
  log.warning() << "warning" << std::endl;
  log.error() << "error" << std::endl;
  log.fatal() << "fatal" << std::endl;

  int lineCount = 0;
  std::string str;
  while (getline(stream, str)) ++lineCount;

  return lineCount;
}

TEST(LoggerTest, LogLevelTest) {
  EXPECT_EQ(testLevel(logLevel::All), 5);
  EXPECT_EQ(testLevel(logLevel::Debug), 5);
  EXPECT_EQ(testLevel(logLevel::Info), 4);
  EXPECT_EQ(testLevel(logLevel::Warning), 3);
  EXPECT_EQ(testLevel(logLevel::Error), 2);
  EXPECT_EQ(testLevel(logLevel::Fatal), 1);
  EXPECT_EQ(testLevel(logLevel::None), 0);
}

TEST(LoggerTest, LogLevelTestDisabled) {
  EXPECT_EQ(testLevel(logLevel::All, false), 0);
  EXPECT_EQ(testLevel(logLevel::Debug, false), 0);
  EXPECT_EQ(testLevel(logLevel::Info, false), 0);
  EXPECT_EQ(testLevel(logLevel::Warning, false), 0);
  EXPECT_EQ(testLevel(logLevel::Error, false), 0);
  EXPECT_EQ(testLevel(logLevel::Fatal, false), 0);
  EXPECT_EQ(testLevel(logLevel::None, false), 0);
}

TEST(LoggerTest, defaultConstructorTest) {
  std::stringstream stream_out;
  std::stringstream stream_err;
  {
    ScopedRedirect redirect(std::cout, stream_out);
    ScopedRedirect redirect2(std::cerr, stream_err);
    Logger log;
    log.fatal() << "test" << std::endl;
    log.error() << "test" << std::endl;
    log.warning() << "test" << std::endl;
  }
  int lineCount = 0;
  std::string str;
  while (getline(stream_out, str)) ++lineCount;
  EXPECT_EQ(lineCount, 0);
  while (getline(stream_err, str)) ++lineCount;
  EXPECT_EQ(lineCount, 2);
}
