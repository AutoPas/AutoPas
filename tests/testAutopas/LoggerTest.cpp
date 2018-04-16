/**
 * @file LoggerTest.cpp
 * @author seckler
 * @date 16.04.18
 */

#include "LoggerTest.h"
#include "autopasIncludes.h"
#include "gtest/gtest.h"

using namespace autopas::log;

int testLevel(logLevel level) {
  std::stringstream stream;
  Logger log(level, &stream);

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

TEST(LoggerTest, defaultConstructorTest) {
  std::stringstream stream;
  {
    ScopedRedirect redirect(std::cout, stream);
    Logger log;
    log.fatal() << "test" << std::endl;
    log.error() << "test" << std::endl;
    log.warning() << "test" << std::endl;
  }
  int lineCount = 0;
  std::string str;
  while (getline(stream, str)) ++lineCount;
  EXPECT_EQ(lineCount, 2);
}
