/**
 * @file LoggerTest.h
 * @author seckler
 * @date 16.04.18
 */

#pragma once

#include <gtest/gtest.h>

#include <ostream>

#include "AutoPasTestBase.h"

class LoggerTest : public AutoPasTestBase {
 public:
  LoggerTest() : AutoPasTestBase(stream) {}
  int testLevel(autopas::Logger::LogLevel level, bool enabled);

 private:
  std::stringstream stream;
};
