/**
 * @file LoggerTest.h
 * @author seckler
 * @date 16.04.18
 */

#pragma once

#include <gtest/gtest.h>

#include <sstream>

#include "AutoPasTestBase.h"
#include "autopas/utils/logging/Logger.h"

class LoggerTest : public AutoPasTestBase {
 public:
  void SetUp() override;
  int testLevel(autopas::Logger::LogLevel level, bool enabled);

 private:
  std::stringstream _stream;
};
