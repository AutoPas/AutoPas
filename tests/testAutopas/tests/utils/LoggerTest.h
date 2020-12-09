/**
 * @file LoggerTest.h
 * @author seckler
 * @date 16.04.18
 */

#pragma once

#include <gtest/gtest.h>
#include <spdlog/common.h>

#include <ostream>

#include "AutoPasTestBase.h"

namespace LoggerTest {

class LoggerTest : public AutoPasTestBase {
 public:
  void SetUp() override;
  void TearDown() override;
  int testLevel(autopas::Logger::LogLevel level, bool enabled);

 private:
  std::stringstream stream;
};

// Helper:

class ScopedRedirect {
 public:
  ScopedRedirect(std::ostream &inOriginal, std::ostream &inRedirect) : mOriginal(inOriginal), mRedirect(inRedirect) {
    mOriginal.rdbuf(mRedirect.rdbuf(mOriginal.rdbuf()));
  }

  ~ScopedRedirect() { mOriginal.rdbuf(mRedirect.rdbuf(mOriginal.rdbuf())); }

 private:
  ScopedRedirect(const ScopedRedirect &);
  ScopedRedirect &operator=(const ScopedRedirect &);

  std::ostream &mOriginal;
  std::ostream &mRedirect;
};
} // end namespace LoggerTest
