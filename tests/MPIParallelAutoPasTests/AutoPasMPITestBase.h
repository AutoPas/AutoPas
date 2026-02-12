/**
 * @file AutoPasMPITestBase.h
 * @author W. Thieme
 * @date 01.05.2020
 */
#pragma once

#include <gtest/gtest.h>

#include "autopas/utils/logging/Logger.h"

#pragma once

class AutoPasMPITestBase : public testing::Test {
 public:
  explicit AutoPasMPITestBase(std::ostream &os = std::cout) : _logger(os) {}

 private:
  autopas::Logger _logger;
};
