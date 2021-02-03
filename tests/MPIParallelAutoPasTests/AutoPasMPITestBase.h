/**
 * @file AutoPasMPITestBase.h
 * @author W. Thieme
 * @date 01.05.2020
 */

#include <gtest/gtest.h>

#include "autopas/utils/logging/Logger.h"

#pragma once

class AutoPasMPITestBase : public testing::Test {
 public:
  AutoPasMPITestBase() { autopas::Logger::create(); }

  virtual ~AutoPasMPITestBase() { autopas::Logger::unregister(); }
};
