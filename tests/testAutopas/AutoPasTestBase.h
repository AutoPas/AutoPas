/**
 * @file AutoPasTestBase.h
 * @author seckler
 * @date 24.04.18
 */

#pragma once

#include <gtest/gtest.h>

#include "autopas/utils/logging/Logger.h"


class AutoPasTestBase : public testing::Test {
 public:
  AutoPasTestBase() { autopas::Logger::create(); }

  virtual ~AutoPasTestBase() { autopas::Logger::unregister(); }
};
