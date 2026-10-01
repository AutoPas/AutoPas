/**
 * @file AutoPasTestBase.h
 * @author seckler
 * @date 24.04.18
 */
#pragma once

#include <gtest/gtest.h>

#include "autopas/utils/WrapOpenMP.h"
#include "NumThreadGuard.h"

class AutoPasTestBase : public testing::Test {
 protected:
  NumThreadGuard _numThreadGuard{autopas::autopas_get_max_threads(), TUNED_THREADS};
};
