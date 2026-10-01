/**
 * @file AutoPasMPITestBase.h
 * @author W. Thieme
 * @date 01.05.2020
 */
#pragma once

#include <gtest/gtest.h>

#include "autopas/utils/WrapOpenMP.h"
#include "testingHelpers/NumThreadGuard.h"

class AutoPasMPITestBase : public testing::Test {
 protected:
  NumThreadGuard _numThreadGuard{autopas::autopas_get_max_threads(), TUNED_THREADS};
};
