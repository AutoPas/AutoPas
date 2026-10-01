/**
 * @file FullParticleCellTest.h
 * @author seckler
 * @date 10.09.19
 */

#pragma once

#include <gtest/gtest.h>

#include "autopas/utils/WrapOpenMP.h"
#include "NumThreadGuard.h"

class FullParticleCellTest : public testing::Test {
 protected:
  NumThreadGuard _numThreadGuard{autopas::autopas_get_max_threads(), TUNED_THREADS};
};