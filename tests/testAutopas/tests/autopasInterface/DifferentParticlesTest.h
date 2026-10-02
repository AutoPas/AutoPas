/**
 * @file DifferentParticlesTest.h
 * @author seckler
 * @date 20.02.2020
 */

#pragma once

#include <gtest/gtest.h>

#include "NumThreadGuard.h"
#include "autopas/utils/WrapOpenMP.h"

class DifferentParticlesTest : public testing::Test {
 protected:
  NumThreadGuard _numThreadGuard{autopas::autopas_get_max_threads(), TUNED_THREADS};
};
