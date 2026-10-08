/**
 * @file RuntimeScheduledLoopCountsTest.h
 * @author S. J. Newcome
 * @date 08/10/2026
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"

/**
 * Tests that the loop counts the traversals report (see autopas::TraversalInterface::getRuntimeScheduledLoopCounts())
 * are the numbers of iterations the OpenMP runtime actually schedules.
 *
 * The OpenMP runtime reports these through the OpenMP tools interface (OMPT). Only libomp supports this, so the test is
 * skipped with other runtimes.
 */
class RuntimeScheduledLoopCountsTest : public AutoPasTestBase {};