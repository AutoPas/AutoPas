/**
 * @file MPIParallelizedStrategyTest.h
 * @author W. Thieme
 * @date 11.06.2020
 */

#pragma once

#include <gtest/gtest.h>
#include <mpi.h>

#include "AutoPasMPITestBase.h"
#include "autopas/selectors/tuningStrategy/ActiveHarmony.h"
#include "autopas/selectors/tuningStrategy/BayesianSearch.h"
#include "autopas/selectors/tuningStrategy/FullSearch.h"
#include "autopas/selectors/tuningStrategy/MPIParallelizedStrategy.h"
#include "autopas/selectors/tuningStrategy/RandomSearch.h"

class MPIParallelizedStrategyTest : public AutoPasMPITestBase {};