/**
 * @file MPIParallelizedStrategyTest.h
 * @author W. Thieme
 * @date 11.06.2020
 */

#pragma once

#include <gtest/gtest.h>
#include <mpi.h>

#include "AutoPasMPITestBase.h"
#include "autopas/tuning/tuningStrategy/ActiveHarmony.h"
#include "autopas/tuning/tuningStrategy/BayesianSearch.h"
#include "autopas/tuning/tuningStrategy/FullSearch.h"
#include "autopas/tuning/tuningStrategy/MPIParallelizedStrategy.h"
#include "autopas/tuning/tuningStrategy/RandomSearch.h"

class MPIParallelizedStrategyTest : public AutoPasMPITestBase {};