/**
 * @file PeriodicBoundariesTest.h
 * @author N. Fottner
 * @date 2/8/19
 */

#include "PeriodicBoundariesTest.h"


TEST_F(PeriodicBoundariesTest,PeriodicVisualization){
    _simulation.simulate();
    ASSERT_TRUE(true);
}