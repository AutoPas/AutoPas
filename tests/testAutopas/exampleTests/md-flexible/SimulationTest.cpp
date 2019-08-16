/**
 * @file PeriodicBoundariesTest.h
 * @author N. Fottner
 * @date 2/8/19
 */

#include "SimulationTest.h"

void SimulationTest::initFillWithParticles(
        autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autopas,
        std::array<unsigned long, 3> particlesPerDim) {
    autopas.setBoxMin(boxmin);
    autopas.setBoxMax(boxmax);
    autopas.init();
    PrintableMolecule dummy;
    GridGenerator::fillWithParticles(autopas, particlesPerDim, 0, 0, dummy, {1, 1, 1}, {0., 0., 0.});
}

TEST_F(){

}


//TEST_F(SimulationTest,PeriodicVisualization){
//    _simulation.simulate();
//    ASSERT_TRUE(true);
//}