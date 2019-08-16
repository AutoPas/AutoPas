/**
 * @file PeriodicBoundariesTest.h
 * @author N. Fottner
 * @date 2/8/19
 */

#include "SimulationTest.h"

void SimulationTest::initFillWithParticles(
        autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autopas,
        std::array<size_t, 3> particlesPerDim,double particleSpacing,double cutoff) {
    double minimalBoxLength = cutoff + 0.2/*default skin value*/;
    std::array<double,3> boxmax = {std::max(particlesPerDim[0]*particleSpacing,minimalBoxLength),std::max(particlesPerDim[1]*particleSpacing,minimalBoxLength),std::max(particlesPerDim[2]*particleSpacing,minimalBoxLength)};
    std::array<double,3> boxmin={0.,0.,0.};
    autopas.setBoxMin(boxmin);
    autopas.setBoxMax(boxmax);
    autopas.init();
    PrintableMolecule dummy;
    GridGenerator::fillWithParticles(autopas, particlesPerDim, 0, 0, dummy, {particleSpacing, particleSpacing, particleSpacing}, {0., 0., 0.});
}

void SimulationTest::smallSzenario(std::array<size_t,3> particlesPerDim,double cutoff,double particleSpacing,double epsilon,double sigma,double mass,int iterations,double delta_t){
    //initializes all members for a smallSzenario to test behaviour of Particles interacting with eachother
    auto autopas = autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>();
    initFillWithParticles(autopas,particlesPerDim,particleSpacing,cutoff);
    ParticlePropertiesLibrary<double, size_t> _particlePropertiesLibrary;
    _particlePropertiesLibrary.addType(0, epsilon, sigma, mass);
    auto functor = autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>>(cutoff, 0.0,_particlePropertiesLibrary);
    TimeDiscretization<decltype(autopas), decltype(_particlePropertiesLibrary)> _timeDiscretization(delta_t, _particlePropertiesLibrary);
    ASSERT_TRUE(autopas.getNumberOfParticles()==particlesPerDim[0]*particlesPerDim[1]*particlesPerDim[2]);
    //starting simulation
    for(int i=0;i< iterations;i++){
        writeVTKFile(autopas,i);
        _timeDiscretization.CalculateX(autopas);
        autopas.iteratePairwise(&functor);
        _timeDiscretization.CalculateV(autopas);
    }
}

TEST_F(SimulationTest,smallSzenarios){
    smallSzenario({8,8,1}, /*particlesPerDim*/
                  1.0, /*cutoff*/
                  1.13, /*particleSpacing*/
                  1.0, /*epsilon*/
                  1.0, /*sigma*/
                  1.0, /*mass*/
                  20, /*number of Iterations*/
                  0.002 /*delta_t*/);

    smallSzenario({8,8,1}, /*particlesPerDim*/
                  1.0, /*cutoff*/
                  1.13, /*particleSpacing*/
                  1.0, /*epsilon*/
                  1.0, /*sigma*/
                  1.0, /*mass*/
                  20, /*number of Iterations*/
                  0.002 /*delta_t*/);



    ASSERT_TRUE(true);
}
