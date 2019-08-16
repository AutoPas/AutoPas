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

void SimulationTest::VisualizeSmallSzenario(std::array<size_t,3> particlesPerDim, double cutoff, double particleSpacing, double epsilon, double sigma, double mass, int iterations, double delta_t, const std::string &filename){
    //initializes all members for a VisualizeSmallSzenario to test behaviour of Particles interacting with eachother
    auto autopas = autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>();
    initFillWithParticles(autopas,particlesPerDim,particleSpacing,cutoff);
    ParticlePropertiesLibrary<double, size_t> _particlePropertiesLibrary;
    _particlePropertiesLibrary.addType(0, epsilon, sigma, mass);
    auto functor = autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>>(cutoff, 0.0,_particlePropertiesLibrary);
    TimeDiscretization<decltype(autopas), decltype(_particlePropertiesLibrary)> _timeDiscretization(delta_t, _particlePropertiesLibrary);
    ASSERT_TRUE(autopas.getNumberOfParticles()==particlesPerDim[0]*particlesPerDim[1]*particlesPerDim[2]);
    //starting simulation
    for(int i=0;i< iterations;i++){
        writeVTKFile(autopas,i,filename);
        _timeDiscretization.CalculateX(autopas);
        autopas.iteratePairwise(&functor);
        _timeDiscretization.CalculateV(autopas);
    }
}

TEST_F(SimulationTest,BasicSzenarios){
    //following Tests are all with cutoff > ParticleSpacing
    //tests Particle movement when particleSpacing = 1.12 * sigma and cutoff = 1.5 > particleSpacing
    auto FirstAutopas = autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>();
    double FirstCutoff = 1.5;
    double FirstEpsilon = 1.0;
    double FirstSigma = 1.0;
    double FirstMass = 1.0;
    double FirstDeltaT = 0.002;
    ParticlePropertiesLibrary<double, size_t> _particlePropertiesLibrary;
    _particlePropertiesLibrary.addType(0, FirstEpsilon, FirstSigma, FirstMass);
    auto FirstFunctor = autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>,true>(FirstCutoff, 0.0,_particlePropertiesLibrary);
    TimeDiscretization<decltype(FirstAutopas), decltype(_particlePropertiesLibrary)> _timeDiscretization(FirstDeltaT, _particlePropertiesLibrary);
    //as particleSpacing=1.12 * sigma and cutoff > particleSpacing, testing if particlesPositions stay the same
    std::array<double,3> boxmax={5.,5.,5.};
    std::array<double,3> boxmin={0.,0.,0.};
    FirstAutopas.setBoxMin(boxmin);
    FirstAutopas.setBoxMax(boxmax);
    FirstAutopas.setCutoff(FirstCutoff);
    FirstAutopas.init();
    PrintableMolecule p1 ({0.,0.,0.},{0.,0.,0.},0);
    PrintableMolecule p2 ({1.12,0.,0.},{0.,0.,0.},1);
    FirstAutopas.addParticle(p1);
    FirstAutopas.addParticle(p2);
    //forces of particle1 are getting lower
    //forces of particle2 are getting higher
    std::array<double,3> particle1ForceBefore={0.,0.,0.};
    std::array<double,3> particle2ForceBefore={0.,0.,0.};
    std::array<double,3> particle1PositionBefore={0.,0.,0.};
    std::array<double,3> particle2PositionBefore={1.12,0.,0.};
    for(int i=0;i < 10000 /*number of iterations*/ ; i++) {
        FirstAutopas.iteratePairwise(&FirstFunctor);
        auto iter=FirstAutopas.begin();
        ASSERT_EQ(particle1PositionBefore,iter->getR());
        particle1PositionBefore = iter->getR();
        //forces of Particle1 are getting lower
        ASSERT_LE(iter->getF(),particle1ForceBefore);
        particle1ForceBefore=iter->getF();
        ++iter;
        ASSERT_EQ(particle2PositionBefore,iter->getR());
        particle2PositionBefore = iter->getR();
        //forces of Particle2 are getting higher
        ASSERT_GE(iter->getF(),particle2ForceBefore);
        particle2ForceBefore = iter->getF();
    }
    //as particleSpacing > 1.12 * sigma and cutoff > particleSpacing, testing if particlesPositions Attrakt each other
    //resets all forces
//    FirstAutopas.deleteAllParticles();
//    std::cout << "RESET" << std::endl;
//    PrintableMolecule p3 ({0.,0.,0.},{0.,0.,0.},0);
//    PrintableMolecule p4 ({1.13,0.,0.},{0.,0.,0.},1);
//    for(int i=0; i<10 /*number of iterations*/ ;i++){
//        FirstAutopas.iteratePairwise(&FirstFunctor);
//        auto iter=FirstAutopas.begin();
//        std::cout << iter->toString() << std::endl;
//        ++iter;
//        std::cout << iter->toString() << std::endl;
//    }
//
//
//
    }

//    for(int i=0;i < 10 /* 10 iterations */;i ++){
//        _timeDiscretization.CalculateX(FirstAutopas);
//
//    }

//build function(particlesPerDim,cutoff,particleSpacing,epsilon,sigma,mass,iterations,delta_t)
//test function(autopas,timeDiskretization)

//system
//reset function

TEST_F(SimulationTest,smallSzenariosVTKVisualization){
    VisualizeSmallSzenario({2, 2, 1}, /*particlesPerDim*/
                           1.0, /*cutoff*/
                           1.12, /*particleSpacing*/
                           1.12, /*epsilon*/
                           1.0, /*sigma*/
                           1.0, /*mass*/
                           20, /*number of Iterations*/
                           0.002, /*delta_t*/
                           "Particles:2*2_Cutoff:1.12_sigma:1.12_PS:1.12" /*filename*/ );

    VisualizeSmallSzenario({2, 2, 1}, /*particlesPerDim*/
                           1.0, /*cutoff*/
                           1.12, /*particleSpacing*/
                           1.0, /*epsilon*/
                           1.0, /*sigma*/
                           1.0, /*mass*/
                           20, /*number of Iterations*/
                           0.002, /*delta_t*/
                           "Particles:2*2_Cutoff:1_sigma:1.12_PS:1.12" /*filename*/ );

    VisualizeSmallSzenario({8, 8, 1}, /*particlesPerDim*/
                           1.0, /*cutoff*/
                           1.13, /*particleSpacing*/
                           1.0, /*epsilon*/
                           1.0, /*sigma*/
                           1.0, /*mass*/
                           20, /*number of Iterations*/
                           0.002, /*delta_t*/
                           "Particles:8*8_Cutoff:1_sigma:1_PS:1.13" /*filename*/ );

    VisualizeSmallSzenario({8, 8, 1}, /*particlesPerDim*/
                           1.0, /*cutoff*/
                           1.13, /*particleSpacing*/
                           1.0, /*epsilon*/
                           1.0, /*sigma*/
                           1.0, /*mass*/
                           20, /*number of Iterations*/
                           0.002, /*delta_t*/
                           "Particles:8*8_Cutoff:1_sigma:1_PS:1.13" /*filename*/ );

    VisualizeSmallSzenario({8, 8, 1}, /*particlesPerDim*/
                           1.0, /*cutoff*/
                           1.13, /*particleSpacing*/
                           1.0, /*epsilon*/
                           1.0, /*sigma*/
                           1.0, /*mass*/
                           20, /*number of Iterations*/
                           0.002, /*delta_t*/
                           "Particles:8*8_Cutoff:1_sigma:1_PS:1.13" /*filename*/ );

    ASSERT_TRUE(true);
}