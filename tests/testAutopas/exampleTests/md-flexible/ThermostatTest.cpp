/**
 * @file ThermostatTest.cpp
 * @author N. Fottner
 * @date 28/08/19.
 */

#include "ThermostatTest.h"

//@todo write real tests that assert something :) WIP

void ThermostatTest::initFillWithParticles(std::array<size_t, 3> particlesPerDim, double particleSpacing, double cutoff) {
    //initializes
    double minimalBoxLength = cutoff + 0.2 /*default skin value*/;
    std::array<double, 3> boxmax = {std::max(particlesPerDim[0] * particleSpacing, minimalBoxLength),
                                    std::max(particlesPerDim[1] * particleSpacing, minimalBoxLength),
                                    std::max(particlesPerDim[2] * particleSpacing, minimalBoxLength)};
    std::array<double, 3> boxmin = {0., 0., 0.};
    _autopas.setBoxMin(boxmin);
    _autopas.setBoxMax(boxmax);
    _autopas.setCutoff(cutoff);
    _autopas.init();
    PrintableMolecule dummy;
    GridGenerator::fillWithParticles(_autopas, particlesPerDim, 0, 0, dummy,
                                     {particleSpacing, particleSpacing, particleSpacing}, {0., 0., 0.});
}


TEST_F(ThermostatTest, BasicTest_BM_ON){
    //initializes 9 particles with velocity 0.
    initFillWithParticles({2,2,2}/*particlesPerdim*/, 1. /*particleSpacing*/,1.5 /*cutoff*/);
    //    for(auto iter= _autopas.begin();iter.isValid();++iter){
    //        std::cout << "Particle Nr: " << iter->getID() << " with velo: " << autopas::ArrayUtils::to_string(iter->getV()) <<std::endl;
    //    }
    double init_temp = 0.01;
    double target_temp = 3.0;
    double delta_temp = 0.001;
    auto _thermostat =Thermostat<decltype(_autopas), std::remove_reference_t<decltype(_particlePropertiesLibrary)>>(init_temp,true,target_temp,delta_temp,_particlePropertiesLibrary);
    std::cout << "TEMPERATURE BEFORE INITIALIZATION:" << std::endl;
    std::cout << _thermostat.temperature(_autopas) << std::endl;
    std::cout << "INITIALIZATION" << std::endl;
    _thermostat.initialize(_autopas);
    std::cout << "current temperature: " << std::endl;  //is output right?
    std::cout << _thermostat.temperature(_autopas) << std::endl;
    std::cout << "APPLICATION" << std::endl;
    for(size_t i=0;i<2000;i++) {
        _thermostat.applyThermo(_autopas);
        if(i % 100==0) {
            std::cout << "Temperature at timeStep: " << i;
            std::cout << "   is: " << _thermostat.temperature(_autopas) << std::endl << std::endl;
        }
    }
    SUCCEED();
}
TEST_F(ThermostatTest, BasicTest_BM_OFF){
    //initializes 9 particles with velocity 0.
    initFillWithParticles({2,2,2}/*particlesPerdim*/, 1. /*particleSpacing*/,1.5 /*cutoff*/);
        for(auto iter= _autopas.begin();iter.isValid();++iter){
                iter->addV({0.5,0.5,0.5});
        }
    double init_temp = 0.01;
    double target_temp = 3.0;
    double delta_temp = 0.001;
    auto _thermostat =Thermostat<decltype(_autopas), std::remove_reference_t<decltype(_particlePropertiesLibrary)>>(init_temp,false,target_temp,delta_temp,_particlePropertiesLibrary);
    std::cout << "TEMPERATURE BEFORE INITIALIZATION:" << std::endl;
    std::cout << _thermostat.temperature(_autopas) << std::endl;
    std::cout << "INITIALIZATION" << std::endl;
    _thermostat.initialize(_autopas);
    std::cout << "current temperature: " << std::endl;  //is output right?
    std::cout << _thermostat.temperature(_autopas) << std::endl;
    std::cout << "APPLICATION" << std::endl;
    for(size_t i=0;i<2000;i++) {
        _thermostat.applyThermo(_autopas);
        if(i % 100==0) {
            std::cout << "Temperature at timeStep: " << i;
            std::cout << "   is: " << _thermostat.temperature(_autopas) << std::endl << std::endl;
        }
    }
    SUCCEED();
}
