/**
 * @file ThermostatTest.cpp
 * @author N. Fottner
 * @date 28/08/19.
 */

#include "ThermostatTest.h"

void ThermostatTest::initFillWithParticles(std::array<size_t, 3> particlesPerDim, double particleSpacing,
                                           double cutoff, autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autopas) {
  double minimalBoxLength = cutoff + 0.2 /*default skin value*/;
  std::array<double, 3> boxmax = {std::max(particlesPerDim[0] * particleSpacing, minimalBoxLength),
                                  std::max(particlesPerDim[1] * particleSpacing, minimalBoxLength),
                                  std::max(particlesPerDim[2] * particleSpacing, minimalBoxLength)};
  std::array<double, 3> boxmin = {0., 0., 0.};
  autopas.setBoxMin(boxmin);
  autopas.setBoxMax(boxmax);
  autopas.setCutoff(cutoff);
  autopas.init();
  PrintableMolecule dummy;
  GridGenerator::fillWithParticles(autopas, particlesPerDim, 0, 0, dummy,
                                   {particleSpacing, particleSpacing, particleSpacing}, {0., 0., 0.});
}

void ThermostatTest::basicApplication(double initT, double targetT, double deltaT,bool initBM,autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autopas){
    double nrApplications = targetT / deltaT;
    auto _thermostat = Thermostat<decltype(autopas), std::remove_reference_t<decltype(_particlePropertiesLibrary)>>(
            initT, initBM, targetT, deltaT, _particlePropertiesLibrary);
    if(initBM){
        _thermostat.initialize(autopas);
    }else {
        //initial velocity value of particles necessary otherwise zero divisions causing error
        for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
            iter->addV({0, 0.1, 0});
        }
        _thermostat.initialize(autopas);
    }
    _thermostat.apply(autopas);
    for (size_t i = 1; i <= nrApplications; i++) {
        EXPECT_NEAR(_thermostat.calcTemperature(autopas), (initT + (deltaT * i)) > targetT ? targetT : initT + (deltaT * i) , absDelta);
        _thermostat.apply(autopas);
    }
    _thermostat.apply(autopas);
    EXPECT_NEAR(_thermostat.calcTemperature(autopas), targetT, absDelta);
}

void ThermostatTest::calcTemperature(size_t particlesPerDimension) {
    auto autopas = autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>();
    initFillWithParticles({particlesPerDimension, particlesPerDimension, particlesPerDimension} /*particlesPerdim*/, 1. /*particleSpacing*/, 1.5 /*cutoff*/,
                          autopas /*autopas Object*/);
    auto _thermostat = Thermostat<decltype(autopas), std::remove_reference_t<decltype(_particlePropertiesLibrary)>>(
            0.1/*initT*/, true /*initBM*/, 5 /*targetT*/, 0.01 /*deltaT*/, _particlePropertiesLibrary);
    for(auto iter=autopas.begin();iter.isValid();++iter){
        iter->setV({0.1,0.1,0.2});
    }
    EXPECT_NEAR(_thermostat.calcTemperature(autopas), 0.02,absDelta);
}

TEST_F(ThermostatTest, Application){
    auto _autopas= autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>();
    initFillWithParticles({2, 2, 2} /*particlesPerdim*/, 1. /*particleSpacing*/, 1.5 /*cutoff*/,_autopas /*autopas Object*/);
    basicApplication(0.01,3.0,0.001,true,_autopas);
    //reset _autopas
    initFillWithParticles({2, 2, 2} /*particlesPerdim*/, 1. /*particleSpacing*/, 1.5 /*cutoff*/,_autopas /*autopas Object*/);
    basicApplication(1.0,7.0,0.01,true,_autopas);
    //dont use current Temp for Brownian Motion initialization
    initFillWithParticles({2, 2, 2} /*particlesPerdim*/, 1. /*particleSpacing*/, 1.5 /*cutoff*/,_autopas /*autopas Object*/);
    basicApplication(0.01,3.0,0.001,false,_autopas);
    //reset _autopas
    initFillWithParticles({2, 2, 2} /*particlesPerdim*/, 1. /*particleSpacing*/, 1.5 /*cutoff*/,_autopas /*autopas Object*/);
    basicApplication(1.0,7.0,0.01,false,_autopas);
}

TEST_F(ThermostatTest, calcTemperature) {
    calcTemperature(1);
    calcTemperature(2);
    calcTemperature(10);
}