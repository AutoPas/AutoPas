#pragma once
#include <gtest/gtest.h>
#include "../../../../examples/md-flexible/PrintableMolecule.h"
#include "../../../../examples/md-flexible/ParticleClassLibrary.h"
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"
using namespace std;

class PCLTest : public AutoPasTestBase{
public:
    PCLTest() : AutoPasTestBase(){
        PrintableMolecule::setMass(1.0);
        PrintableMolecule::setEpsilon(1.0);
        PrintableMolecule::setSigma(1.0);
        epsilon= 1.0;
        sigma =1.0;
        mass = 1.0;
        map<unsigned long, double> universalMap;
        for (unsigned long i = 0; i < 1000; i++) {
            universalMap.emplace(i, 1.0);
        }
        PCL = ParticleClassLibrary(universalMap,universalMap,universalMap);
    }
protected:
    PrintableMolecule dummyParticle;
    double epsilon;
    double sigma;
    double mass;
    ParticleClassLibrary PCL;
};