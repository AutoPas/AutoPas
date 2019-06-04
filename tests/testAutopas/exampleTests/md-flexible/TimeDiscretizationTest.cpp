//
// Created by nicola on 22.05.19.
//


#include <gtest/gtest.h>
#include "autopas/AutoPas.h"
#include "../../testingHelpers/GridGenerator.h"
#include "../../testingHelpers/RandomGenerator.h"
#include "../../../../examples/md-flexible/PrintableMolecule.h"
#include "testingHelpers/commonTypedefs.h"


using namespace std;


TEST(TimeDTEst,CalcX){
    autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> autoPas;
    PrintableMolecule::setEpsilon(1.0);
    PrintableMolecule::setSigma(1.0);
    PrintableMolecule::setMass(1.0);
    autoPas.setBoxMin({0., 0., 0.});
    autoPas.setBoxMax({10., 10., 10.});
    autoPas.setCutoff(1.);
    autoPas.init();
    PrintableMolecule p1({1., 1., 1.}, {0., 0., 0.}, 0);
    autoPas.addParticle(p1);
    PrintableMolecule p2({2., 2., 2.}, {0., 0., 0.}, 1);
    autoPas.addParticle(p2);
    for (auto iter = autoPas.getContainer()->begin() ; iter.isValid(); ++iter) {
        cout << "Particles information" << iter->toString() << "  __END" << endl;
    }
    ASSERT_TRUE(true);
}



//TEST(TimeDTEst,CalcV){}
