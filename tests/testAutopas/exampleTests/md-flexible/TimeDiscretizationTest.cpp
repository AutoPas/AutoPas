//
// Created by nicola on 22.05.19.
//


#include <gtest/gtest.h>
#include "autopas/AutoPas.h"
#include "../../testingHelpers/GridGenerator.h"
#include "../../testingHelpers/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"
#include "../../../../examples/md-flexible/PrintableMolecule.h"
#include "../../../../examples/md-flexible/TimeDiscretization.h"
#include "../../../../examples/md-flexible/MDFlexParser.h"


using namespace std;

string arrayString(array<double,3> input){
    std::ostringstream os;
    for (double i: input){
        os << i;
    }
    std::string str(os.str());
    return str;
}

TEST(TimeDTEst,CalcX){
    //auto container = autopas::LinkedCells<PrintableMolecule,ParticleCell<PrintableMolecule>>({0., 0., 0.},{10., 10., 10.},1);
    autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>* autoPas;
    PrintableMolecule::setEpsilon(1.0);
    PrintableMolecule::setSigma(1.0);
    PrintableMolecule::setMass(1.0);
    autoPas->setBoxMin({0., 0., 0.});
    autoPas->setBoxMax({10., 10., 10.});
    autoPas->setCutoff(1.);
    // erstmal auf linked cells testen
    autoPas->setAllowedContainers({autopas::ContainerOption::linkedCells});
    autoPas->init();
    PrintableMolecule p1({1., 1., 1.}, {0.5, 0.5, 0.5}, 0);
    autoPas->addParticle(p1);
    PrintableMolecule p2({2., 2., 2.}, {0., 0.5, 0.}, 1);
    autoPas->addParticle(p2);
    for (auto iter = autoPas->getContainer()->begin() ; iter.isValid(); ++iter) {
        cout << "Particles information" << iter->toString() << "  __END" << endl;
    }

    double particleD = 0.01;
    int iterations;
    //iterationen beginnend
    TimeDiscretization<decltype(autoPas)> td(particleD);
    auto* functor = new autopas::LJFunctor<PrintableMolecule,ParticleCell<PrintableMolecule>, autopas::FunctorN3Modes::Both, true>(1, 1.0, 1.0, 0.0,{0., 0., 0.},{5., 5., 5.},true);

    //domain vorbeireiten: -Force initialisieren
    autoPas->iteratePairwise(functor);
    //Dokumentation prints
    for (auto iter = autoPas->getContainer()->begin() ; iter.isValid(); ++iter) {
        cout << "Particles information" << iter->toString() << "  __END" << endl;
    }

    //while(iterations < 10){
      //  td.VSCalculateX(autoPas);
        //autoPas->iteratePairwise(functor);
      //  td.VSCalculateV(autopas);
       // iterations++;
   // }




    ASSERT_TRUE(true);
}



//TEST(TimeDTEst,CalcV){}
