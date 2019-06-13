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
#include <math.h>
#include "../../../../src/autopas/utils/ArrayMath.h"

using namespace std;
using namespace autopas;

string arrayString(array<double,3> input){
    std::ostringstream os;
    for (double i: input){
        os << i;
        os << " _ ";
    }
    std::string str(os.str());
    return str;
}

double L2Norm(std::array<double,3> array) {
    double square_sum = 0;
    for (unsigned int i = 0; i < array.size(); i++) {
        square_sum += (array[i] * array[i]);
    }
    return sqrt(square_sum);
}


std::array<double,3> lennardForceCalculation(std::array<double, 3> x1, std::array<double, 3> x2){
    std::array<double, 3> difference = ArrayMath::sub(x1, x2);
    double distance=L2Norm(difference);
    double epsilon=5; double sigma=1;
    return ArrayMath::mulScalar(difference,(24*epsilon)/(distance*distance) * (std::pow(sigma/distance,6) - 2 * std::pow(sigma/distance,12)));
}


TEST(TimeDTest,LennardJonesCalculation){
    //dient zum vergleich zu den Forces die im Functor ausgerechnet werden
    std::array<double,3> x1 = {1,1,1};
    std::array<double,3> x2 = {1.5,1.5,1.5};
    cout << "lennardJonesFunction Calculation:  " << arrayString(lennardForceCalculation(x1,x2)) << endl;
    ASSERT_TRUE(true);
}

//Testet und visualisiert die Kräfte berechnungen und TimeDiscreatization Classe
TEST(TimeDTEst,GeneralForceTest){
    //auto container = autopas::LinkedCells<PrintableMolecule,ParticleCell<PrintableMolecule>>({0., 0., 0.},{10., 10., 10.},1);
    PrintableMolecule::setEpsilon(5.0);
    PrintableMolecule::setSigma(1.0);
    PrintableMolecule::setMass(1.0);
    auto *autoPas = new autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>(std::cout);
    double epsilon = 5.0;double sigma=1.0;double cutoff=2;array<double,3> boxmin={0.,0.,0.}; array<double,3> boxmax={5.,5.,5.};
    PrintableMolecule::setEpsilon(epsilon);
    PrintableMolecule::setSigma(sigma);
    PrintableMolecule::setMass(1.0);
    autoPas->setBoxMax(boxmax);
    autoPas->setCutoff(cutoff);
    // erstmal auf linked cells testen
    autoPas->setAllowedContainers({autopas::ContainerOption::linkedCells});
    autoPas->init();
    PrintableMolecule p1({1., 1., 1.}, {0.5, 0.5, 0.5}, 0);
    autoPas->addParticle(p1);
    PrintableMolecule p2({1.5, 1.5, 1.5}, {0., 0.5, 0.}, 1);
    autoPas->addParticle(p2);
    for (auto iter = autoPas->getContainer()->begin() ; iter.isValid(); ++iter) {
        cout << iter->toString() << endl;
        cout << "ParticleOldF= " << arrayString(iter->getOldf()) << endl;
    }
    double particleD = 0.01;
    int iterations=0;
    //iterationen beginnend
    TimeDiscretization<decltype(autoPas)> td(particleD);
    auto* functor = new autopas::LJFunctor<PrintableMolecule,autopas::ParticleCell<PrintableMolecule>, autopas::FunctorN3Modes::Both, true>(cutoff, epsilon, sigma, 0.0,boxmin,boxmax,true);
    //domain vorbeireiten: -Force initialisieren
    autoPas->iteratePairwise(functor);
    //Dokumentation prints
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "-----AFTER INITIALIZATION----" << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    for (auto iter = autoPas->getContainer()->begin() ; iter.isValid(); ++iter) {
        cout << iter->toString() << "  __END" << endl;
        cout << "ParticleOldF= " << arrayString(iter->getOldf()) << endl;
    }
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "-------ITERATIONS START------" << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;

    while(iterations < 10){
        td.VSCalculateX(autoPas);
        autoPas->iteratePairwise(functor);
        td.VSCalculateV(autoPas);
        iterations++;
        for (auto iter = autoPas->getContainer()->begin() ; iter.isValid(); ++iter) {
            cout << iter->toString() << endl;
            cout << "ParticleOldF= " << arrayString(iter->getOldf()) << endl;
        }
        cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
        cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
        cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;

    }
    ASSERT_TRUE(true);
    delete functor;
    delete autoPas;
}




//TEST(TimeDTEst,CalcV){}