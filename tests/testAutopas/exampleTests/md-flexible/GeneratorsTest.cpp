//
// Created by nicola on 16.06.19.
//

#include <gtest/gtest.h>
#include "autopas/AutoPas.h"
#include "../../testingHelpers/GridGenerator.h"
#include "../../testingHelpers/RandomGenerator.h"
#include "../../testingHelpers/GaussianGenerator.h"
#include "testingHelpers/commonTypedefs.h"
#include "../../../../examples/md-flexible/PrintableMolecule.h"
#include "../../../../examples/md-flexible/TimeDiscretization.h"
#include "../../../../examples/md-flexible/MDFlexParser.h"
#include <math.h>
#include <vector>
#include "../../../../src/autopas/utils/ArrayMath.h"

using namespace std;
using namespace autopas;

//hier die generatoren nur auf PrintableMolecule ausgelegt

void initContainerGrid(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                                                           size_t particlesPerDim, double particelSpacing) {
    std::array<double, 3> boxMin({0., 0., 0.});
    std::array<double, 3> boxMax(
            {(particlesPerDim)*particelSpacing, (particlesPerDim)*particelSpacing, (particlesPerDim)*particelSpacing});

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    PrintableMolecule dummyParticle;
    GridGenerator::fillWithParticles(autopas, {particlesPerDim, particlesPerDim, particlesPerDim}, dummyParticle,
                                     {particelSpacing, particelSpacing, particelSpacing},
                                     {particelSpacing / 2, particelSpacing / 2, particelSpacing / 2});
}


void initContainerGauss(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                                                            double boxLength, size_t numParticles, double distributionMean, double distributionStdDev) {
    std::array<double, 3> boxMin({0., 0., 0.});
    std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    PrintableMolecule dummyParticle;
    GaussianGenerator::fillWithParticles(autopas, numParticles, dummyParticle, distributionMean, distributionStdDev);
}


void initContainerUniform(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                                                              double boxLength, size_t numParticles) {
    std::array<double, 3> boxMin({0., 0., 0.});
    std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    PrintableMolecule dummyParticle;
    RandomGenerator::fillWithParticles(autopas, dummyParticle, numParticles);
}


// FRAGE-FABIO
//wieso kann ich nicht Particle und ParticleCell im konstruktor behalten-> wenn ichs behalte-> no matching constructor fehler
void MolSimTaskGeneration(autopas::AutoPas<PrintableMolecule,FullParticleCell<PrintableMolecule>> &autopas){
    std::array<double, 3> boxMin({0., 0., 0.});
    std::array<double, 3> boxMax({40., 30., 1.});

    std::array<double,3> smallGridBoxMin({15.,15.,0});
    std::array<double,3> bigGridBoxMin({0.,0.,0.});

    std::array<double,3> smallGridBoxMax({23.,23.,0});
    std::array<double,3> bigGridBoxMax({40.,10.,0});

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();
    Particle dummyParticle;
    //small grid:
    //GridGenerator::fillWithParticles(autopas, {8, 8, 8}, dummyParticle,{1., 1., 1.},{0.5,0.5,0.5});

    //GridGenerator::fillWithParticlesOnR(autopas,smallGridBoxMin,{8.,8.,0.},dummyParticle,{1.,1.,1.},{0.5,0.5,0.5});


    //letzte ändereung-- particle zu printablemol. und particlecell zu fullparticlecell
    GridGenerator::fillWithParticlesOnR(autopas,smallGridBoxMin,{8,8,0},dummyParticle,{1.,1.,1.},{0.5,0.5,0.5});

}




TEST(Generator, Behavior){
    PrintableMolecule::setEpsilon(5.0);
    PrintableMolecule::setSigma(1.0);
    PrintableMolecule::setMass(1.0);
    auto *autoPas = new autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>(std::cout);
    double epsilon = 5.0;double sigma=1.0;double cutoff=2;array<double,3> boxmin={0.,0.,0.}; array<double,3> boxmax={5.,5.,5.};
    PrintableMolecule::setEpsilon(epsilon);
    PrintableMolecule::setSigma(sigma);
    PrintableMolecule::setMass(1.0);
    autoPas->setBoxMin(boxmin);
    autoPas->setBoxMax(boxmax);
    autoPas->setCutoff(cutoff);
    autoPas->setAllowedContainers({autopas::ContainerOption::linkedCells});
    autoPas->init();
    //for GRID generator:
    int particlePerDim = 5;
    double particleSpacing= 0.5;
    //initContainerGrid(*autoPas,particlePerDim,particleSpacing);
    cout << "Number of particles generated " << autoPas->getNumberOfParticles() << endl;
    //Uniform generator
    //initContainerUniform(*autoPas,5.,125);
    //Gauß generator
    //initContainerGauss(*autoPas,5.,125,5,2);

    MolSimTaskGeneration(*autoPas);
    for (auto iter = autoPas->getContainer()->begin() ; iter.isValid(); ++iter) {
        cout << iter->toString() << endl;
    }


    //print State -- ugly code , I know
    size_t numParticles = autoPas->getNumberOfParticles();
    string filename = "VtkTestOutput.vtu";
    std::ofstream vtkFile;
    vtkFile.open(filename);
    vtkFile << "# vtk DataFile Version 2.0" << endl;
    vtkFile << "Timestep" << endl;
    vtkFile << "ASCII" << endl;
    vtkFile << "DATASET STRUCTURED_GRID" << endl;
    vtkFile << "DIMENSIONS 1 1 1" << endl;
    vtkFile << "POINTS " << numParticles << " double" << endl;
    for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
        auto pos = iter->getR();
        vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << endl;
    }
    vtkFile.close();

    delete autoPas;
    ASSERT_TRUE(true);
}
