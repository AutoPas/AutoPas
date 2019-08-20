/**
 * @file PeriodicBoundariesTest.h
 * @author N. Fottner
 * @date 2/8/19
 */

#pragma once
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "../../../../examples/md-flexible/BoundaryConditions.h"
#include "../../../../examples/md-flexible/Generator.h"
#include "../../../../examples/md-flexible/Objects.h"
#include "../../../../examples/md-flexible/PrintableMolecule.h"
#include "../../../../examples/md-flexible/TimeDiscretization.h"
#include "../../../../examples/md-flexible/Simulation.h"
#include "../../../../examples/md-flexible/YamlParser.h"
#include "../../../../src/autopas/utils/ArrayMath.h"
#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"

class SimulationTest : public AutoPasTestBase {
 public:
  SimulationTest()
      : AutoPasTestBase(),
        _parser{std::make_shared<YamlParser>()}
      {}

 static void initFillWithParticles(autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autopas,
                           std::array<unsigned long, 3> particlesPerDim,double particleSpacing,double cutoff);

 void VisualizeSmallSzenario(std::array<size_t,3> particlesPerDim, double cutoff, double particleSpacing, double epsilon, double sigma, double mass, int iterations, double delta_t, const std::string &filename);

 static double distanceBetween2Points(const std::array<double, 3>& iPos,const std::array<double, 3>& jPos);

 static void initWithTwoParticlesWithDistance(autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autopas, double distance, double cutoff);

    /**Prints state of current Iteration of Simulation as .vtu file
    * */
  template <class AutoPasTemplate>
  void writeVTKFile(AutoPasTemplate &autopas,size_t iteration,const std::string &filename) {
    using namespace std;
        std::stringstream strstr;
        auto maxNumDigits = 4;
        strstr << filename << "_" << std::setfill('0') << std::setw(maxNumDigits) << iteration << ".vtu";
        std::ofstream vtkFile;
        vtkFile.open(strstr.str());

        vtkFile << "# vtk DataFile Version 2.0" << std::endl;
        vtkFile << "Timestep" << std::endl;
        vtkFile << "ASCII" << std::endl;
        vtkFile << "DATASET STRUCTURED_GRID" << std::endl;
        vtkFile << "DIMENSIONS 1 1 1" << std::endl;
        vtkFile << "POINTS " << autopas.getNumberOfParticles() << " double" << std::endl;

        for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
            auto pos = iter->getR();
            vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
        }
        vtkFile.close();
  }

 protected:
  std::shared_ptr<YamlParser> _parser;
};
