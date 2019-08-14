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
class PeriodicBoundariesTest : public AutoPasTestBase {
 public:
  PeriodicBoundariesTest()
      : AutoPasTestBase(),
        _parser{std::make_shared<YamlParser>()}
        {/* _parser->setFilename("periodic.yaml"); */
        _parser->setFilename("MolSimBlatt2Task3.yaml");
        _parser->parseYamlFile();
        _parser->setFilename("VtkPeriodicOutput");
        _simulation.initialize(_parser);
        }
        /**Prints state of current Iteration of Simulation as .vtu file
         * */
  template <class AutoPasTemplate>
  void writeVTKFile(std::string &vtkFilename, size_t numParticles, AutoPasTemplate &autopas) {
    using namespace std;
    stringstream strstr;
    strstr << vtkFilename;
    // string path = "./vtk";
    std::ofstream vtkFile;
    vtkFile.open(strstr.str());
    vtkFile << "# vtk DataFile Version 2.0" << endl;
    vtkFile << "Timestep" << endl;
    vtkFile << "ASCII" << endl;
    vtkFile << "DATASET STRUCTURED_GRID" << endl;
    vtkFile << "DIMENSIONS 1 1 1" << endl;
    vtkFile << "POINTS " << numParticles << " double" << endl;

    for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
      auto pos = iter->getR();
      vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << endl;
    }
    vtkFile.close();
  }

 protected:
  std::shared_ptr<YamlParser> _parser;
  Simulation<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> _simulation;

};
