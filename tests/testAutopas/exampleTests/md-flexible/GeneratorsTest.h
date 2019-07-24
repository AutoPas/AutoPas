
#pragma once
#include <gtest/gtest.h>
#include <math.h>
#include <vector>
#include "../../../../examples/md-flexible/Generator.h"
#include "../../../../examples/md-flexible/ParticleClassLibrary.h"
#include "../../../../examples/md-flexible/PrintableMolecule.h"
#include "../../../../examples/md-flexible/TimeDiscretization.h"
#include "../../../../src/autopas/utils/ArrayMath.h"
#include "../../../../examples/md-flexible/Objects.h"
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"
#include "../../../../examples/md-flexible/YamlParser.h"
#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
class GeneratorsTest : public AutoPasTestBase {
 public:
  GeneratorsTest()
      : AutoPasTestBase(),
        epsilon{1.0},
        sigma{1.0},
        cutoff{1.},
        boxmin{{0., 0., 0.}},
        boxmax{{5., 5., 5.}},
        PCL{ParticleClassLibrary(epsilon, sigma, 1.0, 800)},
        functor{autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>,
                                   autopas::FunctorN3Modes::Both, true>(cutoff, PCL, 0.0)},
        parser{YamlParser()},
        filename{"testParsing.yaml"}
        { parser.parseYamlFile();}

  void MolSimTaskGeneration(autopas::AutoPas<Particle, FPCell> &autopas);

  template <class AutoPasTemplate>
  void writeVTKFile(string &filename, size_t numParticles, AutoPasTemplate &autopas) {
    stringstream strstr;
    strstr << filename;
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

  double L2Norm(std::array<double, 3> array) {
    double square_sum = 0;
    for (auto e : array) {
      square_sum += (e * e);
    }
    return sqrt(square_sum);
  }

 protected:
  double epsilon;
  double sigma;
  double cutoff;
  array<double, 3> boxmin;
  array<double, 3> boxmax;
  ParticleClassLibrary PCL;
  autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>, autopas::FunctorN3Modes::Both, true>
      functor;
    YamlParser parser;
    std::string filename;
};