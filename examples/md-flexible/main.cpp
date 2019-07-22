/**
 * @file main.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include "Simulation.h"

#include <autopas/utils/MemoryProfiler.h>
#include <yaml-cpp/yaml.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "PrintableMolecule.h"  // includes autopas.h
#include "YamlParser.h"
#include "autopas/AutoPas.h"
#include "autopas/pairwiseFunctors/LJFunctorAVX.h"
#include <vector>
using namespace std;
using namespace autopas;

/**
 * Prints position and forces of all particles in the autopas object.
 * @tparam AutoPasTemplate Template for the templetized autopas type.
 * @param autopas
 */
template <class AutoPasTemplate>
void printMolecules(AutoPasTemplate &autopas) {
  for (auto particleIterator = autopas.begin(); particleIterator.isValid(); ++particleIterator) {
    particleIterator->print();
  }
}

/** Writes a VTK file for the current state of the AutoPas object
 * @tparam AutoPasTemplate Template for the templetized autopas type.
 * @param filename
 * @param numParticles
 * @param autopas
 */
template <class AutoPasTemplate>
void writeVTKFile(string &filename, AutoPasTemplate &autopas) {
  std::ofstream vtkFile;
  vtkFile.open(filename);

  vtkFile << "# vtk DataFile Version 2.0" << endl;
  vtkFile << "Timestep" << endl;
  vtkFile << "ASCII" << endl;
  vtkFile << "DATASET STRUCTURED_GRID" << endl;
  vtkFile << "DIMENSIONS 1 1 1" << endl;
  vtkFile << "POINTS " << autopas.getNumberOfParticles() << " double" << endl;

  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    auto pos = iter->getR();
    vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << endl;
  }

  vtkFile.close();
}

int main(int argc, char **argv) {


    struct Object {
        std::array<double, 3> velocity;
        std::array<size_t, 3> particlesPerDim;
        double particleSpacing;
        std::array<double, 3> boxLength;
        size_t numParticles;
        double distributionMean;
        double distributionStdDev;
        std::array<double, 3> center;
        unsigned long id;
        int radius;
    };

    struct CubeGrid : Object {
//      CubeGrid(): particlesPerDim({20,20,20}), particleSpacing(.4),velocity({0.,0.,0.}){}
        std::array<size_t, 3> particlesPerDim;
        double particleSpacing;
        std::array<double, 3> velocity;
    };
    struct CubeGauss : Object {
        std::array<double, 3> boxLength;
        size_t numParticles;
        double distributionMean;
        double distributionStdDev;
        std::array<double, 3> velocity;
    };
    struct CubeUniform : Object {
        std::array<double, 3> boxLength;
        size_t numParticles;
        std::array<double, 3> velocity;
    };
    struct Sphere : Object {
        std::array<double, 3> center;
        int radius;
        double particleSpacing;
        unsigned long id;
        std::array<double, 3> velocity;
    };

    std::vector<Object> ObjectGenerator = {};

    std::string input = "parsingFile.yaml";
    YAML::Node config = YAML::LoadFile(input);
    //parsed Objecte die Generiert werden sollen

    if (config["Objects"]) {
            for (YAML::const_iterator it = config["Objects"].begin(); it != config["Objects"].end(); ++it) {
                if (it->first.as<std::string>() == "CubeGrid") {
                    CubeGrid C;
                    C.particlesPerDim={it->second["particles-per-Dim"][0].as<unsigned long>(),it->second["particles-per-Dim"][1].as<unsigned long>(),it->second["particles-per-Dim"][2].as<unsigned long>()};
                    C.particleSpacing= it->second["particleSpacing"].as<double>();
                    C.velocity= {it->second["velocity"][0].as<double>(),it->second["velocity"][1].as<double>(),it->second["velocity"][2].as<double>()};
                    ObjectGenerator.emplace_back(C);
                    continue;
                }
                if(it->first.as<std::string>() == "CubeGauss"){
                    CubeGauss C;
                    C.boxLength={it->second["box-length"][0].as<double>(),it->second["box-length"][1].as<double>(),it->second["box-length"][2].as<double>()};
                    C.numParticles = it->second["numberOfParticles"].as<size_t>();
                    C.distributionMean =it->second["distribution-mean"].as<double>();
                    C.distributionStdDev = it->second["distribution-stddeviation"].as<double>();
                    C.velocity= {it->second["velocity"][0].as<double>(),it->second["velocity"][1].as<double>(),it->second["velocity"][2].as<double>()};
                    ObjectGenerator.emplace_back(C);
                    continue;
                }
                if(it->first.as<std::string>() == "CubeUniform"){
                    CubeUniform C;
                    C.boxLength= {it->second["box-length"][0].as<double>(),it->second["box-length"][1].as<double>(),it->second["box-length"][2].as<double>()};
                    C.numParticles = it->second["numberOfParticles"].as<size_t>();
                    C.velocity= {it->second["velocity"][0].as<double>(),it->second["velocity"][1].as<double>(),it->second["velocity"][2].as<double>()};
                    ObjectGenerator.emplace_back(C);                    continue;

                }
                if(it->first.as<std::string>() =="Sphere"){
                    Sphere S;
                    S.center= {it->second["center"][0].as<double>(),it->second["center"][1].as<double>(),it->second["center"][2].as<double>()};
                    S.radius =it->second["radius"].as<int>();
                    S.particleSpacing = it->second["particleSpacing"].as<double>();
                    S.id=it->second["firstId"].as<unsigned long>();
                    S.velocity= {it->second["velocity"][0].as<double>(),it->second["velocity"][1].as<double>(),it->second["velocity"][2].as<double>()};
                    ObjectGenerator.emplace_back(S);
                    continue;
                }

            }
    }
    //@todo komische werte werden hier eingelesen!!!
    cout << ObjectGenerator.front().velocity[0] << endl;


//  Simulation<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> simulation;
//  YamlParser parser;
//  //@todo catch exception and errors for parser
//  //@todo parsing file über die command line übergeben?
//  std::string filename= "parsingFile.yaml";
//  parser.parseInput(filename);
//  parser.printConfig();
//  simulation.initialize(parser);
//  simulation.simulate();
//  simulation.printStatistics();
  // frage FABIO, wenn ich hier manuel den destructor von simlation aufrufe; wieso kriege ich 4 invalid reads(autopas
  // container-traversals usw) und 18 invalid free

  return EXIT_SUCCESS;
}
