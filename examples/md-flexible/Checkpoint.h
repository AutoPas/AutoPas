/**
 * @file Checkpoint.h
 * @author N. Fottner
 * @date 10/9/19
 */
#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include "autopas/AutoPas.h"

/**
 * This class implements the initialization of an AutoPas container from vtk checkpoint files
 */
template <class AutoPasTemplate, class Particle>
class Checkpoint {
 public:
  /**
   * Default constructor.
   */
  Checkpoint() = default;
  /**
   * Default destructor.
   */
  ~Checkpoint() = default;

  /**
   * Reads the Data of all particles from a Vtk file, and adds all Particles
   * with their properties into the AutoPas container
   * @param autopas
   * @param vtkFilename
   */
  static void initDomain(AutoPasTemplate &autopas, const std::string &vtkFilename);

  /**
   * helper function to go to a specific line of a file
   * @param file
   * @param index
   */
  static std::ifstream &GotoLine(std::ifstream &file, int index) {
    file.seekg(std::ios::beg);
    for (int i = 0; i < index - 1; ++i) {
      file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return file;
  }

  /**
   *  Reads from a vtkFile the data of @param numOfP Particles returns a vector of the the data collected
   *  @param file
   *  @param numOfP size of data set to be collected
   */
  static std::vector<std::array<double, 3>> readingData(std::ifstream &file, int numOfP);
};

template <class AutoPasTemplate, class Particle>
void Checkpoint<AutoPasTemplate, Particle>::initDomain(AutoPasTemplate &autopas, const std::string &vtkFilename) {
  // first: reading data from VtkFile:
  std::ifstream infile(vtkFilename);
  std::string extract;
  // offset of vtk header informations
  int dataOffset = 6;
  GotoLine(infile, dataOffset);
  std::string strNumberOfParticles;
  infile >> extract;
  // getting the numberOfParticles
  infile >> strNumberOfParticles;
  int numberOfParticles = std::stoi(strNumberOfParticles);
  // getting all positions for all particles
  GotoLine(infile, dataOffset + 1);
  // position
  std::vector<std::array<double, 3>> positions = readingData(infile, numberOfParticles);
  // skip 3 lines of header informations:
  std::getline(infile, extract);
  std::getline(infile, extract);
  std::getline(infile, extract);
  // velocities
  std::vector<std::array<double, 3>> velocities = readingData(infile, numberOfParticles);
  // skip 2 lines of header informations:
  std::getline(infile, extract);
  std::getline(infile, extract);
  // forces
  std::vector<std::array<double, 3>> forces = readingData(infile, numberOfParticles);

  // creating Particles for checkpoint:
  std::vector<Particle> particles;
  for (auto i = 0; i < numberOfParticles; i++) {
    Particle p;
    p.setR(positions.at(i));
    p.setV(velocities.at(i));
    p.setF(forces.at(i));
    p.setID(i);
    particles.emplace_back(p);
  }
  // second: adding all particle to AutoPas Object
  for (auto it = particles.begin(); it != particles.end(); ++it) {
    autopas.addParticle(*it);
  }
}

template <class AutoPasTemplate, class Particle>
std::vector<std::array<double, 3>> Checkpoint<AutoPasTemplate, Particle>::readingData(std::ifstream &file, int numOfP) {
  std::vector<std::array<double, 3>> data;
  std::string dataString;
  for (int i = 0; i < numOfP; i++) {
    std::getline(file, dataString);
    data.emplace_back(autopas::utils::StringUtils::parseArrayD3(dataString));
    dataString = "";
  }
  return data;
}
