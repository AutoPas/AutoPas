/**
 * @file ParallelVtkWriter.cpp
 * @author J. Körner
 * @date 31.05.2021
 */
#include "ParallelVtkWriter.h"

#include <ctime>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>

#include "mpi.h"

ParallelVtkWriter::ParallelVtkWriter(const std::string &sessionName, const std::string &outputFolder)
  : _sessionName(sessionName) {
  MPI_Comm_rank(MPI_COMM_WORLD, &_mpiRank);

  if (_mpiRank == 0) {
    tryCreateSessionFolder(_sessionName, outputFolder);

    //int numberOfProcesses;
    //MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

    // @todo: create paralle vtk file
  }

  MPI_Bcast(_sessionFolderPath.data(), _sessionFolderPath.size(), MPI_CHAR, 0, MPI_COMM_WORLD);
  
  //_processFolderPath = "Process" + _mpiRank;
  //tryCreateFolder(_processFolderPath, _sessionFolderPath);
}

void ParallelVtkWriter::recordTimestep(const int &currentIteration, const int &maximumNumberOfDigitsInIteration, const autopas::AutoPas<ParticleType> &autoPasContainer){
  std::ostringstream timestepFileName;
  timestepFileName << _sessionFolderPath << "/" << _sessionName << _mpiRank << std::setfill('0')
    << std::setw(maximumNumberOfDigitsInIteration) << currentIteration << ".vtk";

  std::ofstream timestepFile;
  timestepFile.open(timestepFileName.str(), std::ios::out | std::ios::binary);

  if (not timestepFile.is_open()) {
    throw std::runtime_error("Simulation::writeVTKFile(): Failed to open file \"" + timestepFileName.str() + "\"");
  }

  timestepFile << "# vtk DataFile Version 2.0\n"
          << "Timestep\n"
          << "ASCII\n";

  const auto numberOfParticles = autoPasContainer.getNumberOfParticles(autopas::IteratorBehavior::owned);

  // print positions
  timestepFile << "DATASET STRUCTURED_GRID\n"
          << "DIMENSIONS 1 1 1\n"
          << "POINTS " << numberOfParticles << " double\n";

  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    auto pos = particle->getR();
    timestepFile << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
  }
  timestepFile << "\n";

  timestepFile << "POINT_DATA " << numberOfParticles << "\n";
  // print velocities
  timestepFile << "VECTORS velocities double\n";
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    auto v = particle->getV();
    timestepFile << v[0] << " " << v[1] << " " << v[2] << "\n";
  }
  timestepFile << "\n";

  // print Forces
  timestepFile << "VECTORS forces double\n";
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    auto f = particle->getF();
    timestepFile << f[0] << " " << f[1] << " " << f[2] << "\n";
  }
  timestepFile << "\n";

  // print TypeIDs
  timestepFile << "SCALARS typeIds int\n";
  timestepFile << "LOOKUP_TABLE default\n";
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    timestepFile << particle->getTypeId() << "\n";
  }
  timestepFile << "\n";

  // print TypeIDs
  timestepFile << "SCALARS particleIds int\n";
  timestepFile << "LOOKUP_TABLE default\n";
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    timestepFile << particle->getID() << "\n";
  }
  timestepFile << "\n";

  timestepFile.close();
}

void ParallelVtkWriter::tryCreateSessionFolder(const std::string &name, const std::string location){
  time_t rawTime;
  time(&rawTime);

  struct tm* timeInformation;
  timeInformation = localtime(&rawTime);

  char buffer[80];
  strftime(buffer,sizeof(buffer),"%d-%m-%Y %H:%M:%S",timeInformation);
  std::string timeString(buffer);

  std::string _sessionFolderPath = name + "_" + timeString + location;
  tryCreateFolder(name + "_" + timeString, location);
}

void ParallelVtkWriter::tryCreateFolder(const std::string &name, const std::string location){
  try {
    std::filesystem::path newDirectoryPath(location + "/" + name);
    std::filesystem::create_directory(newDirectoryPath);
  }
  catch(std::filesystem::filesystem_error const& ex) {
    throw std::runtime_error("The output location " + location + " passed to ParallelVtkWriter is invalid");
  }
}

