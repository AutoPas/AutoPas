/**
 * @file ParallelVtkWriter.cpp
 * @author J. Körner
 * @date 31.05.2021
 */
#include "ParallelVtkWriter.h"

#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>

ParallelVtkWriter::ParallelVtkWriter(const std::string &sessionName, const std::string &outputFolder)
  : _sessionName(sessionName) {
  MPI_Comm_rank(MPI_COMM_WORLD, &_mpiRank);

  if (_mpiRank == 0) {
    tryCreateSessionFolder(_sessionName, outputFolder);

    int numberOfProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

    // @todo: create paralle vtk file
  }

  MPI_Bcast(_sessionFolderPath.data(), _sessionFolderPath.size(), MPI_CHAR, 0, MPI_COMM_WORLD)
  
  _processFolderPath = "Process" + _mpiRank;
  tryCreateFolder(_processFolderPath, _sessionFolderPath);
}

void ParallelVtkWriter::recordTimestep(const int &currentIteration, const size_t &maximumNumberOfDigitsInIteration, const std::vector<ParticleType> &particles){
  std::string timestepData;

  timestepData.append("# vtk DataFile Version 2.0\n");
  timestepData.append("Timestep\n");
  timestepData.append("ASCII\n");
  timestepData.append("DATASET STRUCTURED_GRID\n");

  // print positions
  timestepData.append("DIMENSIONS 1 1 1\n");
  timestepData.append("POINTS " + particles.size() + " double\n");

  for (const auto &particle : particles) {
    auto pos = particle->getR();
    timestepData.append(pos[0] + " " + pos[1] + " " + pos[2] + "\n");
  }
  timestepData.append("\n");

  // print velocities
  timestepData.append("POINT_DATA " + particles.size() + "\n");
  timestepData.append("VECTORS velocities double\n");
  for (const auto &particle : particles) {
    auto v = particle->getV();
    timestepData.append(v[0] + " " + v[1] + " " + v[2] + "\n");
  }
  timestepData.append("\n");

  // print Forces
  timestepData.append("VECTORS forces double\n");
  for (const auto &particle : particles) {
    auto f = particle->getF();
    timestepData.append(f[0] + " " + f[1] + " " + f[2] + "\n");
  }
  timestepData.append("\n");

  // print TypeIDs
  timestepData.append("SCALARS typeIds int\n";);
  timestepData.append("LOOKUP_TABLE default\n";);
  for (const auto &particle : particles) {
    timestepData.append(particle->getTypeId() + "\n");
  }
  timestepData.append("\n");

  // print TypeIDs
  timestepData.append("SCALARS particleIds int\n");
  timestepData.append("LOOKUP_TABLE default\n");
  for (const auto &particle : particles) {
    timestepData.append(particle->getID() + "\n");
  }
  timestepData.append("\n");

  std::string timestepFileName = _sessionName + "_" + _mpiRank + std::setfill('0')
    + std::setw(maximumNumberOfDigitsInIteration) + _iteration + ".vtk";

  ofstream timestepFile(timestepFileName, ios::out | ios::binary);
  timestepFile.write(timestepData.data(), timestepData.size());
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
    std::filesystem::create_directory(newDirecotryPath);
  }
  catch(std::filesystem::filesystem_error const& ex) {
    throw std::runtime_error("The output location " + location + " passed to ParallelVtkWriter is invalid");
  }
}

