/**
 * @file ParallelVtkWriter.cpp
 * @author J. Körner
 * @date 31.05.2021
 */
#include "ParallelVtkWriter.h"

#include <cstddef>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <utility>

#include "autopas/utils/WrapMPI.h"

ParallelVtkWriter::ParallelVtkWriter(std::string sessionName, const std::string &outputFolder)
    : _sessionName(std::move(sessionName)), _mpiRank(0) {
  // This using directive is necessary, because 'autopas::AUTOPAS_...' variables defined in WrapMPI.h do not exist
  // when compiling with MPI. When compiling without MPI the namespace prefix needs to be used.
  using namespace autopas;

  AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &_mpiRank);

  if (_mpiRank == 0) {
    tryCreateSessionFolder(_sessionName, outputFolder);
  }

  size_t sessionFolderPathLength = _sessionFolderPath.size();
  AutoPas_MPI_Bcast(&sessionFolderPathLength, 1, AUTOPAS_MPI_INT, 0, AUTOPAS_MPI_COMM_WORLD);

  if (_mpiRank != 0) {
    _sessionFolderPath.resize(sessionFolderPathLength);
  }

  AutoPas_MPI_Bcast(&_sessionFolderPath[0], sessionFolderPathLength, AUTOPAS_MPI_CHAR, 0, AUTOPAS_MPI_COMM_WORLD);
}

void ParallelVtkWriter::recordTimestep(const int &currentIteration, const int &maximumNumberOfDigitsInIteration,
                                       const autopas::AutoPas<ParticleType> &autoPasContainer) {
  std::ostringstream timestepFileName;
  timestepFileName << _sessionFolderPath << _sessionName << "_" << _mpiRank << "_" << std::setfill('0')
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

void ParallelVtkWriter::tryCreateSessionFolder(const std::string &name, std::string location) {
  time_t rawTime;
  time(&rawTime);

  struct tm *timeInformation;
  timeInformation = localtime(&rawTime);

  char buffer[80];
  strftime(buffer, sizeof(buffer), "%d%m%Y_%H%M%S", timeInformation);
  std::string timeString(buffer);

  _sessionFolderPath = location + "/" + name + "_" + timeString + "/";
  tryCreateFolder(name + "_" + timeString, location);
}

void ParallelVtkWriter::tryCreateFolder(const std::string &name, const std::string &location) {
  try {
    std::filesystem::path newDirectoryPath(location + "/" + name);
    std::filesystem::create_directory(newDirectoryPath);
  } catch (std::filesystem::filesystem_error const &ex) {
    throw std::runtime_error("The output location " + location + " passed to ParallelVtkWriter is invalid");
  }
}
