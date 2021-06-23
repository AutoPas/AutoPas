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

ParallelVtkWriter::ParallelVtkWriter(std::string sessionName, const std::string &outputFolder, const int &maximumNumberOfDigitsInIteration)
    : _sessionName(std::move(sessionName)), _mpiRank(0), _maximumNumberOfDigitsInIteration(maximumNumberOfDigitsInIteration) {
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &_numberOfRanks);
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &_mpiRank);

  if (_mpiRank == 0) {
    tryCreateSessionAndDataFolders(_sessionName, outputFolder);
  }

  size_t sessionFolderPathLength = _sessionFolderPath.size();
  autopas::AutoPas_MPI_Bcast(&sessionFolderPathLength, 1, AUTOPAS_MPI_INT, 0, AUTOPAS_MPI_COMM_WORLD);
  
  size_t dataFolderPathLength = _dataFolderPath.size();
  autopas::AutoPas_MPI_Bcast(&dataFolderPathLength, 1, AUTOPAS_MPI_INT, 0, AUTOPAS_MPI_COMM_WORLD);

  if (_mpiRank != 0) {
    _sessionFolderPath.resize(sessionFolderPathLength);
    _dataFolderPath.resize(dataFolderPathLength);
  }

  autopas::AutoPas_MPI_Bcast(&_sessionFolderPath[0], sessionFolderPathLength, AUTOPAS_MPI_CHAR, 0,
                             AUTOPAS_MPI_COMM_WORLD);
  autopas::AutoPas_MPI_Bcast(&_dataFolderPath[0], dataFolderPathLength, AUTOPAS_MPI_CHAR, 0,
                             AUTOPAS_MPI_COMM_WORLD);
}

void ParallelVtkWriter::recordTimestep(const int &currentIteration, const autopas::AutoPas<ParticleType> &autoPasContainer) {
  if (_mpiRank == 0) {
    createPvtuFile(currentIteration);
  }

  std::ostringstream timestepFileName;
  timestepFileName << _dataFolderPath << _sessionName << "_" << _mpiRank << "_" << std::setfill('0')
                   << std::setw(_maximumNumberOfDigitsInIteration) << currentIteration << ".vtu";

  std::ofstream timestepFile;
  timestepFile.open(timestepFileName.str(), std::ios::out | std::ios::binary);

  if (not timestepFile.is_open()) {
    throw std::runtime_error("Simulation::writeVTKFile(): Failed to open file \"" + timestepFileName.str() + "\"");
  }

  const auto numberOfParticles = autoPasContainer.getNumberOfParticles(autopas::IteratorBehavior::owned);

  timestepFile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>" << std::endl;
  timestepFile << "<VTKFile byte_order=\"LittleEndian\" type=\"UnstructuredGrid\" version=\"0.1\">" << std::endl;
  timestepFile << "  <UnstructuredGrid>" << std::endl;
  timestepFile << "    <Piece NumberOfCells=\"0\" NumberOfPoints=\"" << numberOfParticles << "\">" << std::endl;
  timestepFile << "      <PointData>" << std::endl;

  // print velocities
  timestepFile << "        <DataArray Name=\"velocities\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">" << std::endl;
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    auto v = particle->getV();
    timestepFile << "        " << v[0] << " " << v[1] << " " << v[2] << std::endl;
  }
  timestepFile << "        </DataArray>" << std::endl;

  // print forces
  timestepFile << "        <DataArray Name=\"forces\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">" << std::endl;
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    auto f = particle->getF();
    timestepFile << "        " << f[0] << " " << f[1] << " " << f[2] << std::endl;
  }
  timestepFile << "        </DataArray>" << std::endl;

  // print type ids
  timestepFile << "        <DataArray Name=\"typeIds\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\">" << std::endl;
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    timestepFile << "        " <<  particle->getTypeId() << std::endl;
  }
  timestepFile << "        </DataArray>" << std::endl;

  // print ids
  timestepFile << "        <DataArray Name=\"ids\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\">" << std::endl;
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    timestepFile << "        " << particle->getID() << std::endl;;
  }
  timestepFile << "        </DataArray>" << std::endl;

  timestepFile << "      </PointData>" << std::endl;
  timestepFile << "      <CellData/>" << std::endl;
  timestepFile << "      <Points>" << std::endl;

  // print positions
  timestepFile << "        <DataArray Name=\"position\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">" << std::endl;
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    auto pos = particle->getR();
    timestepFile << "        " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
  }
  timestepFile << "        </DataArray>" << std::endl;

  timestepFile << "      </Points>" << std::endl;
  timestepFile << "      <Cells>" << std::endl;
  timestepFile << "        <DataArray Name=\"types\" NumberOfComponents=\"0\" format=\"ascii\" type=\"Float32\"/>" << std::endl;
  timestepFile << "      </Cells>" << std::endl;
  timestepFile << "    </Piece>" << std::endl;
  timestepFile << "  </UnstructuredGrid>" << std::endl;
  timestepFile << "</VTKFile>" << std::endl;

  timestepFile.close();
}

void ParallelVtkWriter::tryCreateSessionAndDataFolders(const std::string &name, std::string location) {
  time_t rawTime;
  time(&rawTime);

  struct tm *timeInformation;
  timeInformation = localtime(&rawTime);

  char buffer[80];
  strftime(buffer, sizeof(buffer), "%d%m%Y_%H%M%S", timeInformation);
  std::string timeString(buffer);

  _sessionFolderPath = location + "/" + name + "_" + timeString + "/";
  tryCreateFolder(name + "_" + timeString, location);

  _dataFolderPath = _sessionFolderPath + "data/";
  tryCreateFolder("data", _sessionFolderPath);
}

void ParallelVtkWriter::createPvtuFile(const int &currentIteration){
  std::ostringstream filename;
  filename << _sessionFolderPath << _sessionName << "_" << std::setfill('0') << std::setw(_maximumNumberOfDigitsInIteration)
    << currentIteration << ".pvtu";

  std::ofstream timestepFile;
  timestepFile.open(filename.str(), std::ios::out | std::ios::binary);

  if (not timestepFile.is_open()) {
    throw std::runtime_error("Simulation::writeVTKFile(): Failed to open file \"" + filename.str() + "\"");
  }

  timestepFile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>" << std::endl;
  timestepFile << "<VTKFile byte_order=\"LittleEndian\" type=\"PUnstructuredGrid\" version=\"0.1\">" << std::endl;
  timestepFile << "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
  timestepFile << "    <PPointData>" << std::endl;
  timestepFile << "      <PDataArray Name=\"velocities\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\"/>" << std::endl;
  timestepFile << "      <PDataArray Name=\"forces\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Int32\"/>" << std::endl;
  timestepFile << "      <PDataArray Name=\"typeIds\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\"/>" << std::endl;
  timestepFile << "      <PDataArray Name=\"ids\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\"/>" << std::endl;
  timestepFile << "    </PPointData>" << std::endl;
  timestepFile << "    <PCellData/>" << std::endl;
  timestepFile << "    <PPoints>" << std::endl;
  timestepFile << "      <PDataArray Name=\"points\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\"/>" << std::endl;
  timestepFile << "    </PPoints>" << std::endl;
  timestepFile << "    <PCells>" << std::endl;
  timestepFile << "      <PDataArray Name=\"types\" NumberOfComponents=\"0\" format=\"ascii\" type=\"Float32\"/>" << std::endl;
  timestepFile << "    </PCells>" << std::endl;

  for (int i = 0; i < _numberOfRanks; ++i){
    timestepFile << "    <Piece Source=\"./data/" << _sessionName << "_" << i << "_" << std::setfill('0') << std::setw(_maximumNumberOfDigitsInIteration) << currentIteration << ".vtu\"/>" << std::endl;
  }

  timestepFile << "  </PUnstructuredGrid>" << std::endl;
  timestepFile << "</VTKFile>" << std::endl;

  timestepFile.close();
}

void ParallelVtkWriter::tryCreateFolder(const std::string &name, const std::string &location) {
  try {
    std::filesystem::path newDirectoryPath(location + "/" + name);
    std::filesystem::create_directory(newDirectoryPath);
  } catch (std::filesystem::filesystem_error const &ex) {
    throw std::runtime_error("The output location " + location + " passed to ParallelVtkWriter is invalid");
  }
}
