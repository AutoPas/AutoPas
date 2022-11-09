/**
 * @file ParallelVtkWriter.cpp
 * @author J. Körner
 * @date 31.05.2021
 */
#include "ParallelVtkWriter.h"

#include <cstddef>
#include <fstream>
#include <ios>
#include <iostream>
#include <string>
#include <utility>

#include "autopas/utils/WrapMPI.h"

ParallelVtkWriter::ParallelVtkWriter(std::string sessionName, const std::string &outputFolder,
                                     const int &maximumNumberOfDigitsInIteration)
    : _sessionName(std::move(sessionName)),
      _mpiRank(0),
      _maximumNumberOfDigitsInIteration(maximumNumberOfDigitsInIteration) {
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
  autopas::AutoPas_MPI_Bcast(&_dataFolderPath[0], dataFolderPathLength, AUTOPAS_MPI_CHAR, 0, AUTOPAS_MPI_COMM_WORLD);
}

void ParallelVtkWriter::recordTimestep(size_t currentIteration,
                                       const autopas::AutoPas<ParticleType> &autoPasContainer,
                                       const RegularGridDecomposition &decomposition) {
  recordParticleStates(currentIteration, autoPasContainer);
  recordDomainSubdivision(currentIteration, autoPasContainer.getCurrentConfig(), decomposition);
}

/**
 * @todo: Currently this function runs over all the particles for each property separately.
 * This can be improved by using multiple string streams (one for each property).
 * The streams can be combined to a single output stream after iterating over the particles, once.
 */
void ParallelVtkWriter::recordParticleStates(size_t currentIteration,
                                             const autopas::AutoPas<ParticleType> &autoPasContainer) {
  if (_mpiRank == 0) {
    createPvtuFile(currentIteration);
  }

  std::ostringstream timestepFileName;
  generateFilename("vtu", currentIteration, timestepFileName);

  std::ofstream timestepFile;
  timestepFile.open(timestepFileName.str(), std::ios::out | std::ios::binary);

  if (not timestepFile.is_open()) {
    throw std::runtime_error("Simulation::writeVTKFile(): Failed to open file \"" + timestepFileName.str() + "\"");
  }

  const auto numberOfParticles = autoPasContainer.getNumberOfParticles(autopas::IteratorBehavior::owned);

  timestepFile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n";
  timestepFile << "<VTKFile byte_order=\"LittleEndian\" type=\"UnstructuredGrid\" version=\"0.1\">\n";
  timestepFile << "  <UnstructuredGrid>\n";
  timestepFile << "    <Piece NumberOfCells=\"0\" NumberOfPoints=\"" << numberOfParticles << "\">\n";
  timestepFile << "      <PointData>\n";

  // print velocities
  timestepFile
      << "        <DataArray Name=\"velocities\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">\n";
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    auto v = particle->getV();
    timestepFile << "        " << v[0] << " " << v[1] << " " << v[2] << "\n";
  }
  timestepFile << "        </DataArray>\n";

  // print forces
  timestepFile << "        <DataArray Name=\"forces\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">\n";
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    auto f = particle->getF();
    timestepFile << "        " << f[0] << " " << f[1] << " " << f[2] << "\n";
  }
  timestepFile << "        </DataArray>\n";

  // print type ids
  timestepFile << "        <DataArray Name=\"typeIds\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">\n";
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    timestepFile << "        " << particle->getTypeId() << "\n";
  }
  timestepFile << "        </DataArray>\n";

  // print ids
  timestepFile << "        <DataArray Name=\"ids\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">\n";
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    timestepFile << "        " << particle->getID() << "\n";
    ;
  }
  timestepFile << "        </DataArray>\n";

  timestepFile << "      </PointData>\n";
  timestepFile << "      <CellData/>\n";
  timestepFile << "      <Points>\n";

  // print positions
  timestepFile << "        <DataArray Name=\"positions\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">\n";
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    auto pos = particle->getR();
    timestepFile << "        " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
  }
  timestepFile << "        </DataArray>\n";

  timestepFile << "      </Points>\n";
  timestepFile << "      <Cells>\n";
  timestepFile << "        <DataArray Name=\"types\" NumberOfComponents=\"0\" format=\"ascii\" type=\"Float32\"/>\n";
  timestepFile << "      </Cells>\n";
  timestepFile << "    </Piece>\n";
  timestepFile << "  </UnstructuredGrid>\n";
  timestepFile << "</VTKFile>\n";

  timestepFile.close();
}

void ParallelVtkWriter::recordDomainSubdivision(size_t currentIteration,
                                                const autopas::Configuration &autoPasConfiguration,
                                                const RegularGridDecomposition &decomposition) {
  if (_mpiRank == 0) {
    createPvtsFile(currentIteration, decomposition);
  }

  std::ostringstream timestepFileName;
  generateFilename("vts", currentIteration, timestepFileName);

  std::ofstream timestepFile;
  timestepFile.open(timestepFileName.str(), std::ios::out | std::ios::binary);

  if (not timestepFile.is_open()) {
    throw std::runtime_error("Simulation::writeVTKFile(): Failed to open file \"" + timestepFileName.str() + "\"");
  }

  const std::array<int, 6> wholeExtent = calculateWholeExtent(decomposition);
  const std::array<double, 3> localBoxMin = decomposition.getLocalBoxMin();
  const std::array<double, 3> localBoxMax = decomposition.getLocalBoxMax();

  auto printDataArray = [&](const auto &data, const std::string &type, const std::string name) {
    timestepFile << "        <DataArray type=\"" << type << "\" Name=\"" << name << "\" format=\"ascii\">\n";
    timestepFile << "          " << data << "\n";
    timestepFile << "        </DataArray>\n";
  };

  timestepFile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n";
  timestepFile << "<VTKFile byte_order=\"LittleEndian\" type=\"StructuredGrid\" version=\"0.1\">\n";
  timestepFile << "  <StructuredGrid WholeExtent=\"" << wholeExtent[0] << " " << wholeExtent[1] << " " << wholeExtent[2]
               << " " << wholeExtent[3] << " " << wholeExtent[4] << " " << wholeExtent[5] << "\">\n";
  timestepFile << "    <Piece Extent=\"" << wholeExtent[0] << " " << wholeExtent[1] << " " << wholeExtent[2] << " "
               << wholeExtent[3] << " " << wholeExtent[4] << " " << wholeExtent[5] << "\">\n";
  timestepFile << "      <CellData>\n";
  printDataArray(decomposition.getDomainIndex(), "Int32", "DomainId");
  printDataArray(autoPasConfiguration.cellSizeFactor, "Float32", "CellSizeFactor");
  printDataArray(static_cast<int>(autoPasConfiguration.container), "Int32", "Container");
  printDataArray(static_cast<int>(autoPasConfiguration.dataLayout), "Int32", "DataLayout");
  printDataArray(static_cast<int>(autoPasConfiguration.loadEstimator), "Int32", "LoadEstimator");
  printDataArray(static_cast<int>(autoPasConfiguration.traversal), "Int32", "Traversal");
  printDataArray(static_cast<int>(autoPasConfiguration.newton3), "Int32", "Newton3");
  printDataArray(_mpiRank, "Int32", "Rank");
  timestepFile << "      </CellData>\n";
  timestepFile << "      <Points>\n";
  timestepFile << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  timestepFile << "          " << localBoxMin[0] << " " << localBoxMin[1] << " " << localBoxMin[2] << "\n";
  timestepFile << "          " << localBoxMax[0] << " " << localBoxMin[1] << " " << localBoxMin[2] << "\n";
  timestepFile << "          " << localBoxMin[0] << " " << localBoxMax[1] << " " << localBoxMin[2] << "\n";
  timestepFile << "          " << localBoxMax[0] << " " << localBoxMax[1] << " " << localBoxMin[2] << "\n";
  timestepFile << "          " << localBoxMin[0] << " " << localBoxMin[1] << " " << localBoxMax[2] << "\n";
  timestepFile << "          " << localBoxMax[0] << " " << localBoxMin[1] << " " << localBoxMax[2] << "\n";
  timestepFile << "          " << localBoxMin[0] << " " << localBoxMax[1] << " " << localBoxMax[2] << "\n";
  timestepFile << "          " << localBoxMax[0] << " " << localBoxMax[1] << " " << localBoxMax[2] << "\n";
  timestepFile << "        </DataArray>\n";
  timestepFile << "      </Points>\n";
  timestepFile << "    </Piece>\n";
  timestepFile << "  </StructuredGrid>\n";
  timestepFile << "</VTKFile>\n";

  timestepFile.close();
}

std::array<int, 6> ParallelVtkWriter::calculateWholeExtent(const RegularGridDecomposition &domainDecomposition) {
  std::array<int, 6> wholeExtent;
  std::array<int, 3> domainId = domainDecomposition.getDomainId();
  std::array<int, 3> decomposition = domainDecomposition.getDecomposition();
  for (int i = 0; i < 3; ++i) {
    wholeExtent[2 * i] = domainId[i];
    wholeExtent[2 * i + 1] = std::min(domainId[i] + 1, decomposition[i]);
  }
  return wholeExtent;
}

void ParallelVtkWriter::tryCreateSessionAndDataFolders(const std::string &name, std::string location) {
  time_t rawTime;
  time(&rawTime);

  struct tm timeInformation;
  gmtime_r(&rawTime, &timeInformation);

  char buffer[80];
  strftime(buffer, sizeof(buffer), "%d%m%Y_%H%M%S", &timeInformation);
  std::string timeString(buffer);

  if (not checkFileExists(location)) {
    tryCreateFolder(location, "./");
  }

  _sessionFolderPath = location + "/" + name + "_" + timeString + "/";
  tryCreateFolder(name + "_" + timeString, location);

  _dataFolderPath = _sessionFolderPath + "data/";
  tryCreateFolder("data", _sessionFolderPath);
}

void ParallelVtkWriter::createPvtuFile(size_t currentIteration) {
  std::ostringstream filename;
  filename << _sessionFolderPath << _sessionName << "_" << std::setfill('0')
           << std::setw(_maximumNumberOfDigitsInIteration) << currentIteration << ".pvtu";

  std::ofstream timestepFile;
  timestepFile.open(filename.str(), std::ios::out | std::ios::binary);

  if (not timestepFile.is_open()) {
    throw std::runtime_error("Simulation::writeVTKFile(): Failed to open file \"" + filename.str() + "\"");
  }

  timestepFile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n";
  timestepFile << "<VTKFile byte_order=\"LittleEndian\" type=\"PUnstructuredGrid\" version=\"0.1\">\n";
  timestepFile << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
  timestepFile << "    <PPointData>\n";
  timestepFile
      << "      <PDataArray Name=\"velocities\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\"/>\n";
  timestepFile << "      <PDataArray Name=\"forces\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\"/>\n";
  timestepFile << "      <PDataArray Name=\"typeIds\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\"/>\n";
  timestepFile << "      <PDataArray Name=\"ids\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\"/>\n";
  timestepFile << "    </PPointData>\n";
  timestepFile << "    <PCellData/>\n";
  timestepFile << "    <PPoints>\n";
  timestepFile << "      <PDataArray Name=\"positions\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\"/>\n";
  timestepFile << "    </PPoints>\n";
  timestepFile << "    <PCells>\n";
  timestepFile << "      <PDataArray Name=\"types\" NumberOfComponents=\"0\" format=\"ascii\" type=\"Float32\"/>\n";
  timestepFile << "    </PCells>\n";

  for (int i = 0; i < _numberOfRanks; ++i) {
    timestepFile << "    <Piece Source=\"./data/" << _sessionName << "_" << i << "_" << std::setfill('0')
                 << std::setw(_maximumNumberOfDigitsInIteration) << currentIteration << ".vtu\"/>\n";
  }

  timestepFile << "  </PUnstructuredGrid>\n";
  timestepFile << "</VTKFile>\n";

  timestepFile.close();
}

void ParallelVtkWriter::createPvtsFile(size_t currentIteration, const RegularGridDecomposition &decomposition) {
  std::ostringstream filename;
  filename << _sessionFolderPath << _sessionName << "_" << std::setfill('0')
           << std::setw(_maximumNumberOfDigitsInIteration) << currentIteration << ".pvts";

  std::ofstream timestepFile;
  timestepFile.open(filename.str(), std::ios::out | std::ios::binary);

  if (not timestepFile.is_open()) {
    throw std::runtime_error("Simulation::writeVTKFile(): Failed to open file \"" + filename.str() + "\"");
  }
  const std::array<int, 3> wholeExtent = decomposition.getDecomposition();
  const std::array<double, 3> globalBoxMin = decomposition.getGlobalBoxMin();
  const std::array<double, 3> globalBoxMax = decomposition.getGlobalBoxMax();
  timestepFile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n";
  timestepFile << "<VTKFile byte_order=\"LittleEndian\" type=\"PStructuredGrid\" version=\"0.1\">\n";
  timestepFile << "  <PStructuredGrid WholeExtent=\"0 " << wholeExtent[0] << " 0 " << wholeExtent[1] << " 0 "
               << wholeExtent[2] << "\" GhostLevel=\"0\">\n";
  timestepFile << "    <PPointData/>\n";
  timestepFile << "    <PCellData>\n";
  timestepFile << "      <PDataArray type=\"Int32\" Name=\"DomainId\" />\n";
  timestepFile << "      <PDataArray type=\"Float32\" Name=\"CellSizeFactor\" />\n";
  timestepFile << "      <PDataArray type=\"Int32\" Name=\"Container\" />\n";
  timestepFile << "      <PDataArray type=\"Int32\" Name=\"DataLayout\" />\n";
  timestepFile << "      <PDataArray type=\"Int32\" Name=\"FullConfiguration\" />\n";
  timestepFile << "      <PDataArray type=\"Int32\" Name=\"LoadEstimator\" />\n";
  timestepFile << "      <PDataArray type=\"Int32\" Name=\"Traversal\" />\n";
  timestepFile << "      <PDataArray type=\"Int32\" Name=\"Newton3\" />\n";
  timestepFile << "      <PDataArray type=\"Int32\" Name=\"Rank\" />\n";
  timestepFile << "    </PCellData>\n";
  timestepFile << "    <PPoints>\n";
  timestepFile << "      <DataArray NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">\n";
  timestepFile << "        " << globalBoxMin[0] << " " << globalBoxMin[1] << " " << globalBoxMin[2] << "\n";
  timestepFile << "        " << globalBoxMax[0] << " " << globalBoxMin[1] << " " << globalBoxMin[2] << "\n";
  timestepFile << "        " << globalBoxMin[0] << " " << globalBoxMax[1] << " " << globalBoxMin[2] << "\n";
  timestepFile << "        " << globalBoxMax[0] << " " << globalBoxMax[1] << " " << globalBoxMin[2] << "\n";
  timestepFile << "        " << globalBoxMin[0] << " " << globalBoxMin[1] << " " << globalBoxMax[2] << "\n";
  timestepFile << "        " << globalBoxMax[0] << " " << globalBoxMin[1] << " " << globalBoxMax[2] << "\n";
  timestepFile << "        " << globalBoxMin[0] << " " << globalBoxMax[1] << " " << globalBoxMax[2] << "\n";
  timestepFile << "        " << globalBoxMax[0] << " " << globalBoxMax[1] << " " << globalBoxMax[2] << "\n";
  timestepFile << "      </DataArray>\n";
  timestepFile << "    </PPoints>\n";

  for (int i = 0; i < _numberOfRanks; ++i) {
    std::array<int, 6> pieceExtent = decomposition.getExtentOfSubdomain(i);
    timestepFile << "    <Piece "
                 << "Extent=\"" << pieceExtent[0] << " " << pieceExtent[1] << " " << pieceExtent[2] << " "
                 << pieceExtent[3] << " " << pieceExtent[4] << " " << pieceExtent[5] << "\" "
                 << "Source=\"./data/" << _sessionName << "_" << i << "_" << std::setfill('0')
                 << std::setw(_maximumNumberOfDigitsInIteration) << currentIteration << ".vts\"/>\n";
  }

  timestepFile << "  </PStructuredGrid>\n";
  timestepFile << "</VTKFile>\n";

  timestepFile.close();
}

void ParallelVtkWriter::tryCreateFolder(const std::string &name, const std::string &location) {
  try {
    // filesystem library unfortunately not available on all target systems e.g. Fugaku
    // std::filesystem::path newDirectoryPath(location + "/" + name);
    // std::filesystem::create_directory(newDirectoryPath);
    const auto newDirectoryPath{location + "/" + name};
    mkdir(newDirectoryPath.c_str(), 0777);
  } catch (const std::exception &ex) {
    throw std::runtime_error("The output location " + location +
                             " passed to ParallelVtkWriter is invalid: " + ex.what());
  }
}

void ParallelVtkWriter::generateFilename(const std::string &filetype, size_t currentIteration,
                                         std::ostringstream &filenameStream) {
  filenameStream << _dataFolderPath << _sessionName << "_" << _mpiRank << "_" << std::setfill('0')
                 << std::setw(_maximumNumberOfDigitsInIteration) << currentIteration << "." << filetype;
}
