/**
 * @file ParallelVtkWriter.h
 * @author J. Körner
 * @date 31.05.2021
 */
#pragma once

#include <array>
#include <string>
#include <filesystem>

#include "autopas/AutoPas.h"
#include "autopas/selectors/Configuration.h"
#include "src/TypeDefinitions.h"
#include "src/domainDecomposition/RegularGridDecomposition.h"

/**
 * The ParallelVtkWriter can be used to create vtk-files for MPI parallel processes.
 */
template <class ParticleClass>
class ParallelVtkWriter {
 public:
  /**
   * Constructor.
   * @param sessionName Sets the prefix for every created folder / file.
   * @param outputFolder Sets the folder where the vtk files will be created.
   * @param maximumNumberOfDigitsInIteration The maximum number of digits an iteration index can have.
   */
  ParallelVtkWriter(std::string sessionName, const std::string &outputFolder,
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

  /**
   * Destructor.
   */
  ~ParallelVtkWriter() = default;

  /**
   * Writes the current state of particles and the current domain subdivision into vtk files.
   * @param currentIteration The simulation's current iteration.
   * @param autoPasContainer The AutoPas container whose owned particles will be logged.
   * @param decomposition: The decomposition of the global domain.
   */
  void recordTimestep(const int &currentIteration, const autopas::AutoPas<ParticleClass> &autoPasContainer,
                      const RegularGridDecomposition<ParticleClass> &decomposition) {
    recordParticleStates(currentIteration, autoPasContainer);
    recordDomainSubdivision(currentIteration, autoPasContainer.getCurrentConfig(), decomposition);
  }

 private:
  /**
   * Stores the number of ranks used in the simulation.
   * This information is required when creating the .pvtu file.
   */
  int _numberOfRanks;

  /**
   * Stores the MPI rank of the current process.
   * Every process will write into it's own .vtu file, while the process with rank 0 will
   * create the parallel .pvtu file.
   */
  int _mpiRank;

  /**
   * Stores the session name.
   */
  std::string _sessionName;

  /**
   * Stores the path to the current session's output folder.
   */
  std::string _sessionFolderPath;

  /**
   * Stores the path to the folder where the current session's actual data is stored.
   */
  std::string _dataFolderPath;

  /**
   * Stores the name of output .vtu file for the current process.
   */
  std::string _outputFileName;

  /**
   * Stores the maximum number of digits an iteration can have.
   * This is used to determine the number of leading zeros for each timestep record.
   */
  int _maximumNumberOfDigitsInIteration;

  /**
   * Writes the current state of particles into vtk files.
   * @param currentIteration: The simulations current iteration.
   * @param autoPasContainer The AutoPas container whose owned particles will be logged.
   *
   * @todo: Currently this function runs over all the particles for each property separately.
   * This can be improved by using multiple string streams (one for each property).
   * The streams can be combined to a single output stream after iterating over the particles, once.
   */
  void recordParticleStates(const int &currentIteration, const autopas::AutoPas<ParticleClass> &autoPasContainer) {
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
    timestepFile << "        <DataArray Name=\"typeIds\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\">\n";
    for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
      timestepFile << "        " << particle->getTypeId() << "\n";
    }
    timestepFile << "        </DataArray>\n";

    // print ids
    timestepFile << "        <DataArray Name=\"ids\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\">\n";
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

  /**
   * Writes the current domain subdivision into vtk files.
   * @param currentIteration: The simulations current iteration.
   * @param autoPasConfiguration: The configuration of an autoPasContainer.
   * @param decomposition: The simulations domain decomposition.
   */
  void recordDomainSubdivision(const int &currentIteration, const autopas::Configuration &autoPasConfiguration,
                               const RegularGridDecomposition<ParticleClass> &decomposition) {
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

  /**
   * Calculates the whole extent of the decompositions local domain.
   * The whole extent defines the space this local domain is occupying in the global domain.
   * The layout of the returned array is [ xmin, xmax, ymin, ymax, zmin, zmax ], where x, y and z are coordinates in
   * in the decomposition grid.
   * @param domainDecomposition: The simulations domain decomposition.
   * @return the whole extent of the local domain.
   */
  std::array<int, 6> calculateWholeExtent(const RegularGridDecomposition<ParticleClass> &domainDecomposition) {
    std::array<int, 6> wholeExtent;
    std::array<int, 3> domainId = domainDecomposition.getDomainId();
    std::array<int, 3> decomposition = domainDecomposition.getDecomposition();
    for (int i = 0; i < 3; ++i) {
      wholeExtent[2 * i] = domainId[i];
      wholeExtent[2 * i + 1] = std::min(domainId[i] + 1, decomposition[i]);
    }
    return wholeExtent;
  }

  /**
   * Tries to create a folder for the current writer session and stores it in _sessionFolderPath.
   */
  void tryCreateSessionAndDataFolders(const std::string &name, const std::string location) {
    time_t rawTime;
    time(&rawTime);

    struct tm timeInformation;
    gmtime_r(&rawTime, &timeInformation);

    char buffer[80];
    strftime(buffer, sizeof(buffer), "%d%m%Y_%H%M%S", &timeInformation);
    std::string timeString(buffer);

    if (not std::filesystem::exists(location)) {
      tryCreateFolder(location, "./");
    }

    _sessionFolderPath = location + "/" + name + "_" + timeString + "/";
    tryCreateFolder(name + "_" + timeString, location);

    _dataFolderPath = _sessionFolderPath + "data/";
    tryCreateFolder("data", _sessionFolderPath);
  }

  /**
   * Creates the .pvtu file required to load unstructured grid data from multiple ranks into ParaView.
   * @param currentIteration: The simulation's current iteration.
   */
  void createPvtuFile(const int &currentIteration) {
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
    timestepFile << "      <PDataArray Name=\"forces\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Int32\"/>\n";
    timestepFile << "      <PDataArray Name=\"typeIds\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\"/>\n";
    timestepFile << "      <PDataArray Name=\"ids\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\"/>\n";
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

  /**
   * Creates the .pvts file required to load structured grid data from multiple ranks into ParaView.
   * @param currentIteration: The simulation's current iteration.
   * @param decomposition: The decomposition of the domain.
   */
  void createPvtsFile(const int &currentIteration, const RegularGridDecomposition<ParticleClass> &decomposition) {
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

  /**
   * Tries to create a folder at a location.
   * If the location does not exist this function will throw an error.
   * @param name The name of the new folder.
   * @param location The location where the new folder will be created.
   */
  void tryCreateFolder(const std::string &name, const std::string &location) {
    try {
      std::filesystem::path newDirectoryPath(location + "/" + name);
      std::filesystem::create_directory(newDirectoryPath);
    } catch (std::filesystem::filesystem_error const &ex) {
      throw std::runtime_error("The output location " + location + " passed to ParallelVtkWriter is invalid");
    }
  }

  /**
   * Generates the file name for a given vtk file type.
   * @param currentIteration: The current iteration to record.
   * @param filetype: The vtk file type extension. Pass the extension without the '.'.
   * @param filenameStream: The output string string for the filename.
   */
  void generateFilename(const std::string &filetype, const int &currentIteration, std::ostringstream &filenameStream) {
    filenameStream << _dataFolderPath << _sessionName << "_" << _mpiRank << "_" << std::setfill('0')
                   << std::setw(_maximumNumberOfDigitsInIteration) << currentIteration << "." << filetype;
  }
};
