/**
 * @file ParallelVtkWriter.cpp
 * @author J. Körner
 * @date 31.05.2021
 */
#include "ParallelVtkWriter.h"

#include <algorithm>
#include <cstddef>
#include <ios>
#include <limits>
#include <string>
#include <utility>

#include "autopas/utils/WrapMPI.h"

ParallelVtkWriter::ParallelVtkWriter(std::string sessionName, const std::string &outputFolder,
                                     const int &maximumNumberOfDigitsInIteration)
    : _sessionName(std::move(sessionName)), _maximumNumberOfDigitsInIteration(maximumNumberOfDigitsInIteration) {
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &_numberOfRanks);
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &_mpiRank);

  if (_mpiRank == 0) {
    tryCreateSessionAndDataFolders(_sessionName, outputFolder);
  }

  int sessionFolderPathLength = static_cast<int>(_sessionFolderPath.size());
  autopas::AutoPas_MPI_Bcast(&sessionFolderPathLength, 1, AUTOPAS_MPI_INT, 0, AUTOPAS_MPI_COMM_WORLD);

  int dataFolderPathLength = static_cast<int>(_dataFolderPath.size());
  autopas::AutoPas_MPI_Bcast(&dataFolderPathLength, 1, AUTOPAS_MPI_INT, 0, AUTOPAS_MPI_COMM_WORLD);

  if (_mpiRank != 0) {
    _sessionFolderPath.resize(sessionFolderPathLength);
    _dataFolderPath.resize(dataFolderPathLength);
  }

  autopas::AutoPas_MPI_Bcast(_sessionFolderPath.data(), sessionFolderPathLength, AUTOPAS_MPI_CHAR, 0,
                             AUTOPAS_MPI_COMM_WORLD);
  autopas::AutoPas_MPI_Bcast(_dataFolderPath.data(), dataFolderPathLength, AUTOPAS_MPI_CHAR, 0, AUTOPAS_MPI_COMM_WORLD);
}

void ParallelVtkWriter::recordTimestep(size_t currentIteration, const autopas::AutoPas<ParticleType> &autoPasContainer,
                                       const RegularGridDecomposition &decomposition) const {
  recordParticleStates(currentIteration, autoPasContainer);
  const auto currentConfig = autoPasContainer.getCurrentConfigs();
  recordDomainSubdivision(currentIteration, currentConfig, decomposition);
}

/**
 * @todo: Currently this function runs over all the particles for each property separately.
 * This can be improved by using multiple string streams (one for each property).
 * The streams can be combined to a single output stream after iterating over the particles, once.
 */
void ParallelVtkWriter::recordParticleStates(size_t currentIteration,
                                             const autopas::AutoPas<ParticleType> &autoPasContainer) const {
  if (_mpiRank == 0) {
    createParticlesPvtuFile(currentIteration);
  }

  std::ostringstream timestepFileName;
  generateFilename("Particles", "vtu", currentIteration, timestepFileName);

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
  autoPasContainer.forEach(
      [&](ParticleType const& particle) {
        const auto v = particle.getV();
        timestepFile << "        " << v[0] << " " << v[1] << " " << v[2] << "\n";
      },
      autopas::IteratorBehavior::owned);
  timestepFile << "        </DataArray>\n";

  // print forces
  timestepFile << "        <DataArray Name=\"forces\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">\n";
  autoPasContainer.forEach(
      [&](ParticleType const& particle) {
        const auto f = particle.getF();
        timestepFile << "        " << f[0] << " " << f[1] << " " << f[2] << "\n";
      },
      autopas::IteratorBehavior::owned);
  timestepFile << "        </DataArray>\n";

#if MD_FLEXIBLE_MODE == MULTISITE
  // print quaternions
  timestepFile
      << "        <DataArray Name=\"quaternions\" NumberOfComponents=\"4\" format=\"ascii\" type=\"Float32\">\n";
  autoPasContainer.forEach(
      [&](ParticleType const& particle) {
        const auto q = particle.getQuaternion();
        timestepFile << "        " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << "\n";
      },
      autopas::IteratorBehavior::owned);
  timestepFile << "        </DataArray>\n";

  // print angular velocities
  timestepFile
      << "        <DataArray Name=\"angularVelocities\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">\n";
  autoPasContainer.forEach(
      [&](ParticleType const& particle) {
        const auto angVel = particle.getAngularVel();
        timestepFile << "        " << angVel[0] << " " << angVel[1] << " " << angVel[2] << "\n";
      },
      autopas::IteratorBehavior::owned);
  timestepFile << "        </DataArray>\n";

  // print torques
  timestepFile << "        <DataArray Name=\"torques\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">\n";
  autoPasContainer.forEach(
      [&](ParticleType const& particle) {
        const auto torque = particle.getTorque();
        timestepFile << "        " << torque[0] << " " << torque[1] << " " << torque[2] << "\n";
      },
      autopas::IteratorBehavior::owned);
  timestepFile << "        </DataArray>\n";
#endif

  // print type ids
  timestepFile << "        <DataArray Name=\"typeIds\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">\n";
  autoPasContainer.forEach(
      [&](ParticleType const& particle) { timestepFile << "        " << particle.getTypeId() << "\n"; },
      autopas::IteratorBehavior::owned);
  timestepFile << "        </DataArray>\n";

  // print type ids
  timestepFile << "        <DataArray Name=\"blockIds\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">\n";
  autoPasContainer.forEach(
      [&](ParticleType const& particle) { timestepFile << "        " << particle.getBlockId() << "\n"; },
      autopas::IteratorBehavior::owned);
  timestepFile << "        </DataArray>\n";


  // print ids
  timestepFile << "        <DataArray Name=\"ids\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">\n";
  autoPasContainer.forEach(
      [&](ParticleType const& particle) { timestepFile << "        " << particle.getID() << "\n"; },
      autopas::IteratorBehavior::owned);
  timestepFile << "        </DataArray>\n";

  timestepFile << "      </PointData>\n";
  timestepFile << "      <CellData/>\n";
  timestepFile << "      <Points>\n";

  // print positions
  timestepFile << "        <DataArray Name=\"positions\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">\n";
  const auto boxMax = autoPasContainer.getBoxMax();
  autoPasContainer.forEach(
      [&](ParticleType const& particle) {
        // When we write to the file in ASCII, values are rounded to the precision of the filestream.
        // Since a higher precision results in larger files because more characters are written,
        // and mdflex is not intended as a perfectly precice tool for application scientists,
        // we are fine with the rather low default precision.
        // However, if a particle is very close to the domain border it can happen that the particle position is rounded
        // exactly to the boundary position. This then causes problems when the checkpoint is loaded because boxMax is
        // considered to be not part of the domain, hence such a particle would not be loaded and thus be lost.
        // This function identifies these problematic values and raises the write precision just for this value high
        // enough to be distinguishable from the boundary.
        const auto writeWithDynamicPrecision = [&](double position, double border) {
          const auto initialPrecision = timestepFile.precision();
          // Simple and cheap check if we even need to do anything.
          if (border - position < 0.1) {
            using autopas::utils::Math::roundFloating;
            using autopas::utils::Math::isNearAbs;
            // As long as the used precision results in the two values being indistinguishable increase the precision
            while (isNearAbs(roundFloating(position, timestepFile.precision()), border,
                             std::pow(10, -timestepFile.precision()))) {
              timestepFile << std::setprecision(timestepFile.precision() + 1);
              // Abort if the numbers are indistinguishable beyond machine precision
              constexpr auto machinePrecision = std::numeric_limits<double>::digits10;
              if (timestepFile.precision() > machinePrecision) {
                throw std::runtime_error(
                    "ParallelVtkWriter::writeWithDynamicPrecision(): "
                    "The two given numbers are identical up to " +
                    std::to_string(machinePrecision) +
                    " digits of precision!\n"
                    "Number: " +
                    std::to_string(position) + "\n" + particle.toString());
              }
            }
          }
          // Write with the new precision and then reset it
          timestepFile << position << std::setprecision(initialPrecision);
        };

        const auto pos = particle.getR();
        timestepFile << "        ";
        writeWithDynamicPrecision(pos[0], boxMax[0]);
        timestepFile << " ";
        writeWithDynamicPrecision(pos[1], boxMax[1]);
        timestepFile << " ";
        writeWithDynamicPrecision(pos[2], boxMax[2]);
        timestepFile << "\n";
      },
      autopas::IteratorBehavior::owned);
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

void ParallelVtkWriter::recordDomainSubdivision(
    size_t currentIteration,
    const std::unordered_map<autopas::InteractionTypeOption::Value,
                             std::reference_wrapper<const autopas::Configuration>> &autoPasConfigurations,
    const RegularGridDecomposition &decomposition) const {
  // Extract active interaction types to print them to the .pvtu file.
  std::unordered_set<autopas::InteractionTypeOption::Value> interactionTypes;
  interactionTypes.reserve(autoPasConfigurations.size());
  std::transform(autoPasConfigurations.begin(), autoPasConfigurations.end(),
                 std::inserter(interactionTypes, interactionTypes.end()), [&](const auto &pair) { return pair.first; });
  if (_mpiRank == 0) {
    createRanksPvtuFile(currentIteration, decomposition, interactionTypes);
  }

  std::ostringstream timestepFileName;
  generateFilename("Ranks", "vtu", currentIteration, timestepFileName);

  std::ofstream timestepFile;
  timestepFile.open(timestepFileName.str(), std::ios::out | std::ios::binary);

  if (not timestepFile.is_open()) {
    throw std::runtime_error("Simulation::writeVTKFile(): Failed to open file \"" + timestepFileName.str() + "\"");
  }

  const std::array<double, 3> localBoxMin = decomposition.getLocalBoxMin();
  const std::array<double, 3> localBoxMax = decomposition.getLocalBoxMax();

  auto printDataArray = [&](const auto &data, const std::string &type, const std::string &name) {
    timestepFile << "        <DataArray type=\"" << type << "\" Name=\"" << name << "\" format=\"ascii\">\n";
    timestepFile << "          " << data << "\n";
    timestepFile << "        </DataArray>\n";
  };

  timestepFile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n";
  timestepFile << "<VTKFile byte_order=\"LittleEndian\" type=\"UnstructuredGrid\" version=\"0.1\">\n";
  timestepFile << "  <UnstructuredGrid>\n";
  timestepFile << "    <Piece NumberOfPoints=\"8\" NumberOfCells=\"1\">\n";
  timestepFile << "      <CellData>\n";
  printDataArray(decomposition.getDomainIndex(), "Int32", "DomainId");

  // General Configuration information
  printDataArray(autoPasConfigurations.begin()->second.get().cellSizeFactor, "Float32", "CellSizeFactor");
  printDataArray(static_cast<int>(autoPasConfigurations.begin()->second.get().container), "Int32", "Container");

  // Pairwise Configuration
  if (autoPasConfigurations.find(autopas::InteractionTypeOption::pairwise) != autoPasConfigurations.end()) {
    auto pairwiseConfig = autoPasConfigurations.at(autopas::InteractionTypeOption::pairwise).get();
    printDataArray(static_cast<int>(pairwiseConfig.dataLayout), "Int32", "DataLayout");
    printDataArray(static_cast<int>(pairwiseConfig.loadEstimator), "Int32", "LoadEstimator");
    printDataArray(static_cast<int>(pairwiseConfig.traversal), "Int32", "Traversal");
    printDataArray(static_cast<int>(pairwiseConfig.newton3), "Int32", "Newton3");
  }

  // Triwise Configuration
  if (autoPasConfigurations.find(autopas::InteractionTypeOption::triwise) != autoPasConfigurations.end()) {
    auto triwiseConfig = autoPasConfigurations.at(autopas::InteractionTypeOption::triwise).get();
    printDataArray(static_cast<int>(triwiseConfig.dataLayout), "Int32", "DataLayout-3B");
    printDataArray(static_cast<int>(triwiseConfig.traversal), "Int32", "Traversal-3B");
    printDataArray(static_cast<int>(triwiseConfig.newton3), "Int32", "Newton3-3B");
  }

  printDataArray(_mpiRank, "Int32", "Rank");
  timestepFile << "      </CellData>\n";
  timestepFile << "      <Points>\n";
  timestepFile << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  timestepFile << "          " << localBoxMin[0] << " " << localBoxMin[1] << " " << localBoxMin[2] << "\n";
  timestepFile << "          " << localBoxMin[0] << " " << localBoxMin[1] << " " << localBoxMax[2] << "\n";
  timestepFile << "          " << localBoxMin[0] << " " << localBoxMax[1] << " " << localBoxMin[2] << "\n";
  timestepFile << "          " << localBoxMin[0] << " " << localBoxMax[1] << " " << localBoxMax[2] << "\n";
  timestepFile << "          " << localBoxMax[0] << " " << localBoxMin[1] << " " << localBoxMin[2] << "\n";
  timestepFile << "          " << localBoxMax[0] << " " << localBoxMin[1] << " " << localBoxMax[2] << "\n";
  timestepFile << "          " << localBoxMax[0] << " " << localBoxMax[1] << " " << localBoxMin[2] << "\n";
  timestepFile << "          " << localBoxMax[0] << " " << localBoxMax[1] << " " << localBoxMax[2] << "\n";
  timestepFile << "        </DataArray>\n";
  timestepFile << "      </Points>\n";
  timestepFile << "      <Cells>\n";
  timestepFile << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  timestepFile << "          0 1 2 3 4 5 6 7\n";  // These indices refer to the Points DataArray above.
  timestepFile << "        </DataArray>\n";
  timestepFile << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  timestepFile << "          8\n";  // The cell is defined by 8 points
  timestepFile << "        </DataArray>\n";
  timestepFile << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  timestepFile << "          11\n";  // = VTK_VOXEL
  timestepFile << "        </DataArray>\n";
  timestepFile << "      </Cells>\n";
  timestepFile << "    </Piece>\n";
  timestepFile << "  </UnstructuredGrid>\n";
  timestepFile << "</VTKFile>\n";

  timestepFile.close();
}

void ParallelVtkWriter::tryCreateSessionAndDataFolders(const std::string &name, const std::string &location) {
  if (not checkFileExists(location)) {
    tryCreateFolder(location, "./");
  }

  _sessionFolderPath = location + "/" + name + "/";
  tryCreateFolder(name, location);

  _dataFolderPath = _sessionFolderPath + "data/";
  tryCreateFolder("data", _sessionFolderPath);
}

void ParallelVtkWriter::createParticlesPvtuFile(size_t currentIteration) const {
  std::ostringstream filename;
  filename << _sessionFolderPath << _sessionName << "_Particles_" << std::setfill('0')
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
#if MD_FLEXIBLE_MODE == MULTISITE
  timestepFile
      << "      <PDataArray Name=\"quaternions\" NumberOfComponents=\"4\" format=\"ascii\" type=\"Float32\"/>\n";
  timestepFile
      << "      <PDataArray Name=\"angularVelocities\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\"/>\n";
  timestepFile << "      <PDataArray Name=\"torques\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\"/>\n";
#endif
  timestepFile << "      <PDataArray Name=\"typeIds\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\"/>\n";
  timestepFile << "      <PDataArray Name=\"blockIds\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\"/>\n";
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
    timestepFile << "    <Piece Source=\"./data/" << _sessionName << "_Particles_" << i << "_" << std::setfill('0')
                 << std::setw(_maximumNumberOfDigitsInIteration) << currentIteration << ".vtu\"/>\n";
  }

  timestepFile << "  </PUnstructuredGrid>\n";
  timestepFile << "</VTKFile>\n";

  timestepFile.close();
}

void ParallelVtkWriter::createRanksPvtuFile(
    size_t currentIteration, const RegularGridDecomposition &decomposition,
    const std::unordered_set<autopas::InteractionTypeOption::Value> &interactionTypes) const {
  std::ostringstream filename;
  filename << _sessionFolderPath << _sessionName << "_Ranks_" << std::setfill('0')
           << std::setw(_maximumNumberOfDigitsInIteration) << currentIteration << ".pvtu";

  std::ofstream timestepFile;
  timestepFile.open(filename.str(), std::ios::out | std::ios::binary);

  if (not timestepFile.is_open()) {
    throw std::runtime_error("Simulation::writeVTKFile(): Failed to open file \"" + filename.str() + "\"");
  }
  const auto &globalBoxMin = decomposition.getGlobalBoxMin();
  const auto &globalBoxMax = decomposition.getGlobalBoxMax();
  timestepFile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n";
  timestepFile << "<VTKFile byte_order=\"LittleEndian\" type=\"PUnstructuredGrid\" version=\"0.1\">\n";
  timestepFile << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
  timestepFile << "    <PPointData/>\n";
  timestepFile << "    <PCellData>\n";
  timestepFile << "      <PDataArray type=\"Int32\" Name=\"DomainId\" />\n";

  // General configuration options
  timestepFile << "      <PDataArray type=\"Float32\" Name=\"CellSizeFactor\" />\n";
  timestepFile << "      <PDataArray type=\"Int32\" Name=\"Container\" />\n";

  // Pairwise configuration
  if (interactionTypes.find(autopas::InteractionTypeOption::pairwise) != interactionTypes.end()) {
    timestepFile << "      <PDataArray type=\"Int32\" Name=\"DataLayout\" />\n";
    timestepFile << "      <PDataArray type=\"Int32\" Name=\"LoadEstimator\" />\n";
    timestepFile << "      <PDataArray type=\"Int32\" Name=\"Traversal\" />\n";
    timestepFile << "      <PDataArray type=\"Int32\" Name=\"Newton3\" />\n";
  }

  // Triwise configuration
  if (interactionTypes.find(autopas::InteractionTypeOption::triwise) != interactionTypes.end()) {
    timestepFile << "      <PDataArray type=\"Int32\" Name=\"DataLayout-3B\" />\n";
    timestepFile << "      <PDataArray type=\"Int32\" Name=\"Traversal-3B\" />\n";
    timestepFile << "      <PDataArray type=\"Int32\" Name=\"Newton3-3B\" />\n";
  }

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
    timestepFile << "    <Piece "
                 << "Source=\"./data/" << _sessionName << "_Ranks_" << i << "_" << std::setfill('0')
                 << std::setw(_maximumNumberOfDigitsInIteration) << currentIteration << ".vtu\"/>\n";
  }

  timestepFile << "  </PUnstructuredGrid>\n";
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
    throw std::runtime_error("ParallelVtkWriter::tryCreateFolder(): The output location " + location +
                             " passed to ParallelVtkWriter is invalid: " + ex.what());
  }
}

void ParallelVtkWriter::generateFilename(const std::string &tag, const std::string &fileExtension,
                                         size_t currentIteration, std::ostringstream &filenameStream) const {
  filenameStream << _dataFolderPath << _sessionName << "_" << tag << "_" << _mpiRank << "_" << std::setfill('0')
                 << std::setw(_maximumNumberOfDigitsInIteration) << currentIteration << "." << fileExtension;
}
