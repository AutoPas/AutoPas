/**
 * @file ParallelVtkWriter.h
 * @author J. Körner
 * @date 31.05.2021
 */
#pragma once

#include <array>
#include <string>

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
                    const int &maximumNumberOfDigitsInIteration);

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
                      const RegularGridDecomposition<ParticleClass> &decomposition);

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
   */
  void recordParticleStates(const int &currentIteration, const autopas::AutoPas<ParticleClass> &autoPasContainer);

  /**
   * Writes the current domain subdivision into vtk files.
   * @param currentIteration: The simulations current iteration.
   * @param autoPasConfiguration: The configuration of an autoPasContainer.
   * @param decomposition: The simulations domain decomposition.
   */
  void recordDomainSubdivision(const int &currentIteration, const autopas::Configuration &autoPasConfiguration,
                               const RegularGridDecomposition<ParticleClass> &decomposition);

  /**
   * Calculates the whole extent of the decompositions local domain.
   * The whole extent defines the space this local domain is occupying in the global domain.
   * The layout of the returned array is [ xmin, xmax, ymin, ymax, zmin, zmax ], where x, y and z are coordinates in
   * in the decomposition grid.
   * @param domainDecomposition: The simulations domain decomposition.
   * @return the whole extent of the local domain.
   */
  std::array<int, 6> calculateWholeExtent(const RegularGridDecomposition<ParticleClass> &domainDecomposition);

  /**
   * Tries to create a folder for the current writer session and stores it in _sessionFolderPath.
   */
  void tryCreateSessionAndDataFolders(const std::string &name, const std::string location);

  /**
   * Creates the .pvtu file required to load unstructured grid data from multiple ranks into ParaView.
   * @param currentIteration: The simulation's current iteration.
   */
  void createPvtuFile(const int &currentIteration);

  /**
   * Creates the .pvts file required to load structured grid data from multiple ranks into ParaView.
   * @param currentIteration: The simulation's current iteration.
   * @param decomposition: The decomposition of the domain.
   */
  void createPvtsFile(const int &currentIteration, const RegularGridDecomposition<ParticleClass> &decomposition);

  /**
   * Tries to create a folder at a location.
   * If the location does not exist this function will throw an error.
   * @param name The name of the new folder.
   * @param location The location where the new folder will be created.
   */
  void tryCreateFolder(const std::string &name, const std::string &location);

  /**
   * Generates the file name for a given vtk file type.
   * @param currentIteration: The current iteration to record.
   * @param filetype: The vtk file type extension. Pass the extension without the '.'.
   * @param filenameStream: The output string string for the filename.
   */
  void generateFilename(const std::string &filetype, const int &currentIteration, std::ostringstream &filenameStream);
};
