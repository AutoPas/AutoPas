/**
 * @file ParallelVtkWriter.h
 * @author J. Körner
 * @date 31.05.2021
 */
#pragma once

#include <sys/stat.h>

#include <array>
#include <string>
#include <unordered_set>

#include "autopas/AutoPas.h"
#include "autopas/tuning/Configuration.h"
#include "src/TypeDefinitions.h"
#include "src/domainDecomposition/RegularGridDecomposition.h"

/**
 * The ParallelVtkWriter can be used to create vtk-files for MPI parallel processes.
 */
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
  void recordTimestep(size_t currentIteration, const autopas::AutoPas<ParticleType> &autoPasContainer,
                      const RegularGridDecomposition &decomposition) const;

 private:
  /**
   * Checks if a file with the given path exists.
   * @param filename
   * @return True iff the file exists.
   */
  static bool checkFileExists(const std::string &filename) {
    struct stat buffer {};
    return (stat(filename.c_str(), &buffer) == 0);
  }

  /**
   * Stores the number of ranks used in the simulation.
   * This information is required when creating the .pvtu file.
   */
  int _numberOfRanks{};

  /**
   * Stores the MPI rank of the current process.
   * Every process will write into it's own .vtu file, while the process with rank 0 will
   * create the parallel .pvtu file.
   */
  int _mpiRank{};

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
  void recordParticleStates(size_t currentIteration, const autopas::AutoPas<ParticleType> &autoPasContainer) const;

  /**
   * Writes the current domain subdivision into vtk files.
   * @param currentIteration: The simulations current iteration.
   * @param autoPasConfigurations: All current configuration of an autopas container (pairwise, triwise).
   * @param decomposition: The simulations domain decomposition.
   */
  void recordDomainSubdivision(
      size_t currentIteration,
      const std::unordered_map<autopas::InteractionTypeOption::Value,
                               std::reference_wrapper<const autopas::Configuration>> &autoPasConfigurations,
      const RegularGridDecomposition &decomposition) const;

  /**
   * Tries to create a folder for the current writer session and stores it in _sessionFolderPath.
   */
  void tryCreateSessionAndDataFolders(const std::string &name, const std::string &location);

  /**
   * Creates the .pvtu file for particle data that references all particle data vtu files from this timestep.
   *
   * @note For visualization in ParaView the .pvtu files need to be loaded.
   *
   * @param currentIteration: The simulation's current iteration.
   */
  void createParticlesPvtuFile(size_t currentIteration) const;

  /**
   * Creates the .pvtu file for rank data that references all rank data vtu files from this timestep
   *
   * @note For visualization in ParaView the .pvtu files need to be loaded.
   *
   * @param currentIteration: The simulation's current iteration.
   * @param decomposition: The decomposition of the domain.
   * @param interactionTypes: Interaction types that are considered in the current simulation.
   */
  void createRanksPvtuFile(size_t currentIteration, const RegularGridDecomposition &decomposition,
                           const std::unordered_set<autopas::InteractionTypeOption::Value> &interactionTypes) const;

  /**
   * Tries to create a folder at a location.
   * If the location does not exist this function will throw an error.
   * @param name The name of the new folder.
   * @param location The location where the new folder will be created.
   */
  static void tryCreateFolder(const std::string &name, const std::string &location);

  /**
   * Generates the file name for a given vtk file type.
   * @param tag String tag that is inserted in the file name
   * @param currentIteration: The current iteration to record.
   * @param fileExtension: The vtk file type extension. Pass the extension without the '.'.
   * @param filenameStream: The output string for the filename.
   */
  void generateFilename(const std::string &tag, const std::string &fileExtension, size_t currentIteration,
                        std::ostringstream &filenameStream) const;
};
