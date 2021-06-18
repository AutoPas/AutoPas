/**
 * @file ParallelVtkWriter.h
 * @author J. Körner
 * @date 31.05.2021
 */
#pragma once

#include <string>

#include "autopas/AutoPas.h"
#include "src/TypeDefinitions.h"

/**
 * The ParallelVtkWriter can be used to create vtk-files for MPI parallel processes.
 */
class ParallelVtkWriter {
 public:
  /**
   * Constructor.
   * @param sessionName Sets the prefix for every created folder / file.
   * @param outputFolder Sets the folder where the vtk files will be created.
   */
  ParallelVtkWriter(std::string sessionName, const std::string &outputFolder);

  /**
   * Destructor.
   */
  ~ParallelVtkWriter() = default;

  /**
   * Writes the current state of particles into a vtk file.
   * @param currentIteration The simulation's current iteration.
   * @param maximumNumberOfDigitsInIteration The maximum number of digits an iteration index may have.
   * @param autoPasContainer The AutoPas container whose owned particles will be logged.
   */
  void recordTimestep(const int &currentIteration, const int &maximumNumberOfDigitsInIteration,
                      const autopas::AutoPas<ParticleType> &autoPasContainer);

 private:
  /**
   * Stores the MPI rank of the current process.
   * Every process will write into it's own vtk file, while the process with rank 0 will
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
   * Stores the name of output .vtu file for the current process.
   */
  std::string _outputFileName;

  /**
   * Tries to create a folder for the current writer session and stores it in _sessionFolderPath.
   */
  void tryCreateSessionFolder(const std::string &name, const std::string location);

  /**
   * Tries to create a folder at a location.
   * If the location does not exist this function will throw an error.
   * @param name The name of the new folder.
   * @param location The location where the new folder will be created.
   */
  void tryCreateFolder(const std::string &name, const std::string &location);
};
