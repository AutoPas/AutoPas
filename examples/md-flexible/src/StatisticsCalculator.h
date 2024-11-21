/**
 * @file StatisticsCalculator.h
 * @author Joon Kim
 * @date 21.11.2024
 */
#pragma once

#include <sys/stat.h>

#include <tuple>
#include <array>
#include <string>

#include "autopas/Autopas.h"
#include "src/TypeDefinitions.h"

/**
 * The StatisticsCalculator can be used to calculate statistics and to save them into a file.
 */
class StatisticsCalculator {
 public:
  /**
   * Constructor.
   * @param sessionName Sets the prefix for every created folder / file.
   * @param outputFolder Sets the folder where the statistics files will be created.
   */
  StatisticsCalculator(std::string sessionName, const std::string &outputFolder);

  /**
   * Destructor.
   */
  ~StatisticsCalculator() = default;

  /**
   * Writes the statistics of the current state of the simulation into a file.
   */
  void recordStatistics(size_t currentIteration, const autopas::AutoPas<ParticleType> &autopasContainer,
                        ParticlePropertiesLibraryType &particlePropertiesLib);

 private:
  /**
   * Calculates the statistics of the current state of the simulation.
   * @param autopasContainer
   * @param particlePropertiesLib
   * @return tuple of doubles containing the statistics.
   */
  static std::tuple<double, double, double> calculateStatistics(const autopas::AutoPas<ParticleType> &autopasContainer,
                                         ParticlePropertiesLibraryType &particlePropertiesLib) ;

  /**
   * Generates the output file (.csv) for the statistics.
   */
  void generateOutputFile();

  /**
   * Tries to create a folder for the current writer session and stores it in _sessionFolderPath.
   */
  void tryCreateStatisticsFolders(const std::string &name, const std::string &location);

  /**
   * Tries to create a folder at a location.
   * If the location does not exist this function will throw an error.
   * @param name The name of the new folder.
   * @param location The location where the new folder will be created.
   */
  void tryCreateFolder(const std::string &name, const std::string &location);

  /**
   * Checks if a file with the given path exists.
   * @param filename
   * @return True iff the file exists.
   */
  static bool checkFileExists(const std::string &filename) {
    struct stat buffer{};
    return (stat(filename.c_str(), &buffer) == 0);
  }

  /**
   * Writes a value of statistics into the output file.
   * @tparam T
   * @param os
   * @param value
   */
  template <typename T>
  void writeValue(std::ostream &os, const T &value) {
    os << value;
  }

  /**
   * Writes multiple values of statistics recursively into the output file.
   * @tparam T
   * @tparam Args
   * @param os
   * @param value
   * @param args
   */
  template <typename T, typename... Args>
  void writeValue(std::ostream &os, const T &value, const Args &...args) {
    os << value << ",";
    writeValue(os, args...);
  }

  /**
   * Writes a tuple of statistics into the output file.
   * @tparam Args
   * @param os
   * @param tuple
   */
  template <typename... Args>
  void writeValue(std::ostream &os, const std::tuple<Args...> &tuple) {
    writeTupleElements(os, tuple, std::index_sequence_for<Args...>{});
  }

  /**
   * Writes the elements of a tuple into the output file.
   * @tparam Tuple
   * @tparam Indices
   * @param os
   * @param tuple
   */
  template <typename Tuple, std::size_t... Indices>
  void writeTupleElements(std::ostream &os, const Tuple &tuple, std::index_sequence<Indices...>) {
    ((os << (Indices == 0 ? "" : ",") << std::get<Indices>(tuple)), ...);
  }

  template <typename... Args>
  void writeRow(const Args &...args) {
    if (outputFile.is_open()) {
      writeValue(outputFile, args...);
      outputFile << "\n";
    } else {
      throw std::runtime_error("StatisticsCalculator::writeRow(): Could not open file " + _outputFileName);
    }
  }


  /**
   * The file stream to write the statistics to.
   */
  std::ofstream outputFile;

  /**
   * Stores the session name.
   */
  std::string _sessionName;

  /**
   * Stores the path to the current session's output folder.
   */
  std::string _sessionFolderPath;

  /**
   * Stores the path to the folder where the current session's actual statistics are stored.
   */
  std::string _statisticsFolderPath;

  /**
   * Stores the name of output .vtu file for the current process.
   */
  std::string _outputFileName;
};