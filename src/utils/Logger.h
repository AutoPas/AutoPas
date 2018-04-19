/**
 * LoggerInterface class
 *
 * @file Logger.h
 * @author seckler
 * @date 16.04.18
 */

#pragma once

#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>

namespace autopas {
/* we use a separate namespace because we have some global definitions for
 * the log level */
namespace log {
class Logger;
}
/**
 * the global logger instance of autopas
 */
extern std::unique_ptr<log::Logger> logger;
namespace log {

/**
 * list of available log levels
 * For each level a name has to be specified in the constructor of the Logger()
 * class.
 * This name will be prepended later on to the log message. */
typedef enum {
  None = 0,    /* supress output */
  Fatal = 1,   /* program exit */
  Error = 2,   /* program corrected */
  Warning = 4, /* perhaps wrong */
  Info = 8,    /* user info */
  Debug = 16,  /* detailed info for debugging */
  All
} logLevel;

/** @brief The Logger class provides a simple interface to handle log messages.
 *
 * Provides easy interface to handle log messages. Initialize either with
 * output level and stream or output level and filename or use default
 * constructor values (Error, &(std::cout)).
 * For writing log messages use fatal(), error(), warning(), info() or debug()
 * as with normal streams, e.g.
 * log.error() << "Wrong parameter." << std::endl;
 * Please include std::endl statements at the end of output as they will flush
 * the stream buffers.
 */
class Logger {
 public:
  /**
   * Constructor of the Logger class.
   * It provides different levels of output, which are defined using a logLevel.
   * @param level the level of the output
   * @param os the normal output stream, default is std::cout
   * @param es error stream, default is std::cerr
   * @param identifier you can mark this stream with a special identifier, use
   * this e.g. to mark mpi ranks
   */
  Logger(logLevel level = logLevel::Error, std::ostream *os = &(std::cout),
         std::ostream *es = &(std::cerr), std::string identifier = "0")
      : _log_level(level),
        _msg_log_level(logLevel::Error),
        _enabled(true),
        _filename(""),
        _log_stream(os),
        _error_stream(es),
        logLevelNames(),
        _starttime(),
        _identifier(identifier) {
    init_starting_time();
    this->init_log_levels();
  }

  /**
   * Constructor to print values to a file.
   * The filename is setup as "${prefix}${identifier}.log"
   * @overload
   * @param level the output level of the logger. See logLevel
   * @param prefix the prefix for the file name
   * @param identifier you can mark this stream with a special identifier, use
   * this e.g. to mark mpi ranks
   */
  Logger(logLevel level, std::string prefix, std::string identifier = "0")
      : _log_level(level),
        _msg_log_level(logLevel::Error),
        _enabled(true),
        _filename(""),
        _log_stream(nullptr),
        logLevelNames(),
        _starttime(),
        _identifier(identifier) {
    init_starting_time();
    std::stringstream filenamestream;
    filenamestream << prefix;
    filenamestream << _identifier;
    filenamestream << ".log";
    _filename = filenamestream.str();
    _log_stream = new std::ofstream(_filename.c_str());
    _error_stream = _log_stream;
  }
  /**
   * Destructor, flushes stream
   */
  ~Logger() {
    *_log_stream << std::flush;
    *_error_stream << std::flush;
    if (_filename != "") {
      (static_cast<std::ofstream *>(_log_stream))->close();
      // don't close the error_stream as it is the same as the log stream!
    }
  }

  // don't allow copy-construction
  Logger(const Logger &) = delete;

  // don't allow assignment
  Logger &operator=(const Logger &) = delete;

  /**
   * enable or disable the logger.
   * @param enabled
   */
  void setEnabled(bool enabled) { _enabled = enabled; }

  /**
   * General output template for variables, strings, etc.
   * @tparam T
   * @param t
   * @return
   */
  template <typename T>
  Logger &operator<<(T const &t) {
    if (_msg_log_level <= _log_level && _enabled) *_current_stream << t;
    return *this;
  }

  /**
   * @overload
   * Specialized version for manipulators. e.g. std::endl
   * @param f
   * @return
   */
  Logger &operator<<(std::ostream &(*f)(std::ostream &)) {
    if (_msg_log_level <= _log_level && _enabled) *_current_stream << f;
    return *this;
  }

  /**
   * Specialized version for ios_base, e.g. for std::hex.
   * @overload
   * @param f
   * @return
   */
  Logger &operator<<(std::ios_base &(*f)(std::ios_base &)) {
    if (_msg_log_level <= _log_level && _enabled) f(*_current_stream);
    return *this;
  }
  /**
   * Specialized version for basic_ios
   * @tparam Ch
   * @tparam Tr
   * @param f
   * @return
   */
  template <class Ch, class Tr>
  Logger &operator<<(std::basic_ios<Ch, Tr> &(*f)(std::basic_ios<Ch, Tr> &)) {
    if (_msg_log_level <= _log_level && _enabled) f(*_current_stream);
    return *this;
  }

  /// Add log info in front of messages
  /**
   * Get a logger to the specific logLevel.
   * @param level
   * @return Reference to the logger
   */
  Logger &msg_level(logLevel level) {
    _msg_log_level = level;
    if (_msg_log_level <= _log_level && _enabled) {
      if (_msg_log_level > 4) {
        _current_stream = _log_stream;
      } else {
        _current_stream = _error_stream;
      }
      // Include timestamp
      auto t_sys = std::chrono::system_clock::now();
      std::time_t t_t = std::chrono::system_clock::to_time_t(t_sys);
      auto lt = localtime(&t_t);
      std::stringstream timestampstream;
      // maybe sprintf is easier here...
      timestampstream << std::setfill('0') << std::setw(4)
                      << (1900 + lt->tm_year) << std::setw(2)
                      << (1 + lt->tm_mon) << std::setw(2) << lt->tm_mday << "T"
                      << std::setw(2) << lt->tm_hour << std::setw(2)
                      << lt->tm_min << std::setw(2) << lt->tm_sec;
      *_current_stream << logLevelNames[level] << ":\t" << timestampstream.str()
                       << " ";
      auto t_hi_res = std::chrono::high_resolution_clock::now();
      *_current_stream << std::fixed << std::setw(8)
                       << std::chrono::duration_cast<std::chrono::microseconds>(
                              t_hi_res - _starttime)
                                  .count() /
                              1.E6
                       << "\t";

      *_current_stream << "[" << _identifier << "]\t";
    }
    return *this;
  }

  /**
   * Get a logger to write fatal error messages.
   * @return reference to the logger to print things.
   */
  Logger &fatal() { return msg_level(Fatal); }
  /**
   * Get a logger to write error messages.
   * @return reference to the logger to print things.
   */
  Logger &error() { return msg_level(Error); }
  /**
   * Get a logger to write warnings.
   * @return reference to the logger to print things.
   */
  Logger &warning() { return msg_level(Warning); }
  /**
   * Get a logger to write informative messages.
   * @return reference to the logger to print things.
   */
  Logger &info() { return msg_level(Info); }
  /**
   * Get a logger to write debug messages.
   * @return reference to the logger to print things.
   */
  Logger &debug() { return msg_level(Debug); }

  /// set log level
  logLevel set_log_level(logLevel l) {
    _log_level = l;
    return _log_level;
  }
  /// return log level
  logLevel get_log_level() { return _log_level; }

 private:
  logLevel _log_level;
  logLevel _msg_log_level;
  bool _enabled;
  std::string _filename;
  std::ostream *_log_stream;
  std::ostream *_error_stream;
  std::ostream *_current_stream;
  std::map<logLevel, std::string> logLevelNames;
  std::chrono::time_point<std::chrono::high_resolution_clock> _starttime;

  std::string _identifier;

  /// initilaize the list of log levels with the corresponding short names
  void init_log_levels() {
    logLevelNames.insert(
        std::pair<logLevel, std::string>(Fatal, "FATAL ERROR"));
    logLevelNames.insert(std::pair<logLevel, std::string>(Error, "ERROR"));
    logLevelNames.insert(std::pair<logLevel, std::string>(Warning, "WARNING"));
    logLevelNames.insert(std::pair<logLevel, std::string>(Info, "INFO"));
    logLevelNames.insert(std::pair<logLevel, std::string>(Debug, "DEBUG"));
  }

  /// initialize starting time
  void init_starting_time() {
    _starttime = std::chrono::high_resolution_clock::now();
  }

}; /* end of class Logger */
}  // namespace log
}  // namespace autopas
