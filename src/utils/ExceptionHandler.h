/**
 * @file ExceptionHandler.h
 * @author seckler
 * @date 17.04.18
 */

#pragma once

#include <cstdlib>
#include <functional>
#include "Logger.h"

namespace autopas {
namespace utils {
/**
   * enum that defines the behavior of the expectionhandling
   * please check the enum values for a more detailed description
   */
enum ExceptionBehavior {
  ignore,                    /// ignore all exceptions
  throwException,            /// throw the exception
  printAbort,                /// print the exception and
  printCustomAbortFunction,  /// print the exception and call a custom abort
  /// function
};

/**
 * Defines and handles the throwing and printing of exceptions.
 * This class defines what should happen if an error occurs within AutoPas.
 * For a detailed list please check the enum ExceptionBehavior
 */
class ExceptionHandler {
 public:

  /**
   * Set the behavior of the handler
   * @param behavior the behavior
   */
  static void setBehavior(ExceptionBehavior behavior) {
    std::lock_guard<std::mutex> guard(exceptionMutex);
    _behavior = behavior;
  }

  /**
   * Handle an exception derived by std::exception
   * @param e the exception to be handled
   */
  static void exception(const std::exception& e) {
    std::lock_guard<std::mutex> guard(exceptionMutex);
    switch (_behavior) {
      case ignore:
        // do nothing
        break;
      case throwException:
        throw e;
      case printAbort:
        AutoPasLogger->error("{}\naborting", e.what());
        AutoPasLogger->flush();
        std::abort();
      case printCustomAbortFunction:
        spdlog::get("AutoPasLog");
        AutoPasLogger->error("{}\nusing custom abort function", e.what());
        AutoPasLogger->flush();
        _customAbortFunction();
        break;
      default:
        break;
    }
  }

  /**
   * Handles an exception that is defined using the input string
   * @param exceptionString the string to describe the exception
   */
  static void exception(const std::string& exceptionString) {
    // no lock here, as a different public function is called!!!
    AutoPasException autoPasException(exceptionString);
    exception(autoPasException);
  }

  /**
   * Set a custom abort function
   * @param function the custom abort function
   */
  static void setCustomAbortFunction(std::function<void()> function) {
    std::lock_guard<std::mutex> guard(exceptionMutex);
    _customAbortFunction = function;
  }

 private:
  static std::mutex exceptionMutex;
  static ExceptionBehavior _behavior;
  static std::function<void()> _customAbortFunction;

  class AutoPasException : public std::exception {
   public:
    explicit AutoPasException(const std::string& description)
        : _description(description){};

    virtual const char* what() const throw() override {
      return _description.c_str();
    }

   private:
    std::string _description;
  };
};
}
}