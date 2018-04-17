/**
 * @file ExceptionHandler.h
 * @author seckler
 * @date 17.04.18
 */

#pragma once

#include <functional>
#include "utils/Logger.h"

namespace autopas {
namespace utils {
/**
 * Defines and handles the throwing and printing of exceptions.
 * This class defines what should happen if an error occurs within AutoPas.
 * For a detailed list please check the enum ExceptionBehavior
 */
class ExceptionHandler {
 public:
  /**
   * enum that defines the behavior of the expectionhandling
   * please check the enum values for a more detailed description
   */
  enum ExceptionBehavior {
    ignore,                    /// ignore all exceptions, t
    throwException,            /// throw the exception
    printAbort,                /// print the exception and
    printCustomAbortFunction,  /// print the exception and call a custom abort
                               /// function
  };

  /**
   * Constructor of the ExceptionHandler class
   * @param behavior the behavior the handler should use
   */
  explicit ExceptionHandler(ExceptionBehavior behavior)
      : _behavior(behavior), _customAbortFunction(abort) {}

  /**
   * Set the behavior of the handler
   * @param behavior the behavior
   */
  void setBehavior(ExceptionBehavior behavior) { _behavior = behavior; }

  /**
   * Handle an exception derived by std::exception
   * @param e the exception to be handled
   */
  void exception(const std::exception& e) {
    switch (_behavior) {
      case ignore:
        // do nothing
        break;
      case throwException:
        throw e;
      case printAbort:
        autopas::logger->fatal() << e.what() << std::endl;
        autopas::logger->fatal() << "aborting" << std::endl;
        abort();
      case printCustomAbortFunction:
        autopas::logger->fatal() << e.what() << std::endl;
        autopas::logger->fatal() << "using custom abort function" << std::endl;
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
  void exception(const std::string& exceptionString) {
    AutoPasException autoPasException(exceptionString);
    exception(autoPasException);
  }

  /**
   * Set a custom abort function
   * @param function the custom abort function
   */
  void setCustomAbortFunction(std::function<void()> function){
    _customAbortFunction = function;
  }

 private:
  ExceptionBehavior _behavior;
  std::function<void()> _customAbortFunction;

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