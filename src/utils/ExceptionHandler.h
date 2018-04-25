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
   * Handle an exception derived by std::exception.
   * If the behavior is set to throw and this function is called in a catch
   * clause, instead of using the passed exception the underlying error is
   * rethrown.
   * @param e the exception to be handled
   * @tparam Exception the type of the exception, needed as throw only uses the
   * static type of e
   */
  template <class Exception>
  static void exception(const Exception e) {
    std::lock_guard<std::mutex> guard(exceptionMutex);
    std::exception_ptr p;
    switch (_behavior) {
      case ignore:
        // do nothing
        break;
      case throwException:
        p = std::current_exception();
        if (p == std::exception_ptr()) {
          throw e;
        } else {
          std::rethrow_exception(p);
        }
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

 public:
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

/**
 * Handles an exception that is defined using the input string
 * @param exceptionString the string to describe the exception
 */
template <>
void ExceptionHandler::exception(const std::string exceptionString);

/**
 * Handles an exception that is defined using the input string
 * @param exceptionString the string to describe the exception
 */
template <>
void ExceptionHandler::exception(const char* const exceptionString);
}
}