/**
 * @file ExceptionHandler.h
 * @author seckler
 * @date 17.04.18
 */

#pragma once

#include <cstdlib>
#include <functional>
#include <mutex>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>

#include "autopas/utils/logging/Logger.h"

namespace autopas::utils {

/**
 * Enum that defines the behavior of the expection handling.
 */
enum ExceptionBehavior {
  /**
   * Ignore all exceptions
   */
  ignore,
  /**
   * Throw the exception
   */
  throwException,
  /**
   * Print the exception and
   */
  printAbort,
  /**
   * Print the exception and call a custom abort function
   */
  printCustomAbortFunction,
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
  static void setBehavior(ExceptionBehavior behavior);

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
    switch (_behavior) {
      case throwException:
        throw e;  // NOLINT
      default:
        nonThrowException(e);
    }
  }

  /**
   * Handles an exception that is defined using exceptionString as well as multiple arguments.
   * It uses the same fmt library that is also used in spdlog and AutoPasLogger. Call it using e.g.:
   * exception("failure, because {} is not less than {}", 4, 3);
   *
   * @tparam First the template type of the first argument
   * @tparam Args types of multiple
   * @param exceptionString the basic exception string
   * @param first the first argument
   * @param args more arguments
   * @note First is needed to differentiate this from the template with exception(const Exception e)
   * @note this is a variadic function, and can thus incorporate an arbitrary amount of arguments
   */
  template <typename First, typename... Args>
  static void exception(std::string exceptionString, First first, Args... args);  // recursive variadic function

  /**
   * Rethrows the current exception or prints it.
   * Depending on the set behavior the currently active exception is either
   * rethrown, printed or otherwise handled.
   * @note Use this only inside a catch clause.
   */
  static void rethrow();

  /**
   * Set a custom abort function
   * @param function the custom abort function
   */
  static void setCustomAbortFunction(std::function<void()> function);

 private:
  static std::mutex exceptionMutex;
  static ExceptionBehavior _behavior;
  static std::function<void()> _customAbortFunction;

  static void nonThrowException(const std::exception &e);

  template <typename T>
  static auto formatArg(T &&arg) {
    using DecayedT = std::decay_t<T>;
    if constexpr (requires { std::forward<T>(arg).to_string(); }) {
      return std::forward<T>(arg).to_string();
    } else if constexpr (std::is_arithmetic_v<DecayedT> or std::is_convertible_v<T, std::string>) {
      return std::forward<T>(arg);
    } else if constexpr (requires(std::ostream &os) { os << std::forward<T>(arg); }) {
      std::ostringstream ss;
      ss << std::forward<T>(arg);
      return ss.str();
    } else if constexpr (std::is_enum_v<DecayedT>) {
      return static_cast<std::underlying_type_t<DecayedT>>(arg);
    } else {
      return std::forward<T>(arg);
    }
  }

 public:
  /**
   * Default exception class for autopas exceptions.
   * @note normally generated using ExceptionHandler::exception("some string")
   */
  class AutoPasException : public std::exception {
   public:
    /**
     * constructor
     * @param description a descriptive string
     */
    explicit AutoPasException(std::string description);

    /**
     * Copy constructor.
     * @param exception
     */
    AutoPasException(const AutoPasException &exception);

    ~AutoPasException() override;

    /**
     * returns the description
     * @return
     */
    [[nodiscard]] const char *what() const noexcept override;

   private:
    std::string _description;
  };
};

/**
 * Handles an exception that is defined using the input string
 * @param e the string to describe the exception
 */
template <>
void ExceptionHandler::exception(const std::string e);  // NOLINT

/**
 * Handles an exception that is defined using the input string
 * @param e the string to describe the exception
 */
template <>
void ExceptionHandler::exception(const char *const e);  // NOLINT

template <typename First, typename... Args>
void ExceptionHandler::exception(std::string exceptionString, First first, Args... args) {
  std::string s = fmt::format(fmt::runtime(exceptionString), formatArg(first), formatArg(args)...);
  exception(s);
}

}  // namespace autopas::utils
