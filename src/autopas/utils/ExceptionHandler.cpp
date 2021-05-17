/**
 * @file ExceptionHandler.cpp
 * @author seckler
 * @date 24.04.18
 */

#include "autopas/utils/ExceptionHandler.h"

std::mutex autopas::utils::ExceptionHandler::exceptionMutex;
autopas::utils::ExceptionBehavior autopas::utils::ExceptionHandler::_behavior = ExceptionBehavior::throwException;
std::function<void()> autopas::utils::ExceptionHandler::_customAbortFunction = abort;  // NOLINT

template <>
void autopas::utils::ExceptionHandler::exception(const std::string e) {  // NOLINT
  // no lock here, as a different public function is called!!!
  AutoPasException autoPasException(e);
  exception(autoPasException);
}

template <>
void autopas::utils::ExceptionHandler::exception(const char *const e) {
  exception(std::string(e));
}

void autopas::utils::ExceptionHandler::rethrow() {
  std::lock_guard<std::mutex> guard(exceptionMutex);
  std::exception_ptr p = std::current_exception();
  if (p == std::exception_ptr()) {
    exception("ExceptionHandler::rethrow() called outside of a catch block");
  }
  switch (_behavior) {
    case throwException:
      std::rethrow_exception(p);
    default:
      try {
        std::rethrow_exception(p);
      } catch (std::exception &e) {
        nonThrowException(e);
      }
  }
}
