/**
 * @file ExceptionHandler.cpp
 * @author seckler
 * @date 24.04.18
 */

#include "ExceptionHandler.h"

std::mutex autopas::utils::ExceptionHandler::exceptionMutex;
autopas::utils::ExceptionBehavior autopas::utils::ExceptionHandler::_behavior =
    ExceptionBehavior::throwException;
std::function<void()> autopas::utils::ExceptionHandler::_customAbortFunction =
    abort;

template <>
void autopas::utils::ExceptionHandler::exception(
    const std::string exceptionString) {
  // no lock here, as a different public function is called!!!
  AutoPasException autoPasException(exceptionString);
  exception(autoPasException);
}

template <>
void autopas::utils::ExceptionHandler::exception(
    const char* const exceptionString) {
  exception(std::string(exceptionString));
}