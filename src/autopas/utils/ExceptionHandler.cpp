/**
 * @file ExceptionHandler.cpp
 * @author seckler
 * @date 24.04.18
 */

#include "autopas/utils/ExceptionHandler.h"

#include "autopas/utils/WrapMPI.h"

std::mutex autopas::utils::ExceptionHandler::exceptionMutex;
autopas::utils::ExceptionBehavior autopas::utils::ExceptionHandler::_behavior = ExceptionBehavior::throwException;
std::function<void()> autopas::utils::ExceptionHandler::_customAbortFunction = abort;  // NOLINT

template <>
void autopas::utils::ExceptionHandler::exception(const std::string e) {  // NOLINT
  // no lock here, as a different public function is called!!!
  int myRank{};
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &myRank);
  AutoPasException autoPasException("Rank " + std::to_string(myRank) + " : " + e);
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

void autopas::utils::ExceptionHandler::setBehavior(autopas::utils::ExceptionBehavior behavior) {
  std::lock_guard<std::mutex> guard(exceptionMutex);
  _behavior = behavior;
}

void autopas::utils::ExceptionHandler::setCustomAbortFunction(std::function<void()> function) {
  std::lock_guard<std::mutex> guard(exceptionMutex);
  _customAbortFunction = std::move(function);
}

void autopas::utils::ExceptionHandler::nonThrowException(const std::exception &e) {
  switch (_behavior) {
    case ignore:
      // do nothing
      break;
    case printAbort:
      AutoPasLog(error, "{}\naborting", e.what());
      std::abort();
    case printCustomAbortFunction:
      spdlog::get("AutoPasLog");
      AutoPasLog(error, "{}\nusing custom abort function", e.what());
      _customAbortFunction();
      break;
    default:
      break;
  }
}

autopas::utils::ExceptionHandler::AutoPasException::AutoPasException(std::string description)
    : _description(std::move(description)) {}

autopas::utils::ExceptionHandler::AutoPasException::AutoPasException(
    const autopas::utils::ExceptionHandler::AutoPasException &exception) = default;

autopas::utils::ExceptionHandler::AutoPasException::~AutoPasException() = default;

const char *autopas::utils::ExceptionHandler::AutoPasException::what() const noexcept { return _description.c_str(); }
