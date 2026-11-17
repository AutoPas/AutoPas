/**
 * @file PythonInterpreter.cpp
 * @author S. Newcome
 * @date 21.10.2026
 */

#include "autopas/utils/PythonInterpreter.h"

#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
#include <pybind11/embed.h>
#endif

#include <mutex>

#include "autopas/utils/ExceptionHandler.h"

namespace autopas::utils {

namespace {
/**
 * Whether a PythonInterpreter already finalized the interpreter in this process.
 */
[[maybe_unused]] bool pythonWasFinalized{false};
}  // namespace

PythonInterpreter::PythonInterpreter() {
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  // There can only be one Python interpreter per process. If one is already running (e.g. started by a parent
  // project), pybind11 uses it automatically, so only start one if none is running.
  if (not Py_IsInitialized()) {
    if (pythonWasFinalized) {
      utils::ExceptionHandler::exception(
          "PythonInterpreter: Python was already finalized in this process, and libraries like numpy cannot be loaded "
          "again in a restarted interpreter. To use Python in several AutoPas objects one after another, keep the "
          "interpreter alive in between by holding the shared pointer returned by "
          "autopas::utils::PythonInterpreter::getSharedInterpreter().");
    }
    pybind11::initialize_interpreter(false);
    _thisClassStartedInterpreter = true;
  }
#else
  utils::ExceptionHandler::exception(
      "PythonInterpreter constructed but AUTOPAS_ENABLE_PYTHON_BASED_TUNING=OFF! "
      "Set this CMake variable to ON to use Python-based tuning strategies.");
#endif
}

PythonInterpreter::~PythonInterpreter() {
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  if (_thisClassStartedInterpreter) {
    pybind11::finalize_interpreter();
    pythonWasFinalized = true;
  }
#endif
}

std::shared_ptr<PythonInterpreter> PythonInterpreter::getSharedInterpreter() {
  // We use a static weak_ptr to ensure each caller points to the same interpreter, if one exists. As weak_ptrs do not
  // own that which they point to, we need not worry about the static weak_ptr keeping the interpreter alive after the
  // last caller's destruction, the same way a static shared_ptr would.

  // We use a locked mutex to avoid multiple threads attempting this at once (not currently possible but provides us
  // with some future-proofing). This is static so all threads have the same mutex.
  static std::mutex interpreterMutex;
  static std::weak_ptr<PythonInterpreter> interpreterWeak;
  // Lock the mutex to make this thread-safe.
  const std::scoped_lock lock(interpreterMutex);
  // weak_ptr::lock() returns a shared_ptr that is null if no interpreter has yet been created or non-null if it has
  // been created (set below).
  auto interpreter = interpreterWeak.lock();
  if (not interpreter) {
    // weak_ptr returned a null pointer => create a shared interpreter. std::make_shared cannot call the private
    // constructor, so use new. Otherwise, this is the same as std::make_shared: the shared_ptr still deletes the
    // interpreter when the last owner releases it.
    interpreter = std::shared_ptr<PythonInterpreter>(new PythonInterpreter());
    // Point the static weak reference at the new interpreter so the next caller's .lock() finds and shares it.
    interpreterWeak = interpreter;
  }
  return interpreter;
}

}  // namespace autopas::utils