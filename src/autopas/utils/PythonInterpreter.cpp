/**
 * @file PythonInterpreter.cpp
 * @author S. Newcome
 * @date 21.10.2026
 */

#include "autopas/utils/PythonInterpreter.h"

#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
#include <pybind11/embed.h>
#endif

#include "autopas/utils/ExceptionHandler.h"

namespace autopas::utils {

PythonInterpreter::PythonInterpreter() {
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  // There can only be one Python interpreter per process. If one is already running (e.g. started by a parent
  // project), pybind11 uses it automatically, so only start one if none is running.
  if (not Py_IsInitialized()) {
    pybind11::initialize_interpreter();
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
  }
#endif
}

}  // namespace autopas::utils