/**
 * @file AutoPasPython.cpp
 * @date 03.11.2026
 *
 * Entry point of the Python module autopas._core.
 */
#include "AutoPasPython.h"

#include <pybind11/pybind11.h>

#include "autopas/utils/WrapOpenMP.h"

namespace py = pybind11;

PYBIND11_MODULE(_core, module) {
  module.doc() = "Python interface of the particle simulation library AutoPas.";
  module.attr("__version__") = AUTOPAS_PYTHON_VERSION;
#if defined(MD_FLEXIBLE_CALC_GLOBALS)
  module.attr("__calculates_globals__") = true;
#else
  module.attr("__calculates_globals__") = false;
#endif
  module.def("max_threads", &autopas::autopas_get_max_threads, "The number of threads used. Set via OMP_NUM_THREADS.");

  autopasPython::bindNativeFunctors(module);
  autopasPython::bindPythonFunctors(module);
  autopasPython::bindParticleProperties(module);
  autopasPython::bindGenerators(module);
  autopasPython::bindSimulation(module);
}
