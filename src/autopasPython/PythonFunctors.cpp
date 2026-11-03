/**
 * @file PythonFunctors.cpp
 * @date 03.11.2026
 */
#include "PythonFunctors.h"

#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <memory>
#include <optional>
#include <stdexcept>

#include "FunctorHandle.h"
#include "PythonFunctorAdapters.h"
#include "autopas/utils/WrapOpenMP.h"
#include "templateInstantiations/computeInteractionsPythonPairwise.h"
#include "templateInstantiations/computeInteractionsPythonTriwise.h"

//! @cond Doxygen_Suppress
extern template class autopas::AutoPas<ParticleType>;
//! @endcond

namespace py = pybind11;

namespace autopasPython {
namespace {

/**
 * A functor implemented in Python. The Python subclass derives from this, without a trampoline, because C++ looks up
 * the methods by name once per traversal and not via a virtual call per pair.
 * @tparam Adapter_T The AutoPas functor which forwards to Python.
 * @tparam interactionType_v Its interaction type.
 */
template <class Adapter_T, autopas::InteractionTypeOption::Value interactionType_v>
class PythonFunctor final : public FunctorHandle {
 public:
  /**
   * Constructor.
   * @param cutoff Only (sets of) particles within this distance reach Python. Defaults to the one of the simulation.
   * @param allowsNewton3
   * @param allowsNonNewton3
   */
  PythonFunctor(std::optional<double> cutoff, bool allowsNewton3, bool allowsNonNewton3)
      : _cutoff(cutoff), _allowsNewton3(allowsNewton3), _allowsNonNewton3(allowsNonNewton3) {
    if (not(allowsNewton3 or allowsNonNewton3)) {
      throw std::invalid_argument("A functor has to allow at least one of Newton 3 on and off.");
    }
  }

  [[nodiscard]] autopas::InteractionTypeOption::Value interactionType() const override { return interactionType_v; }
  [[nodiscard]] bool isPython() const override { return true; }
  [[nodiscard]] bool allowsNewton3() const override { return _allowsNewton3; }
  [[nodiscard]] bool allowsNonNewton3() const override { return _allowsNonNewton3; }

  [[nodiscard]] std::string name() const override {
    const py::gil_scoped_acquire gil{};
    return py::str(py::type::of(self()).attr("__name__"));
  }

  void bind(double cutoff, ParticlePropertiesLibraryType & /*particlePropertiesLibrary*/,
            bool relevantForTuning) override {
    if (_adapter) {
      throw std::invalid_argument(name() + " is already used by another AutoPas object.");
    }
    if (_cutoff.value_or(cutoff) > cutoff) {
      throw std::invalid_argument("The cutoff of " + name() + " exceeds the cutoff of the simulation.");
    }
    _adapter = std::make_unique<Adapter_T>(_cutoff.value_or(cutoff),
                                           PythonFunctorTraits{name(), _allowsNewton3, _allowsNonNewton3,
                                                               relevantForTuning},
                                           _dispatcher);
  }

  bool compute(AutoPasType &autoPas) override {
    {
      const py::gil_scoped_acquire gil{};
      _dispatcher.acquire(self());
    }
    /// Restores the state also if AutoPas itself throws.
    struct TraversalScope {
      /// All kernels are serialized by the GIL, so more threads would only compete for it. Measured with 512 particles
      /// on 10 threads, this makes an iteration 3.6 times slower than on a single thread.
      explicit TraversalScope(PythonDispatcher &dispatcher)
          : _dispatcher(dispatcher), _maxThreads(autopas::autopas_get_max_threads()) {
        autopas::autopas_set_num_threads(1);
      }
      ~TraversalScope() {
        autopas::autopas_set_num_threads(_maxThreads);
        _dispatcher.release();
      }
      PythonDispatcher &_dispatcher;
      int _maxThreads;
    };
    bool tuning{};
    {
      const TraversalScope scope{_dispatcher};
      tuning = autoPas.computeInteractions(_adapter.get());
    }
    _dispatcher.rethrow();
    return tuning;
  }

  [[nodiscard]] double potentialEnergy() const override { return callDouble("potential_energy"); }

  [[nodiscard]] double virial() const override { return callDouble("virial"); }

  /**
   * The Python object of this functor.
   * @return The instance of the Python subclass.
   */
  [[nodiscard]] py::object self() const { return py::cast(this, py::return_value_policy::reference); }

 private:
  [[nodiscard]] double callDouble(const char *method) const {
    const py::gil_scoped_acquire gil{};
    return self().attr(method)().template cast<double>();
  }

  std::optional<double> _cutoff;
  bool _allowsNewton3;
  bool _allowsNonNewton3;
  PythonDispatcher _dispatcher{};
  std::unique_ptr<Adapter_T> _adapter{};
};

/**
 * Exposes the base class for functors written in Python.
 */
template <class Handle>
void bindBase(py::module_ &module, const char *name, const char *doc, const char *aosSignature) {
  py::class_<Handle, FunctorHandle>(module, name, doc)
      .def(py::init<std::optional<double>, bool, bool>(), py::arg("cutoff") = py::none(),
           py::arg("allows_newton3") = true, py::arg("allows_non_newton3") = true)
      .def(
          "aos",
          [aosSignature](const py::object &, const py::args &, const py::kwargs &) {
            PyErr_SetString(PyExc_NotImplementedError,
                            (std::string("Functors have to override ") + aosSignature).c_str());
            throw py::error_already_set();
          },
          aosSignature)
      .def(
          "init_traversal", [](const py::object &) {}, "Called once before each traversal.")
      .def(
          "end_traversal", [](const py::object &, bool) {}, py::arg("newton3"), "Called once after each traversal.")
      // These shadow Functor's versions, which would otherwise call back into the Python method and never terminate.
      .def(
          "potential_energy", [](const py::object &) { return std::numeric_limits<double>::quiet_NaN(); },
          "Override this to feed the energy history.")
      .def(
          "virial", [](const py::object &) { return std::numeric_limits<double>::quiet_NaN(); },
          "Override this to feed the virial history.");
}

/**
 * Copies a three dimensional array to NumPy.
 */
py::array_t<double> toArray(const std::array<double, 3> &values) { return py::array_t<double>(3, values.data()); }

}  // namespace

void bindPythonFunctors(py::module_ &module) {
  py::class_<ParticleType>(module, "Particle",
                           "A molecule. The instances handed to a functor are only valid during the call.")
      .def(py::init<const std::array<double, 3> &, const std::array<double, 3> &, unsigned long, unsigned long>(),
           py::arg("r"), py::arg("v") = std::array<double, 3>{0., 0., 0.}, py::arg("id") = 0ul,
           py::arg("type_id") = 0ul)
      .def_property_readonly(
          "r", [](const ParticleType &p) { return toArray(p.getR()); }, "Position (copy).")
      .def_property_readonly(
          "v", [](const ParticleType &p) { return toArray(p.getV()); }, "Velocity (copy).")
      .def_property(
          "f", [](const ParticleType &p) { return toArray(p.getF()); }, &ParticleType::setF, "Force (copy).")
      .def_property_readonly(
          "old_f", [](const ParticleType &p) { return toArray(p.getOldF()); }, "Force of the previous iteration (copy).")
      .def_property_readonly("id", &ParticleType::getID)
      .def_property_readonly("type_id", &ParticleType::getTypeId)
      .def_property_readonly("is_owned", &ParticleType::isOwned)
      .def_property_readonly("is_halo", &ParticleType::isHalo)
      .def("add_f", &ParticleType::addF, py::arg("f"), "Adds to the force.")
      .def("sub_f", &ParticleType::subF, py::arg("f"), "Subtracts from the force.")
      .def("__repr__", [](const ParticleType &p) {
        return "Particle(id=" + std::to_string(p.getID()) + ", type_id=" + std::to_string(p.getTypeId()) + ")";
      });

  bindBase<PythonFunctor<PythonPairwiseAdapter, autopas::InteractionTypeOption::pairwise>>(
      module, "PairwiseFunctor",
      "Base class for pairwise functors written in Python. Runs with the AoS data layout on one thread only, so this is "
      "orders of magnitude slower than native functors.",
      "aos(self, i, j, newton3)");
  bindBase<PythonFunctor<PythonTriwiseAdapter, autopas::InteractionTypeOption::triwise>>(
      module, "TriwiseFunctor",
      "Base class for triwise functors written in Python. Runs with the AoS data layout on one thread only, so this is "
      "orders of magnitude slower than native functors.",
      "aos(self, i, j, k, newton3)");
}

}  // namespace autopasPython
