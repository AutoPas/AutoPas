/**
 * @file NativeFunctors.cpp
 * @date 03.11.2026
 */
#include "NativeFunctors.h"

#include <memory>
#include <stdexcept>
#include <type_traits>

#include "FunctorHandle.h"

// All of these are instantiated in md-flexible's templateInstantiations.
//! @cond Doxygen_Suppress
extern template class autopas::AutoPas<ParticleType>;
#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
extern template bool autopas::AutoPas<ParticleType>::computeInteractions(LJFunctorTypeAutovec *);
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_AVX) && defined(__AVX__)
extern template bool autopas::AutoPas<ParticleType>::computeInteractions(LJFunctorTypeAVX *);
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_SVE) && defined(__ARM_FEATURE_SVE)
extern template bool autopas::AutoPas<ParticleType>::computeInteractions(LJFunctorTypeSVE *);
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_HWY)
extern template bool autopas::AutoPas<ParticleType>::computeInteractions(LJFunctorTypeHWY *);
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC)
extern template bool autopas::AutoPas<ParticleType>::computeInteractions(ATMFunctor *);
#endif
//! @endcond

namespace py = pybind11;

namespace autopasPython {
namespace {

/**
 * A functor implemented in C++.
 * @tparam Functor_T The AutoPas functor.
 * @tparam interactionType_v Its interaction type.
 */
template <class Functor_T, autopas::InteractionTypeOption::Value interactionType_v>
class NativeFunctor final : public FunctorHandle {
 public:
  /**
   * Constructor.
   * @param name Name of the Python class.
   */
  explicit NativeFunctor(std::string name) : _name(std::move(name)) {}

  [[nodiscard]] autopas::InteractionTypeOption::Value interactionType() const override { return interactionType_v; }
  [[nodiscard]] bool isPython() const override { return false; }
  [[nodiscard]] bool allowsNewton3() const override { return true; }
  [[nodiscard]] bool allowsNonNewton3() const override { return true; }
  [[nodiscard]] std::string name() const override { return _name; }

  void bind(double cutoff, ParticlePropertiesLibraryType &particlePropertiesLibrary, bool relevantForTuning) override {
    if (_functor) {
      throw std::invalid_argument(_name + " is already used by another AutoPas object.");
    }
    if (not relevantForTuning) {
      throw std::invalid_argument("Only one native functor per interaction type can be used.");
    }
    _functor = std::make_unique<Functor_T>(cutoff, std::ref(particlePropertiesLibrary));
  }

  bool compute(AutoPasType &autoPas) override {
    const auto tuning = autoPas.computeInteractions(_functor.get());
    if constexpr (mdFlexibleTypeDefs::calcGlobals) {
      _potentialEnergy = _functor->getPotentialEnergy();
      _virial = _functor->getVirial();
    }
    return tuning;
  }

  [[nodiscard]] double potentialEnergy() const override { return _potentialEnergy; }
  [[nodiscard]] double virial() const override { return _virial; }

 private:
  std::string _name;
  std::unique_ptr<Functor_T> _functor{};
  double _potentialEnergy{std::numeric_limits<double>::quiet_NaN()};
  double _virial{std::numeric_limits<double>::quiet_NaN()};
};

/**
 * Exposes one native functor as Python class without any further configuration: The cutoff and the particle
 * properties are given by the simulation it is handed over to.
 */
template <class Functor_T, autopas::InteractionTypeOption::Value interactionType_v>
void bindNative(py::module_ &module, const char *name, const char *doc) {
  using Handle = NativeFunctor<Functor_T, interactionType_v>;
  py::class_<Handle, FunctorHandle>(module, name, doc).def(py::init([name]() { return new Handle(name); }));
}

}  // namespace

void bindNativeFunctors(py::module_ &module) {
  py::class_<FunctorHandle>(module, "Functor", "Base class of everything which can be handed over to AutoPas.")
      .def_property_readonly("name", &FunctorHandle::name, "Name of the functor.")
      .def_property_readonly(
          "interaction_type",
          [](const FunctorHandle &functor) {
            return autopas::InteractionTypeOption{functor.interactionType()}.to_string();
          },
          "Either 'pairwise' or 'triwise'.")
      .def("potential_energy", &FunctorHandle::potentialEnergy, "Potential energy of the latest iteration, or NaN.")
      .def("virial", &FunctorHandle::virial, "Virial of the latest iteration, or NaN.");

  using autopas::InteractionTypeOption;
#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
  bindNative<LJFunctorTypeAutovec, InteractionTypeOption::pairwise>(module, "LennardJones",
                                                                    "12-6 Lennard-Jones, auto vectorized.");
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_AVX) && defined(__AVX__)
  bindNative<LJFunctorTypeAVX, InteractionTypeOption::pairwise>(module, "LennardJonesAVX",
                                                                "12-6 Lennard-Jones, AVX intrinsics.");
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_SVE) && defined(__ARM_FEATURE_SVE)
  bindNative<LJFunctorTypeSVE, InteractionTypeOption::pairwise>(module, "LennardJonesSVE",
                                                                "12-6 Lennard-Jones, SVE intrinsics.");
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_HWY)
  bindNative<LJFunctorTypeHWY, InteractionTypeOption::pairwise>(module, "LennardJonesHWY",
                                                                "12-6 Lennard-Jones, Google Highway SIMD. With the SoA data layout, potential "
                                                                "energy and virial are wrong as soon as there are halo particles.");
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC)
  bindNative<ATMFunctor, InteractionTypeOption::triwise>(module, "AxilrodTellerMuto",
                                                         "Axilrod-Teller-Muto three-body potential, auto vectorized.");
#endif
}

}  // namespace autopasPython
