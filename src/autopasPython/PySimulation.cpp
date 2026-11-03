/**
 * @file PySimulation.cpp
 * @date 03.11.2026
 */
#include "PySimulation.h"

#include <pybind11/stl.h>
#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

#include "ConfigConversion.h"
#include "Generators.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"
#include "src/Thermostat.h"
#include "src/TimeDiscretization.h"
#include "src/configuration/YamlParser.h"
#include "src/configureAutoPas.h"

//! @cond Doxygen_Suppress
extern template class autopas::AutoPas<ParticleType>;
//! @endcond

namespace py = pybind11;
using py::literals::operator""_a;

namespace autopasPython {
namespace {

/**
 * Loads the configuration. Generators in the dictionary are moved to the objects of the configuration directly, as
 * they have no representation in YAML.
 */
void parseConfig(const py::object &config, MDFlexConfig &result) {
  YAML::Node node{};
  if (py::isinstance<py::dict>(config)) {
    auto dict = py::dict(config.attr("copy")());
    if (dict.contains("Objects") and not py::isinstance<py::dict>(dict["Objects"])) {
      for (const auto &object : py::iterable(dict["Objects"])) {
        result.particleObjects.push_back(object.cast<std::shared_ptr<Object>>());
      }
      PyDict_DelItemString(dict.ptr(), "Objects");
    }
    node = toYaml(dict);
  } else {
    const auto path = py::str(py::module_::import("os").attr("fspath")(config)).cast<std::string>();
    if (not std::filesystem::exists(path)) {
      PyErr_SetString(PyExc_FileNotFoundError, ("The configuration file does not exist: " + path).c_str());
      throw py::error_already_set();
    }
    result.yamlFilename.value = path;
    node = YAML::LoadFile(path);
  }
  const auto errors = MDFlexParser::YamlParser::parseYamlNode(node, result);
  if (not errors.empty()) {
    std::ostringstream message;
    message << "Invalid configuration:";
    for (const auto &error : errors) {
      message << "\n" << error;
    }
    throw std::invalid_argument(message.str());
  }
}

/**
 * Removes what is not supported from a set of options.
 */
template <class Option_T, class Predicate_T>
void restrictOptions(std::set<Option_T> &options, Predicate_T keep, const std::string &functorName,
                     const std::string &optionName) {
  std::erase_if(options, [&](const auto &option) { return not keep(option); });
  if (options.empty()) {
    throw std::invalid_argument("The configuration leaves no " + optionName + " which " + functorName + " supports.");
  }
}

}  // namespace

PySimulation::PySimulation(py::object config, py::object functors,
                           std::shared_ptr<ParticlePropertiesLibraryType> particleProperties, py::object objects,
                           size_t historyInterval)
    : _particleProperties(std::move(particleProperties)), _historyInterval(historyInterval) {
  parseConfig(config, _config);
  if (not objects.is_none()) {
    for (const auto &object : py::iterable(objects)) {
      _config.particleObjects.push_back(object.cast<std::shared_ptr<Object>>());
    }
  }
  _config.finalizeConfiguration(/*initPropertiesLibrary*/ _particleProperties == nullptr);

  if (_particleProperties) {
    if (_particleProperties->getCutoff() != _config.cutoff.value) {
      throw std::invalid_argument("The cutoff of the ParticlePropertiesLibrary differs from the one of the simulation.");
    }
    _particleProperties->calculateMixingCoefficients();
  } else {
    _particleProperties = _config.getParticlePropertiesLibrary();
  }

  setupFunctors(std::move(functors));

  _domainDecomposition = std::make_shared<RegularGridDecomposition>(_config);

  if (not _config.vtkFileName.value.empty()) {
    _vtkWriter.emplace(_config.vtkFileName.value, _config.vtkOutputFolder.value,
                       std::to_string(_config.iterations.value).size());
  }

  _autoPas = std::make_unique<AutoPasType>(std::cout);
  configureAutoPas(*_autoPas, _config, _domainDecomposition->getLocalBoxMin(), _domainDecomposition->getLocalBoxMax(),
                   "Rank0_" + _config.outputSuffix.value);
  // The functors take precedence over what the configuration states.
  std::set<autopas::InteractionTypeOption> interactionTypes{};
  for (const auto *functor : _functors) {
    interactionTypes.insert(functor->interactionType());
  }
  _autoPas->setAllowedInteractionTypeOptions(interactionTypes);
  _autoPas->init();

  if (_autoPas->searchSpaceIsTrivial() and _config.tuningPhases.value > 0) {
    throw std::invalid_argument(
        "The search space must not be trivial if the simulation time is limited by the number tuning phases.");
  }

  if (not _config.checkpointfile.value.empty()) {
    _config.flushParticles();
    _config.loadParticlesFromCheckpoint(0, 1);
  }
  for (auto &particle : _config.particles) {
    validateParticle(particle);
    _autoPas->addParticle(particle);
  }
  _config.flushParticles();

  if (_config.useThermostat.value and _config.deltaT.value != 0) {
    if (_config.addBrownianMotion.value) {
      Thermostat::addBrownianMotion(*_autoPas, *_particleProperties, _config.initTemperature.value);
    }
    Thermostat::apply(*_autoPas, *_particleProperties, _config.initTemperature.value,
                      std::numeric_limits<double>::max());
  }
}

PySimulation::~PySimulation() = default;

void PySimulation::setupFunctors(py::object functors) {
  // Without functors, make the same choice as md-flexible.
  if (functors.is_none()) {
    const auto module = py::module_::import("autopas._core");
    auto chosen = py::list{};
    const auto add = [&](const char *name) {
      if (not py::hasattr(module, name)) {
        throw std::invalid_argument(std::string("The configuration requests a functor which was not compiled: ") +
                                    name);
      }
      chosen.append(module.attr(name)());
    };
    switch (_config.functorOption.value) {
      case MDFlexConfig::FunctorOption::lj12_6:
        add("LennardJones");
        break;
      case MDFlexConfig::FunctorOption::lj12_6_AVX:
        add("LennardJonesAVX");
        break;
      case MDFlexConfig::FunctorOption::lj12_6_SVE:
        add("LennardJonesSVE");
        break;
      case MDFlexConfig::FunctorOption::lj12_6_HWY:
        add("LennardJonesHWY");
        break;
      case MDFlexConfig::FunctorOption::none:
        break;
    }
    if (_config.functorOption3B.value == MDFlexConfig::FunctorOption3B::at) {
      add("AxilrodTellerMuto");
    }
    if (chosen.empty()) {
      for (const auto *fallback : {"LennardJonesHWY", "LennardJonesAVX", "LennardJonesSVE", "LennardJones"}) {
        if (py::hasattr(module, fallback)) {
          add(fallback);
          break;
        }
      }
      PyErr_WarnEx(PyExc_UserWarning,
                   ("No functor was specified. Defaulting to " + py::str(chosen[0]).cast<std::string>()).c_str(), 1);
    }
    functors = chosen;
  } else if (py::isinstance<FunctorHandle>(functors)) {
    py::list single{};
    single.append(functors);
    functors = single;
  }

  std::vector<py::object> objects{};
  for (const auto &functor : py::iterable(functors)) {
    if (not py::isinstance<FunctorHandle>(functor)) {
      throw py::type_error("The functors have to be instances of autopas.Functor.");
    }
    objects.push_back(py::reinterpret_borrow<py::object>(functor));
  }
  if (objects.empty()) {
    throw std::invalid_argument("At least one functor is required.");
  }
  // Native functors are compiled as relevant for tuning. Hence, they have to be the ones in charge.
  std::ranges::stable_sort(objects, {},
                           [](const py::object &object) { return object.cast<FunctorHandle *>()->isPython(); });

  std::map<autopas::InteractionTypeOption::Value, size_t> numPerType{};
  bool anyPython = false;
  for (const auto &object : objects) {
    auto *functor = object.cast<FunctorHandle *>();
    const auto type = functor->interactionType();
    const auto relevantForTuning = numPerType[type]++ == 0;
    if (not functor->isPython() and not relevantForTuning) {
      throw std::invalid_argument("Only one native functor per interaction type can be used, but got a second one: " +
                                  functor->name());
    }
    anyPython |= functor->isPython();

    // One configuration is shared by all functors of a type, so all have to support it.
    const auto pairwise = type == autopas::InteractionTypeOption::pairwise;
    auto &newton3 = pairwise ? _config.newton3Options.value : _config.newton3Options3B.value;
    restrictOptions(
        newton3,
        [&](const autopas::Newton3Option &option) {
          return option == autopas::Newton3Option::enabled ? functor->allowsNewton3() : functor->allowsNonNewton3();
        },
        functor->name(), "Newton 3 option");
    if (functor->isPython()) {
      auto &dataLayouts = pairwise ? _config.dataLayoutOptions.value : _config.dataLayoutOptions3B.value;
      restrictOptions(
          dataLayouts, [](const autopas::DataLayoutOption &option) { return option == autopas::DataLayoutOption::aos; },
          functor->name(), "data layout (Python functors need AoS)");
    }

    functor->bind(_config.cutoff.value, *_particleProperties, relevantForTuning);
    _functors.push_back(functor);
  }
  _functorObjects = std::move(objects);

  // @todo Remove this as soon as LJFunctorHWY treats halo particles correctly: Its SoA kernels use a mask for
  // "not a dummy" as if it meant "owned". Forces are not affected, but halo pairs enter energy and virial with full
  // weight, so the history of a simulation depends on the configuration selected by the tuner.
  if constexpr (mdFlexibleTypeDefs::calcGlobals) {
    const auto usesHWY = std::ranges::any_of(_functors, [](const auto *f) { return f->name() == "LennardJonesHWY"; });
    if (usesHWY and _config.dataLayoutOptions.value.contains(autopas::DataLayoutOption::soa)) {
      if (PyErr_WarnEx(PyExc_UserWarning,
                       "LennardJonesHWY reports a wrong potential energy and virial for the SoA data layout if there "
                       "are halo particles. Forces are correct. Use LennardJones or restrict data-layout to AoS if you "
                       "need the energies.",
                       1) != 0) {
        throw py::error_already_set();
      }
    }
  }

  // The benchmark would call the functors with made up particles.
  if (anyPython) {
    _config.useSortingThresholdBenchmark.value = false;
  }
}

void PySimulation::validateParticle(const ParticleType &particle) const {
  if (static_cast<long>(particle.getTypeId()) >= _particleProperties->getNumberRegisteredSiteTypes()) {
    throw std::invalid_argument("Particle " + std::to_string(particle.getID()) + " has the unknown type " +
                                std::to_string(particle.getTypeId()) + ".");
  }
  if (not autopas::utils::inBox(particle.getR(), _config.boxMin.value, _config.boxMax.value)) {
    throw std::invalid_argument("Particle " + std::to_string(particle.getID()) + " is outside of the domain.");
  }
}

void PySimulation::ensureAccessible() const {
  if (_stepping) {
    throw std::runtime_error("This can not be accessed while an iteration is being computed.");
  }
  if (_failed) {
    throw std::runtime_error(
        "An earlier iteration failed half way through, so the particles are in an inconsistent state. Create a new "
        "AutoPas object.");
  }
}

void PySimulation::run(const py::object &self, std::optional<size_t> iterations, const py::object &callback,
                       size_t callbackInterval) {
  ensureAccessible();
  if (_running) {
    throw std::runtime_error("The simulation is already running. It can not be started from within a callback.");
  }
  if (callbackInterval == 0) {
    throw std::invalid_argument("callback_interval has to be positive.");
  }

  /// Resets the flag no matter how the loop is left. Everything up to the last completed iteration is valid, so a
  /// KeyboardInterrupt or a raising callback do not prevent resuming.
  struct RunningGuard {
    explicit RunningGuard(bool &flag) : _flag(flag) { _flag = true; }
    ~RunningGuard() { _flag = false; }
    bool &_flag;
  } guard{_running};

  const auto target = iterations ? _iteration + *iterations : 0ul;
  const auto needsMoreIterations = [&]() {
    if (iterations) {
      return _iteration < target;
    }
    return _iteration < _config.iterations.value or _numTuningPhasesCompleted < _config.tuningPhases.value;
  };

  while (needsMoreIterations()) {
    _stepping = true;
    try {
      const py::gil_scoped_release noGil{};
      step();
    } catch (...) {
      _stepping = false;
      _failed = true;
      throw;
    }
    _stepping = false;

    if (PyErr_CheckSignals() != 0) {
      throw py::error_already_set();
    }
    if (not callback.is_none() and not _paused and _iteration % callbackInterval == 0) {
      const auto result = callback(self);
      if (result.is(py::bool_(false))) {
        break;
      }
    }
  }

  if (_vtkWriter.has_value()) {
    _vtkWriter->recordTimestep(_iteration, *_autoPas, *_domainDecomposition);
  }
}

void PySimulation::step() {
  const auto deltaT = _config.deltaT.value;

  if (_vtkWriter.has_value() and _iteration % _config.vtkWriteFrequency.value == 0) {
    _vtkWriter->recordTimestep(_iteration, *_autoPas, *_domainDecomposition);
  }

  if (deltaT != 0 and not _paused) {
    TimeDiscretization::calculatePositionsAndResetForces(*_autoPas, *_particleProperties, deltaT,
                                                         _config.globalForce.value, _config.fastParticlesThrow.value);
  }

  // Also when not moving, because this advances the iteration counter of AutoPas.
  auto emigrants = _autoPas->updateContainer();

  if (deltaT != 0 and not _paused) {
    _domainDecomposition->exchangeMigratingParticles(*_autoPas, emigrants);
    _domainDecomposition->reflectParticlesAtBoundaries(*_autoPas, *_particleProperties);
    _domainDecomposition->exchangeHaloParticles(*_autoPas);
  }

  computeInteractions();

  if (_config.pauseSimulationDuringTuning.value) {
    if (_currentIterationIsTuningIteration and not _previousIterationWasTuningIteration) {
      _paused = true;
    }
    if (_previousIterationWasTuningIteration and not _currentIterationIsTuningIteration) {
      // Reset the forces accumulated during tuning and compute those of the latest iteration again.
      for (auto particle = _autoPas->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
        particle->setF(_config.globalForce.value);
      }
      computeInteractions();
      _paused = false;
    }
  }

  if (deltaT != 0 and not _paused) {
    TimeDiscretization::calculateVelocities(*_autoPas, *_particleProperties, deltaT);
    if (_config.useThermostat.value and _iteration % _config.thermostatInterval.value == 0) {
      Thermostat::apply(*_autoPas, *_particleProperties, _config.targetTemperature.value, _config.deltaTemp.value);
    }
  }

  if (not _paused) {
    if (_historyInterval != 0 and _iteration % _historyInterval == 0) {
      record();
    }
    ++_iteration;
  }
}

void PySimulation::computeInteractions() {
  _previousIterationWasTuningIteration = _currentIterationIsTuningIteration;
  _currentIterationIsTuningIteration = false;
  for (auto *functor : _functors) {
    _currentIterationIsTuningIteration |= functor->compute(*_autoPas);
  }
  if (_currentIterationIsTuningIteration) {
    ++_numTuningIterations;
  } else if (_previousIterationWasTuningIteration) {
    ++_numTuningPhasesCompleted;
  }
}

void PySimulation::record() {
  double potentialEnergySum = 0.;
  double virialSum = 0.;
  for (const auto *functor : _functors) {
    potentialEnergySum += functor->potentialEnergy();
    virialSum += functor->virial();
  }
  const auto kinetic = kineticEnergy();
  const auto count = _autoPas->getNumberOfParticles(autopas::IteratorBehavior::owned);
  _historyIteration.push_back(_iteration);
  _historyPotentialEnergy.push_back(potentialEnergySum);
  _historyVirial.push_back(virialSum);
  _historyKineticEnergy.push_back(kinetic);
  _historyTemperature.push_back(count == 0 ? std::numeric_limits<double>::quiet_NaN()
                                           : kinetic / (1.5 * static_cast<double>(count)));
  _historyTuning.push_back(_currentIterationIsTuningIteration);
}

double PySimulation::kineticEnergy() const {
  double energy = 0.;
  AUTOPAS_OPENMP(parallel reduction(+ : energy))
  for (auto particle = _autoPas->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    const auto &v = particle->getV();
    energy += 0.5 * _particleProperties->getSiteMass(particle->getTypeId()) *
              (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  }
  return energy;
}

double PySimulation::temperature() const {
  ensureAccessible();
  const auto count = numParticles();
  return count == 0 ? std::numeric_limits<double>::quiet_NaN() : kineticEnergy() / (1.5 * static_cast<double>(count));
}

double PySimulation::potentialEnergy() const {
  double sum = 0.;
  for (const auto *functor : _functors) {
    sum += functor->potentialEnergy();
  }
  return sum;
}

double PySimulation::virial() const {
  double sum = 0.;
  for (const auto *functor : _functors) {
    sum += functor->virial();
  }
  return sum;
}

size_t PySimulation::numParticles() const {
  ensureAccessible();
  return _autoPas->getNumberOfParticles(autopas::IteratorBehavior::owned);
}

void PySimulation::addParticles(const py::array_t<double> &positions, const py::object &velocities,
                                const py::object &typeIds, const py::object &ids) {
  ensureAccessible();
  if (positions.ndim() != 2 or positions.shape(1) != 3) {
    throw std::invalid_argument("positions has to be of shape (n, 3).");
  }
  const auto n = positions.shape(0);

  const auto numpy = py::module_::import("numpy");
  const auto broadcast = [&](const py::object &values, const py::object &dtype, std::vector<py::ssize_t> shape) {
    return py::array(numpy.attr("broadcast_to")(numpy.attr("asarray")(values, "dtype"_a = dtype), py::cast(shape)));
  };
  const auto velocitiesArray = py::array_t<double>(
      broadcast(velocities.is_none() ? py::object(py::float_(0.)) : velocities, py::dtype::of<double>(), {n, 3}));
  const auto typeIdsArray = py::array_t<unsigned long>(
      broadcast(typeIds.is_none() ? py::object(py::int_(0)) : typeIds, py::dtype::of<unsigned long>(), {n}));

  unsigned long nextId = 0;
  for (auto particle = _autoPas->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    nextId = std::max(nextId, particle->getID() + 1);
  }
  std::optional<py::array_t<unsigned long>> idsArray{};
  if (not ids.is_none()) {
    idsArray = py::array_t<unsigned long>(broadcast(ids, py::dtype::of<unsigned long>(), {n}));
  }

  const auto r = positions.unchecked<2>();
  const auto v = velocitiesArray.unchecked<2>();
  const auto t = typeIdsArray.unchecked<1>();

  // Validate everything first, to add either all or none.
  std::vector<ParticleType> particles{};
  particles.reserve(n);
  for (py::ssize_t i = 0; i < n; ++i) {
    particles.emplace_back(std::array<double, 3>{r(i, 0), r(i, 1), r(i, 2)},
                           std::array<double, 3>{v(i, 0), v(i, 1), v(i, 2)},
                           idsArray ? idsArray->at(i) : nextId + static_cast<unsigned long>(i), t(i));
    validateParticle(particles.back());
  }
  for (const auto &particle : particles) {
    _autoPas->addParticle(particle);
  }
}

py::dict PySimulation::particles() const {
  ensureAccessible();
  std::vector<ParticleType> particles{};
  particles.reserve(_autoPas->getNumberOfParticles(autopas::IteratorBehavior::owned));
  for (auto particle = _autoPas->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    particles.push_back(*particle);
  }
  return particlesToArrays(std::move(particles));
}

py::object PySimulation::particleAttribute(const char *key) const { return particles()[key]; }

py::dict PySimulation::history() const {
  std::vector<uint8_t> tuning(_historyTuning.begin(), _historyTuning.end());
  return py::dict("iteration"_a = py::array(py::cast(_historyIteration)),
                  "potential_energy"_a = py::array(py::cast(_historyPotentialEnergy)),
                  "kinetic_energy"_a = py::array(py::cast(_historyKineticEnergy)),
                  "temperature"_a = py::array(py::cast(_historyTemperature)),
                  "virial"_a = py::array(py::cast(_historyVirial)),
                  "tuning"_a = py::array(py::cast(tuning)).attr("astype")("bool"));
}

py::dict PySimulation::currentConfiguration() const {
  ensureAccessible();
  py::dict result{};
  for (const auto &[type, configuration] : _autoPas->getCurrentConfigs()) {
    result[py::str(autopas::InteractionTypeOption{type}.to_string())] = configuration.get().toString();
  }
  return result;
}

py::list PySimulation::functors() const {
  py::list result{};
  for (const auto &object : _functorObjects) {
    result.append(object);
  }
  return result;
}

std::string PySimulation::repr() const {
  std::ostringstream os;
  os << "AutoPas(particles=" << _autoPas->getNumberOfParticles(autopas::IteratorBehavior::owned)
     << ", iteration=" << _iteration << ", functors=[";
  for (size_t i = 0; i < _functors.size(); ++i) {
    os << (i == 0 ? "" : ", ") << _functors[i]->name();
  }
  os << "])";
  return os.str();
}

void bindSimulation(py::module_ &module) {
  py::class_<PySimulation>(module, "AutoPas", R"doc(
A particle simulation on top of an AutoPas container.

Args:
    config: Dictionary with the keys of md-flexible's YAML files, or the path to such a file.
        In a dictionary, "Objects" may also be a list of generators.
    functors: A functor or a list of functors. If omitted, the keys functor and functor-3b of the
        configuration decide. At most one native functor per interaction type can be used.
    particle_properties: A ParticlePropertiesLibrary, which replaces the key Sites of the configuration.
    objects: Further generators of particles.
    history_interval: Every how many iterations thermodynamic values are recorded. 0 disables this.

Call the object to start the simulation or to resume it from its current state.
)doc")
      .def(py::init<py::object, py::object, std::shared_ptr<ParticlePropertiesLibraryType>, py::object, size_t>(),
           py::arg("config"), py::kw_only(), py::arg("functors") = py::none(),
           py::arg("particle_properties") = nullptr, py::arg("objects") = py::none(), py::arg("history_interval") = 1)
      .def(
          "__call__",
          [](py::object self, std::optional<size_t> iterations, const py::object &callback, size_t callbackInterval) {
            self.cast<PySimulation &>().run(self, iterations, callback, callbackInterval);
            return self;
          },
          py::arg("iterations") = py::none(), py::kw_only(), py::arg("callback") = py::none(),
          py::arg("callback_interval") = 1,
          R"doc(
Starts the simulation or resumes it from its current state.

Args:
    iterations: The number of further iterations. If omitted, runs until the end defined by the
        configuration, which does nothing if that was reached already.
    callback: Called with this object in between iterations. Return False to stop early.
    callback_interval: Every how many iterations the callback is called.

Returns:
    This object.
)doc")
      .def("add_particles", &PySimulation::addParticles, py::arg("positions"), py::kw_only(),
           py::arg("velocities") = py::none(), py::arg("type_ids") = py::none(), py::arg("ids") = py::none(),
           "Adds particles. All of them have to be inside the domain and of a known type.")
      .def("particles", &PySimulation::particles, "Copies of all particle attributes, as dictionary of arrays.")
      .def_property_readonly(
          "positions", [](const PySimulation &s) { return s.particleAttribute("positions"); }, "Sorted by ID.")
      .def_property_readonly(
          "velocities", [](const PySimulation &s) { return s.particleAttribute("velocities"); }, "Sorted by ID.")
      .def_property_readonly(
          "forces", [](const PySimulation &s) { return s.particleAttribute("forces"); }, "Sorted by ID.")
      .def_property_readonly(
          "ids", [](const PySimulation &s) { return s.particleAttribute("ids"); }, "Sorted by ID.")
      .def_property_readonly(
          "type_ids", [](const PySimulation &s) { return s.particleAttribute("type_ids"); }, "Sorted by ID.")
      .def_property_readonly("history", &PySimulation::history, "Thermodynamic values per recorded iteration.")
      .def_property_readonly("current_configuration", &PySimulation::currentConfiguration,
                             "The algorithm configurations AutoPas currently uses.")
      .def_property_readonly("iteration", &PySimulation::iteration, "Number of completed iterations.")
      .def_property_readonly("num_tuning_phases_completed", &PySimulation::numTuningPhasesCompleted)
      .def_property_readonly("num_particles", &PySimulation::numParticles, "Number of owned particles.")
      .def_property_readonly("potential_energy", &PySimulation::potentialEnergy,
                             "Of the latest iteration. NaN if a functor does not provide it.")
      .def_property_readonly("virial", &PySimulation::virial,
                             "Of the latest iteration. NaN if a functor does not provide it.")
      .def_property_readonly("kinetic_energy", &PySimulation::kineticEnergy)
      .def_property_readonly("temperature", &PySimulation::temperature, "In units where k_B = 1.")
      .def_property_readonly("functors", &PySimulation::functors)
      .def_property_readonly("particle_properties", &PySimulation::particleProperties)
      .def_property_readonly(
          "box_min", [](const PySimulation &s) { return s.config().boxMin.value; }, "Lower corner of the domain.")
      .def_property_readonly(
          "box_max", [](const PySimulation &s) { return s.config().boxMax.value; }, "Upper corner of the domain.")
      .def_property_readonly(
          "cutoff", [](const PySimulation &s) { return s.config().cutoff.value; })
      .def_property_readonly(
          "delta_t", [](const PySimulation &s) { return s.config().deltaT.value; })
      .def(
          "config_summary", [](const PySimulation &s) { return s.config().to_string(); },
          "The full configuration in md-flexible's notation.")
      .def("__repr__", &PySimulation::repr);
}

}  // namespace autopasPython
