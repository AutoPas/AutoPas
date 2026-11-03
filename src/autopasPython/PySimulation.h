/**
 * @file PySimulation.h
 * @date 03.11.2026
 */
#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "FunctorHandle.h"
#include "src/ParallelVtkWriter.h"
#include "src/configuration/MDFlexConfig.h"
#include "src/domainDecomposition/RegularGridDecomposition.h"

namespace autopasPython {

/**
 * What Python calls AutoPas: An AutoPas container together with the state of a simulation.
 *
 * This reuses md-flexible's configuration, integrator, thermostat, and boundary handling, but implements its own main
 * loop instead of using md-flexible's Simulation. That class selects the functor from a compile time list, has no
 * notion of pausing the loop and handing back control, and is tailored to MPI. The Python interface runs on one rank.
 */
class PySimulation {
 public:
  /**
   * Constructor.
   * @param config Dictionary, or the path to a YAML file.
   * @param functors One functor or a sequence thereof. If none, the functors follow from the configuration.
   * @param particleProperties Replaces the properties given in the configuration.
   * @param objects Particle generators, in addition to the ones in the configuration.
   * @param historyInterval Every how many iterations thermodynamic values are recorded. 0 turns this off.
   */
  PySimulation(pybind11::object config, pybind11::object functors,
               std::shared_ptr<ParticlePropertiesLibraryType> particleProperties, pybind11::object objects,
               size_t historyInterval);

  ~PySimulation();

  /**
   * Starts the simulation or resumes it from the current state.
   * @param self The Python object of this instance. It is handed to the callback.
   * @param iterations The number of further iterations. If none, runs until the end defined by the configuration.
   * @param callback Called with self. Return False to stop early.
   * @param callbackInterval Every how many iterations the callback is called.
   */
  void run(const pybind11::object &self, std::optional<size_t> iterations, const pybind11::object &callback,
           size_t callbackInterval);

  /**
   * Adds particles to the container.
   * @param positions Array of shape (n, 3).
   * @param velocities Array of shape (n, 3). Defaults to zero.
   * @param typeIds Either one type for all, or an array of shape (n). Defaults to type 0.
   * @param ids Array of shape (n). Defaults to continuing after the largest id present.
   */
  void addParticles(const pybind11::array_t<double> &positions, const pybind11::object &velocities,
                    const pybind11::object &typeIds, const pybind11::object &ids);

  /**
   * All owned particles as dictionary of NumPy arrays, sorted by ID.
   * @return dictionary with the keys positions, velocities, forces, ids, and type_ids.
   */
  [[nodiscard]] pybind11::dict particles() const;

  /**
   * One of the arrays of particles().
   * @param key
   * @return Array.
   */
  [[nodiscard]] pybind11::object particleAttribute(const char *key) const;

  /**
   * The recorded thermodynamic values.
   * @return dictionary of NumPy arrays of equal length.
   */
  [[nodiscard]] pybind11::dict history() const;

  /**
   * The algorithm configurations AutoPas currently uses.
   * @return Dictionary from interaction type to the string representation of the configuration.
   */
  [[nodiscard]] pybind11::dict currentConfiguration() const;

  /**
   * Sum over the potential energies of all functors in the latest iteration.
   * @return Energy, NaN if at least one functor does not provide it or no iteration was done.
   */
  [[nodiscard]] double potentialEnergy() const;

  /**
   * Sum over the virials of all functors in the latest iteration.
   * @return Virial, NaN if at least one functor does not provide it or no iteration was done.
   */
  [[nodiscard]] double virial() const;

  /**
   * Current kinetic energy of the owned particles.
   * @return kinetic energy.
   */
  [[nodiscard]] double kineticEnergy() const;

  /**
   * Current temperature, in units where the Boltzmann constant is one.
   * @return Temperature, NaN if there are no particles.
   */
  [[nodiscard]] double temperature() const;

  /**
   * Number of owned particles.
   * @return number of owned particles.
   */
  [[nodiscard]] size_t numParticles() const;

  /**
   * Number of completed iterations.
   * @return iteration.
   */
  [[nodiscard]] size_t iteration() const { return _iteration; }

  /**
   * Number of completed tuning phases.
   * @return tuning phases.
   */
  [[nodiscard]] size_t numTuningPhasesCompleted() const { return _numTuningPhasesCompleted; }

  /**
   * The Python objects of the functors in the order in which they are evaluated.
   * @return list of functors.
   */
  [[nodiscard]] pybind11::list functors() const;

  /**
   * The particle properties library.
   * @return The library. Changes take effect immediately, after calculate_mixing_coefficients().
   */
  [[nodiscard]] std::shared_ptr<ParticlePropertiesLibraryType> particleProperties() const { return _particleProperties; }

  /**
   * The configuration of this simulation.
   * @return configuration
   */
  [[nodiscard]] const MDFlexConfig &config() const { return _config; }

  /**
   * Human readable summary.
   * @return summary.
   */
  [[nodiscard]] std::string repr() const;

 private:
  /**
   * Raises if the container is being worked on, which can only be observed from another Python thread or from a
   * functor. Between iterations, e.g. in a callback, everything is accessible.
   */
  void ensureAccessible() const;

  /**
   * Resolves the functors, restricts the search space to what they support, and binds them.
   * @param functors
   */
  void setupFunctors(pybind11::object functors);

  /**
   * Sanity checks a particle before it enters the container.
   * @param particle
   */
  void validateParticle(const ParticleType &particle) const;

  /**
   * One iteration of velocity Verlet. Same order of operations as Simulation::run() of md-flexible, minus MPI.
   * @note The GIL is not held when this is called.
   */
  void step();

  /**
   * Evaluates all functors.
   */
  void computeInteractions();

  /**
   * Stores the current thermodynamic values.
   */
  void record();

  MDFlexConfig _config;
  std::shared_ptr<RegularGridDecomposition> _domainDecomposition;
  std::shared_ptr<ParticlePropertiesLibraryType> _particleProperties;
  std::unique_ptr<AutoPasType> _autoPas;
  std::optional<ParallelVtkWriter> _vtkWriter{};

  /**
   * Python owns the functors, which also keeps the instances of Python subclasses alive. Sorted such that the functor
   * relevant for tuning comes first per interaction type, as the others reuse the configuration it selects.
   */
  std::vector<pybind11::object> _functorObjects{};
  std::vector<FunctorHandle *> _functors{};

  size_t _iteration{0};
  size_t _numTuningPhasesCompleted{0};
  size_t _numTuningIterations{0};
  bool _currentIterationIsTuningIteration{false};
  bool _previousIterationWasTuningIteration{false};
  bool _paused{false};

  bool _running{false};
  bool _stepping{false};
  bool _failed{false};

  size_t _historyInterval;
  std::vector<size_t> _historyIteration{};
  std::vector<double> _historyPotentialEnergy{};
  std::vector<double> _historyKineticEnergy{};
  std::vector<double> _historyTemperature{};
  std::vector<double> _historyVirial{};
  std::vector<bool> _historyTuning{};
};

/**
 * Binds the simulation, which is the AutoPas class in Python.
 * @param module
 */
void bindSimulation(pybind11::module_ &module);

}  // namespace autopasPython
