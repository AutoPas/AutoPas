/**
 * @file AutoPasConfigEndToEndTest.h
 * @author S. Newcome
 * @date 24.09.2026
 */

#pragma once

#include <gtest/gtest.h>

#include <array>
#include <compare>
#include <map>
#include <optional>
#include <ostream>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "AutoPasTestBase.h"
#include "autopas/options/InteractionTypeOption.h"
#include "autopas/tuning/Configuration.h"

/**
 * Types that parameterize the AutoPasConfigEndToEndTest.
 */
namespace AutoPasConfigEndToEndTestHelper {

/**
 * Enum for what we do to the particle count every iteration.
 */
enum ParticleChangeMode {
  // We have chosen the values explicitly, s.t., this enum can be used using bit manipulation, i.e., deleteAndAdd
  // enables both bits for deleteParticles and addParticles.
  noChanges = 0b00,
  deleteParticles = 0b01,
  addParticles = 0b10,
  deleteAndAdd = 0b11
};

/**
 * How the surroundings of the simulation box are treated.
 */
enum HaloMode {
  /**
   * There are no halo particles. Particles that leave the box are deleted.
   */
  noHalo,
  /**
   * The box is a subdomain surrounded by other subdomains. The halo region is filled with independent particles of the
   * same distribution, which move independently of the owned particles. Particles that leave the box become halo
   * particles and halo particles that enter the box become owned particles.
   */
  subdomainHalo,
  /**
   * All boundaries are periodic. Particles that leave the box re-enter at the opposite side and the halo particles are
   * periodic images of owned particles.
   */
  periodicHalo,
};

/**
 * Describes the initial particle distribution of a test case.
 * Particles are placed on an FCC lattice with the given density in the box and, for HaloMode::subdomainHalo, in the
 * halo region around it. If denseSphere is set, a sphere of density AutoPasConfigEndToEndTest::_sphereDensity and
 * radius AutoPasConfigEndToEndTest::_sphereRadius is placed at the center of the box and the background lattice is
 * removed there.
 */
struct Scenario {
  /**
   * Upper corner of the simulation box. The lower corner is AutoPasConfigEndToEndTest::_boxMin.
   */
  std::array<double, 3> boxMax;
  /**
   * Number of particles per unit volume of the background FCC lattice.
   */
  double density;
  /**
   * If true, a dense sphere is placed at the center of the box.
   */
  bool denseSphere;
  /**
   * If true, no DirectSum configurations are tested for this scenario because they would be too expensive.
   */
  bool skipDirectSum;

  /**
   * Compares all members that define the simulated system.
   * @param rhs
   * @return
   */
  auto operator<=>(const Scenario &rhs) const {
    return std::tie(boxMax, density, denseSphere) <=> std::tie(rhs.boxMax, rhs.density, rhs.denseSphere);
  }

  /**
   * Equality of all members that define the simulated system.
   * @param rhs
   * @return
   */
  bool operator==(const Scenario &rhs) const { return (*this <=> rhs) == 0; }

  /**
   * String representation that can be used in a parameterized test name.
   * @return
   */
  [[nodiscard]] std::string toString() const;
};

/**
 * Stream operator for Scenario, used by gtest to print the parameter.
 * @param os
 * @param scenario
 * @return
 */
inline std::ostream &operator<<(std::ostream &os, const Scenario &scenario) { return os << scenario.toString(); }

/**
 * Parameters of one test: everything that defines a simulation, except the configuration.
 */
using TestingTuple = std::tuple<Scenario, HaloMode, ParticleChangeMode, autopas::InteractionTypeOption>;

}  // namespace AutoPasConfigEndToEndTestHelper

/**
 * The tests in this class compare the forces and global values of all configurations with a reference configuration.
 * Every configuration is simulated end-to-end through the AutoPas interface, i.e. with updateContainer(), particle
 * exchange, halo particles, and computeInteractions(), over _numTimesteps iterations. Positions are not updated by
 * forces but by deterministic random perturbations, so the tested configuration and the reference use exactly the same
 * particle positions.
 *
 * Each test covers one simulation (scenario, halo mode, particle change mode, interaction type) and loops over all
 * valid configurations, so the reference is only calculated once per test.
 */
class AutoPasConfigEndToEndTest : public AutoPasTestBase,
                                  public ::testing::WithParamInterface<AutoPasConfigEndToEndTestHelper::TestingTuple> {
 public:
  using ParticleChangeMode = AutoPasConfigEndToEndTestHelper::ParticleChangeMode;
  using HaloMode = AutoPasConfigEndToEndTestHelper::HaloMode;
  using Scenario = AutoPasConfigEndToEndTestHelper::Scenario;

  /**
   * Key that specifies a simulation, i.e. all test parameters. Tests with the same key share the same reference.
   */
  using mykey_t = std::tuple<Scenario,                                // scenario
                             HaloMode,                                // haloMode
                             ParticleChangeMode,                      // particleChangeMode
                             autopas::InteractionTypeOption::Value>;  // interaction type

  /**
   * Struct to hold global values
   */
  struct Globals {
    /**
     * Potential energy of all owned particles.
     */
    double potentialEnergy{};
    /**
     * Virial of all owned particles.
     */
    double virial{};
  };

  /**
   * Results of one timestep.
   */
  struct StepResult {
    /**
     * Forces of all owned particles, indexed by particle id. Entries of ids that are not owned are zero.
     */
    std::vector<std::array<double, 3>> forces;
    /**
     * Number of owned particles during the force calculation.
     */
    size_t numOwned{};
    /**
     * Global values of the force calculation.
     */
    Globals globals;
  };

  /**
   * Struct to hold parameters that might differ for different interaction types
   */
  struct TraversalTestParams {
    /**
     * Percentage of particles that are marked as deleted in every iteration, directly after the perturbation.
     */
    double deletionPercentage;
    /**
     * Number of particles, in percent of the initial number of owned particles, that are added in every iteration.
     */
    double additionPercentage;
    /**
     * Initial particle distributions that are tested for this interaction type.
     */
    std::vector<Scenario> scenarios;
  };

  /**
   * Generates the reference simulation results for a simulation that is specified by the given key, if they do not
   * exist yet. For the reference a linked cells algorithm and c08 (pairwise) or c01 (triwise) traversal without sorting
   * particles is used.
   * @param key The key that specifies the simulation.
   * @return The reference results of every step. Empty if the reference configuration was rejected.
   */
  static const std::vector<StepResult> &generateReference(const mykey_t &key);

  /**
   * Generates the parameters of all tests, i.e. the cartesian product of all scenarios, halo modes, and particle change
   * modes for each interaction type.
   * @return Vector of test parameters.
   */
  static std::vector<AutoPasConfigEndToEndTestHelper::TestingTuple> getTestParams();

  /**
   * Checks whether the given configuration is expected to be rejected in the given scenario.
   * @param config
   * @param scenario
   * @return True if the configuration is in _expectedSkips for the scenario.
   */
  static bool isExpectedSkip(const autopas::Configuration &config, const Scenario &scenario);

  /**
   * A map of the test scenarios run for each interaction type.
   */
  static std::unordered_map<autopas::InteractionTypeOption::Value, TraversalTestParams> params;

  /**
   * Lower corner of the simulation box.
   */
  static constexpr std::array<double, 3> _boxMin{0, 0, 0};
  /**
   * Cutoff radius of the interactions.
   */
  static constexpr double _cutoff{2.5};
  /**
   * Verlet skin.
   */
  static constexpr double _skin{_cutoff * 0.1};
  /**
   * Number of simulated timesteps per test.
   */
  static constexpr size_t _numTimesteps{3};
  /**
   * Verlet rebuild frequency. Together with _numTimesteps, this ensures that there is at least one iteration in which
   * the neighbor lists are reused.
   */
  static constexpr unsigned int _rebuildFrequency{2};
  /**
   * Particles are moved by a uniform random vector from [-_perturbationWidth/2, _perturbationWidth/2]^3 before every
   * force calculation. A single perturbation must be shorter than _skin/2, so neighbor lists stay valid for one step
   * without a rebuild.
   */
  static constexpr double _perturbationWidth{0.14};
  // The longest possible perturbation is the half diagonal of the perturbation cube: sqrt(3) * _perturbationWidth / 2.
  static_assert(3. * (_perturbationWidth / 2.) * (_perturbationWidth / 2.) < (_skin / 2.) * (_skin / 2.),
                "A single perturbation must be shorter than half the skin to keep neighbor lists valid for one step.");
  /**
   * Radius of the dense sphere of scenarios with Scenario::denseSphere.
   */
  static constexpr double _sphereRadius{5.};
  /**
   * Density of the dense sphere of scenarios with Scenario::denseSphere.
   */
  static constexpr double _sphereDensity{1.};

 protected:
  /**
   * Simulates _numTimesteps iterations with the given configuration through the AutoPas interface.
   *
   * Every iteration:
   *  1. All owned and halo particles are moved by a deterministic random perturbation.
   *  2. Optionally, particles are deleted.
   *  3. Optionally, new particles are added.
   *  4. The container is updated via AutoPas::updateContainer().
   *  5. Particles are exchanged according to the HaloMode, like a user would do with neighboring subdomains.
   *  6. Forces are calculated via AutoPas::computeInteractions().
   * Since no position depends on calculated forces, floating point differences can not accumulate over the steps.
   *
   * @param config The configuration to use. It is the only configuration in the search space.
   * @param key The key that specifies the simulation.
   * @param useSorting If the CellFunctor should apply sorting of particles.
   * @param rejectionMessage If AutoPas rejects the configuration, the exception message is written here.
   * @return Results of every step, or std::nullopt if AutoPas rejected the configuration.
   */
  static std::optional<std::vector<StepResult>> simulate(const autopas::Configuration &config, const mykey_t &key,
                                                         bool useSorting, std::string &rejectionMessage);

  /**
   * Implementation of simulate() for a specific functor type.
   * @tparam Functor_T
   * @param functor
   * @param config
   * @param key
   * @param useSorting
   * @return Results of every step.
   */
  template <class Functor_T>
  static std::vector<StepResult> simulateImpl(Functor_T &functor, const autopas::Configuration &config,
                                              const mykey_t &key, bool useSorting);

  /**
   * Lennard-Jones epsilon.
   */
  static constexpr double _epsilon{1.};
  /**
   * Lennard-Jones sigma.
   */
  static constexpr double _sigma{1.};
  /**
   * Axilrod-Teller-Muto nu.
   */
  static constexpr double _nu{1.};

  /**
   * Configurations that are expected to be rejected by AutoPas as not applicable, per scenario.
   * Every rejection that is not listed here is a test failure, as is a listed configuration that turns out to be
   * applicable.
   */
  static std::map<Scenario, std::set<autopas::Configuration>> _expectedSkips;

  /**
   * Cache of the reference results of all keys that have been simulated in this process.
   */
  static inline std::map<mykey_t, std::vector<StepResult>> _reference{};
};
