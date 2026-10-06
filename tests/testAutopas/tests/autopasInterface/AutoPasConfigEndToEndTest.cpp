/**
 * @file AutoPasConfigEndToEndTest.cpp
 * @author S. Newcome
 * @date 24.09.2026
 */

#include "AutoPasConfigEndToEndTest.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iterator>
#include <random>
#include <sstream>

#include "autopas/AutoPasDecl.h"
#include "autopas/LogicHandler.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/NumberSetFinite.h"
#include "autopas/utils/generators/ClosestPackedGenerator.h"
#include "autopas/utils/generators/PseudoContainer.h"
#include "autopas/utils/generators/SphereGenerator.h"
#include "autopas/utils/generators/UniformGenerator.h"
#include "autopas/utils/inBox.h"
#include "testingHelpers/GenerateValidConfigurations.h"
#include "testingHelpers/commonTypedefs.h"

extern template class autopas::AutoPas<Molecule>;
extern template bool autopas::AutoPas<Molecule>::computeInteractions(LJFunctorGlobals *);
extern template bool autopas::AutoPas<Molecule>::computeInteractions(ATMFunctorGlobals *);

using AutoPasConfigEndToEndTestHelper::HaloMode;
using AutoPasConfigEndToEndTestHelper::ParticleChangeMode;
using AutoPasConfigEndToEndTestHelper::Scenario;

namespace {

/**
 * Wraps a position that left the box back into it, as if all boundaries were periodic.
 * @param pos
 * @param boxMin
 * @param boxMax
 * @return position inside [boxMin, boxMax).
 */
std::array<double, 3> wrapIntoBox(std::array<double, 3> pos, const std::array<double, 3> &boxMin,
                                  const std::array<double, 3> &boxMax) {
  for (size_t d = 0; d < 3; ++d) {
    const double boxLength = boxMax[d] - boxMin[d];
    if (pos[d] < boxMin[d]) {
      pos[d] += boxLength;
    } else if (pos[d] >= boxMax[d]) {
      pos[d] = pos[d] - boxLength;
    }
  }
  return pos;
}

/**
 * Deterministic pseudo random numbers in [0, 1) that only depend on the particle id and the timestep.
 * Because they do not depend on the order in which particles are visited, all configurations and the reference make
 * the exact same random decisions. std::seed_seq is used as a hash function.
 * @param id Particle id.
 * @param timestep
 * @return Four independent numbers: the perturbation in x, y, and z, and the number for the deletion decision.
 */
std::array<double, 4> randomNumbers(const size_t id, const size_t timestep) {
  std::seed_seq seedSequence{id, timestep};
  std::array<std::uint32_t, 4> values{};
  seedSequence.generate(values.begin(), values.end());
  std::array<double, 4> numbers{};
  for (size_t i = 0; i < numbers.size(); ++i) {
    // Dividing by 2^32 maps the 32-bit values to [0, 1).
    numbers[i] = values[i] / 4294967296.;
  }
  return numbers;
}

/**
 * Generates the particles of a scenario on an FCC lattice in [boxMin - haloWidth, boxMax + haloWidth). The lattice
 * spacing is chosen such that the number of particles matches the density of the scenario, see
 * autopas::generators::ClosestPackedGenerator::optimizeSpacingForDensity(). The dense sphere is generated with
 * autopas::generators::SphereGenerator.
 * Ids are sequential starting at 0.
 * @param scenario
 * @param haloWidth Width of the region around the box that is also filled.
 * @return
 */
std::vector<Molecule> generateInitialParticles(const Scenario &scenario, const double haloWidth) {
  using namespace autopas::utils::ArrayMath::literals;
  using autopas::generators::ClosestPackedGenerator::LatticeStructure;
  constexpr auto boxMin = AutoPasConfigEndToEndTest::_boxMin;

  const auto backgroundMin = boxMin - haloWidth;
  const auto backgroundMax = scenario.boxMax + haloWidth;
  const double backgroundSpacing = autopas::generators::ClosestPackedGenerator::optimizeSpacingForDensity(
      backgroundMin, backgroundMax, scenario.density, LatticeStructure::FCC, /*centered*/ false);
  std::vector<Molecule> fccMolecules;
  autopas::generators::PseudoContainer fccMoleculesWrapper(fccMolecules);
  autopas::generators::ClosestPackedGenerator::fillWithParticles(fccMoleculesWrapper, backgroundMin, backgroundMax,
                                                                 Molecule{}, backgroundSpacing, LatticeStructure::FCC,
                                                                 /*centered*/ false);
  if (not scenario.denseSphere) {
    return fccMolecules;
  }

  const auto center = (boxMin + scenario.boxMax) * 0.5;
  constexpr double radius = AutoPasConfigEndToEndTest::_sphereRadius;
  // The sphere is a simple cubic grid, which has one particle per spacing^3.
  const double sphereSpacing = std::cbrt(1. / AutoPasConfigEndToEndTest::_sphereDensity);

  const auto distanceToCenter = [&](const Molecule &p) { return autopas::utils::ArrayMath::L2Norm(p.getR() - center); };
  std::vector<Molecule> allMolecules;
  // Remove the background around the sphere, with a margin to avoid particles overlapping.
  std::ranges::copy_if(fccMolecules, std::back_inserter(allMolecules),
                       [&](const Molecule &p) { return distanceToCenter(p) >= radius + sphereSpacing / 2.; });

  const auto radiusInParticles = static_cast<int>(std::round(radius / sphereSpacing)) - 1;
  autopas::generators::SphereGenerator::iteratePositions(center, radiusInParticles, sphereSpacing,
                                                         [&](const std::array<double, 3> &pos) {
                                                           Molecule particle{};
                                                           particle.setR(pos);
                                                           allMolecules.push_back(particle);
                                                         });
  for (size_t id = 0; id < allMolecules.size(); ++id) {
    allMolecules[id].setID(id);
  }
  return allMolecules;
}

/**
 * Generates periodic images of all owned particles within the interaction length of the boundaries.
 * Images keep the id of the particle they are copied from.
 * @param autoPas
 * @param interactionLength
 * @return The images, which all lie outside the box.
 */
std::vector<Molecule> generatePeriodicImages(autopas::AutoPas<Molecule> &autoPas, const double interactionLength) {
  const auto &boxMin = autoPas.getBoxMin();
  const auto &boxMax = autoPas.getBoxMax();

  std::vector<Molecule> images;
  for (auto iter = autoPas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    const auto &pos = iter->getR();
    // Check every face, edge, and corner of the box, given as a direction from the box center. The particle has an
    // image on the opposite side if it is within the interaction length of all boundaries of this direction.
    for (const int x : {-1, 0, 1}) {
      for (const int y : {-1, 0, 1}) {
        for (const int z : {-1, 0, 1}) {
          if (x == 0 and y == 0 and z == 0) {
            continue;
          }
          const std::array<int, 3> direction{x, y, z};
          bool nearBoundaries = true;
          std::array<double, 3> shift{};
          for (size_t d = 0; d < 3; ++d) {
            const double boxLength = boxMax[d] - boxMin[d];
            if (direction[d] == -1) {
              nearBoundaries = nearBoundaries and pos[d] < boxMin[d] + interactionLength;
              shift[d] = boxLength;
            } else if (direction[d] == 1) {
              nearBoundaries = nearBoundaries and pos[d] >= boxMax[d] - interactionLength;
              shift[d] = -boxLength;
            }
          }
          if (nearBoundaries) {
            auto image = *iter;
            image.addR(shift);
            images.push_back(image);
          }
        }
      }
    }
  }
  return images;
}

/**
 * Human-readable name of a halo mode for test names.
 * @param haloMode
 * @return
 */
std::string toString(HaloMode haloMode) {
  switch (haloMode) {
    case HaloMode::noHalo:
      return "NoHalo";
    case HaloMode::subdomainHalo:
      return "SubdomainHalo";
    case HaloMode::periodicHalo:
      return "PeriodicHalo";
  }
  return "UnknownHaloMode";
}

/**
 * Human-readable name of a particle change mode for test names.
 * @param particleChangeMode
 * @return
 */
std::string toString(ParticleChangeMode particleChangeMode) {
  switch (particleChangeMode) {
    case ParticleChangeMode::noChanges:
      return "NoChanges";
    case ParticleChangeMode::deleteParticles:
      return "Delete";
    case ParticleChangeMode::addParticles:
      return "Add";
    case ParticleChangeMode::deleteAndAdd:
      return "DeleteAndAdd";
  }
  return "UnknownParticleChangeMode";
}

}  // namespace

std::string Scenario::toString() const {
  std::stringstream resStream;
  resStream << "Box" << boxMax[0] << "x" << boxMax[1] << "x" << boxMax[2] << "_Rho" << density
            << (denseSphere ? "_Sphere" : "");
  return resStream.str();
}

std::optional<std::vector<AutoPasConfigEndToEndTest::StepResult>> AutoPasConfigEndToEndTest::simulate(
    const autopas::Configuration &config, const mykey_t &key, const bool useSorting, std::string &rejectionMessage) {
  const auto interactionType = std::get<autopas::InteractionTypeOption::Value>(key);
  try {
    if (interactionType == autopas::InteractionTypeOption::pairwise) {
      LJFunctorGlobals functor{_cutoff};
      functor.setParticleProperties(_epsilon * 24, _sigma * _sigma);
      return simulateImpl(functor, config, key, useSorting);
    } else if (interactionType == autopas::InteractionTypeOption::triwise) {
      ATMFunctorGlobals functor{_cutoff};
      functor.setParticleProperties(_nu);
      return simulateImpl(functor, config, key, useSorting);
    }
  } catch (const autopas::utils::ExceptionHandler::AutoPasException &autoPasException) {
    const std::string message = autoPasException.what();
    // We expect some configurations to be inapplicable to the domain. We handle these silently and check them in our
    // list of expected inapplicable configs.
    if (message.find("Rejected the only configuration in the search space!") != std::string::npos) {
      rejectionMessage = message;
      return std::nullopt;
    }
    // else, something actually bad happened
    throw;
  }
  return std::vector<StepResult>{};
}

template <class Functor_T>
std::vector<AutoPasConfigEndToEndTest::StepResult> AutoPasConfigEndToEndTest::simulateImpl(
    Functor_T &functor, const autopas::Configuration &config, const mykey_t &key, bool useSorting) {
  using namespace autopas::utils::ArrayMath::literals;
  const auto &[scenario, haloMode, particleChangeMode, interactionType] = key;

  autopas::AutoPas<Molecule> autoPas;
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::warn);
  autoPas.setBoxMin(_boxMin);
  autoPas.setBoxMax(scenario.boxMax);
  autoPas.setCutoff(_cutoff);
  autoPas.setVerletSkin(_skin);
  autoPas.setVerletRebuildFrequency(_rebuildFrequency);
  autoPas.setAllowedInteractionTypeOptions({interactionType});
  autoPas.setAllowedContainers({config.container});
  autoPas.setAllowedTraversals({config.traversal}, interactionType);
  autoPas.setAllowedLoadEstimators({config.loadEstimator});
  autoPas.setAllowedDataLayouts({config.dataLayout}, interactionType);
  autoPas.setAllowedNewton3Options({config.newton3}, interactionType);
  autoPas.setAllowedCellSizeFactors(autopas::NumberSetFinite<double>(std::set<double>{config.cellSizeFactor}));
  autoPas.setAllowedVecPatterns({config.vecPattern}, interactionType);
  const size_t sortingThreshold = useSorting ? 5 : std::numeric_limits<size_t>::max();
  autoPas.setAoSSortingThreshold(sortingThreshold);
  autoPas.setSoASortingThreshold(sortingThreshold);
  autoPas.init();

  const auto &boxMin = autoPas.getBoxMin();
  const auto &boxMax = autoPas.getBoxMax();
  constexpr double interactionLength = _cutoff + _skin;
  const auto haloBoxMin = boxMin - interactionLength;
  const auto haloBoxMax = boxMax + interactionLength;

  // Adds a particle as owned if it is in the box, as halo if it is in the halo region, or drops it otherwise.
  auto addByPosition = [&](const Molecule &p) {
    if (autopas::utils::inBox(p.getR(), boxMin, boxMax)) {
      autoPas.addParticle(p);
    } else if (autopas::utils::inBox(p.getR(), haloBoxMin, haloBoxMax)) {
      autoPas.addHaloParticle(p);
    }
  };

  // Initial particles. Only the subdomain mode fills the halo region.
  const auto initialParticles =
      generateInitialParticles(scenario, haloMode == HaloMode::subdomainHalo ? interactionLength : 0.);
  for (const auto &p : initialParticles) {
    addByPosition(p);
  }
  const auto numInitialOwned = autoPas.getNumberOfParticles(autopas::IteratorBehavior::owned);
  const auto numParticlesToAdd = static_cast<size_t>(
      std::round(static_cast<double>(numInitialOwned) * params[interactionType].additionPercentage / 100.));
  size_t nextId = initialParticles.size();
  const double deletionFraction = params[interactionType].deletionPercentage / 100.;

  // The perturbation of a particle: maps the first three random numbers from [0, 1) to [-width/2, width/2).
  auto perturbation = [](const std::array<double, 4> &numbers) -> std::array<double, 3> {
    return {(numbers[0] - 0.5) * _perturbationWidth, (numbers[1] - 0.5) * _perturbationWidth,
            (numbers[2] - 0.5) * _perturbationWidth};
  };
  const bool deleteParticles = particleChangeMode & ParticleChangeMode::deleteParticles;

  std::vector<StepResult> results;
  for (size_t timestep = 0; timestep < _numTimesteps; ++timestep) {
    // 1. Perturb all owned particles and reset the forces.
    for (auto iter = autoPas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
      iter->addR(perturbation(randomNumbers(iter->getID(), timestep)));
      iter->setF({0., 0., 0.});
    }

    // 2. Delete owned particles.
    if (deleteParticles) {
      for (auto iter = autoPas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
        if (randomNumbers(iter->getID(), timestep)[3] < deletionFraction) {
          autoPas.deleteParticle(iter);
        }
      }
    }

    // 3. Update the container.
    // This deletes the halo particles. In a real multi-rank simulation, these would communicated anew from the
    // neighboring MPI ranks, but as we merely mimic a multi-rank simulation in the subdomain halo mode, we just save
    // the halo particles and later add them back in. The saved copies are then perturbed and potentially deleted like
    // owned particles, as the neighboring ranks would do with their particles.
    std::vector<Molecule> haloParticles;
    if (haloMode == HaloMode::subdomainHalo) {
      for (auto iter = autoPas.begin(autopas::IteratorBehavior::halo); iter.isValid(); ++iter) {
        const auto numbers = randomNumbers(iter->getID(), timestep);
        if (deleteParticles and numbers[3] < deletionFraction) {
          continue;
        }
        auto haloParticle = *iter;
        haloParticle.addR(perturbation(numbers));
        haloParticle.setF({0., 0., 0.});
        haloParticles.push_back(haloParticle);
      }
    }
    const auto leavingParticles = autoPas.updateContainer();

    // 4. Exchange particles with the surroundings.
    switch (haloMode) {
      case HaloMode::noHalo:
        // Leaving particles are deleted.
        break;
      case HaloMode::subdomainHalo:
        // Leaving particles become halo particles, halo particles that entered the box become owned.
        for (const auto &p : leavingParticles) {
          addByPosition(p);
        }
        for (const auto &p : haloParticles) {
          addByPosition(p);
        }
        break;
      case HaloMode::periodicHalo:
        for (auto p : leavingParticles) {
          p.setR(wrapIntoBox(p.getR(), boxMin, boxMax));
          autoPas.addParticle(p);
        }
        for (const auto &image : generatePeriodicImages(autoPas, interactionLength)) {
          autoPas.addHaloParticle(image);
        }
        break;
    }

    // 5. Add new particles. The positions only depend on the timestep, so all configurations add the same particles.
    if (particleChangeMode & ParticleChangeMode::addParticles) {
      std::mt19937 generator(timestep);
      for (size_t i = 0; i < numParticlesToAdd; ++i) {
        const auto pos = autopas::generators::UniformGenerator::randomPosition(generator, boxMin, boxMax);
        autoPas.addParticle(Molecule(pos, {0., 0., 0.}, nextId++));
      }
    }

    // 6. Calculate forces.
    autoPas.computeInteractions(&functor);

    StepResult result{.forces = std::vector<std::array<double, 3>>(nextId), .numOwned = 0, .globals = {}};
    for (auto iter = autoPas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
      result.forces.at(iter->getID()) = iter->getF();
      ++result.numOwned;
    }
    result.globals = {.potentialEnergy = functor.getPotentialEnergy(), .virial = functor.getVirial()};
    results.push_back(std::move(result));
  }

#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
  // This test is intended to have at least one iteration of reusing neighbor lists previously built, so error if this
  // is not the case.
  EXPECT_GT(autoPas.getMeanRebuildFrequency(), 1.) << "The neighbor lists were rebuilt in every iteration.";
#endif

  return results;
}

const std::vector<AutoPasConfigEndToEndTest::StepResult> &AutoPasConfigEndToEndTest::generateReference(
    const mykey_t &key) {
  const auto interactionType = std::get<autopas::InteractionTypeOption::Value>(key);
  // Calculate reference forces. For the reference forces we switch off sorting.
  if (not _reference.contains(key)) {
    const auto referenceConfig =
        interactionType == autopas::InteractionTypeOption::pairwise
            ? autopas::Configuration(autopas::ContainerOption::linkedCells, 1.0, autopas::TraversalOption::lc_c08,
                                     autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                                     autopas::Newton3Option::enabled, autopas::InteractionTypeOption::pairwise,
                                     autopas::VectorizationPatternOption::NA)
            : autopas::Configuration(autopas::ContainerOption::linkedCells, 1.0, autopas::TraversalOption::lc_c01,
                                     autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                                     autopas::Newton3Option::disabled, autopas::InteractionTypeOption::triwise,
                                     autopas::VectorizationPatternOption::NA);
    std::string rejectionMessage;
    const auto reference = simulate(referenceConfig, key, false, rejectionMessage);
    if (not reference) {
      // An empty reference would let every comparison pass trivially.
      ADD_FAILURE() << "The reference configuration was rejected: " << rejectionMessage;
    }
    _reference[key] = reference.value_or(std::vector<StepResult>{});
  }
  return _reference[key];
}

bool AutoPasConfigEndToEndTest::isExpectedSkip(const autopas::Configuration &config, const Scenario &scenario) {
  const auto skippedConfigs = _expectedSkips.find(scenario);
  return skippedConfigs != _expectedSkips.end() and skippedConfigs->second.contains(config);
}

/**
 * This tests all valid configurations against a reference configuration.
 */
TEST_P(AutoPasConfigEndToEndTest, configTest) {
  const auto &[scenario, haloMode, particleChangeMode, interactionType] = GetParam();
  const mykey_t key{scenario, haloMode, particleChangeMode, interactionType};

  constexpr double tolerance = 1.0e-10;

  const auto &reference = generateReference(key);
  ASSERT_EQ(reference.size(), _numTimesteps) << "No valid reference.";

  for (const auto &config : generateAllValidConfigurations(interactionType)) {
    if (scenario.skipDirectSum and config.container == autopas::ContainerOption::directSum) {
      continue;
    }
    // Todo: Remove this when the AxilrodTeller functor implements SoA
    if (interactionType == autopas::InteractionTypeOption::triwise and
        config.dataLayout == autopas::DataLayoutOption::soa) {
      continue;
    }
    SCOPED_TRACE(config.toShortString(false));

    const auto startTime = std::chrono::steady_clock::now();
    std::string rejectionMessage;
    std::optional<std::vector<StepResult>> calculated;
    try {
      calculated = simulate(config, key, true, rejectionMessage);
    } catch (const std::exception &e) {
      // Report the exception as a failure of this configuration and continue with the remaining configurations.
      ADD_FAILURE() << "Exception thrown for configuration " << config.toShortString(false) << " in scenario "
                    << scenario.toString() << ":\n"
                    << e.what();
      continue;
    }
    const auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count();
    // Timings per configuration, used to balance the cost of tests.
    RecordProperty(config.toShortString(false, true), std::to_string(duration));

    const bool expectedSkip = isExpectedSkip(config, scenario);
    if (not calculated) {
      if (not expectedSkip) {
        ADD_FAILURE() << "Unexpectedly rejected configuration " << config.toShortString(false) << " in scenario "
                      << scenario.toString() << ":\n"
                      << rejectionMessage;
      }
      continue;
    }
    if (expectedSkip) {
      ADD_FAILURE() << "Configuration " << config.toShortString(false) << " is expected to be rejected in scenario "
                    << scenario.toString() << " but was applicable.";
    }

    ASSERT_EQ(calculated->size(), reference.size());
    for (size_t timestep = 0; timestep < calculated->size(); ++timestep) {
      const auto &[calculatedForces, calculatedNumOwned, calculatedGlobals] = (*calculated)[timestep];
      const auto &[referenceForces, referenceNumOwned, referenceGlobals] = reference[timestep];

      EXPECT_EQ(calculatedNumOwned, referenceNumOwned) << "Step: " << timestep;
      ASSERT_EQ(calculatedForces.size(), referenceForces.size()) << "Step: " << timestep;

      for (size_t id = 0; id < calculatedForces.size(); ++id) {
        // The tolerance is relative to the magnitude of the whole force vector, as single components can be arbitrarily
        // close to zero due to cancellation.
        const double forceTolerance = autopas::utils::ArrayMath::L2Norm(referenceForces[id]) * tolerance;
        for (unsigned int d = 0; d < 3; ++d) {
          if (not(std::abs(calculatedForces[id][d] - referenceForces[id][d]) <= forceTolerance)) {
            ADD_FAILURE() << "Force mismatch. Step: " << timestep << " Dim: " << d << " Particle id: " << id
                          << " Calculated: " << calculatedForces[id][d] << " Reference: " << referenceForces[id][d]
                          << " Tolerance: " << forceTolerance;
          }
        }
      }

      EXPECT_NEAR(calculatedGlobals.potentialEnergy, referenceGlobals.potentialEnergy,
                  std::abs(tolerance * referenceGlobals.potentialEnergy))
          << "Step: " << timestep;

      EXPECT_NEAR(calculatedGlobals.virial, referenceGlobals.virial, std::abs(tolerance * referenceGlobals.virial))
          << "Step: " << timestep;
    }
  }
}

/**
 * Lambda to generate a readable string out of the parameters of this test.
 */
static auto testName = [](const auto &info) {
  const auto &[scenario, haloMode, particleChangeMode, interactionType] = info.param;
  std::string res = interactionType.to_string() + "_" + scenario.toString() + "_" + toString(haloMode) + "_" +
                    toString(particleChangeMode);
  std::ranges::replace(res, '-', '_');
  std::ranges::replace(res, '.', '_');
  return res;
};

std::vector<AutoPasConfigEndToEndTestHelper::TestingTuple> AutoPasConfigEndToEndTest::getTestParams() {
  std::vector<AutoPasConfigEndToEndTestHelper::TestingTuple> testParams{};
  for (auto interactionType : autopas::InteractionTypeOption::getMostOptions()) {
    for (const auto &scenario : params[interactionType].scenarios) {
      for (HaloMode haloMode : {HaloMode::noHalo, HaloMode::subdomainHalo, HaloMode::periodicHalo}) {
        for (ParticleChangeMode particleChangeMode :
             {ParticleChangeMode::noChanges, ParticleChangeMode::deleteParticles, ParticleChangeMode::addParticles,
              ParticleChangeMode::deleteAndAdd}) {
          testParams.emplace_back(scenario, haloMode, particleChangeMode, interactionType);
        }
      }
    }
  }
  return testParams;
}

/**
 * Generates the cartesian product of the given boxes and densities as scenarios without dense spheres.
 * @param boxes
 * @param densities
 * @return
 */
static std::vector<Scenario> uniformScenarios(const std::vector<std::array<double, 3>> &boxes,
                                              const std::vector<double> &densities) {
  std::vector<Scenario> scenarios;
  for (const auto &boxMax : boxes) {
    for (const auto density : densities) {
      scenarios.push_back({.boxMax = boxMax,
                           .density = density,
                           .denseSphere = false /*denseSphere*/,
                           .skipDirectSum = false /*skipDirectSum*/});
    }
  }
  return scenarios;
}

std::unordered_map<autopas::InteractionTypeOption::Value, AutoPasConfigEndToEndTest::TraversalTestParams>
    AutoPasConfigEndToEndTest::params = [] {
      // Add small box scenarios at a range of densities.
      const std::vector<double> densities{0.05, 0.3, 1.};
      const std::vector<std::array<double, 3>> smallBoxes{{3., 3., 3.}, {10., 10., 10.}, {3., 50., 3.}};
      auto pairwiseScenarios = uniformScenarios(smallBoxes, densities);

      // Add additionally a couple large box scenarios, but only for lower densities.
      pairwiseScenarios.push_back({{50., 50., 50.}, 0.05, false, false});
      pairwiseScenarios.push_back({{50., 50., 50.}, 0.3, false, false});
      pairwiseScenarios.push_back({{50., 50., 50.}, 0.05, true /*denseSphere*/, false});

      // For triwise, start with the same small box scenarios
      auto triwiseScenarios = uniformScenarios(smallBoxes, densities);
      // Only add one large box scenario at low density, and skip direct sum.
      triwiseScenarios.push_back({{50., 50., 50.}, 0.05, false, true /*skipDirectSum*/});

      return std::unordered_map<autopas::InteractionTypeOption::Value, TraversalTestParams>{
          {autopas::InteractionTypeOption::pairwise,
           {.deletionPercentage = 10., .additionPercentage = 10., .scenarios = pairwiseScenarios}},
          {autopas::InteractionTypeOption::triwise,
           {.deletionPercentage = 10., .additionPercentage = 10., .scenarios = triwiseScenarios}},
      };
    }();

std::map<AutoPasConfigEndToEndTest::Scenario, std::set<autopas::Configuration>>
    AutoPasConfigEndToEndTest::_expectedSkips = [] {
      const auto allContainers = autopas::ContainerOption::getAllOptions();
      const auto allLoadEstimators = autopas::LoadEstimatorOption::getAllOptions();
      const auto allDataLayouts = autopas::DataLayoutOption::getAllOptions();
      const auto allNewton3Options = autopas::Newton3Option::getAllOptions();

      // The c04 traversals are not applicable with cell size factor 0.5.
      const std::set<autopas::TraversalOption> c04Traversals{autopas::TraversalOption::lc_c04,
                                                             autopas::TraversalOption::lc_c04_HCP,
                                                             autopas::TraversalOption::lc_c04_combined_SoA};
      const auto c04SmallCells =
          generateAllValidConfigurations(autopas::InteractionTypeOption::all, allContainers, c04Traversals,
                                         allLoadEstimators, allDataLayouts, allNewton3Options, {0.5});

      // lc_c04 is not applicable with any cell size factor if some dimension is shorter than two interaction lengths.
      const auto c04ThinBoxes = generateAllValidConfigurations(autopas::InteractionTypeOption::all, allContainers,
                                                               {autopas::TraversalOption::lc_c04}, allLoadEstimators,
                                                               allDataLayouts, allNewton3Options, {1.0, 1.5});

      std::map<Scenario, std::set<autopas::Configuration>> expectedSkips;
      for (const auto &[interactionType, testParams] : params) {
        for (const auto &scenario : testParams.scenarios) {
          auto &skips = expectedSkips[scenario];
          skips.insert(c04SmallCells.begin(), c04SmallCells.end());
          if (std::ranges::any_of(scenario.boxMax,
                                  [](double boxLength) { return boxLength < 2 * (_cutoff + _skin); })) {
            skips.insert(c04ThinBoxes.begin(), c04ThinBoxes.end());
          }
        }
      }
      return expectedSkips;
    }();

INSTANTIATE_TEST_SUITE_P(Generated, AutoPasConfigEndToEndTest,
                         ::testing::ValuesIn(AutoPasConfigEndToEndTest::getTestParams()), testName);
