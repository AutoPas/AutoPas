/**
 * @file TraversalComparison.cpp
 * @author humig
 * @date 12.07.19
 */

#include "TraversalComparison.h"

#include <string>

#include "autopas/tuning/selectors/ContainerSelector.h"
#include "autopas/tuning/selectors/TraversalSelector.h"
#include "autopas/utils/StringUtils.h"
#include "autopasTools/generators/UniformGenerator.h"

/**
 * Generates a random 3d shift with the given magnitude. The shift is uniformly distributed on a sphere with radius
 * magnitude.
 * @param magnitude
 * @param generator
 * @return shift vector
 */
std::array<double, 3> randomShift(double magnitude, std::mt19937 &generator) {
  std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  double theta = 2 * M_PI * uniform01(generator);
  double phi = acos(1 - 2 * uniform01(generator));
  double x = sin(phi) * cos(theta) * magnitude;
  double y = sin(phi) * sin(theta) * magnitude;
  double z = cos(phi) * magnitude;
  return {x, y, z};
}

/**
 * Adds shifts with the given magnitude to all particles.
 * The shifts are generated in order of the particle id and with a fixed seed to ensure a reproducible behavior.
 * @tparam ContainerType
& * @param container
 * @param magnitude
 * @param numTotalParticles
 */
template <class ContainerType>
void TraversalComparison::executeShift(ContainerType &container, double magnitude, size_t numTotalParticles) {
  std::vector<std::array<double, 3>> shiftVectorByID(numTotalParticles);
  constexpr unsigned seed = 42;
  std::mt19937 generator(seed);
  for (auto &elem : shiftVectorByID) {
    elem = randomShift(magnitude, generator);
  }
  size_t numIteratedParticles = 0;
  for (auto iter = container.begin(autopas::IteratorBehavior::ownedOrHaloOrDummy); iter != container.end(); ++iter) {
    if (not iter->isDummy()) {
      iter->addR(shiftVectorByID[iter->getID()]);
    }
    ++numIteratedParticles;
  }
  EXPECT_EQ(numIteratedParticles, numTotalParticles);
}

/**
 * Marks a certain percentage of all particles as deleted.
 * The particles are selected randomly with the given seed.
 * @tparam ContainerT
 * @param container
 * @param numTotalParticles
 * @param seed
 * @param interactionT
 */
template <typename ContainerT>
void TraversalComparison::markSomeParticlesAsDeleted(ContainerT &container, size_t numTotalParticles, unsigned seed,
                                                     autopas::InteractionTypeOption interactionT) {
  // Here, we delete about deletionPercentage % of all particles.

  double deletionPercentage = params[interactionT].deletionPercentage;
  std::mt19937 generator(seed);
  std::uniform_real_distribution<double> uniform0_100(0.0, 100.0);
  // We create a vector of numMolecules + numHaloMolecules size.
  // The N-th entry of this vector indicates whether we delete the molecule with id N
  std::vector<bool> doDelete(numTotalParticles);
  // Generate the sequence of random numbers.
  std::generate(std::begin(doDelete), std::end(doDelete), [deletionPercentage, &uniform0_100, &generator]() {
    // Set to true if we are within deletionPercentage
    return uniform0_100(generator) < deletionPercentage;
  });
  for (auto &mol : container) {
    if (doDelete[mol.getID()]) {
      autopas::internal::markParticleAsDeleted(mol);
    }
  }
}

template <bool globals>
std::tuple<std::vector<std::array<double, 3>>, TraversalComparison::Globals> TraversalComparison::calculateForces(
    autopas::ContainerOption containerOption, autopas::TraversalOption traversalOption,
    autopas::DataLayoutOption dataLayoutOption, autopas::Newton3Option newton3Option, double cellSizeFactor,
    mykey_t key, bool useSorting) {
  auto [numParticles, numHaloParticles, boxMax, doSlightShift, particleDeletionPosition, _ /*globals*/,
        interactionType] = key;
  std::vector<std::array<double, 3>> calculatedForces;
  Globals calculatedGlobals;

  if (interactionType == autopas::InteractionTypeOption::pairwise) {
    mdLib::LJFunctor<Molecule, true /*applyShift*/, false /*useMixing*/, autopas::FunctorN3Modes::Both,
                     globals /*calculateGlobals*/>
        functor{_cutoff};
    functor.setParticleProperties(_eps * 24, _sig * _sig);
    std::tie(calculatedForces, calculatedGlobals) = calculateForcesImpl<decltype(functor), globals>(
        functor, containerOption, traversalOption, dataLayoutOption, newton3Option, cellSizeFactor, key, useSorting);
  } else if (interactionType == autopas::InteractionTypeOption::triwise) {
    mdLib::AxilrodTellerMutoFunctor<Molecule, false /*useMixing*/, autopas::FunctorN3Modes::Both,
                                    globals /*calculateGlobals*/>
        functor{_cutoff};
    functor.setParticleProperties(_nu);
    std::tie(calculatedForces, calculatedGlobals) = calculateForcesImpl<decltype(functor), globals>(
        functor, containerOption, traversalOption, dataLayoutOption, newton3Option, cellSizeFactor, key, useSorting);
  }

  return {calculatedForces, calculatedGlobals};
}
/**
 * Calculates the forces for a given configuration.
 * @param containerOption Specifies the container.
 * @param traversalOption Specifies the traversal.
 * @param dataLayoutOption Specifies the data layout.
 * @param newton3Option Specifies whether the newton3 optimization should be used or not.
 * @param numMolecules The number of molecules.
 * @param numHaloMolecules The number of halo molecules.
 * @param boxMax The maximum of the simulation box. The minimum is {0.,0.,0.}
 * @param cellSizeFactor The cell size factor.
 * @param doSlightShift Specifies whether to add random shifts of size skin/2 to all particles after the neighbor list
 * generation.
 * @param useSorting For traversals that use the CellFunctor: if the CellFunctor should apply sorting of particles
 * @return Tuple of forces for all particles, ordered by particle id, and global values.
 */
template <typename Functor, bool globals>
std::tuple<std::vector<std::array<double, 3>>, TraversalComparison::Globals> TraversalComparison::calculateForcesImpl(
    Functor functor, autopas::ContainerOption containerOption, autopas::TraversalOption traversalOption,
    autopas::DataLayoutOption dataLayoutOption, autopas::Newton3Option newton3Option, double cellSizeFactor,
    mykey_t key, bool useSorting) {
  auto [numParticles, numHaloParticles, boxMax, doSlightShift, particleDeletionPosition, _ /*globals*/,
        interactionType] = key;

  // Construct container
  constexpr double skin = _cutoff * 0.1;
  constexpr unsigned int rebuildFrequency = 1;
  const size_t sortingThreshold = useSorting ? 5 : std::numeric_limits<size_t>::max();
  const auto containerInfo = autopas::ContainerSelectorInfo{_boxMin,
                                                            boxMax,
                                                            _cutoff,
                                                            cellSizeFactor,
                                                            skin,
                                                            32,
                                                            sortingThreshold,
                                                            autopas::LoadEstimatorOption::none,
                                                            _orderCellsByMortonIndex,
                                                            _useOptimizedLJFunctor,
                                                            _useCompactAoS,
                                                            _reserveVLSizes,
                                                            _bucketSortParticles,
                                                            _sortVerletLists,
                                                            _sortingFrequency};
  auto container = autopas::ContainerSelector<Molecule>::generateContainer(containerOption, containerInfo);

  autopasTools::generators::UniformGenerator::fillWithParticles(*container, Molecule({0., 0., 0.}, {0., 0., 0.}, 0),
                                                                container->getBoxMin(), container->getBoxMax(),
                                                                numParticles);
  EXPECT_EQ(container->size(), numParticles) << "Wrong number of molecules inserted!";
  autopasTools::generators::UniformGenerator::fillWithHaloParticles(
      *container, Molecule({0., 0., 0.}, {0., 0., 0.}, numParticles /*initial ID*/), container->getCutoff(),
      numHaloParticles);
  EXPECT_EQ(container->size(), numParticles + numHaloParticles) << "Wrong number of halo molecules inserted!";
  const auto config =
      autopas::Configuration{containerOption,  cellSizeFactor, traversalOption, autopas::LoadEstimatorOption::none,
                             dataLayoutOption, newton3Option,  interactionType};

  auto traversal = autopas::TraversalSelector::generateTraversalFromConfig<Molecule, decltype(functor)>(
      config, functor, container->getTraversalSelectorInfo());

  if (not traversal) {
    return {};
  }

  if (particleDeletionPosition & DeletionPosition::beforeLists) {
    markSomeParticlesAsDeleted(*container, numParticles + numHaloParticles, 19, interactionType);
  }

  container->rebuildNeighborLists(traversal.get());

  if (doSlightShift) {
    executeShift(*container, skin / 2, numParticles + numHaloParticles);
  }

  if (particleDeletionPosition & DeletionPosition::afterLists) {
    markSomeParticlesAsDeleted(*container, numParticles + numHaloParticles, 99, interactionType);
  }

  functor.initTraversal();
  container->computeInteractions(traversal.get());

  functor.endTraversal(newton3Option);

  std::vector<std::array<double, 3>> forces(numParticles);
  for (auto it = container->begin(autopas::IteratorBehavior::owned); it.isValid(); ++it) {
    EXPECT_TRUE(it->isOwned());
    forces.at(it->getID()) = it->getF();
  }

  if (globals) {
    return {forces, {functor.getPotentialEnergy(), functor.getVirial()}};
  } else {
    return {forces, {0., 0.}};
  }
}

/**
 * Generates the reference for a simulation configuration that is specified by the given key.
 * For the reference a linked cells algorithm and c08 traversal without sorting particles is used.
 * @param key The key that specifies the simulation.
 */
template <bool globals>
void TraversalComparison::generateReference(mykey_t key) {
  auto [numParticles, numHaloParticles, boxMax, doSlightShift, particleDeletionPosition, _, interactionType] = key;
  std::vector<std::array<double, 3>> calculatedForces;
  Globals calculatedGlobals;
  // Calculate reference forces. For the reference forces we switch off sorting. For the forces that are calculated in
  // tests and compared against the reference, sorting is enabled.
  if (_forcesReference.count(key) == 0) {
    if (interactionType == autopas::InteractionTypeOption::pairwise) {
      std::tie(calculatedForces, calculatedGlobals) =
          calculateForces<globals>(autopas::ContainerOption::linkedCells, autopas::TraversalOption::lc_c08,
                                   autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled, 1., key, false);
    } else if (interactionType == autopas::InteractionTypeOption::triwise) {
      std::tie(calculatedForces, calculatedGlobals) =
          calculateForces<globals>(autopas::ContainerOption::linkedCells, autopas::TraversalOption::lc_c01,
                                   autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, 1., key, false);
    }
    _forcesReference[key] = calculatedForces;
    _globalValuesReference[key] = calculatedGlobals;
  }
}

/**
 * This tests a given configuration against a reference configuration.
 */
TEST_P(TraversalComparison, traversalTest) {
  auto [containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles, numHaloParticles, boxMax,
        cellSizeFactor, doSlightShift, particleDeletionPosition, globals, interactionType] = GetParam();
  // Todo: Remove this when the AxilrodTeller functor implements SoA
  if (interactionType == autopas::InteractionTypeOption::triwise and
      dataLayoutOption == autopas::DataLayoutOption::soa) {
    GTEST_SKIP_("SoAs are not yet implemented for triwise traversals/functors.");
  }

  TraversalComparison::mykey_t key{numParticles, numHaloParticles, boxMax, doSlightShift, particleDeletionPosition,
                                   globals,      interactionType};

  // empirically determined and set near the minimal possible value for the given number of particles
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  constexpr double rel_err_tolerance = 1.0e-10;
  constexpr double rel_err_tolerance_globals = 1.0e-10;

  std::vector<std::array<double, 3>> calculatedForces;
  Globals calculatedGlobals;
  if (globals) {
    std::tie(calculatedForces, calculatedGlobals) = calculateForces<true>(
        containerOption, traversalOption, dataLayoutOption, newton3Option, cellSizeFactor, key, true);
    generateReference<true>(key);
  } else {
    std::tie(calculatedForces, calculatedGlobals) = calculateForces<false>(
        containerOption, traversalOption, dataLayoutOption, newton3Option, cellSizeFactor, key, true);
    generateReference<false>(key);
  }

  if (calculatedForces.empty()) {
    GTEST_SKIP_("Not applicable!");
  }

  for (size_t i = 0; i < numParticles; ++i) {
    for (unsigned int d = 0; d < 3; ++d) {
      const double calculatedForce = calculatedForces[i][d];
      const double referenceForce = _forcesReference[key][i][d];
      EXPECT_NEAR(calculatedForce, referenceForce, std::fabs(calculatedForce * rel_err_tolerance))
          << "Dim: " << d << " Particle id: " << i;
    }
  }

  auto &globalValuesReferenceRef = _globalValuesReference[key];
  if (globals) {
    EXPECT_NE(calculatedGlobals.upot, 0);
    EXPECT_NEAR(calculatedGlobals.upot, globalValuesReferenceRef.upot,
                std::abs(rel_err_tolerance_globals * globalValuesReferenceRef.upot));

    EXPECT_NE(calculatedGlobals.virial, 0);
    EXPECT_NEAR(calculatedGlobals.virial, globalValuesReferenceRef.virial,
                std::abs(rel_err_tolerance_globals * globalValuesReferenceRef.virial));
  }
}

/**
 * Lambda to generate a readable string out of the parameters of this test.
 */
static auto toString = [](const auto &info) {
  auto [containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles, numHaloParticles, boxMax,
        cellSizeFactor, doSlightShift, particleDeletionPosition, globals, interactionType] = info.param;
  std::stringstream resStream;
  resStream << containerOption.to_string() << "_" << traversalOption.to_string()
            << (interactionType == autopas::InteractionTypeOption::triwise ? "_3B" : "") << "_"
            << dataLayoutOption.to_string() << "_"
            << (newton3Option == autopas::Newton3Option::enabled ? "_N3" : "_noN3") << "_NP" << numParticles << "_NH"
            << numHaloParticles << "_" << boxMax[0] << "_" << boxMax[1] << "_" << boxMax[2] << "_CSF_" << cellSizeFactor
            << "_" << (doSlightShift ? "withShift" : "noshift")
            << (particleDeletionPosition == DeletionPosition::never ? "_NoDeletions" : "")
            << (particleDeletionPosition & DeletionPosition::beforeLists ? "_DeletionsBeforeLists" : "")
            << (particleDeletionPosition & DeletionPosition::afterLists ? "_DeletionsAfterLists" : "")
            << (globals ? "_globals" : "_noGlobals");
  std::string res = resStream.str();
  std::replace(res.begin(), res.end(), '-', '_');
  std::replace(res.begin(), res.end(), '.', '_');
  return res;
};

/**
 * Function to generate all possible configurations.
 * @return
 */
auto TraversalComparison::getTestParams() {
  std::vector<TestingTuple> testParams{};
  auto containerOption = autopas::ContainerOption::verletListsSoA;
  // for (auto containerOption : autopas::ContainerOption::getAllOptions()) {
    for (auto interactionType : autopas::InteractionTypeOption::getMostOptions()) {
      for (auto traversalOption :
           autopas::compatibleTraversals::allCompatibleTraversals(containerOption, interactionType)) {
        // for (auto dataLayoutOption : autopas::DataLayoutOption::getAllOptions()) {
        auto dataLayoutOption = autopas::DataLayoutOption::soa;
          for (auto newton3Option : autopas::Newton3Option::getAllOptions()) {
            for (auto numParticles : params[interactionType].numParticles) {
              for (auto boxMax : params[interactionType].boxMax) {
                for (double cellSizeFactor : params[interactionType].cellSizeFactors) {
                  for (auto numHalo : params[interactionType].numHaloParticles) {
                    for (bool slightMove : {true, false}) {
                      for (bool globals : {true, /*false*/}) {
                        for (DeletionPosition particleDeletionPosition :
                             {DeletionPosition::never, /*DeletionPosition::beforeLists, DeletionPosition::afterLists,*/
                              DeletionPosition::beforeAndAfterLists}) {
                          testParams.emplace_back(containerOption, traversalOption, dataLayoutOption, newton3Option,
                                                  numParticles, numHalo, boxMax, cellSizeFactor, slightMove,
                                                  particleDeletionPosition, globals, interactionType);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        // }
      }
    }
  // }
  return testParams;
}

std::unordered_map<autopas::InteractionTypeOption::Value, TraversalComparison::TraversalTestParams>
    TraversalComparison::params = {{autopas::InteractionTypeOption::pairwise,
                                    {30.,                              // deletionPercentage
                                     {100, 2000},                      // numParticles
                                     {200},                            // numHaloParticles
                                     {{3., 3., 3.}, {10., 10., 10.}},  // boxMax
                                     {0.5, 1., 2.}}},                  // cellSizeFactor
                                   {autopas::InteractionTypeOption::triwise,
                                    {10.,                           // deletionPercentage
                                     {100, 400},                    // numParticles
                                     {200},                         // numHaloParticles
                                     {{3., 3., 3.}, {6., 6., 6.}},  // boxMax
                                     {0.5, 1., 1.5}}}};             // cellSizeFactor

INSTANTIATE_TEST_SUITE_P(Generated, TraversalComparison, ::testing::ValuesIn(TraversalComparison::getTestParams()),
                         toString);
