/**
 * @file HGridTraversalComparison.cpp
 * @author atacann
 * @date 15.01.2024
 */
// Mostly copied from TraversalComparison.cpp, minor changes to use different molecule types and test HGrid only

#include "HGridTraversalComparison.h"

#include <string>

#include "autopas/tuning/selectors/ContainerSelector.h"
#include "autopas/tuning/selectors/TraversalSelector.h"
#include "autopas/utils/StaticCellSelector.h"
#include "autopas/utils/StringUtils.h"
#include "autopasTools/generators/UniformGenerator.h"
#include "molecularDynamicsLibrary/LJFunctorAVX.h"

std::shared_ptr<ParticlePropertiesLibrary<>> HGridTraversalComparison::_particlePropertiesLibrary =
    std::make_shared<ParticlePropertiesLibrary<>>(1.);

/**
 * Generates a random 3d shift with the given magnitude. The shift is uniformly distributed on a sphere with radius
 * magnitude.
 * @param magnitude
 * @param generator
 * @return shift vector
 */
std::array<double, 3> HGridTraversalComparison::randomShift(double magnitude, std::mt19937 &generator) {
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
void HGridTraversalComparison::executeShift(ContainerType &container, double magnitude, size_t numTotalParticles) {
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
void HGridTraversalComparison::markSomeParticlesAsDeleted(ContainerT &container, size_t numTotalParticles,
                                                          unsigned seed, autopas::InteractionTypeOption interactionT) {
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
std::tuple<std::vector<std::array<double, 3>>, HGridTraversalComparison::Globals>
HGridTraversalComparison::calculateForces(autopas::ContainerOption containerOption,
                                          autopas::TraversalOption traversalOption,
                                          autopas::DataLayoutOption dataLayoutOption,
                                          autopas::Newton3Option newton3Option, double cellSizeFactor, mykey_t key,
                                          bool useSorting) {
  auto [numParticles, numHaloParticles, boxMax, doSlightShift, particleDeletionPosition, _ /*globals*/,
        interactionType] = key;
  std::vector<std::array<double, 3>> calculatedForces;
  Globals calculatedGlobals;
  if (interactionType == autopas::InteractionTypeOption::pairwise) {
    mdLib::LJFunctor<Molecule, true /*applyShift*/, true /*useMixing*/, autopas::FunctorN3Modes::Both,
                     globals /*calculateGlobals*/, false /* countFLOPs */, true /* relevantForTuning */,
                     true /* scalingCutoff */>
        functor{_cutoff, *_particlePropertiesLibrary};
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
std::tuple<std::vector<std::array<double, 3>>, HGridTraversalComparison::Globals>
HGridTraversalComparison::calculateForcesImpl(Functor functor, autopas::ContainerOption containerOption,
                                              autopas::TraversalOption traversalOption,
                                              autopas::DataLayoutOption dataLayoutOption,
                                              autopas::Newton3Option newton3Option, double cellSizeFactor, mykey_t key,
                                              bool useSorting) {
  auto [numParticles, numHaloParticles, boxMax, doSlightShift, particleDeletionPosition, _ /*globals*/,
        interactionType] = key;

  // Construct container
  autopas::ContainerSelector<Molecule> selector{_boxMin, boxMax, _cutoff, {0.5, 1., 1.5}};
  constexpr double skin = _cutoff * 0.1;
  constexpr unsigned int rebuildFrequency = 1;
  selector.selectContainer(containerOption, autopas::ContainerSelectorInfo{cellSizeFactor, skin, rebuildFrequency, 32,
                                                                           autopas::LoadEstimatorOption::none});
  auto &container = selector.getCurrentContainer();
  std::array<size_t, 3> numParticlesArr = {numParticles / 3, numParticles / 3, numParticles - 2 * (numParticles / 3)};
  std::array<size_t, 3> numHaloParticlesArr = {numHaloParticles / 3, numHaloParticles / 3,
                                               numHaloParticles - 2 * (numHaloParticles / 3)};
  unsigned long moleculeId = 0;
  for (size_t i = 0; i < 3; ++i) {
    autopasTools::generators::UniformGenerator::fillWithParticles(
        container, Molecule({0., 0., 0.}, {0., 0., 0.}, moleculeId, i), container.getBoxMin(), container.getBoxMax(),
        numParticlesArr[i], i);
    moleculeId += numParticlesArr[i];
  }

  EXPECT_EQ(container.size(), numParticles) << "Wrong number of molecules inserted!";

  for (size_t i = 0; i < 3; ++i) {
    autopasTools::generators::UniformGenerator::fillWithHaloParticles(
        container, Molecule({0., 0., 0.}, {0., 0., 0.}, moleculeId, i), container.getCutoff(), numHaloParticlesArr[i],
        static_cast<unsigned int>(i));
    moleculeId += numHaloParticlesArr[i];
  }

  EXPECT_EQ(container.size(), numParticles + numHaloParticles) << "Wrong number of halo molecules inserted!";
  auto traversal =
      autopas::utils::withStaticCellType<Molecule>(container.getParticleCellTypeEnum(), [&](auto particleCellDummy) {
        auto traversalUniquePtr =
            autopas::TraversalSelector<decltype(particleCellDummy)>::template generateTraversal<decltype(functor)>(
                traversalOption, functor, container.getTraversalSelectorInfo(), dataLayoutOption, newton3Option);

        // set useSorting of the traversal if it can be cast to a CellTraversal and uses the CellFunctor
        if (auto *cellTraversalPtr =
                dynamic_cast<autopas::CellTraversal<decltype(particleCellDummy)> *>(traversalUniquePtr.get())) {
          cellTraversalPtr->setSortingThreshold(useSorting ? 5 : std::numeric_limits<size_t>::max());
        }

        return traversalUniquePtr;
      });

  if (not traversal->isApplicable()) {
    return {};
  }

  if (particleDeletionPosition & DeletionPosition::beforeLists) {
    markSomeParticlesAsDeleted(container, numParticles + numHaloParticles, 19, interactionType);
  }

  container.rebuildNeighborLists(traversal.get());

  if (doSlightShift) {
    // need to set max displacement because particles will get shifted
    container.setMaxDisplacement(skin / 2);
    executeShift(container, skin / 2, numParticles + numHaloParticles);
  }

  if (particleDeletionPosition & DeletionPosition::afterLists) {
    markSomeParticlesAsDeleted(container, numParticles + numHaloParticles, 99, interactionType);
  }

  functor.initTraversal();
  container.computeInteractions(traversal.get());

  functor.endTraversal(newton3Option);

  std::vector<std::array<double, 3>> forces(numParticles);
  for (auto it = container.begin(autopas::IteratorBehavior::owned); it.isValid(); ++it) {
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
void HGridTraversalComparison::generateReference(mykey_t key) {
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
    }
    _forcesReference[key] = calculatedForces;
    _globalValuesReference[key] = calculatedGlobals;
  }
}

/**
 * This tests a given configuration against a reference configuration.
 */
TEST_P(HGridTraversalComparison, traversalTest) {
  auto [containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles, numHaloParticles, boxMax,
        cellSizeFactor, doSlightShift, particleDeletionPosition, globals, interactionType] = GetParam();

  HGridTraversalComparison::mykey_t key{numParticles, numHaloParticles, boxMax, doSlightShift, particleDeletionPosition,
                                        globals,      interactionType};

  // empirically determined and set near the minimal possible value for the given number of particles
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  constexpr double rel_err_tolerance = 1.0e-10;
  constexpr double rel_err_tolerance_globals = 1.0e-10;
  Molecule::setParticlePropertiesLibrary(_particlePropertiesLibrary);

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
auto HGridTraversalComparison::getTestParams() {
  // initialize site Ids with mandatory mass parameter
  for (auto [siteTypeId, mass] : std::vector<std::pair<int, double>>{{0, 1.}, {1, 2.}, {2, 1.}}) {
    _particlePropertiesLibrary->addSiteType(siteTypeId, mass);
  }
  // initialize LJ parameters
  for (auto [siteTypeId, eps_sigma] :
       std::vector<std::pair<int, std::pair<double, double>>>{{0, {1., 1.5}}, {1, {1., 1.}}, {2, {1., 0.5}}}) {
    auto [epsilon, sigma] = eps_sigma;
    _particlePropertiesLibrary->addLJParametersToSite(siteTypeId, epsilon, sigma);
  }
  // initialize AT parameters
  for (auto [siteTypeId, nu] : std::vector<std::pair<int, double>>{{0, 1.}, {1, 1.}, {2, 1.}}) {
    _particlePropertiesLibrary->addATParametersToSite(siteTypeId, nu);
  }
  _particlePropertiesLibrary->calculateMixingCoefficients(true);
  Molecule::setParticlePropertiesLibrary(_particlePropertiesLibrary);
  Molecule::setCutoffMultiplier(1.0);
  std::vector<TestingTuple> testParams{};
  for (auto containerOption : {autopas::ContainerOption::hierarchicalGrid}) {
    for (auto interactionType : {autopas::InteractionTypeOption::pairwise}) {
      for (auto traversalOption :
           autopas::compatibleTraversals::allCompatibleTraversals(containerOption, interactionType)) {
        for (auto dataLayoutOption : autopas::DataLayoutOption::getAllOptions()) {
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
        }
      }
    }
  }
  return testParams;
}

std::unordered_map<autopas::InteractionTypeOption::Value, HGridTraversalComparison::TraversalTestParams>
    HGridTraversalComparison::params = {
        {autopas::InteractionTypeOption::pairwise,
         {30.,                               // deletionPercentage
          {100, 2000},                       // numParticles
          {200},                             // numHaloParticles
          {{3., 3.5, 4.}, {10., 10., 10.}},  // boxMax
          {0.5, 1., 2.}}},                   // cellSizeFactor
};

INSTANTIATE_TEST_SUITE_P(Generated, HGridTraversalComparison,
                         ::testing::ValuesIn(HGridTraversalComparison::getTestParams()), toString);
