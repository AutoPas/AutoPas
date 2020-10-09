/**
 * @file TraversalComparison.cpp
 * @author humig
 * @date 12.07.19
 */

#include "TraversalComparison.h"

#include <string>

#include "autopas/selectors/ContainerSelector.h"
#include "autopas/selectors/TraversalSelector.h"
#include "autopas/utils/StringUtils.h"
#include "autopasTools/generators/RandomGenerator.h"

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
 * @tparam ContainerPtrType
 * @param containerPtr
 * @param magnitude
 * @param totalNumParticles
 */
template <class ContainerPtrType>
void TraversalComparison::executeShift(ContainerPtrType containerPtr, double magnitude, size_t totalNumParticles) {
  std::vector<std::array<double, 3>> shiftVectorByID(totalNumParticles);
  constexpr unsigned seed = 42;
  std::mt19937 generator(seed);
  for (auto &elem : shiftVectorByID) {
    elem = randomShift(magnitude, generator);
  }
  size_t numIteratedParticles = 0;
  for (auto iter = containerPtr->begin(autopas::IteratorBehavior::haloOwnedAndDummy); iter != containerPtr->end();
       ++iter) {
    if (not iter->isDummy()) {
      iter->addR(shiftVectorByID[iter->getID()]);
    }
    ++numIteratedParticles;
  }
  EXPECT_EQ(numIteratedParticles, totalNumParticles);
}

template <typename ContainerT>
void markSomeParticlesAsDeleted(ContainerT &container, size_t numTotalParticles, unsigned seed) {
  // Here, we delete about deletionPercentage % of all particles.
  double deletionPercentage = 30;
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
  for (auto &mol : *container) {
    if (doDelete[mol.getID()]) {
      autopas::internal::markParticleAsDeleted(mol);
    }
  }
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
 * @return Tuple of forces for all particles, ordered by particle id, and global values.
 */
template <bool globals>
std::tuple<std::vector<std::array<double, 3>>, TraversalComparison::Globals> TraversalComparison::calculateForces(
    autopas::ContainerOption containerOption, autopas::TraversalOption traversalOption,
    autopas::DataLayoutOption dataLayoutOption, autopas::Newton3Option newton3Option, size_t numMolecules,
    size_t numHaloMolecules, std::array<double, 3> boxMax, double cellSizeFactor, bool doSlightShift,
    DeletionPosition particleDeletionPosition) {
  // Construct container
  autopas::ContainerSelector<Molecule, FMCell> selector{_boxMin, boxMax, _cutoff};
  constexpr double skin = _cutoff * 0.1;
  selector.selectContainer(
      containerOption, autopas::ContainerSelectorInfo{cellSizeFactor, skin, 32, autopas::LoadEstimatorOption::none});
  auto container = selector.getCurrentContainer();
  autopas::LJFunctor<Molecule, FMCell, true /*applyShift*/, false /*useMixing*/, autopas::FunctorN3Modes::Both,
                     globals /*calculateGlobals*/>
      functor{_cutoff};
  functor.setParticleProperties(_eps * 24, _sig * _sig);

  autopasTools::generators::RandomGenerator::fillWithParticles(*container, Molecule({0., 0., 0.}, {0., 0., 0.}, 0),
                                                               container->getBoxMin(), container->getBoxMax(),
                                                               numMolecules);

  autopasTools::generators::RandomGenerator::fillWithHaloParticles(
      *container, Molecule({0., 0., 0.}, {0., 0., 0.}, numMolecules /*initial ID*/), container->getCutoff(),
      numHaloMolecules);

  auto traversal = autopas::TraversalSelector<FMCell>::generateTraversal(
      traversalOption, functor, container->getTraversalSelectorInfo(), dataLayoutOption, newton3Option);
  if (not traversal->isApplicable()) {
    return {};
  }

  if (particleDeletionPosition & DeletionPosition::beforeLists) {
    markSomeParticlesAsDeleted(container, numMolecules + numHaloMolecules, 19);
  }

  container->rebuildNeighborLists(traversal.get());

  if (doSlightShift) {
    executeShift(container, skin / 2, numMolecules + numHaloMolecules);
  }

  if (particleDeletionPosition & DeletionPosition::afterLists) {
    markSomeParticlesAsDeleted(container, numMolecules + numHaloMolecules, 99);
  }

  functor.initTraversal();
  container->iteratePairwise(traversal.get());
  functor.endTraversal(newton3Option);

  std::vector<std::array<double, 3>> forces(numMolecules);
  for (auto it = container->begin(autopas::IteratorBehavior::ownedOnly); it.isValid(); ++it) {
    EXPECT_TRUE(it->isOwned());
    forces.at(it->getID()) = it->getF();
  }

  if (globals) {
    return {forces, {functor.getUpot(), functor.getVirial()}};
  } else {
    return {forces, {0., 0.}};
  }
}

/**
 * Generates the reference for a simulation configuration that is specified by the given key.
 * For the reference a linked cells algorithm is used.
 * @param key The key that specifies the simulation.
 */
void TraversalComparison::generateReference(mykey_t key) {
  auto [numParticles, numHaloParticles, boxMax, doSlightShift, particleDeletionPosition, globals] = key;
  // Calculate reference forces
  if (globals) {
    auto [calculatedForces, calculatedGlobals] =
        calculateForces<true>(autopas::ContainerOption::linkedCells, autopas::TraversalOption::lc_c08,
                              autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled, numParticles,
                              numHaloParticles, boxMax, 1., doSlightShift, particleDeletionPosition);
    _forcesReference[key] = calculatedForces;
    _globalValuesReference[key] = calculatedGlobals;
  } else {
    auto [calculatedForces, calculatedGlobals] =
        calculateForces<false>(autopas::ContainerOption::linkedCells, autopas::TraversalOption::lc_c08,
                               autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled, numParticles,
                               numHaloParticles, boxMax, 1., doSlightShift, particleDeletionPosition);
    _forcesReference[key] = calculatedForces;
    _globalValuesReference[key] = calculatedGlobals;
  }
}

/**
 * This tests a given configuration against a reference configuration.
 */
TEST_P(TraversalComparison, traversalTest) {
  auto [containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles, numHaloParticles, boxMax,
        cellSizeFactor, doSlightShift, particleDeletionPosition, globals] = GetParam();

  // empirically determined and set near the minimal possible value for 2000 particles
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  constexpr double rel_err_tolerance = 1.0e-10;
  constexpr double rel_err_tolerance_globals = 1.0e-12;

  std::vector<std::array<double, 3>> calculatedForces;
  Globals calculatedGlobals;
  if (globals) {
    std::tie(calculatedForces, calculatedGlobals) =
        calculateForces<true>(containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles,
                              numHaloParticles, boxMax, cellSizeFactor, doSlightShift, particleDeletionPosition);
  } else {
    std::tie(calculatedForces, calculatedGlobals) =
        calculateForces<false>(containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles,
                               numHaloParticles, boxMax, cellSizeFactor, doSlightShift, particleDeletionPosition);
  }
  if (calculatedForces.empty()) {
    GTEST_SKIP_("Not applicable!");
  }

  TraversalComparison::mykey_t key{numParticles,  numHaloParticles,         boxMax,
                                   doSlightShift, particleDeletionPosition, globals};
  if (_forcesReference.count(key) == 0) {
    generateReference(key);
  }

  for (size_t i = 0; i < numParticles; ++i) {
    for (unsigned int d = 0; d < 3; ++d) {
      double calculatedForce = calculatedForces[i][d];
      double referenceForce = _forcesReference[key][i][d];
      EXPECT_NEAR(calculatedForce, referenceForce, std::fabs(calculatedForce * rel_err_tolerance))
          << "Particle id: " << i;
    }
  }

  auto &globalValuesReferenceRef = _globalValuesReference[key];
  if (globals) {
    EXPECT_NE(calculatedGlobals.upot, 0);
    EXPECT_NEAR(calculatedGlobals.upot, globalValuesReferenceRef.upot,
                rel_err_tolerance_globals * globalValuesReferenceRef.upot);

    EXPECT_NE(calculatedGlobals.virial, 0);
    EXPECT_NEAR(calculatedGlobals.virial, globalValuesReferenceRef.virial,
                rel_err_tolerance_globals * globalValuesReferenceRef.virial);
  }
}

/**
 * Lambda to generate a readable string out of the parameters of this test.
 */
static auto toString = [](const auto &info) {
  auto [containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles, numHaloParticles, boxMax,
        cellSizeFactor, doSlightShift, particleDeletionPosition, globals] = info.param;
  std::stringstream resStream;
  resStream << containerOption.to_string() << "_" << traversalOption.to_string() << "_" << dataLayoutOption.to_string()
            << "_" << (newton3Option == autopas::Newton3Option::enabled ? "_N3" : "_noN3") << "_NP" << numParticles
            << "_NH" << numHaloParticles << "_" << boxMax[0] << "_" << boxMax[1] << "_" << boxMax[2] << "_CSF_"
            << cellSizeFactor << "_" << (doSlightShift ? "withShift" : "noshift")
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
  std::vector<TestingTuple> params{};
  for (auto containerOption : autopas::ContainerOption::getAllOptions()) {
    for (auto traversalOption : autopas::compatibleTraversals::allCompatibleTraversals(containerOption)) {
      for (auto dataLayoutOption : autopas::DataLayoutOption::getAllOptions()) {
        for (auto newton3Option : autopas::Newton3Option::getAllOptions()) {
          for (auto numParticles : {100ul, 2000ul}) {
            for (auto boxMax : std::vector<std::array<double, 3>>{{3., 3., 3.}, {10., 10., 10.}}) {
              for (double cellSizeFactor : {0.5, 1., 2.}) {
                for (auto numHalo : {/*0ul,*/ 200ul}) {
                  for (bool slightMove : {true, false}) {
                    for (bool globals : {true, false}) {
                      for (DeletionPosition particleDeletionPosition :
                           {DeletionPosition::never, /*DeletionPosition::beforeLists, DeletionPosition::afterLists,*/
                            DeletionPosition::beforeAndAfterLists}) {
                        if (dataLayoutOption == autopas::DataLayoutOption::Value::cuda and
                            traversalOption == autopas::TraversalOption::Value::lc_c01_cuda and (boxMax[0] < 5.) and
                            (numParticles > 500)) {
                          // LJFunctor for cuda doesn't support this, yet: see
                          // https://github.com/AutoPas/AutoPas/issues/419
                          /// @todo reenable
                          continue;
                        }
                        params.emplace_back(containerOption, traversalOption, dataLayoutOption, newton3Option,
                                            numParticles, numHalo, boxMax, cellSizeFactor, slightMove,
                                            particleDeletionPosition, globals);
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
  return params;
}

INSTANTIATE_TEST_SUITE_P(Generated, TraversalComparison, ::testing::ValuesIn(TraversalComparison::getTestParams()),
                         toString);
