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

using ::testing::_;  // anything is ok
using ::testing::Bool;
using ::testing::Combine;
using ::testing::Values;
using ::testing::ValuesIn;

std::array<double, 3> randomShift(double magnitude, std::mt19937 &generator) {
  std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  double theta = 2 * M_PI * uniform01(generator);
  double phi = acos(1 - 2 * uniform01(generator));
  double x = sin(phi) * cos(theta) * magnitude;
  double y = sin(phi) * sin(theta) * magnitude;
  double z = cos(phi) * magnitude;
  return {x, y, z};
}

template <class Container>
void executeSlightShift(Container &container, double magnitude, unsigned long totalNumParticles) {
  std::vector<std::array<double, 3>> shiftVectorByID(totalNumParticles);
  unsigned seed = 42;
  std::mt19937 generator(seed);
  for (auto &elem : shiftVectorByID) {
    elem = randomShift(magnitude, generator);
  }
  for (auto it = container->begin(autopas::IteratorBehavior::ownedOnly); it.isValid(); ++it) {
    it->addR(shiftVectorByID[it->getID()]);
  }
  // assumes that halo particles have other IDs than owned particles!
  for (auto it = container->begin(autopas::IteratorBehavior::haloOnly); it.isValid(); ++it) {
    it->addR(shiftVectorByID[it->getID()]);
  }
}

std::tuple<std::vector<std::array<double, 3>>, std::array<double, 2>> TraversalComparison::calculateForces(
    autopas::ContainerOption containerOption, autopas::TraversalOption traversalOption,
    autopas::DataLayoutOption dataLayoutOption, autopas::Newton3Option newton3Option, unsigned long numMolecules,
    unsigned long numHaloMolecules, std::array<double, 3> boxMax, double cellSizeFactor, bool doSlightShift) {
  // Construct container
  autopas::ContainerSelector<Molecule, FMCell> selector{_boxMin, boxMax, _cutoff};
  double skin = _cutoff * 0.1;
  selector.selectContainer(containerOption, autopas::ContainerSelectorInfo{cellSizeFactor, skin, 32});
  auto container = selector.getCurrentContainer();
  autopas::LJFunctor<Molecule, FMCell, true /*applyShift*/, false /*useMixing*/, autopas::FunctorN3Modes::Both,
                     true /*calculateGlobals*/>
      functor{_cutoff};
  functor.setParticleProperties(_eps * 24, _sig * _sig);

  auto traversal = autopas::TraversalSelector<FMCell>::generateTraversal(
      traversalOption, functor, container->getTraversalSelectorInfo(), dataLayoutOption, newton3Option);
  if (not traversal->isApplicable()) {
    return {};
  }

  autopasTools::generators::RandomGenerator::fillWithParticles(
      *container, autopas::MoleculeLJ({0., 0., 0.}, {0., 0., 0.}, 0), container->getBoxMin(), container->getBoxMax(),
      numMolecules);

  autopasTools::generators::RandomGenerator::fillWithHaloParticles(
      *container, autopas::MoleculeLJ({0., 0., 0.}, {0., 0., 0.}, numMolecules), container->getCutoff(),
      numHaloMolecules);

  container->rebuildNeighborLists(traversal.get());

  if (doSlightShift) {
    executeSlightShift(container, skin / 2, numMolecules + numHaloMolecules);
  }

  functor.initTraversal();
  container->iteratePairwise(traversal.get());
  functor.endTraversal(newton3Option);

  std::vector<std::array<double, 3>> forces(numMolecules);
  for (auto it = container->begin(autopas::IteratorBehavior::ownedOnly); it.isValid(); ++it) {
    forces.at(it->getID()) = it->getF();
  }

  return {forces, {functor.getUpot(), functor.getVirial()}};
}

void TraversalComparison::generateReference(mykey_t key) {
  auto [numParticles, numHaloParticles, boxMax, doSlightShift] = key;
  // Calculate reference forces
  auto [calculatedForces, calculatedGlobals] = calculateForces(
      autopas::ContainerOption::linkedCells, autopas::TraversalOption::c08, autopas::DataLayoutOption::aos,
      autopas::Newton3Option::enabled, numParticles, numHaloParticles, boxMax, 1., doSlightShift);
  _forcesReference[key] = calculatedForces;
  _globalValuesReference[key] = calculatedGlobals;
}

TEST_P(TraversalComparison, traversalTest) {
  auto [containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles, numHaloParticles, boxMax,
        cellSizeFactor, doSlightShift] = GetParam();

  // empirically determined and set near the minimal possible value for 2000 particles
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1.0e-11;
  double rel_err_tolerance_globals = 1.0e-12;

  auto [calculatedForces, calculatedGlobals] =
      calculateForces(containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles, numHaloParticles,
                      boxMax, cellSizeFactor, doSlightShift);
  if (calculatedForces.empty()) {
    GTEST_SKIP_("Not applicable!");
  }

  TraversalComparison::mykey_t key{numParticles, numHaloParticles, boxMax, doSlightShift};
  if (_forcesReference.count(key) == 0) {
    generateReference(key);
  }

  for (size_t i = 0; i < numParticles; ++i) {
    for (unsigned int d = 0; d < 3; ++d) {
      double calculatedForce = calculatedForces[i][d];
      double referenceForce = _forcesReference[key][i][d];
      EXPECT_NEAR(calculatedForce, referenceForce, std::fabs(calculatedForce * rel_err_tolerance));
    }
  }

  auto &globalValuesReferenceRef = _globalValuesReference[key];
  for (size_t index : {0, 1}) {
    EXPECT_NE(calculatedGlobals[index], 0);
    EXPECT_NEAR(calculatedGlobals[index], globalValuesReferenceRef[index],
                rel_err_tolerance_globals * globalValuesReferenceRef[index]);
  }
}

static auto toString = [](const auto &info) {
  auto [containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles, numHaloParticles, boxMax,
        cellSizeFactor, doSlightShift] = info.param;
  std::stringstream resStream;
  resStream << containerOption.to_string() << traversalOption.to_string() << dataLayoutOption.to_string()
            << newton3Option.to_string() << "_NP" << numParticles << "_NH" << numHaloParticles << "_" << boxMax[0]
            << "_" << boxMax[1] << "_" << boxMax[2] << "_CSF_" << cellSizeFactor << "_"
            << (doSlightShift ? "shift" : "noshift");
  std::string res = resStream.str();
  std::replace(res.begin(), res.end(), '-', '_');
  std::replace(res.begin(), res.end(), '.', '_');
  return res;
};

auto TraversalComparison::getTestParams() {
  std::vector<TestingTuple> params{};
  for (auto containerOption : autopas::ContainerOption::getAllOptions()) {
    for (auto traversalOption : autopas::compatibleTraversals::allCompatibleTraversals(containerOption)) {
      for (auto dataLayoutOption : autopas::DataLayoutOption::getAllOptions()) {
        for (auto newton3Option : autopas::Newton3Option::getAllOptions()) {
          for (auto numParticles : _numParticlesVector) {
            for (auto boxMax : _boxMaxVector) {
              for (double cellSizeFactor : {0.5, 1., 2.}) {
                for (double numHalo : _numHaloVector) {
                  for (bool slightMove : {true, false}) {
                    if (dataLayoutOption == autopas::DataLayoutOption::Value::cuda and
                        traversalOption == autopas::TraversalOption::Value::c01Cuda and boxMax[0] < 5. and
                        numParticles > 500) {
                      // LJFunctor for cuda doesn't support this, yet: see https://github.com/AutoPas/AutoPas/issues/419
                      continue;
                    }
                    params.emplace_back(containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles,
                                        numHalo, boxMax, cellSizeFactor, slightMove);
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