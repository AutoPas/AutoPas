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

std::tuple<std::vector<std::array<double, 3>>, std::array<double, 2>> TraversalComparison::calculateForces(
    autopas::ContainerOption containerOption, autopas::TraversalOption traversalOption,
    autopas::DataLayoutOption dataLayoutOption, autopas::Newton3Option newton3Option, unsigned long numMolecules,
    std::array<double, 3> boxMax, double cellSizeFactor) {
  // Construct container
  autopas::ContainerSelector<Molecule, FMCell> selector{_boxMin, boxMax, _cutoff};
  selector.selectContainer(containerOption, autopas::ContainerSelectorInfo{cellSizeFactor, _cutoff * 0.1, 32});
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

  container->rebuildNeighborLists(traversal.get());

  functor.initTraversal();
  container->iteratePairwise(traversal.get());
  functor.endTraversal(newton3Option);

  std::vector<std::array<double, 3>> forces(numMolecules);
  for (auto it = container->begin(); it.isValid(); ++it) {
    forces.at(it->getID()) = it->getF();
  }

  return {forces, {functor.getUpot(), functor.getVirial()}};
}

void TraversalComparison::SetUpTestSuite() {
  autopas::Logger::create();

  // Calculate reference forces
  for (auto numParticles : _numParticlesVector) {
    for (auto boxMax : _boxMaxVector) {
      auto [calculatedForces, calculatedGlobals] =
          calculateForces(autopas::ContainerOption::linkedCells, autopas::TraversalOption::c08,
                          autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled, numParticles, boxMax, 1.);
      _forcesReference[{numParticles, boxMax}] = calculatedForces;
      _globalValuesReference[{numParticles, boxMax}] = calculatedGlobals;
    }
  }

  autopas::Logger::unregister();
}

TEST_P(TraversalComparison, traversalTest) {
  auto [containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles, boxMax, cellSizeFactor] =
      GetParam();

  // empirically determined and set near the minimal possible value for 2000 particles
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1.0e-10;
  double rel_err_tolerance_globals = 1.0e-10;

  auto [calculatedForces, calculatedGlobals] = calculateForces(containerOption, traversalOption, dataLayoutOption,
                                                               newton3Option, numParticles, boxMax, cellSizeFactor);
  if (calculatedForces.empty()) {
    GTEST_SKIP_("Not applicable!");
  }

  for (size_t i = 0; i < numParticles; ++i) {
    for (unsigned int d = 0; d < 3; ++d) {
      double calculatedForce = calculatedForces[i][d];
      double referenceForce = _forcesReference[{numParticles, boxMax}][i][d];
      EXPECT_NEAR(calculatedForce, referenceForce, std::fabs(calculatedForce * rel_err_tolerance));
    }
  }

  auto &globalValuesReferenceRef = _globalValuesReference[{numParticles, boxMax}];
  for (size_t index : {0, 1}) {
    EXPECT_NE(calculatedGlobals[index], 0);
    EXPECT_NEAR(calculatedGlobals[index], globalValuesReferenceRef[index],
                rel_err_tolerance_globals * globalValuesReferenceRef[index]);
  }
}

static auto toString = [](const auto &info) {
  auto [containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles, boxMax, cellSizeFactor] =
      info.param;
  std::stringstream resStream;
  resStream << containerOption.to_string() << traversalOption.to_string() << dataLayoutOption.to_string()
            << newton3Option.to_string() << "_" << numParticles << "_" << boxMax[0] << "_" << boxMax[1] << "_"
            << boxMax[2] << "_CSF_" << cellSizeFactor;
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
                if (dataLayoutOption == autopas::DataLayoutOption::Value::cuda and
                    traversalOption == autopas::TraversalOption::Value::c01Cuda and boxMax[0] < 5. and
                    numParticles > 500) {
                  // LJFunctor for cuda doesn't support this, yet: see https://github.com/AutoPas/AutoPas/issues/419
                  continue;
                }
                params.emplace_back(containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles,
                                    boxMax, cellSizeFactor);
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