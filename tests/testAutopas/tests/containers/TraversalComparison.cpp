/**
 * @file TraversalComparison.cpp
 * @author humig
 * @date 12.07.19
 */

#include "TraversalComparison.h"

#include <autopas/selectors/ContainerSelector.h>
#include <autopas/selectors/TraversalSelector.h>
#include "autopas/utils/StringUtils.h"

#include <string>

using ::testing::_;  // anything is ok
using ::testing::Bool;
using ::testing::Combine;
using ::testing::Values;
using ::testing::ValuesIn;

std::vector<std::array<double, 3>> TraversalComparison::calculateForces(autopas::ContainerOption containerOption,
                                                                        autopas::TraversalOption traversalOption,
                                                                        autopas::DataLayoutOption dataLayoutOption,
                                                                        autopas::Newton3Option newton3Option,
                                                                        unsigned long numMolecules,
                                                                        std::array<double, 3> boxMax) {
  // Construct container
  autopas::ContainerSelector<Molecule, FMCell> selector{_boxMin, boxMax, _cutoff};
  selector.selectContainer(containerOption, autopas::ContainerSelectorInfo{1.0, _cutoff * 0.1});
  auto container = selector.getCurrentContainer();

  auto traversal = autopas::TraversalSelector<FMCell>::generateTraversal(
      traversalOption, _functor, container->getTraversalSelectorInfo(), dataLayoutOption, newton3Option);
  if (not traversal->isApplicable()) {
    return {};
  }

  RandomGenerator::fillWithParticles(*container, autopas::MoleculeLJ({0., 0., 0.}, {0., 0., 0.}, 0), numMolecules);

  container->rebuildNeighborLists(traversal.get());
  container->iteratePairwise(traversal.get());

  std::vector<std::array<double, 3>> forces(numMolecules);
  for (auto it = container->begin(); it.isValid(); ++it) {
    autopas::MoleculeLJ &m = *it;
    forces.at(m.getID()) = m.getF();
  }

  return forces;
}
void TraversalComparison::SetUpTestSuite() {
  autopas::Logger::create();

  autopas::MoleculeLJ::setEpsilon(_eps);
  autopas::MoleculeLJ::setSigma(_sig);

  // Calculate reference forces
  for (auto numParticles : _numParticlesVector) {
    for (auto boxMax : _boxMaxVector) {
      auto calculatedForces =
          calculateForces(autopas::ContainerOption::linkedCells, autopas::TraversalOption::c08,
                          autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled, numParticles, boxMax);
      _forcesReference[{numParticles, boxMax}] = calculatedForces;
    }
  }

  autopas::Logger::unregister();
}

TEST_P(TraversalComparison, traversalTest) {
  auto [containerOption, traversalOption, dataLayoutOption, newton3Option] = GetParam();

  // empirically determined and set near the minimal possible value for 2000 particles
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1.0e-10;

  for (auto numParticles : _numParticlesVector) {
    for (auto boxMax : _boxMaxVector) {
      auto calculatedForces =
          calculateForces(containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles, boxMax);
      if (calculatedForces.empty()) {
        GTEST_SKIP_("Not applicable!");
      }

      for (int i = 0; i < numParticles; ++i) {
        for (int d = 0; d < 3; ++d) {
          double calculatedForce = calculatedForces[i][d];
          double referenceForce = _forcesReference[{numParticles, boxMax}][i][d];
          EXPECT_NEAR(calculatedForce, referenceForce, std::fabs(calculatedForce * rel_err_tolerance));
        }
      }
    }
  }
}

static auto toString = [](testing::TestParamInfo<std::tuple<autopas::ContainerOption, autopas::TraversalOption,
                                                            autopas::DataLayoutOption, autopas::Newton3Option>>
                              param) {
  auto [containerOption, traversalOption, dataLayoutOption, newton3Option] = param.param;
  std::string res{};
  res += autopas::utils::StringUtils::to_string(containerOption);
  res += autopas::utils::StringUtils::to_string(traversalOption);
  res += dataLayoutOption == autopas::DataLayoutOption::aos ? "AoS" : "SoA";
  res += newton3Option == autopas::Newton3Option::enabled ? "Newton3" : "NoNewton3";
  res.erase(std::remove(res.begin(), res.end(), '-'), res.end());
  return res;
};

static auto getTestParams() {
  std::vector<
      std::tuple<autopas::ContainerOption, autopas::TraversalOption, autopas::DataLayoutOption, autopas::Newton3Option>>
      params{};
  auto containerOptions = autopas::allContainerOptions;
  // @todo: Readd verlet cluster lists as soon as the iterator works without dummy particles.
  containerOptions.erase(autopas::ContainerOption::verletClusterLists);
  for (auto containerOption : containerOptions) {
    for (auto traversalOption : autopas::compatibleTraversals::allCompatibleTraversals(containerOption)) {
      for (auto dataLayoutOption : {autopas::DataLayoutOption::aos, autopas::DataLayoutOption::soa}) {
        for (auto newton3Option : {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled}) {
          params.emplace_back(containerOption, traversalOption, dataLayoutOption, newton3Option);
        }
      }
    }
  }
  return params;
}

INSTANTIATE_TEST_SUITE_P(Generated, TraversalComparison, ::testing::ValuesIn(getTestParams()), toString);