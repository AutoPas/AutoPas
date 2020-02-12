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

  // autopasTools::generators::RandomGenerator::fillWithParticles(
  //    *container, autopas::MoleculeLJ({0., 0., 0.}, {0., 0., 0.}, 0), container->getBoxMin(), container->getBoxMax(),
  //    numMolecules, 42);
  size_t id = 0;
  container->addParticle(autopas::MoleculeLJ({0.10041, 0.989893, 2.07191}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.26746, 0.618795, 0.750385}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.90968, 2.59087, 0.904967}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.0747718, 1.09498, 2.29614}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.953581, 0.407185, 0.320236}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.28201, 0.245017, 1.65287}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.69396, 2.22981, 2.94528}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.656549, 1.36319, 1.55596}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.72153, 1.66569, 2.18564}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.88662, 2.38838, 1.75115}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.824299, 2.48879, 2.74104}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.89621, 0.756254, 0.359835}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.646592, 2.66593, 2.9507}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.55156, 2.7407, 1.04568}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.847702, 0.694282, 1.45287}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.16794, 2.9763, 1.69788}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.82081, 1.67025, 0.927694}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.76609, 2.3268, 2.29088}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.32205, 1.04833, 0.956575}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.507687, 2.93495, 0.34496}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.25883, 0.759249, 2.83375}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.99987, 0.655455, 0.590008}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.35971, 1.30205, 0.255938}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.31041, 2.85361, 2.99664}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.356094, 0.701307, 0.690922}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.80896, 1.86924, 0.667219}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.506842, 1.69005, 2.33747}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.43454, 1.45614, 1.66427}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.725419, 2.77819, 2.7126}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.68199, 0.285873, 2.64755}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.02695, 2.54471, 0.406799}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.86071, 1.54458, 1.06225}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.45072, 0.904293, 2.3643}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.70665, 0.214706, 2.21791}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.70329, 0.570799, 2.91921}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.394215, 2.37976, 1.78846}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.06143, 2.8866, 0.478508}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.398906, 1.32114, 1.93465}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.06318, 2.04656, 1.71283}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.77578, 0.728551, 1.9987}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.42333, 2.7555, 1.54341}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.83013, 1.61621, 0.0879934}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.89239, 1.06693, 0.992286}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.25669, 0.773581, 1.20699}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.47459, 0.476874, 1.77779}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.39381, 0.87109, 1.15755}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.182266, 1.93252, 1.04415}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.660775, 2.33143, 2.36529}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.59542, 1.39461, 1.41185}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.30825, 0.170393, 2.1404}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.306954, 1.59373, 1.8959}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.85037, 0.423859, 0.512115}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.93836, 0.316246, 1.57904}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.93064, 2.57293, 2.35262}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.13764, 1.04753, 2.8295}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.91543, 2.44134, 0.700588}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.07298, 2.6236, 2.63311}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.11713, 0.284378, 1.96454}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.48242, 2.8798, 0.359153}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.89427, 1.18805, 0.529547}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.03466, 1.495, 2.12327}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.930567, 0.345367, 2.54713}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.44268, 2.28373, 2.86338}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.0217254, 2.21437, 2.43631}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.37435, 0.352007, 0.483842}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.20385, 0.267434, 2.92518}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.90444, 1.34041, 2.54878}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.53755, 0.45754, 2.83316}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.50209, 1.93996, 2.71296}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.86124, 1.83422, 0.901005}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.39079, 0.868889, 2.39601}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.51406, 1.79946, 2.74137}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.0612, 0.242138, 2.0251}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.924576, 0.263864, 1.23947}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.360889, 2.63821, 1.59148}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.844731, 1.84206, 1.85891}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.76991, 1.7465, 0.199322}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.318691, 1.28405, 0.656862}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.151851, 2.78614, 2.59682}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.86481, 1.64738, 1.43105}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.765812, 1.03818, 2.29993}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.161819, 2.55224, 1.09939}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.90319, 0.613438, 1.34153}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.92829, 1.53801, 1.60539}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.167764, 1.8989, 1.24361}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.75924, 2.74363, 0.0856671}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.618151, 0.513543, 1.83217}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.817473, 0.832235, 0.116213}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.47433, 0.984085, 2.90235}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.07116, 0.848893, 1.54974}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.5022, 1.6147, 2.58791}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.80213, 1.77652, 2.14015}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.90152, 1.67972, 2.75359}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.24305, 0.608011, 1.2916}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({2.84844, 0.775774, 0.190505}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.09205, 2.53501, 2.93414}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.17772, 0.153166, 0.447682}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.00988284, 0.970639, 1.27992}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.126096, 2.44497, 2.264}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.0284479, 0.516129, 0.112895}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.57818, 0.0183292, 1.7276}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({1.16609, 1.82046, 0.504124}, {0., 0., 0.}, id++));
  container->addParticle(autopas::MoleculeLJ({0.306245, 1.72199, 2.18384}, {0., 0., 0.}, id++));

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
      EXPECT_NEAR(calculatedForce, referenceForce, std::fabs(calculatedForce * rel_err_tolerance))
          << "particle id: " << i;
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
              // for (double cellSizeFactor : {0.5, 1., 2.}) {
              for (double cellSizeFactor : {0.5}) {
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