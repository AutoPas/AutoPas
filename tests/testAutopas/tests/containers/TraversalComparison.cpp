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
static double rand_doub(double fmin, double fmax){
  double f = (double) rand()/((double)RAND_MAX);
  return fmin + f*(fmax - fmin);
}
template <int nmol>
static std::array<Particle, nmol> gen_rand_particles(int bmvindex){
  std::array<Particle, nmol> res;
  for(int i = 0; i < nmol; ++i){
    res[i] = autopas::MoleculeLJ({rand_doub(0.,static_cast<double>(TraversalComparison::_boxMaxVector[bmvindex][0])),
                                  rand_doub(0.,static_cast<double>(TraversalComparison::_boxMaxVector[bmvindex][1])),
                                  rand_doub(0.,static_cast<double>(TraversalComparison::_boxMaxVector[bmvindex][2]))},
                                 {0.,0.,0.}, i);
  }
  return res;
}
static int counter_debug = 0;
template <int nmol>
static std::array<Particle, nmol> make_particles(){
  std::array<std::array<double, 3>, nmol> pos;
  int counter = 0;
  pos[counter++] = {0.209266, 2.84798, 1.57799}; //important
  pos[counter++] = {0.424808, 1.82091, 0.0489017}; //important
  pos[counter++] = {0.326426, 2.99677, 0.654771}; //important
  counter_debug = counter;
  std::array<Particle, nmol> res;
  for (int i = 0; i < counter; ++i) {
    res[i] = autopas::MoleculeLJ(pos[i],{0.,0.,0.},i);
  }
  return res;
}
static std::array<Particle, 2000> rand_particles1 = gen_rand_particles<2000>(0);
static std::array<Particle, 2000> rand_particles2 = gen_rand_particles<2000>(1);
static std::array<Particle, 2000> rand_particles3 = gen_rand_particles<2000>(2);
static std::array<Particle, 2000> rand_particles4 = gen_rand_particles<2000>(3);
static std::array<Particle, 2000> minimal_particles = make_particles<2000>();
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
  for(; id < counter_debug;) {
//    double x = rand_doub(0.,boxMax[0]);
//    double y = rand_doub(0.,boxMax[1]);
//    double z = rand_doub(0.,boxMax[2]);
    auto pos = minimal_particles[id].getR();
    container->addParticle(autopas::MoleculeLJ(pos, {0., 0., 0.}, id++));

  }
//  container->addParticle(autopas::MoleculeLJ({0.506842, 1.69005, 2.33747}, {0., 0., 0.}, id++));
//  container->addParticle(autopas::MoleculeLJ({0.398906, 1.32114, 1.93465}, {0., 0., 0.}, id++));
//  //container->addParticle(autopas::MoleculeLJ({2.0,0.5,0.5}, {0., 0., 0.}, id++));
//  //container->addParticle(autopas::MoleculeLJ({1.5,0.5,0.5}, {0., 0., 0.}, id++));
//  //container->addParticle(autopas::MoleculeLJ({0.346245, 1.72199, 2.18384}, {0., 0., 0.}, id++));
//  //container->addParticle(autopas::MoleculeLJ({0.306954, 1.59373, 1.8959}, {0., 0., 0.}, id++));
//  container->addParticle(autopas::MoleculeLJ({0.824299, 2.48879, 2.74104}, {0., 0., 0.}, id++));
//  //container->addParticle(autopas::MoleculeLJ({2.75, 0.5, 0.5}, {0., 0., 0.}, id++));
  //container->addParticle(autopas::MoleculeLJ({0.646592, 2.66593, 2.9507}, {0., 0., 0.}, id++));
  //container->addParticle(autopas::MoleculeLJ({0.725419, 2.77819, 2.7126}, {0., 0., 0.}, id++));
  //container->addParticle(autopas::MoleculeLJ({0.97298, 2.6236, 2.63311}, {0., 0., 0.}, id++));
  //container->addParticle(autopas::MoleculeLJ({1.09205, 2.53501, 2.93414}, {0., 0., 0.}, id++));

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
          calculateForces(autopas::ContainerOption::verletClusterLists, autopas::TraversalOption::verletClusters,
                          autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, numParticles, boxMax, 1.);
      _forcesReference[{numParticles, boxMax}] = calculatedForces;
      _globalValuesReference[{numParticles, boxMax}] = calculatedGlobals;
    }
  }

  //autopas::Logger::unregister();
}

TEST_P(TraversalComparison, traversalTest) {
  auto [containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles, boxMax, cellSizeFactor] =
      GetParam();
  //SetUpTestSuite();

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

  std::vector<std::array<double, 3>> ref = _forcesReference[{numParticles, boxMax}];
  for (size_t i = 0; i < numParticles; ++i) {
    for (unsigned int d = 0; d < 3; ++d) {
      double calculatedForce = calculatedForces[i][d];
      double referenceForce = ref[i][d];
      EXPECT_NEAR(calculatedForce, referenceForce, std::fabs(calculatedForce * rel_err_tolerance))
          << "particle id: " << i;
    }
  }

  auto &globalValuesReferenceRef = _globalValuesReference[{numParticles, boxMax}];
  for (size_t index : {0, 1}) {
//    EXPECT_NE(calculatedGlobals[index], 0);
//    EXPECT_NEAR(calculatedGlobals[index], globalValuesReferenceRef[index],
//                rel_err_tolerance_globals * globalValuesReferenceRef[index]);
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
  auto containerOption = autopas::ContainerOption::linkedCells; {
    for (auto traversalOption : {autopas::TraversalOption::c01}) {
      for (auto dataLayoutOption : {autopas::DataLayoutOption::aos}) {
        for (auto newton3Option : {autopas::Newton3Option::disabled}) {
          for (auto numParticles : _numParticlesVector) {
            for (auto boxMax : _boxMaxVector) {
              // for (double cellSizeFactor : {0.5, 1., 2.}) {
              for (double cellSizeFactor : {0.5, 1., 2.}) {
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