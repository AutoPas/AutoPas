/**
 * @file C04SoATraversalTest.cpp
 * @author C. Menges
 * @date 06.07.2019
 */

#include "C04SoATraversalTest.h"
#include "testingHelpers/GridGenerator.h"
#include "testingHelpers/NumThreadGuard.h"

using ::testing::_;

/**
 * Tests whether the interaction between two cells is correct
 */
TEST_F(C04SoATraversalTest, testTraversal) {
  std::array<size_t, 3> edgeLength = {3, 3, 3};

  autopas::LJFunctor<autopas::Particle, FPCell> functor(1., 1., 1., 1.);
  std::vector<FPCell> cells;
  cells.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);

  autopas::Particle defaultParticle;
  defaultParticle.setR({1.75, 2.1, 1.75});
  cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 2ul, 1ul, edgeLength)].addParticle(defaultParticle);
  defaultParticle.setR({1.75, 1.6, 1.75});
  cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 1ul, 1ul, edgeLength)].addParticle(defaultParticle);

  NumThreadGuard numThreadGuard(1);

  autopas::C04SoATraversal<FPCell, autopas::LJFunctor<autopas::Particle, FPCell>, autopas::DataLayoutOption::soa, true>
      c04SoATraversal(edgeLength, &functor, 1, {1., 1., 1.});
  c04SoATraversal.setCellsToTraverse(cells);
  c04SoATraversal.initTraversal();
  c04SoATraversal.traverseParticlePairs();
  c04SoATraversal.endTraversal();

  size_t num = 0;
  for (auto &cell : cells) {
    num += cell.numParticles();
  }
  EXPECT_EQ(num, 2);
  auto iter = cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 1ul, 1ul, edgeLength)].begin();
  auto iter2 = cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 2ul, 1ul, edgeLength)].begin();
  EXPECT_DOUBLE_EQ(iter->getF()[1], -390144.0);
  EXPECT_EQ(-iter->getF()[1], iter2->getF()[1]);
}

constexpr double cutoff = 1.;
constexpr double skin = 0.2;
constexpr std::array<double, 3> boxMin{0., 0., 0.};
constexpr std::array<double, 3> boxMax{10., 10., 10.};
constexpr double eps = 1.;
constexpr double sigma = 1.;
constexpr double shift = 0.1;
constexpr std::array<double, 3> zeroArr = {0., 0., 0.};

/**
 * Tests whether the interaction between two cells is correct
 */
TEST_F(C04SoATraversalTest, testHaloTraversal) {
    // create AutoPas object
    autopas::AutoPas<Molecule, FMCell> autoPas;

    autoPas.setAllowedContainers({autopas::ContainerOption::linkedCells});
    autoPas.setAllowedTraversals({autopas::TraversalOption::c04SoA});
    autoPas.setAllowedDataLayouts({autopas::DataLayoutOption::soa});
    autoPas.setAllowedNewton3Options({autopas::Newton3Option::enabled});
    autoPas.setAllowedCellSizeFactors(autopas::NumberSetFinite<double>(std::set<double>({0.5})));

    autoPas.setBoxMin(boxMin);
    autoPas.setBoxMax(boxMax);
    autoPas.setCutoff(cutoff);
    autoPas.setVerletSkin(skin);
    autoPas.setVerletRebuildFrequency(2);
    autoPas.setNumSamples(2);
    // init autopas
    autoPas.init();

    double mul = 5.;
    double mid = 5.;
    double x_diff=0.0, y_diff=1.0, z_diff=0.0;
    std::array<double, 3> edge{x_diff * mul + 5.0, y_diff * mul + mid, z_diff * mul + mid};

    std::array<double, 3> diff = {x_diff * 1., y_diff * 1., z_diff * 1.};
    double distance = .5;
    diff = autopas::ArrayMath::mulScalar(autopas::ArrayMath::normalize(diff), distance / 2.);

    auto pos1 = autopas::ArrayMath::sub(edge, diff);
    auto pos2 = autopas::ArrayMath::add(edge, diff);
    long id = 0;
    Molecule particle1(pos1, {0., 0., 0.}, id++);
    autoPas.addParticle(particle1);
    Molecule particle2(pos2, {0., 0., 0.}, id++);
    autoPas.addOrUpdateHaloParticle(particle2);
    autopas::LJFunctor<Molecule, FMCell, autopas::FunctorN3Modes::Both, true /*calculate globals*/> functor(cutoff, eps,
                                                                                                            sigma, shift);
    autoPas.iteratePairwise(&functor);

    std::cout << "Par1: " << autoPas.begin()->toString() << std::endl;
    std::cout << "Par2: " << (++autoPas.begin())->toString() << std::endl;
}