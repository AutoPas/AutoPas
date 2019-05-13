/**
 * @file AutoPasVerletLikeInterfaceTest.cpp
 * @author seckler
 * @date 13.05.19
 */

#include "AutoPasInterfaceTest.h"
#include "testingHelpers/commonTypedefs.h"

constexpr double cutoff = 1.;
constexpr double skin = 0.2;
constexpr std::array<double, 3> boxMin{0., 0., 0.};
constexpr std::array<double, 3> boxMax{10., 10., 10.};
constexpr double eps = 1.;
constexpr double sigma = 1.;
constexpr double shift = 0.1;

template <typename AutoPasT>
void defaultInit(AutoPasT& autoPas) {
  autoPas.setBoxMin(boxMin);
  autoPas.setBoxMax(boxMax);
  autoPas.setCutoff(cutoff);
  autoPas.setVerletSkin(skin);
  // init autopas
  autoPas.init();
}
auto identifyAndSendLeavingParticles(autopas::AutoPas<Molecule, FMCell>& autoPas) { return std::vector<Molecule>{}; }

auto identifyAndSendHaloParticles(autopas::AutoPas<Molecule, FMCell>& autoPas) { return std::vector<Molecule>{}; }

void addEnteringParticles(autopas::AutoPas<Molecule, FMCell>& autoPas, std::vector<Molecule> enteringParticles) {}

void addHaloParticles(autopas::AutoPas<Molecule, FMCell>& autoPas, std::vector<Molecule> haloParticles) {}

template <typename Functor>
void doSimulationLoop(autopas::AutoPas<Molecule, FMCell>& autoPas, Functor* functor) {
  // 1. update Container
  autoPas.updateContainer();

  // 2. leaving and entering particles
  // 2a. identify and send leaving particles (iteration over halo)
  auto leavingParticles = identifyAndSendLeavingParticles(autoPas);

  // 2b. get+add entering particles (addParticle)
  const auto& enteringParticles = leavingParticles;
  addEnteringParticles(autoPas, enteringParticles);

  // 3. halo particles
  // 3a. identify and send inner particles that are in the halo of other autopas instances or itself.
  auto sendHaloParticles = identifyAndSendHaloParticles(autoPas);

  // 3b. get halo particles
  const auto& recvHaloParticles = sendHaloParticles;
  addHaloParticles(autoPas, recvHaloParticles);

  // 4. iteratePairwise
  autoPas.iteratePairwise(functor);
}

void testVerletInterface(autopas::ContainerOption containerOption) {
  // create AutoPas object
  autopas::AutoPas<Molecule, FMCell> autoPas;
  autoPas.setAllowedContainers(std::vector<autopas::ContainerOption>{containerOption});

  defaultInit(autoPas);

  // create two particles with distance .5
  double distance = .5;
  std::array<double, 3> pos1{9.99, 5., 5.};
  std::array<double, 3> distVec{0., distance, 0.};
  std::array<double, 3> pos2 = autopas::ArrayMath::sub(pos1, distVec);

  Molecule particle1(pos1, {0., 0., 0.}, 0);
  Molecule particle2(pos2, {0., 0., 0.}, 1);

  // add the two particles!
  autoPas.addParticle(particle1);
  autoPas.addParticle(particle2);

  autopas::LJFunctor<Molecule, FMCell> functor(cutoff, eps, sigma, shift, boxMin, boxMax, true);
  // do first simulation loop
  doSimulationLoop(autoPas, &functor);

  // update positions a bit and do loop again
  std::array<double, 3> moveVec{skin / 3., 0., 0.};
  for (auto iter = autoPas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    iter->setR(autopas::ArrayMath::add(iter->getR(), moveVec));
  }

  // do second simulation loop
  doSimulationLoop(autoPas, &functor);
}

TEST_F(AutoPasInterfaceTest, InterfaceTest) {
  // this test checks the correct behavior of the autopas interface.
  testVerletInterface(autopas::ContainerOption::linkedCells);
}