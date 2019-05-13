/**
 * @file AutoPasVerletLikeInterfaceTest.cpp
 * @author seckler
 * @date 13.05.19
 */

#include "AutoPasInterfaceTest.h"
#include "testingHelpers/commonTypedefs.h"

constexpr double skin = 0.2;

template <typename AutoPasT>
void defaultInit(AutoPasT& autoPas) {
  std::array<double, 3> boxMin{0., 0., 0.};
  std::array<double, 3> boxMax{10., 10., 10.};
  autoPas.setBoxMin(boxMin);
  autoPas.setBoxMax(boxMax);
  autoPas.setCutoff(1.);
  autoPas.setVerletSkin(skin);
  // init autopas
  autoPas.init();
}

void doSimulationLoop(autopas::AutoPas<Particle, FPCell>& autoPas) {
  // 1. update Container
  autoPas.updateContainer();

  // 2. leaving and entering particles
  // 2a. identify and send leaving particles (iteration over halo)

  // 2b. get entering particles (addParticle)

  // 3. halo particles
  // 3a. identify and send inner particles that are in the halo of other autopas instances or itself.

  // 3b. get halo particles

  // 4. iteratePairwise
}

void testVerletInterface(autopas::ContainerOption containerOption) {
  // create AutoPas object
  autopas::AutoPas<Particle, FPCell> autoPas;
  autoPas.setAllowedContainers(std::vector<autopas::ContainerOption>{containerOption});

  defaultInit(autoPas);

  // create two particles with distance .5
  double distance = .5;
  std::array<double, 3> pos1{9.99, 5., 5.};
  std::array<double, 3> distVec{0., distance, 0.};
  std::array<double, 3> pos2 = autopas::ArrayMath::sub(pos1, distVec);

  Particle particle1(pos1, {0., 0., 0.}, 0);
  Particle particle2(pos2, {0., 0., 0.}, 1);

  // add the two particles!
  autoPas.addParticle(particle1);
  autoPas.addParticle(particle2);

  // do first simulation loop
  doSimulationLoop(autoPas);

  // update positions a bit and do loop again
  std::array<double, 3> moveVec{skin/3.,0.,0.};
  for (auto iter = autoPas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    iter->setR(autopas::ArrayMath::add(iter->getR(), moveVec));
  }

  // do second simulation loop
  doSimulationLoop(autoPas);

}

TEST_F(AutoPasInterfaceTest, InterfaceTest) {
  // this test checks the correct behavior of the autopas interface.
  testVerletInterface(autopas::ContainerOption::linkedCells);
}