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

std::vector<Molecule> convertToEnteringParticles(const std::vector<Molecule>& leavingParticles) {
  std::vector<Molecule> enteringParticles{leavingParticles};
  for (auto& p : enteringParticles) {
    auto pos = p.getR();
    for (auto dim = 0; dim < 3; dim++) {
      if (pos[dim] < boxMin[dim]) {
        // should at most be boxMax
        pos[dim] = std::min(std::nextafter(boxMax[dim], -1), pos[dim] + (boxMax[dim] - boxMin[dim]));
      } else if (pos[dim] >= boxMax[dim]) {
        // should at least be boxMin
        pos[dim] = std::max(boxMin[dim], pos[dim] - (boxMax[dim] - boxMin[dim]));
      }
    }
    p.setR(pos);
  }
  return enteringParticles;
}

auto identifyAndSendHaloParticles(autopas::AutoPas<Molecule, FMCell>& autoPas) {
  std::vector<Molecule> haloParticles;

  for (short x : {-1, 0, 1}) {
    for (short y : {-1, 0, 1}) {
      for (short z : {-1, 0, 1}) {
        std::array<short, 3> direction{x, y, z};
        std::array<double, 3> min{}, max{}, shiftVec{};
        for (size_t dim = 0; dim < 3; ++dim) {
          // The search domain has to be enlarged as the position of the particles is not certain.
          if (direction[dim] == -1) {
            min[dim] = boxMin[dim] - skin;
            max[dim] = boxMin[dim] + cutoff + skin;
          } else if (direction[dim] == 1) {
            min[dim] = boxMax[dim] - cutoff - skin;
            max[dim] = boxMax[dim] + skin;
          } else {  // 0
            min[dim] = boxMin[dim] - skin;
            min[dim] = boxMax[dim] + skin;
          }
          shiftVec[dim] = -(boxMax[dim] - boxMin[dim]) * direction[dim];
        }
        // here it is important to only iterate over the owned particles!
        for (auto iter = autoPas.getRegionIterator(min, max, autopas::IteratorBehavior::ownedOnly); iter.isValid();
             ++iter) {
          auto particleCopy = *iter;
          particleCopy.addR(shiftVec);
          haloParticles.push_back(particleCopy);
        }
      }
    }
  }

  return haloParticles;
}

void addEnteringParticles(autopas::AutoPas<Molecule, FMCell>& autoPas, std::vector<Molecule> enteringParticles) {
  for (auto& p : enteringParticles) {
    autoPas.addParticle(p);
  }
}

void addHaloParticles(autopas::AutoPas<Molecule, FMCell>& autoPas, std::vector<Molecule> haloParticles) {
  for (auto& p : haloParticles) {
    autoPas.addHaloParticle(p);
  }
}

template <typename Functor>
void doSimulationLoop(autopas::AutoPas<Molecule, FMCell>& autoPas, Functor* functor) {
  // 1. update Container; return value is vector of invalid = leaving particles!
  auto invalidParticles = autoPas.updateContainer();

  // 2. leaving and entering particles
  const auto& sendLeavingParticles = invalidParticles;
  // 2b. get+add entering particles (addParticle)
  const auto& enteringParticles = convertToEnteringParticles(sendLeavingParticles);
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