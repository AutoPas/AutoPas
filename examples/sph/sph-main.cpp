//
// Created by seckler on 19.01.18.
//

#include <array>
#include <cassert>
#include <iostream>
#include "autopas.h"
#include "sph/autopassph.h"

typedef autopas::LinkedCells<
    autopas::sph::SPHParticle,
    autopas::FullParticleCell<autopas::sph::SPHParticle>>
    LCContainer;

template <class Container>
void SetupIC(Container& sphSystem, double* end_time,
             const std::array<double, 3>& bBoxMax) {
  // Place SPH particles

  const double dx = 1.0 / 128.0;
  unsigned int i = 0;
  for (double x = 0; x < bBoxMax[0] * 0.5; x += dx) {
    for (double y = 0; y < bBoxMax[1]; y += dx) {
      for (double z = 0; z < bBoxMax[2]; z += dx) {
        // std::array<double, 3> r, std::array<double, 3> v, unsigned long id,
        // double mass, double smth, double snds
        autopas::sph::SPHParticle ith({x, y, z}, {0, 0, 0}, i++, 0.75, 0.012,
                                      0.);
        ith.setDensity(1.0);
        // ith.pos.x = x;
        // ith.pos.y = y;
        // ith.pos.z = z;
        // ith.dens = 1.0;
        // ith.mass = 0.75;
        // ith.eng = 2.5;
        // ith.id = i++;
        // ith.smth = 0.012;
        sphSystem.addParticle(ith);
      }
    }
  }
  for (double x = bBoxMax[0] * 0.5; x < bBoxMax[0] * 1.; x += dx * 2.0) {
    for (double y = 0; y < bBoxMax[1]; y += dx) {
      for (double z = 0; z < bBoxMax[2]; z += dx) {
        // std::array<double, 3> r, std::array<double, 3> v, unsigned long id,
        // double mass, double smth, double snds
        autopas::sph::SPHParticle ith({x, y, z}, {0, 0, 0}, i++, 0.75, 0.012,
                                      0.);
        ith.setDensity(0.5);
        // ith.pos.x = x;
        // ith.pos.y = y;
        // ith.pos.z = z;
        // ith.dens = 0.5;
        // ith.mass = 0.75;
        // ith.eng = 2.5;
        // ith.id = i++;
        // ith.smth = 0.012;
        sphSystem.addParticle(ith);
      }
    }
  }
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->setMass(part->getMass() * bBoxMax[0] * bBoxMax[1] * bBoxMax[2] /
                  (double)(i));
  }
  std::cout << "# of particles is... " << i << std::endl;

  // Set the end time
  *end_time = 0.12;
  // Fin.
  std::cout << "setup..." << std::endl;
}

template <class Container>
void Initialize(Container& sphSystem) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->calcPressure();
  }
}

template <class Container>
double getTimeStepGlobal(Container& sphSystem) {
  double dt = 1.0e+30;  // set VERY LARGE VALUE
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    dt = std::min(dt, part->getDt());
  }
  return dt;
}

template <class Container>
void InitialKick(Container& sphSystem, const double dt) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    // TODO: start here again; plus array math
    part->setVel_half(autopas::arrayMath::add(
        part->getV(),
        autopas::arrayMath::mulScalar(part->getAcceleration(), 0.5 * dt)));
    part->setEng_half(part->getEnergy() + 0.5 * dt * part->getEngDot());
  }
}

template <class Container>
void FullDrift(Container& sphSystem, const double dt) {
  // time becomes t + dt;
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->addR(autopas::arrayMath::mulScalar(part->getVel_half(), dt));
  }
}

template <class Container>
void Predict(LCContainer& sphSystem, const double dt) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->addV(autopas::arrayMath::addScalar(part->getAcceleration(),dt));
    part->addEnergy(part->getEngDot()*dt);
  }
}

template <class Container>
void FinalKick(Container& sphSystem, const double dt) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->setV(autopas::arrayMath::add(
        part->getVel_half(),
        autopas::arrayMath::mulScalar(part->getAcceleration(), 0.5 * dt)));
    part->setEnergy(part->getEng_half() + 0.5 * dt * part->getEngDot());
  }
}

template <class Container>
void setPressure(Container& sphSystem) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->calcPressure();
  }
}

int main() {
  std::array<double, 3> boxMin({0., 0., 0.}), boxMax{};
  boxMax[0] = 1.0;
  boxMax[1] = boxMax[2] = boxMax[0] / 8.0;
  double cutoff = 3.;

  LCContainer sphSystem(boxMin, boxMax, cutoff);
  double t_end;
  SetupIC(sphSystem, &t_end, boxMax);
  Initialize(sphSystem);
}