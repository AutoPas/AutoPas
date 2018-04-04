//
// Created by seckler on 19.01.18.
//

#include <array>
#include <iostream>
#include "autopas.h"
#include "sph/autopassph.h"

typedef autopas::LinkedCells<
    autopas::sph::SPHParticle,
    autopas::FullParticleCell<autopas::sph::SPHParticle>>
    LCContainer;

void SetupIC(LCContainer& sphSystem, double* end_time,
             const std::array<double, 3>& bBoxMax) {
  // Place SPH particles
  std::cout << "setup... started" << std::endl;
  const double dx = 1.0 / 128.0;
  unsigned int i = 0;
  for (double x = 0; x < bBoxMax[0] * 0.5; x += dx) {  // NOLINT
    for (double y = 0; y < bBoxMax[1]; y += dx) {      // NOLINT
      for (double z = 0; z < bBoxMax[2]; z += dx) {    // NOLINT
        // std::array<double, 3> r, std::array<double, 3> v, unsigned long id,
        // double mass, double smth, double snds
        autopas::sph::SPHParticle ith({x, y, z}, {0, 0, 0}, i++, 0.75, 0.012,
                                      0.);
        ith.setDensity(1.0);
        ith.setEnergy(2.5);
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
  for (double x = bBoxMax[0] * 0.5; x < bBoxMax[0] * 1.;
       x += dx * 2.0) {                              // NOLINT
    for (double y = 0; y < bBoxMax[1]; y += dx) {    // NOLINT
      for (double z = 0; z < bBoxMax[2]; z += dx) {  // NOLINT
        // std::array<double, 3> r, std::array<double, 3> v, unsigned long id,
        // double mass, double smth, double snds
        autopas::sph::SPHParticle ith({x, y, z}, {0, 0, 0}, i++, 0.75, 0.012,
                                      0.);
        ith.setDensity(0.5);
        ith.setEnergy(2.5);
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
  std::cout << "setup... completed" << std::endl;
}

void Initialize(LCContainer& sphSystem) {
  std::cout << "initialize... started" << std::endl;
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->calcPressure();
  }
  std::cout << "initialize... completed" << std::endl;
}

double getTimeStepGlobal(LCContainer& sphSystem) {
  double dt = 1.0e+30;  // set VERY LARGE VALUE
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->calcDt();
    dt = std::min(dt, part->getDt());
  }
  std::cout << "the time step dt is..." << dt << std::endl;
  return dt;
}

void leapfrogInitialKick(LCContainer& sphSystem, const double dt) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    // TODO: start here again; plus array math
    part->setVel_half(autopas::arrayMath::add(
        part->getV(),
        autopas::arrayMath::mulScalar(part->getAcceleration(), 0.5 * dt)));
    part->setEng_half(part->getEnergy() + 0.5 * dt * part->getEngDot());
  }
}

void leapfrogFullDrift(LCContainer& sphSystem, const double dt) {
  // time becomes t + dt;
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->addR(autopas::arrayMath::mulScalar(part->getVel_half(), dt));
  }
}

void leapfrogPredict(LCContainer& sphSystem, const double dt) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->addV(autopas::arrayMath::addScalar(part->getAcceleration(), dt));
    part->addEnergy(part->getEngDot() * dt);
  }
}

void leapfrogFinalKick(LCContainer& sphSystem, const double dt) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->setV(autopas::arrayMath::add(
        part->getVel_half(),
        autopas::arrayMath::mulScalar(part->getAcceleration(), 0.5 * dt)));
    part->setEnergy(part->getEng_half() + 0.5 * dt * part->getEngDot());
  }
}

void setPressure(LCContainer& sphSystem) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->calcPressure();
  }
}

void periodicBoundaryUpdate(LCContainer& sphSystem,
                            std::array<double, 3> boxMin,
                            std::array<double, 3> boxMax) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    auto posVec = part->getR();
    for (unsigned int dim = 0; dim < 2; dim++) {
      auto& pos = posVec[dim];
      while (pos < boxMin[dim]) {
        pos += boxMax[dim] - boxMin[dim];
      }
      while (pos > boxMax[dim]) {
        pos -= boxMax[dim] - boxMin[dim];
      }
      if (pos == boxMax[dim]) {
        pos = boxMin[dim];
      }
    }
    part->setR(posVec);
  }
}

/**
 * updates the halo particles
 * this is done by copying the boundary particles to the halo particles
 */
void updateHaloParticles() {
  // TODO: implement
}

/**
 * deletes the halo particles
 */
void deleteHaloParticles() {
  // TODO: implement
}

void densityPressureHydroForce(LCContainer& sphSystem) {
  // declare the used functors
  autopas::sph::SPHCalcDensityFunctor densityFunctor;
  autopas::sph::SPHCalcHydroForceFunctor hydroForceFunctor;

  // 1.first calculate density
  // 1.1 to calculate the density we need the halo particles
  updateHaloParticles();
  // 1.2 then calculate density
  std::cout << "calculation of density... started" << std::endl;
  sphSystem.iteratePairwiseAoS2(&densityFunctor);
  std::cout << "calculation of density... completed" << std::endl;
  // 1.3 delete halo particles, as their values are no longer valid
  deleteHaloParticles();

  // 2. then update pressure
  std::cout << "calculation of pressure... started" << std::endl;
  setPressure(sphSystem);
  std::cout << "calculation of pressure... completed" << std::endl;

  // 0.3 then calculate hydro force
  // 0.3.1 to calculate the density we need the halo particles
  updateHaloParticles();
  // 0.3.2 then calculate hydro force
  std::cout << "calculation of hydroforces... started" << std::endl;
  sphSystem.iteratePairwiseAoS2(&hydroForceFunctor);
  std::cout << "calculation of hydroforces... completed" << std::endl;
  // 0.3.3 delete halo particles, as their values are no longer valid
  deleteHaloParticles();
}

void printConservativeVariables(LCContainer& sphSystem) {
  std::array<double, 3> momSum = {0., 0., 0.};  // total momentum
  double energySum = 0.0;  // total enegry
  for (auto it = sphSystem.begin(); it.isValid(); ++it) {
    momSum = autopas::arrayMath::add(momSum,autopas::arrayMath::mulScalar(it->getV(), it->getMass()));
    energySum += (it->getEnergy() + 0.5 * autopas::arrayMath::dot(it->getV(), it->getV())) *
           it->getMass();
  }
  printf("%.16e\n", energySum);
  printf("%.16e\n", momSum[0]);
  printf("%.16e\n", momSum[1]);
  printf("%.16e\n", momSum[2]);
}

int main() {
  std::array<double, 3> boxMin({0., 0., 0.}), boxMax{};
  boxMax[0] = 1.;
  boxMax[1] = boxMax[2] = boxMax[0] / 8.0;
  double cutoff = .012;

  LCContainer sphSystem(boxMin, boxMax, cutoff);
  double dt;
  double t_end;
  SetupIC(sphSystem, &t_end, boxMax);
  Initialize(sphSystem);

  // 0.1 ---- GET INITIAL FORCES OF SYSTEM ----
  densityPressureHydroForce(sphSystem);

  // 0.2 get time step
  dt = getTimeStepGlobal(sphSystem);
  //---- INITIAL FORCES ARE NOW CALCULATED ----

  printConservativeVariables(sphSystem);

  // 1 ---- START MAIN LOOP ----
  size_t step = 0;
  for (double time = 0.; time < t_end; time += dt, ++step) {
    std::cout << "time step " << step << "(t = " << time << ")... started"
              << std::endl;
    // 1.1 Leap frog: Initial Kick & Full Drift
    leapfrogInitialKick(sphSystem, dt);
    leapfrogFullDrift(sphSystem, dt);
    // 1.2 adjust positions based on boundary conditions (here: periodic)
    periodicBoundaryUpdate(sphSystem, boxMin, boxMax);
    // 1.3 Leap frog: predict
    leapfrogPredict(sphSystem, dt);
    // 1.4 Calculate density, pressure and hydrodynamic forces
    densityPressureHydroForce(sphSystem);
    // 1.5 get time step
    dt = getTimeStepGlobal(sphSystem);
    // 1.6 Leap frog: final Kick
    leapfrogFinalKick(sphSystem, dt);
    std::cout << "time step " << step << "(t = " << time << ")... completed"
              << std::endl;
    printConservativeVariables(sphSystem);
  }
}