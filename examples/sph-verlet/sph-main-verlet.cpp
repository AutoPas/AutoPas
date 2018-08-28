/**
 * @file sph-main-verlet.cpp
 * @date 10.05.2018
 * @author seckler
 */

#include <array>
#include <cmath>
#include <iostream>

#include "autopas/autopasIncludes.h"
#include "autopas/sph/autopassph.h"

typedef autopas::VerletLists<autopas::sph::SPHParticle> Container;
typedef autopas::C08Traversal<autopas::FullParticleCell<autopas::sph::SPHParticle>,
                              autopas::sph::SPHCalcHydroForceFunctor, false, false>
    DummyTraversal;

std::map<std::array<int, 3>, std::vector<autopas::sph::SPHParticle*>> sph_verlet_particle_list;

void SetupIC(Container& sphSystem, double* end_time, const std::array<double, 3>& bBoxMax) {
  // Place SPH particles
  std::cout << "setup... started" << std::endl;
  const double dx = 1.0 / 128.0;
  unsigned int i = 0;
  for (double x = 0; x < bBoxMax[0] * 0.5; x += dx) {  // NOLINT
    for (double y = 0; y < bBoxMax[1]; y += dx) {      // NOLINT
      for (double z = 0; z < bBoxMax[2]; z += dx) {    // NOLINT
        // std::array<double, 3> r, std::array<double, 3> v, unsigned long id,
        // double mass, double smth, double snds
        autopas::sph::SPHParticle ith({x, y, z}, {0, 0, 0}, i++, 0.75, 0.012, 0.);
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
  for (double x = bBoxMax[0] * 0.5; x < bBoxMax[0] * 1.; x += dx * 2.0) {  // NOLINT
    for (double y = 0; y < bBoxMax[1]; y += dx) {                          // NOLINT
      for (double z = 0; z < bBoxMax[2]; z += dx) {                        // NOLINT
        // std::array<double, 3> r, std::array<double, 3> v, unsigned long id,
        // double mass, double smth, double snds
        autopas::sph::SPHParticle ith({x, y, z}, {0, 0, 0}, i++, 0.75, 0.012, 0.);
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
    part->setMass(part->getMass() * bBoxMax[0] * bBoxMax[1] * bBoxMax[2] / (double)(i));
  }
  std::cout << "# of particles is... " << i << std::endl;

  // Set the end time
  *end_time = .12;
  // Fin.
  std::cout << "setup... completed" << std::endl;
}

void Initialize(Container& sphSystem) {
  std::cout << "initialize... started" << std::endl;
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->calcPressure();
  }
  std::cout << "initialize... completed" << std::endl;
}

double getTimeStepGlobal(Container& sphSystem) {
  double dt = 1.0e+30;  // set VERY LARGE VALUE
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->calcDt();
    if (part->getDt() < 0.002) {
      std::cout << "small time step for particle " << part->getID() << " at [" << part->getR()[0] << ", "
                << part->getR()[1] << ", " << part->getR()[2] << "]" << std::endl;
    }
    dt = std::min(dt, part->getDt());
  }
  std::cout << "the time step dt is..." << dt << std::endl;
  return dt;
}

void leapfrogInitialKick(Container& sphSystem, const double dt) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->setVel_half(
        autopas::ArrayMath::add(part->getV(), autopas::ArrayMath::mulScalar(part->getAcceleration(), 0.5 * dt)));
    part->setEng_half(part->getEnergy() + 0.5 * dt * part->getEngDot());
  }
}

void leapfrogFullDrift(Container& sphSystem, const double dt) {
  // time becomes t + dt;
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::ownedOnly); part.isValid(); ++part) {
    part->addR(autopas::ArrayMath::mulScalar(part->getVel_half(), dt));
  }
}

void leapfrogPredict(Container& sphSystem, const double dt) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->addV(autopas::ArrayMath::mulScalar(part->getAcceleration(), dt));
    part->addEnergy(part->getEngDot() * dt);
  }
}

void leapfrogFinalKick(Container& sphSystem, const double dt) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->setV(
        autopas::ArrayMath::add(part->getVel_half(), autopas::ArrayMath::mulScalar(part->getAcceleration(), 0.5 * dt)));
    part->setEnergy(part->getEng_half() + 0.5 * dt * part->getEngDot());
  }
}

void setPressure(Container& sphSystem) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->calcPressure();
  }
}

void periodicBoundaryUpdate(Container& sphSystem, std::array<double, 3> boxMin, std::array<double, 3> boxMax) {
  std::vector<autopas::sph::SPHParticle> invalidParticles;
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    auto posVec = part->getR();
    bool modified = false;
    for (unsigned int dim = 0; dim < 3; dim++) {
      auto& pos = posVec[dim];
      while (pos < boxMin[dim]) {
        pos += boxMax[dim] - boxMin[dim];
        modified = true;
      }
      while (pos > boxMax[dim]) {
        pos -= boxMax[dim] - boxMin[dim];
        modified = true;
      }
      if (pos == boxMax[dim]) {
#if 0
        pos = boxMin[dim];
#else
        // next smaller double before pos/boxMax
        pos = nextafter(pos, 0.);
#endif
        modified = true;
      }
    }
    if (modified) {
      part->setR(posVec);
      invalidParticles.push_back(*part);
      part.deleteCurrentParticle();
      // we have moved particles very far, so we have to delete them at
      // this position and store them again!
      // we can not add particles while an iterator is active, so we have to do
      // that later.
    }
  }
  for (auto p : invalidParticles) {
    // add moved particles again
    sphSystem.addParticle(p);
  }
}

/**
 * get the required region for the regionparticleiterator.
 * @param boxMin
 * @param boxMax
 * @param diff
 * @param reqMin
 * @param reqMax
 * @param cutoff
 * @param shift
 */
void getRequiredHalo(double boxMin, double boxMax, int diff, double& reqMin, double& reqMax, double cutoff,
                     double& shift) {
  if (diff == 0) {
    reqMin = boxMin;
    reqMax = boxMax;
    shift = 0;
  } else if (diff == -1) {
    reqMin = boxMax - cutoff;
    reqMax = boxMax;
    shift = boxMin - boxMax;
  } else if (diff == 1) {
    reqMin = boxMin;
    reqMax = boxMin + cutoff;
    shift = boxMax - boxMin;
  }
}

/**
 *
 * updates the halo particles
 * this is done by copying the boundary particles to the halo particles plus
 * adding appropriate shifts
 * @param sphSystem
 * @param addParticles decides whether the particles should be added (true) or
 * whether the particle should be updated (needed for verlet lists)
 */
void updateHaloParticles(Container& sphSystem, bool addParticles) {
  std::array<double, 3> boxMin = sphSystem.getBoxMin();
  std::array<double, 3> boxMax = sphSystem.getBoxMax();
  std::array<double, 3> requiredHaloMin{0., 0., 0.}, requiredHaloMax{0., 0., 0.};
  std::array<int, 3> diff{0, 0, 0};
  std::array<double, 3> shift{0., 0., 0.};
  double cutoff = sphSystem.getCutoff();
  for (diff[0] = -1; diff[0] < 2; diff[0]++) {
    for (diff[1] = -1; diff[1] < 2; diff[1]++) {
      for (diff[2] = -1; diff[2] < 2; diff[2]++) {
        if (not diff[0] and not diff[1] and not diff[2]) {
          // at least one dimension has to be non-zero
          std::cout << "skipping diff: " << diff[0] << ", " << diff[1] << ", " << diff[2] << std::endl;
          continue;
        }
        // figure out from where we get our halo particles
        for (int i = 0; i < 3; ++i) {
          getRequiredHalo(boxMin[i], boxMax[i], diff[i], requiredHaloMin[i], requiredHaloMax[i], cutoff, shift[i]);
        }
        if (addParticles) {
          sph_verlet_particle_list[diff].clear();
          for (auto iterator = sphSystem.getRegionIterator(requiredHaloMin, requiredHaloMax); iterator.isValid();
               ++iterator) {
            sph_verlet_particle_list[diff].push_back(&*iterator);
            autopas::sph::SPHParticle p = *iterator;
            p.addR(shift);
            // add the halo particle
            sphSystem.addHaloParticle(p);
          }
        } else {
          for (auto particleptr : sph_verlet_particle_list.at(diff)) {
            autopas::sph::SPHParticle p = *particleptr;
            p.addR(shift);
            // update halo particle
            sphSystem.updateHaloParticle(p);
          }
        }
      }
    }
  }
}

/**
 * deletes the halo particles
 * @param sphSystem
 */
void deleteHaloParticles(Container& sphSystem) { sphSystem.deleteHaloParticles(); }

void densityPressureHydroForce(Container& sphSystem) {
  // declare the used functors
  autopas::sph::SPHCalcDensityFunctor densityFunctor;
  autopas::sph::SPHCalcHydroForceFunctor hydroForceFunctor;

  std::cout << "\nhaloupdate\n" << std::endl;

  // 1.first calculate density
  // 1.1 to calculate the density we need the halo particles
  updateHaloParticles(sphSystem, sphSystem.needsRebuild());

  // 1.2 then calculate density
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->setDensity(0.);
    densityFunctor.AoSFunctor(*part, *part);
    part->setDensity(part->getDensity() / 2);
  }
  std::cout << "calculation of density... started" << std::endl;
  DummyTraversal dummyTraversal({0, 0, 0}, &hydroForceFunctor);
  sphSystem.iteratePairwiseAoS(&densityFunctor, &dummyTraversal);
  std::cout << "calculation of density... completed" << std::endl;
  // 1.3 delete halo particles, as their values are no longer valid
  // deleteHaloParticles(sphSystem);

  // 2. then update pressure
  std::cout << "calculation of pressure... started" << std::endl;
  setPressure(sphSystem);
  std::cout << "calculation of pressure... completed" << std::endl;

  // 0.3 then calculate hydro force
  // 0.3.1 to calculate the density we need the halo particles
  updateHaloParticles(sphSystem, false);

  // 0.3.2 then calculate hydro force
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    // self interaction leeds to:
    // 1) vsigmax = 2*part->getSoundSpeed()
    // 2) no change in acceleration
    part->setVSigMax(2 * part->getSoundSpeed());
    part->setAcceleration(std::array<double, 3>{0., 0., 0.});
    part->setEngDot(0.);
  }
  std::cout << "calculation of hydroforces... started" << std::endl;
  sphSystem.iteratePairwiseAoS(&hydroForceFunctor, &dummyTraversal);
  std::cout << "calculation of hydroforces... completed" << std::endl;
  // 0.3.3 delete halo particles, as their values are no longer valid
  // deleteHaloParticles(sphSystem);
}

void printConservativeVariables(Container& sphSystem) {
  std::array<double, 3> momSum = {0., 0., 0.};  // total momentum
  double energySum = 0.0;                       // total energy
  for (auto it = sphSystem.begin(autopas::IteratorBehavior::ownedOnly); it.isValid(); ++it) {
    momSum = autopas::ArrayMath::add(momSum, autopas::ArrayMath::mulScalar(it->getV(), it->getMass()));
    energySum += (it->getEnergy() + 0.5 * autopas::ArrayMath::dot(it->getV(), it->getV())) * it->getMass();
  }
  printf("%.16e\n", energySum);
  printf("%.16e\n", momSum[0]);
  printf("%.16e\n", momSum[1]);
  printf("%.16e\n", momSum[2]);
}

int main() {
  autopas::Logger::create();

  unsigned int rebuildFrequency = 6;
  double skintocutoff = 0.04;

  AutoPasLogger->set_level(spdlog::level::level_enum::debug);
  AutoPasLogger->set_pattern("[%n] [%l] %v");
  std::array<double, 3> boxMin({0., 0., 0.}), boxMax{};
  boxMax[0] = 1.;
  boxMax[1] = boxMax[2] = boxMax[0] / 8.0;
  double cutoff = 0.03;  // 0.012*2.5=0.03; where 2.5 = kernel support radius

  // Container sphSystem(boxMin, boxMax, cutoff);
  Container sphSystem(boxMin, boxMax, cutoff, skintocutoff * cutoff /*skin*/,
                      rebuildFrequency /*every second time step*/ /*rebuild frequency*/);
  double dt;
  double t_end;
  SetupIC(sphSystem, &t_end, boxMax);
  Initialize(sphSystem);

  // 0.1 ---- GET INITIAL FORCES OF SYSTEM ----
  densityPressureHydroForce(sphSystem);

  std::cout << "\n----------------------------" << std::endl;

  // 0.2 get time step
  dt = getTimeStepGlobal(sphSystem);
  //---- INITIAL FORCES ARE NOW CALCULATED ----

  printConservativeVariables(sphSystem);

  // 1 ---- START MAIN LOOP ----
  size_t step = 0;
  autopas::utils::Timer perlooptimer;
  for (double time = 0.; time < t_end; time += dt, ++step) {
    perlooptimer.start();
    std::cout << "\n-------------------------\ntime step " << step << "(t = " << time << ")..." << std::endl;
    // 1.1 Leap frog: Initial Kick & Full Drift
    leapfrogInitialKick(sphSystem, dt);
    leapfrogFullDrift(sphSystem, dt);

    // for verlet-lists this only needs to be done, if particles moved too far.
    if (sphSystem.needsRebuild() or sphSystem.isContainerUpdateNeeded()) {
      // ensure that there are no halo particles if we need to update the
      // container
      deleteHaloParticles(sphSystem);

      // 1.2.1 positions have changed, so the container needs to be updated!
      sphSystem.updateContainer();

      // 1.2.2 adjust positions based on boundary conditions (here: periodic)
      periodicBoundaryUpdate(sphSystem, boxMin, boxMax);
    }

    // 1.3 Leap frog: predict
    leapfrogPredict(sphSystem, dt);

    // 1.4 Calculate density, pressure and hydrodynamic forces
    densityPressureHydroForce(sphSystem);

    // 1.5 get time step
    dt = getTimeStepGlobal(sphSystem);

    // 1.6 Leap frog: final Kick
    leapfrogFinalKick(sphSystem, dt);

    //    for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    //      printf(
    //          "%lu\t%lf\t%lf\t%lf\t%lf\t%lf\t"
    //          "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
    //          part->getID(), part->getMass(), part->getR()[0],
    //          part->getR()[1],
    //          part->getR()[2], part->getV()[0], part->getV()[1],
    //          part->getV()[2],
    //          part->getDensity(), part->getEnergy(), part->getPressure(),
    //          part->getAcceleration()[0], part->getAcceleration()[1],
    //          part->getAcceleration()[2], part->getEngDot(), part->getDt());
    //    }

    printConservativeVariables(sphSystem);
    std::cout << "time in iteration " << step << ": " << perlooptimer.stop() << std::endl;
  }
}
