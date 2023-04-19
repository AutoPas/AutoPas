/**
 * @file sph-main.cpp
 * @date 19.01.2018
 * @author seckler
 */

#include <array>
#include <cmath>
#include <iostream>

#include "autopas/AutoPas.h"
#include "autopas/sph/autopassph.h"
#include "autopas/utils/ArrayUtils.h"

using Particle = autopas::sph::SPHParticle;
using AutoPasContainer = autopas::AutoPas<Particle>;

void SetupIC(AutoPasContainer &sphSystem, double *end_time, const std::array<double, 3> &bBoxMax) {
  // Place SPH particles
  std::cout << "setup... started" << std::endl;
  const double dx = 1.0 / 128.0;
  unsigned int i = 0;
  for (double x = 0; x < bBoxMax[0] * 0.5; x += dx) {  // NOLINT
    for (double y = 0; y < bBoxMax[1]; y += dx) {      // NOLINT
      for (double z = 0; z < bBoxMax[2]; z += dx) {    // NOLINT
        Particle ith({x, y, z}, {0, 0, 0}, i++, 0.75, 0.012, 0.);
        ith.setDensity(1.0);
        ith.setEnergy(2.5);
        sphSystem.addParticle(ith);
      }
    }
  }
  for (double x = bBoxMax[0] * 0.5; x < bBoxMax[0] * 1.; x += dx * 2.0) {  // NOLINT
    for (double y = 0; y < bBoxMax[1]; y += dx) {                          // NOLINT
      for (double z = 0; z < bBoxMax[2]; z += dx) {                        // NOLINT
        Particle ith({x, y, z}, {0, 0, 0}, i++, 0.75, 0.012, 0.);
        ith.setDensity(0.5);
        ith.setEnergy(2.5);
        sphSystem.addParticle(ith);
      }
    }
  }
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->setMass(part->getMass() * bBoxMax[0] * bBoxMax[1] * bBoxMax[2] / (double)(i));
  }
  std::cout << "# of particles is... " << i << std::endl;

  // Set the end time
  *end_time = .018;
  // end_time for original example: (expects 50 timesteps)
  // *end_time = .12;
  // Fin.
  std::cout << "setup... completed" << std::endl;
}

void Initialize(AutoPasContainer &sphSystem) {
  std::cout << "initialize... started" << std::endl;
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->calcPressure();
  }
  std::cout << "initialize... completed" << std::endl;
}

double getTimeStepGlobal(AutoPasContainer &sphSystem) {
  double dt = 1.0e+30;  // set VERY LARGE VALUE
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->calcDt();
    if (part->getDt() < 0.002) {
      std::cout << "small time step for particle " << part->getID() << " at ["
                << autopas::utils::ArrayUtils::to_string(part->getR()) << "]" << std::endl;
    }
    dt = std::min(dt, part->getDt());
  }
  std::cout << "the time step dt is..." << dt << std::endl;
  return dt;
}

void leapfrogInitialKick(AutoPasContainer &sphSystem, const double dt) {
  using namespace autopas::utils::ArrayMath::literals;

  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->setVel_half(part->getV() + (part->getAcceleration() * (0.5 * dt)));
    part->setEng_half(part->getEnergy() + 0.5 * dt * part->getEngDot());
  }
}

void leapfrogFullDrift(AutoPasContainer &sphSystem, const double dt) {
  using namespace autopas::utils::ArrayMath::literals;

  // time becomes t + dt;
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->addR(part->getVel_half() * dt);
  }
}

void leapfrogPredict(AutoPasContainer &sphSystem, const double dt) {
  using namespace autopas::utils::ArrayMath::literals;

  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->addV(part->getAcceleration() * dt);
    part->addEnergy(part->getEngDot() * dt);
  }
}

void leapfrogFinalKick(AutoPasContainer &sphSystem, const double dt) {
  using namespace autopas::utils::ArrayMath::literals;

  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->setV(part->getVel_half() + (part->getAcceleration() * (0.5 * dt)));
    part->setEnergy(part->getEng_half() + 0.5 * dt * part->getEngDot());
  }
}

void setPressure(AutoPasContainer &sphSystem) {
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->calcPressure();
  }
}

void addEnteringParticles(AutoPasContainer &sphSystem, std::vector<Particle> &invalidParticles) {
  std::array<double, 3> boxMin = sphSystem.getBoxMin();
  std::array<double, 3> boxMax = sphSystem.getBoxMax();

  for (auto &p : invalidParticles) {
    // first we have to correct the position of the particles, s.t. they lie inside of the box.
    auto pos = p.getR();
    for (auto dim = 0; dim < 3; dim++) {
      if (pos[dim] < boxMin[dim]) {
        // has to be smaller than boxMax
        pos[dim] = std::min(std::nextafter(boxMax[dim], -1), pos[dim] + (boxMax[dim] - boxMin[dim]));
      } else if (pos[dim] >= boxMax[dim]) {
        // should at least be boxMin
        pos[dim] = std::max(boxMin[dim], pos[dim] - (boxMax[dim] - boxMin[dim]));
      }
    }
    p.setR(pos);
    // add moved particles again
    sphSystem.addParticle(p);
  }
}

/**
 * Get the required region for the regionparticleiterator.
 * @param boxMin
 * @param boxMax
 * @param diff
 * @param reqMin
 * @param reqMax
 * @param cutoff
 * @param shift
 */
void getRequiredHalo(double boxMin, double boxMax, int diff, double &reqMin, double &reqMax, double cutoff,
                     double &shift) {
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
 * Updates the halo particles.
 * This is done by copying the boundary particles to the halo particles plus
 * adding appropriate shifts
 * @param sphSystem
 */
void updateHaloParticles(AutoPasContainer &sphSystem) {
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
          std::cout << "skipping diff: " << autopas::utils::ArrayUtils::to_string(diff) << std::endl;
          continue;
        }
        // figure out from where we get our halo particles
        for (int i = 0; i < 3; ++i) {
          getRequiredHalo(boxMin[i], boxMax[i], diff[i], requiredHaloMin[i], requiredHaloMax[i], cutoff, shift[i]);
        }
        for (auto iterator =
                 sphSystem.getRegionIterator(requiredHaloMin, requiredHaloMax, autopas::IteratorBehavior::owned);
             iterator.isValid(); ++iterator) {
          Particle p = *iterator;
          p.addR(shift);
          sphSystem.addHaloParticle(p);
        }
      }
    }
  }
}

void densityPressureHydroForce(AutoPasContainer &sphSystem) {
  // declare the used functors
  autopas::sph::SPHCalcDensityFunctor<Particle> densityFunctor;
  autopas::sph::SPHCalcHydroForceFunctor<Particle> hydroForceFunctor;

  std::cout << "\nhaloupdate\n" << std::endl;

  // 1.first calculate density
  // 1.1 to calculate the density we need the halo particles
  updateHaloParticles(sphSystem);

  std::cout << "haloparticles... ";
  int haloparts = 0, innerparts = 0;
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    if (part->isHalo()) {
      haloparts++;
    } else {
      innerparts++;
    }
  }
  std::cout << haloparts << std::endl;
  std::cout << "particles... " << innerparts << std::endl;

  // 1.2 then calculate density
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->setDensity(0.);
    densityFunctor.AoSFunctor(*part, *part);
    part->setDensity(part->getDensity() / 2);
  }

  std::cout << "calculation of density... started" << std::endl;
  sphSystem.iteratePairwise(&densityFunctor);
  std::cout << "calculation of density... completed" << std::endl;
  // 1.3 delete halo particles, as their values are no longer valid
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::halo); part.isValid(); ++part) {
    sphSystem.deleteParticle(part);
  }

  // 2. then update pressure
  std::cout << "calculation of pressure... started" << std::endl;
  setPressure(sphSystem);
  std::cout << "calculation of pressure... completed" << std::endl;

  // 3 then calculate hydro force
  // 3.1 to calculate the density we need the halo particles
  updateHaloParticles(sphSystem);

  std::cout << "haloparticles... ";
  haloparts = 0, innerparts = 0;
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    if (part->isHalo()) {
      haloparts++;
    } else {
      innerparts++;
    }
  }
  std::cout << haloparts << std::endl;
  std::cout << "particles... " << innerparts << std::endl;

  // 3.2 then calculate hydro force
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    // self interaction leeds to:
    // 1) vsigmax = 2*part->getSoundSpeed()
    // 2) no change in acceleration
    part->setVSigMax(2 * part->getSoundSpeed());
    part->setAcceleration(std::array<double, 3>{0., 0., 0.});
    part->setEngDot(0.);
  }

  std::cout << "calculation of hydroforces... started" << std::endl;

  sphSystem.iteratePairwise(&hydroForceFunctor);
  std::cout << "calculation of hydroforces... completed" << std::endl;
}

void printConservativeVariables(AutoPasContainer &sphSystem) {
  using namespace autopas::utils::ArrayMath::literals;

  std::array<double, 3> momSum = {0., 0., 0.};  // total momentum
  double energySum = 0.0;                       // total energy
  for (auto it = sphSystem.begin(autopas::IteratorBehavior::owned); it.isValid(); ++it) {
    momSum = momSum + it->getV() * it->getMass();
    energySum += (it->getEnergy() + 0.5 * autopas::utils::ArrayMath::dot(it->getV(), it->getV())) * it->getMass();
  }
  printf("Energy     : %.16e\n", energySum);
  for (int i = 0; i < 3; ++i) {
    printf("Momentum[%d]: %.16e\n", i, momSum[i]);
    if (std::abs(momSum[i]) > 1.e-15) {
      std::stringstream ss;
      ss << std::setprecision(15) << "ERROR: The total momentum should cancel out (should be <1e-15 but is "
         << std::abs(momSum[i]) << ")!";
      throw std::runtime_error(ss.str());
    }
  }
}

int main() {
  std::array<double, 3> boxMin({0., 0., 0.}), boxMax{};
  boxMax[0] = 1.;
  boxMax[1] = boxMax[2] = boxMax[0] / 8.0;
  double cutoff = 0.03;               // 0.012*2.5=0.03; where 2.5 = kernel support radius
  unsigned int rebuildFrequency = 6;  // has to be multiple of two, as there are two functor calls per iteration.
  double skinToCutoffRatio = 0.15;

  AutoPasContainer sphSystem;
  sphSystem.setNumSamples(
      6);  // has to be multiple of 2, should also be multiple of rebuildFrequency (but this is not necessary).
  sphSystem.setBoxMin(boxMin);
  sphSystem.setBoxMax(boxMax);
  sphSystem.setCutoff(cutoff);
  sphSystem.setVerletSkinPerTimestep(skinToCutoffRatio * cutoff / rebuildFrequency);
  sphSystem.setVerletRebuildFrequency(rebuildFrequency);

  // In case you want to use another tuning strategy, you can do that using:
  // sphSystem.setTuningStrategyOption(autopas::TuningStrategyOption::activeHarmony);

  // Debug output of AutoPas can be enabled using:
  // autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);

  std::set<autopas::ContainerOption> allowedContainers{autopas::ContainerOption::linkedCells,
                                                       autopas::ContainerOption::verletLists,
                                                       autopas::ContainerOption::verletListsCells};
  sphSystem.setAllowedContainers(allowedContainers);

  sphSystem.init();

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

    // 1.2.1 positions have changed, so the container needs to be updated!
    auto invalidParticles = sphSystem.updateContainer();

    // 1.2.2 adjust positions based on boundary conditions (here: periodic)
    addEnteringParticles(sphSystem, invalidParticles);

    // 1.3 Leap frog: predict
    leapfrogPredict(sphSystem, dt);
    // 1.4 Calculate density, pressure and hydrodynamic forces
    densityPressureHydroForce(sphSystem);
    // 1.5 get time step
    dt = getTimeStepGlobal(sphSystem);
    // 1.6 Leap frog: final Kick
    leapfrogFinalKick(sphSystem, dt);

    printConservativeVariables(sphSystem);
    std::cout << "time in iteration " << step << ": " << perlooptimer.stop() << std::endl;
  }
}
