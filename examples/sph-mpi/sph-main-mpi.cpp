//
// Created by seckler on 19.01.18.
//

#include <mpi.h>
#include <array>
#include <iostream>
#include "autopasIncludes.h"
#include "sph/autopassph.h"

typedef autopas::LinkedCells<
    autopas::sph::SPHParticle,
    autopas::FullParticleCell<autopas::sph::SPHParticle>>
    Container;

// typedef autopas::DirectSum<
//    autopas::sph::SPHParticle,
//    autopas::FullParticleCell<autopas::sph::SPHParticle>>
//    Container;

void SetupIC(Container& sphSystem, double* end_time,
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
        if (autopas::inBox(ith.getR(), sphSystem.getBoxMin(),
                           sphSystem.getBoxMax())) {
          sphSystem.addParticle(ith);
        }
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
        if (autopas::inBox(ith.getR(), sphSystem.getBoxMin(),
                           sphSystem.getBoxMax())) {
          sphSystem.addParticle(ith);
        }
      }
    }
  }
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->setMass(part->getMass() * bBoxMax[0] * bBoxMax[1] * bBoxMax[2] /
                  (double)(i));
  }

  // we are incrementing i independent of whether we add a particle to the local
  // sphsystem, or not. So it is still the global number of particles.
  std::cout << "total # of particles is... " << i << std::endl;

  // Set the end time
  *end_time = 0.12;
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

double getTimeStepGlobal(Container& sphSystem, MPI::Comm& comm) {
  double dt = 1.0e+30;  // set VERY LARGE VALUE
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->calcDt();
    if (part->getDt() < 0.002) {
      std::cout << "small time step for particle " << part->getID() << " at ["
                << part->getR()[0] << ", " << part->getR()[1] << ", "
                << part->getR()[2] << "]" << std::endl;
    }
    dt = std::min(dt, part->getDt());
  }

  // MPI, global reduction of minimum
  comm.Allreduce(MPI::IN_PLACE, &dt, 1, MPI::DOUBLE, MPI::MIN);

  std::cout << "the time step dt is..." << dt << std::endl;
  return dt;
}

void leapfrogInitialKick(Container& sphSystem, const double dt) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->setVel_half(autopas::arrayMath::add(
        part->getV(),
        autopas::arrayMath::mulScalar(part->getAcceleration(), 0.5 * dt)));
    part->setEng_half(part->getEnergy() + 0.5 * dt * part->getEngDot());
  }
}

void leapfrogFullDrift(Container& sphSystem, const double dt) {
  // time becomes t + dt;
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->addR(autopas::arrayMath::mulScalar(part->getVel_half(), dt));
  }
}

void leapfrogPredict(Container& sphSystem, const double dt) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->addV(autopas::arrayMath::mulScalar(part->getAcceleration(), dt));
    part->addEnergy(part->getEngDot() * dt);
  }
}

void leapfrogFinalKick(Container& sphSystem, const double dt) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->setV(autopas::arrayMath::add(
        part->getVel_half(),
        autopas::arrayMath::mulScalar(part->getAcceleration(), 0.5 * dt)));
    part->setEnergy(part->getEng_half() + 0.5 * dt * part->getEngDot());
  }
}

void setPressure(Container& sphSystem) {
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->calcPressure();
  }
}

void periodicBoundaryUpdate(Container& sphSystem, std::array<double, 3> boxMin,
                            std::array<double, 3> boxMax) {
  std::vector<autopas::sph::SPHParticle> invalidParticles;
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    auto posVec = part->getR();
    for (unsigned int dim = 0; dim < 3; dim++) {
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
    invalidParticles.push_back(*part);
    part.deleteCurrentParticle();
  }
  for (auto p: invalidParticles) {
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
void getRequiredHalo(double boxMin, double boxMax, int diff, double& reqMin,
                     double& reqMax, double cutoff, double& shift) {
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
 */
void updateHaloParticles(Container& sphSystem) {
  // TODO: mpi exchanges!
  std::array<double, 3> boxMin = sphSystem.getBoxMin();
  std::array<double, 3> boxMax = sphSystem.getBoxMax();
  std::array<double, 3> requiredHaloMin, requiredHaloMax;
  std::array<int, 3> diff;
  std::array<double, 3> shift;
  double cutoff = sphSystem.getCutoff();
  for (diff[0] = -1; diff[0] < 2; diff[0]++) {
    for (diff[1] = -1; diff[1] < 2; diff[1]++) {
      for (diff[2] = -1; diff[2] < 2; diff[2]++) {
        if (not diff[0] and not diff[1] and not diff[2]) {
          // at least one dimension has to be non-zero
          std::cout << "skipping diff: " << diff[0] << ", " << diff[1] << ", "
                    << diff[2] << std::endl;
          continue;
        }
        // figure out from where we get our halo particles
        for (int i = 0; i < 3; ++i) {
          getRequiredHalo(boxMin[i], boxMax[i], diff[i], requiredHaloMin[i],
                          requiredHaloMax[i], cutoff, shift[i]);
        }
        for (auto iterator =
                 sphSystem.getRegionIterator(requiredHaloMin, requiredHaloMax);
             iterator.isValid(); ++iterator) {
          autopas::sph::SPHParticle p = *iterator;
          p.addR(shift);
          sphSystem.addHaloParticle(p);
        }
      }
    }
  }
}

/**
 * deletes the halo particles
 * @param sphSystem
 */
void deleteHaloParticles(Container& sphSystem) {
  sphSystem.deleteHaloParticles();
}

void densityPressureHydroForce(Container& sphSystem) {
  // declare the used functors
  autopas::sph::SPHCalcDensityFunctor densityFunctor;
  autopas::sph::SPHCalcHydroForceFunctor hydroForceFunctor;

  std::cout << "\nhaloupdate\n" << std::endl;

  // 1.first calculate density
  // 1.1 to calculate the density we need the halo particles
  updateHaloParticles(sphSystem);

  std::cout << "haloparticles... ";
  int haloparts = 0, innerparts = 0;
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    if (not autopas::inBox(part->getR(), sphSystem.getBoxMin(),
                           sphSystem.getBoxMax())) {
      haloparts++;
    } else {
      innerparts++;
    }
  }
  std::cout << haloparts << std::endl;
  std::cout << "particles... " << innerparts << std::endl;

  // 1.2 then calculate density
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    part->setDensity(0.);
    densityFunctor.AoSFunctor(*part, *part);
    part->setDensity(part->getDensity() / 2);
  }
  std::cout << "calculation of density... started" << std::endl;
  sphSystem.iteratePairwiseAoS2(&densityFunctor);
  std::cout << "calculation of density... completed" << std::endl;
  // 1.3 delete halo particles, as their values are no longer valid
  deleteHaloParticles(sphSystem);

  // 2. then update pressure
  std::cout << "calculation of pressure... started" << std::endl;
  setPressure(sphSystem);
  std::cout << "calculation of pressure... completed" << std::endl;

  // 0.3 then calculate hydro force
  // 0.3.1 to calculate the density we need the halo particles
  updateHaloParticles(sphSystem);

  std::cout << "haloparticles... ";
  haloparts = 0, innerparts = 0;
  for (auto part = sphSystem.begin(); part.isValid(); ++part) {
    if (not autopas::inBox(part->getR(), sphSystem.getBoxMin(),
                           sphSystem.getBoxMax())) {
      haloparts++;
    } else {
      innerparts++;
    }
  }
  std::cout << haloparts << std::endl;
  std::cout << "particles... " << innerparts << std::endl;

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
  sphSystem.iteratePairwiseAoS2(&hydroForceFunctor);
  std::cout << "calculation of hydroforces... completed" << std::endl;
  // 0.3.3 delete halo particles, as their values are no longer valid
  deleteHaloParticles(sphSystem);
}

void printConservativeVariables(Container& sphSystem, MPI::Comm& comm) {
  std::array<double, 3> momSum = {0., 0., 0.};  // total momentum
  double energySum = 0.;                        // total energy
  for (auto it = sphSystem.begin(); it.isValid(); ++it) {
    momSum = autopas::arrayMath::add(
        momSum, autopas::arrayMath::mulScalar(it->getV(), it->getMass()));
    energySum += (it->getEnergy() +
                  0.5 * autopas::arrayMath::dot(it->getV(), it->getV())) *
                 it->getMass();
  }

  // MPI: global reduction
  if (comm.Get_rank() == 0) {
    comm.Reduce(MPI::IN_PLACE, &energySum, 1, MPI::DOUBLE, MPI::SUM, 0);
    comm.Reduce(MPI::IN_PLACE, momSum.data(), 3, MPI::DOUBLE, MPI::SUM, 0);
    printf("%.16e\n", energySum);
    printf("%.16e\n", momSum[0]);
    printf("%.16e\n", momSum[1]);
    printf("%.16e\n", momSum[2]);
  } else {
    comm.Reduce(&energySum, &energySum, 1, MPI::DOUBLE, MPI::SUM, 0);
    comm.Reduce(momSum.data(), momSum.data(), 3, MPI::DOUBLE, MPI::SUM, 0);
  }
}

MPI::Cartcomm getDecomposition(const std::array<double, 3> globalMin,
                               const std::array<double, 3> globalMax,
                               std::array<double, 3>& localMin,
                               std::array<double, 3>& localMax) {
  int numProcs = MPI::COMM_WORLD.Get_size();
  std::array<int, 3> gridSize{0, 0, 0};
  MPI::Compute_dims(numProcs, 3, gridSize.data());
  std::array<bool, 3> period = {true, true, true};
  auto cart =
      MPI::COMM_WORLD.Create_cart(3, gridSize.data(), period.data(), false);
  std::cout << "MPI grid dimensions: " << gridSize[0] << ", " << gridSize[1]
            << ", " << gridSize[2] << std::endl;

  int rank = cart.Get_rank();
  std::array<int, 3> coords;
  cart.Get_coords(rank, 3, coords.data());

  for (int i = 0; i < 3; ++i) {
    localMin[i] =
        coords[i] * (globalMax[i] - globalMin[i]) / gridSize[i] + globalMin[i];
    localMax[i] =
        (coords[i] + 1) * (globalMax[i] - globalMin[i]) / gridSize[i] +
        globalMin[i];
    if (coords[i] == 0) {
      localMin[i] = globalMin[i];
    } else if (coords[i] == gridSize[i] - 1) {
      localMax[i] = globalMax[i];
    }
  }

  std::cout << "MPI coordinate of current process: " << coords[0] << ", "
            << coords[1] << ", " << coords[2] << std::endl;
  return cart;
}

int main(int argc, char* argv[]) {
  MPI::Init(argc, argv);

  std::array<double, 3> globalBoxMin({0., 0., 0.}), globalBoxMax{};
  globalBoxMax[0] = 1.;
  globalBoxMax[1] = globalBoxMax[2] = globalBoxMax[0] / 8.0;
  double cutoff = 0.03;  // 0.012*2.5=0.03; where 2.5 = kernel support radius

  std::array<double, 3> localBoxMin{}, localBoxMax{};

  // get the decomposition -- get the local box of the current process from the
  // global box
  MPI::Cartcomm comm =
      getDecomposition(globalBoxMin, globalBoxMax, localBoxMin, localBoxMax);

  Container sphSystem(localBoxMin, localBoxMax, cutoff);
  double dt;
  double t_end;

  // only adds particles to the current process, due to the right localBoxMin,
  // localBoxMax
  SetupIC(sphSystem, &t_end, globalBoxMax);
  Initialize(sphSystem);

  // 0.1 ---- GET INITIAL FORCES OF SYSTEM ----
  densityPressureHydroForce(sphSystem);

  std::cout << "\n----------------------------" << std::endl;

  // 0.2 get time step
  dt = getTimeStepGlobal(sphSystem, comm);
  //---- INITIAL FORCES ARE NOW CALCULATED ----

  printConservativeVariables(sphSystem, comm);

  // 1 ---- START MAIN LOOP ----
  size_t step = 0;
  for (double time = 0.; time < t_end && step < 2; time += dt, ++step) {
    std::cout << "\n-------------------------\ntime step " << step
              << "(t = " << time << ")..." << std::endl;
    // 1.1 Leap frog: Initial Kick & Full Drift
    leapfrogInitialKick(sphSystem, dt);
    leapfrogFullDrift(sphSystem, dt);  // changes position

    // 1.2.1 positions have changed, so the container needs to be updated!
    sphSystem.updateContainer();

    // 1.2.2 adjust positions based on boundary conditions (here: periodic)
    periodicBoundaryUpdate(sphSystem, globalBoxMin, globalBoxMax);

    // 1.3 Leap frog: predict
    leapfrogPredict(sphSystem, dt);
    // 1.4 Calculate density, pressure and hydrodynamic forces
    densityPressureHydroForce(sphSystem);
    // 1.5 get time step
    dt = getTimeStepGlobal(sphSystem, comm);
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

    printConservativeVariables(sphSystem, comm);
  }

  MPI::Finalize();
}
