/**
 * @file sph-main-mpi.cpp
 * @date 10.04.2018
 * @author seckler
 */

#include <mpi.h>

#include <array>
#include <cmath>
#include <iostream>

#include "autopas/AutoPas.h"
#include "autopas/sph/autopassph.h"

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
        if (autopas::utils::inBox(ith.getR(), sphSystem.getBoxMin(), sphSystem.getBoxMax())) {
          sphSystem.addParticle(ith);
        }
      }
    }
  }
  for (double x = bBoxMax[0] * 0.5; x < bBoxMax[0] * 1.; x += dx * 2.0) {  // NOLINT
    for (double y = 0; y < bBoxMax[1]; y += dx) {                          // NOLINT
      for (double z = 0; z < bBoxMax[2]; z += dx) {                        // NOLINT
        Particle ith({x, y, z}, {0, 0, 0}, i++, 0.75, 0.012, 0.);
        ith.setDensity(0.5);
        ith.setEnergy(2.5);
        if (autopas::utils::inBox(ith.getR(), sphSystem.getBoxMin(), sphSystem.getBoxMax())) {
          sphSystem.addParticle(ith);
        }
      }
    }
  }
  if (i == 0) {
    throw std::runtime_error("No particle added in sph-main-mpi::SetupIC.");
  }

  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->setMass(part->getMass() * bBoxMax[0] * bBoxMax[1] * bBoxMax[2] / (double)(i));
  }

  // we are incrementing i independent of whether we add a particle to the local
  // sphsystem, or not. So it is still the global number of particles.
  std::cout << "total # of particles is... " << i << std::endl;

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

double getTimeStepGlobal(AutoPasContainer &sphSystem, MPI_Comm &comm) {
  double dt = 1.0e+30;  // set VERY LARGE VALUE
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->calcDt();
    if (part->getDt() < 0.002) {
      std::cout << "small time step for particle " << part->getID() << " at [" << part->getR()[0] << ", "
                << part->getR()[1] << ", " << part->getR()[2] << "]" << std::endl;
    }
    dt = std::min(dt, part->getDt());
  }

  // MPI, global reduction of minimum
  MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, comm);

  std::cout << "the time step dt is..." << dt << std::endl;
  return dt;
}

void leapfrogInitialKick(AutoPasContainer &sphSystem, const double dt) {
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->setVel_half(autopas::utils::ArrayMath::add(
        part->getV(), autopas::utils::ArrayMath::mulScalar(part->getAcceleration(), 0.5 * dt)));
    part->setEng_half(part->getEnergy() + 0.5 * dt * part->getEngDot());
  }
}

void leapfrogFullDrift(AutoPasContainer &sphSystem, const double dt) {
  // time becomes t + dt;
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->addR(autopas::utils::ArrayMath::mulScalar(part->getVel_half(), dt));
  }
}

void leapfrogPredict(AutoPasContainer &sphSystem, const double dt) {
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->addV(autopas::utils::ArrayMath::mulScalar(part->getAcceleration(), dt));
    part->addEnergy(part->getEngDot() * dt);
  }
}

void leapfrogFinalKick(AutoPasContainer &sphSystem, const double dt) {
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->setV(autopas::utils::ArrayMath::add(part->getVel_half(),
                                              autopas::utils::ArrayMath::mulScalar(part->getAcceleration(), 0.5 * dt)));
    part->setEnergy(part->getEng_half() + 0.5 * dt * part->getEngDot());
  }
}

void setPressure(AutoPasContainer &sphSystem) {
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->calcPressure();
  }
}

int getSendRecvPartner(const std::array<int, 3> diff, MPI_Comm &comm, bool recvPartner) {
  int neighbor;
  std::array<int, 3> mycoords{0, 0, 0};
  std::array<int, 3> neighborcoords{0, 0, 0};
  std::array<int, 3> dimSize{0, 0, 0}, dummy{0, 0, 0};
  MPI_Cart_get(comm, 3, dimSize.data(), dummy.data(), mycoords.data());
  for (int i = 0; i < 3; ++i) {
    if (recvPartner) {
      neighborcoords[i] = (mycoords[i] + dimSize[i] - diff[i]) % dimSize[i];
    } else {
      neighborcoords[i] = (mycoords[i] + diff[i]) % dimSize[i];
    }
  }
  MPI_Cart_rank(comm, neighborcoords.data(), &neighbor);
  return neighbor;
}

int getReceivePartner(const std::array<int, 3> diff, MPI_Comm &comm) { return getSendRecvPartner(diff, comm, true); }

void issueSend(std::vector<Particle> &sendParticles, const std::array<int, 3> diff, MPI_Comm &comm,
               MPI_Request &sendRequest, std::vector<double> &buffer) {
  int neighbor = getSendRecvPartner(diff, comm, false);

  for (auto &p : sendParticles) {
    std::vector<double> serialized = p.serialize();
    buffer.insert(std::end(buffer), std::begin(serialized), std::end(serialized));
  }
  MPI_Isend(buffer.data(), buffer.size(), MPI_DOUBLE, neighbor, 3, comm, &sendRequest);
}

void receive(std::vector<Particle> &receiveParticles, const std::array<int, 3> diff, MPI_Comm &comm) {
  int neighbor = getReceivePartner(diff, comm);

  int length;
  MPI_Status status;
  MPI_Probe(neighbor, 3, comm, &status);
  MPI_Get_count(&status, MPI_DOUBLE, &length);
  std::vector<double> recvBuffer(length);
  MPI_Recv(recvBuffer.data(), length, MPI_DOUBLE, neighbor, 3, comm, MPI_STATUS_IGNORE);

  for (size_t i = 0; i < (size_t)length;) {
    auto p = Particle::deserialize(recvBuffer.data(), i);
    receiveParticles.push_back(p);
  }
}

void waitSend(MPI_Request &sendRequest) { MPI_Wait(&sendRequest, MPI_STATUS_IGNORE); }

/**
 * Get the haloregion of particles we need to send to the neighbour with the
 * specific diff.
 * Needed for the RegionParticleIterator.
 * @param boxMin
 * @param boxMax
 * @param diff
 * @param sendMin
 * @param sendMax
 * @param cutoff
 * @param shift
 * @param globalBoxMin
 * @param globalBoxMax
 */
void getSendHalo(double boxMin, double boxMax, int diff, double &sendMin, double &sendMax, double cutoff, double &shift,
                 const double globalBoxMin, const double globalBoxMax) {
  if (diff == 0) {
    sendMin = boxMin;
    sendMax = boxMax;
    shift = 0;
  } else if (diff == -1) {
    sendMin = boxMin;
    sendMax = boxMin + cutoff;
    if (boxMin == globalBoxMin) {
      shift = globalBoxMax - globalBoxMin;
    } else {
      shift = 0;
    }
  } else if (diff == 1) {
    sendMin = boxMax - cutoff;
    sendMax = boxMax;
    if (boxMax == globalBoxMax) {
      shift = globalBoxMin - globalBoxMax;
    } else {
      shift = 0;
    }
  }
}

/**
 * Get the region of leaving particles we need to send to the neighbour with the
 * specific diff.
 * Needed for the RegionParticleIterator.
 * @param boxMin
 * @param boxMax
 * @param diff
 * @param sendMin output
 * @param sendMax output
 * @param cutoff
 * @param shift output
 * @param globalBoxMin
 * @param globalBoxMax
 */
void getSendLeaving(double boxMin, double boxMax, int diff, double &sendMin, double &sendMax, double cutoff,
                    double &shift, const double globalBoxMin, const double globalBoxMax) {
  if (diff == 0) {
    sendMin = boxMin;
    sendMax = boxMax;
    shift = 0;
  } else if (diff == -1) {
    sendMin = boxMin - cutoff;
    sendMax = boxMin;
    if (boxMin == globalBoxMin) {
      shift = globalBoxMax - globalBoxMin;
    } else {
      shift = 0;
    }
  } else if (diff == 1) {
    sendMin = boxMax;
    sendMax = boxMax + cutoff;
    if (boxMax == globalBoxMax) {
      shift = globalBoxMin - globalBoxMax;
    } else {
      shift = 0;
    }
  }
}

/**
 * updates the halo particles
 * this is done by sending the boundary particles to othere processes + shifting
 * them around for periodic boundaries
 * @param sphSystem
 * @param comm
 * @param globalBoxMin
 * @param globalBoxMax
 */
void updateHaloParticles(AutoPasContainer &sphSystem, MPI_Comm &comm, const std::array<double, 3> &globalBoxMin,
                         const std::array<double, 3> &globalBoxMax) {
  std::array<double, 3> boxMin = sphSystem.getBoxMin();
  std::array<double, 3> boxMax = sphSystem.getBoxMax();
  std::array<double, 3> requiredHaloMin{0, 0, 0}, requiredHaloMax{0, 0, 0};
  std::array<int, 3> diff{0, 0, 0};
  std::array<double, 3> shift{0, 0, 0};
  double cutoff = sphSystem.getCutoff();
  double skin= sphSystem.getVerletSkin();
  std::vector<double> buffer;
  for (diff[0] = -1; diff[0] < 2; diff[0]++) {
    for (diff[1] = -1; diff[1] < 2; diff[1]++) {
      for (diff[2] = -1; diff[2] < 2; diff[2]++) {
        std::vector<Particle> sendParticles;
        if (not diff[0] and not diff[1] and not diff[2]) {
          // at least one dimension has to be non-zero
          continue;
        }
        // figure out which particles we send
        for (int i = 0; i < 3; ++i) {
          getSendHalo(boxMin[i], boxMax[i], diff[i], requiredHaloMin[i], requiredHaloMax[i], cutoff, shift[i],
                      globalBoxMin[i], globalBoxMax[i]);
        }

        for (auto iterator =
                 sphSystem.getRegionIterator(requiredHaloMin, requiredHaloMax, autopas::IteratorBehavior::owned);
             iterator.isValid(); ++iterator) {
          Particle p = *iterator;  // copies Particle
          p.addR(shift);
          sendParticles.push_back(p);
        }
        MPI_Request sendRequest;
        issueSend(sendParticles, diff, comm, sendRequest, buffer);

        std::vector<Particle> receiveParticles;
        receive(receiveParticles, diff, comm);
        for (auto &particle : receiveParticles) {
          sphSystem.addHaloParticle(particle);
        }
        waitSend(sendRequest);
        buffer.clear();
      }
    }
  }
}

void periodicBoundaryUpdate(AutoPasContainer &sphSystem, MPI_Comm &comm, const std::vector<Particle> &invalidParticles,
                            std::array<double, 3> globalBoxMin, std::array<double, 3> globalBoxMax) {
  std::array<double, 3> boxMin = sphSystem.getBoxMin();
  std::array<double, 3> boxMax = sphSystem.getBoxMax();
  std::array<double, 3> requiredHaloMin{0, 0, 0}, requiredHaloMax{0, 0, 0};
  std::array<int, 3> diff{0, 0, 0};
  std::array<double, 3> shift{0, 0, 0};
  double cutoff = sphSystem.getCutoff();

  std::vector<double> buffer;
  for (diff[0] = -1; diff[0] < 2; diff[0]++) {
    for (diff[1] = -1; diff[1] < 2; diff[1]++) {
      for (diff[2] = -1; diff[2] < 2; diff[2]++) {
        std::vector<Particle> sendParticles;
        if (not diff[0] and not diff[1] and not diff[2]) {
          // at least one dimension has to be non-zero
          continue;
        }
        // figure out which particles we send
        for (int i = 0; i < 3; ++i) {
          getSendLeaving(boxMin[i], boxMax[i], diff[i], requiredHaloMin[i], requiredHaloMax[i], cutoff, shift[i],
                         globalBoxMin[i], globalBoxMax[i]);
        }
        for (const auto &p : invalidParticles) {
          if (autopas::utils::inBox(p.getR(), requiredHaloMin, requiredHaloMax)) {
            // std::cout << "sending particle at (" << p.getR()[0] << ", "
            //          << p.getR()[1] << ", " << p.getR()[2] << ")" << std::endl;
            auto pCopy = p;
            pCopy.addR(shift);
            sendParticles.push_back(pCopy);
          }
        }
        MPI_Request sendRequest;
        issueSend(sendParticles, diff, comm, sendRequest, buffer);

        std::vector<Particle> receiveParticles;
        receive(receiveParticles, diff, comm);
        for (auto &particle : receiveParticles) {
          for (int i = 0; i < 3; i++) {
            auto r = particle.getR();
            if (r[i] == sphSystem.getBoxMax()[i]) {  // required if -1e20 + 3 =
                                                     // 3 and not < 3 //
                                                     // floating point errors
              r[i] = nextafter(r[i], 0);
              particle.setR(r);
            }
          }
          sphSystem.addParticle(particle);
        }
        waitSend(sendRequest);
        buffer.clear();
      }
    }
  }
}

void densityPressureHydroForce(AutoPasContainer &sphSystem, MPI_Comm &comm, const std::array<double, 3> &globalBoxMin,
                               const std::array<double, 3> &globalBoxMax) {
  // declare the used functors
  autopas::sph::SPHCalcDensityFunctor<Particle> densityFunctor;
  autopas::sph::SPHCalcHydroForceFunctor<Particle> hydroForceFunctor;

  // 1.first calculate density
  // 1.1 to calculate the density we need the halo particles
  updateHaloParticles(sphSystem, comm, globalBoxMin, globalBoxMax);

  // 1.2 then calculate density
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->setDensity(0.);
    densityFunctor.AoSFunctor(*part, *part);
    part->setDensity(part->getDensity() / 2);
  }

  sphSystem.iteratePairwise(&densityFunctor);

  // 1.3 delete halo particles, as their values are no longer valid
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::halo); part.isValid(); ++part) {
    sphSystem.deleteParticle(part);
  }

  // 2. then update pressure
  setPressure(sphSystem);

  // 3 then calculate hydro force
  // 3.1 to calculate the density we need the halo particles
  updateHaloParticles(sphSystem, comm, globalBoxMin, globalBoxMax);

  // 3.2 then calculate hydro force
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    // self interaction leeds to:
    // 1) vsigmax = 2*part->getSoundSpeed()
    // 2) no change in acceleration
    part->setVSigMax(2 * part->getSoundSpeed());
    part->setAcceleration(std::array<double, 3>{0., 0., 0.});
    part->setEngDot(0.);
  }

  sphSystem.iteratePairwise(&hydroForceFunctor);
}

void printConservativeVariables(AutoPasContainer &sphSystem, MPI_Comm &comm) {
  std::array<double, 3> momSum = {0., 0., 0.};  // total momentum
  double energySum = 0.;                        // total energy
  for (auto it = sphSystem.begin(autopas::IteratorBehavior::owned); it.isValid(); ++it) {
    momSum = autopas::utils::ArrayMath::add(momSum, autopas::utils::ArrayMath::mulScalar(it->getV(), it->getMass()));
    energySum += (it->getEnergy() + 0.5 * autopas::utils::ArrayMath::dot(it->getV(), it->getV())) * it->getMass();
  }

  // MPI: global reduction
  int myrank;
  MPI_Comm_rank(comm, &myrank);
  if (myrank == 0) {
    MPI_Reduce(MPI_IN_PLACE, &energySum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(MPI_IN_PLACE, momSum.data(), 3, MPI_DOUBLE, MPI_SUM, 0, comm);
    printf("%.16e\n", energySum);
    for (int i = 0; i < 3; ++i) {
      printf("%.16e\n", momSum[i]);
      if (std::abs(momSum[i]) > 1.e-15) {
        throw std::runtime_error("ERROR: bad moment sum detected (should be small, but isn't!");
      }
    }
  } else {
    MPI_Reduce(&energySum, &energySum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(momSum.data(), momSum.data(), 3, MPI_DOUBLE, MPI_SUM, 0, comm);
  }
}

MPI_Comm getDecomposition(const std::array<double, 3> globalMin, const std::array<double, 3> globalMax,
                          std::array<double, 3> &localMin, std::array<double, 3> &localMax) {
  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  std::array<int, 3> gridSize{0, 0, 0};
  MPI_Dims_create(numProcs, 3, gridSize.data());
  std::array<int, 3> period = {1, 1, 1};
  MPI_Comm cart;
  MPI_Cart_create(MPI_COMM_WORLD, 3, gridSize.data(), period.data(), false, &cart);
  std::cout << "MPI grid dimensions: " << gridSize[0] << ", " << gridSize[1] << ", " << gridSize[2] << std::endl;

  int rank;
  MPI_Comm_rank(cart, &rank);
  std::array<int, 3> coords{};
  MPI_Cart_coords(cart, rank, 3, coords.data());

  for (int i = 0; i < 3; ++i) {
    localMin[i] = coords[i] * (globalMax[i] - globalMin[i]) / gridSize[i] + globalMin[i];
    localMax[i] = (coords[i] + 1) * (globalMax[i] - globalMin[i]) / gridSize[i] + globalMin[i];
    if (coords[i] == 0) {
      localMin[i] = globalMin[i];
    } else if (coords[i] == gridSize[i] - 1) {
      localMax[i] = globalMax[i];
    }
  }

  std::cout << "MPI coordinate of current process: " << coords[0] << ", " << coords[1] << ", " << coords[2]
            << std::endl;
  return cart;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  autopas::Logger::create();

  std::array<double, 3> globalBoxMin({0., 0., 0.}), globalBoxMax{};
  globalBoxMax[0] = 1.;
  globalBoxMax[1] = globalBoxMax[2] = globalBoxMax[0] / 8.0;
  double cutoff = 0.03;               // 0.012*2.5=0.03; where 2.5 = kernel support radius
  unsigned int rebuildFrequency = 6;  // has to be multiple of 2
  double skinToCutoffRatio = 0.1;

  std::array<double, 3> localBoxMin{}, localBoxMax{};

  // get the decomposition -- get the local box of the current process from the
  // global box
  MPI_Comm comm = getDecomposition(globalBoxMin, globalBoxMax, localBoxMin, localBoxMax);

  AutoPasContainer sphSystem;
  sphSystem.setNumSamples(
      6);  // has to be multiple of 2, should also be multiple of rebuildFrequency (but this is not necessary)
  sphSystem.setBoxMin(localBoxMin);
  sphSystem.setBoxMax(localBoxMax);
  sphSystem.setCutoff(cutoff);
  sphSystem.setVerletSkinPerTimestep(skinToCutoffRatio * cutoff);
  sphSystem.setVerletRebuildFrequency(rebuildFrequency);

  std::set<autopas::ContainerOption> allowedContainers{autopas::ContainerOption::linkedCells,
                                                       autopas::ContainerOption::verletLists,
                                                       autopas::ContainerOption::verletListsCells};
  sphSystem.setAllowedContainers(allowedContainers);

  int rank;
  MPI_Comm_rank(comm, &rank);

  sphSystem.setOutputSuffix("Rank" + std::to_string(rank) + "_");
  sphSystem.init();

  double dt;
  double t_end;

  // only adds particles to the current process, due to the right localBoxMin,
  // localBoxMax
  SetupIC(sphSystem, &t_end, globalBoxMax);
  Initialize(sphSystem);

  // 0.1 ---- GET INITIAL FORCES OF SYSTEM ----
  densityPressureHydroForce(sphSystem, comm, globalBoxMin, globalBoxMax);

  std::cout << "\n----------------------------" << std::endl;

  // 0.2 get time step
  dt = getTimeStepGlobal(sphSystem, comm);
  //---- INITIAL FORCES ARE NOW CALCULATED ----

  printConservativeVariables(sphSystem, comm);

  // 1 ---- START MAIN LOOP ----
  size_t step = 0;
  for (double time = 0.; time < t_end && step < 55; time += dt, ++step) {
    if (rank == 0) {
      std::cout << "\n-------------------------\ntime step " << step << "(t = " << time << ")..." << std::endl;
    }
    // 1.1 Leap frog: Initial Kick & Full Drift
    leapfrogInitialKick(sphSystem, dt);
    leapfrogFullDrift(sphSystem, dt);  // changes position

    // 1.2.1 positions have changed, so the container needs to be updated!
    auto invalidParticles = sphSystem.updateContainer();

    // 1.2.2 adjust positions based on boundary conditions (here: periodic)
    periodicBoundaryUpdate(sphSystem, comm, invalidParticles, globalBoxMin, globalBoxMax);

    // 1.3 Leap frog: predict
    leapfrogPredict(sphSystem, dt);
    // 1.4 Calculate density, pressure and hydrodynamic forces
    densityPressureHydroForce(sphSystem, comm, globalBoxMin, globalBoxMax);
    // 1.5 get time step
    dt = getTimeStepGlobal(sphSystem, comm);
    // 1.6 Leap frog: final Kick
    leapfrogFinalKick(sphSystem, dt);

    //    for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
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
  std::cout << "-----------------\nfinished" << std::endl;
  sphSystem.finalize();
  MPI_Comm_free(&comm);
  MPI_Finalize();
}
