/**
 * @file main.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include "mpi.h"

#include <chrono>
#include <thread>
#include <sys/types.h>
#include <unistd.h>

#include "simulations/MDFlexMPI.h"
#include "simulations/MDFlexSingleRank.h"

/**
 * The main function for md-flexible.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv) {
  int gdbBlocker = 0;

  pid_t pid = getpid();
  std::cout << "ProcessId: " << pid << std::endl;

  using namespace std::chrono_literals;

  while (gdbBlocker == 0) {
    std::this_thread::sleep_for(5000ms);
  }

  MPI_Init(&argc, &argv);
  MDFlexMPI simulation(3, argc, argv);
  simulation.run();

  MPI_Finalize();
  return EXIT_SUCCESS;
}
