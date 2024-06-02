/**
 * @file Auto4OMPTest.cpp
 * @author MehdiHachicha
 * @date 02.06.2024
 * Simple loop to test Auto4OMP.
 * When compiling: -I/path/to/auto4omp/omp.h
 */

#include "omp.h"

#include <chrono>
#include <iostream>
#include <string>

int main() {
  // Plain loop:
  std::cout << "Plain loop:" << std::endl;

  // Initialize the timer [1].
  auto plainLoopStart = std::chrono::steady_clock::now();

  // Loop.
  for (int i = 1; i <= 1000000; i++) {
    std::cout << "Iteration " << std::to_string(i) << " \r" << std::flush;
  }

  // Print the time [1, 2].
  auto plainLoopEnd = std::chrono::steady_clock::now();
  long plainLoopTime = std::chrono::duration_cast<std::chrono::milliseconds>(plainLoopEnd - plainLoopStart).count();

  std::cout << "Time: " << std::to_string(plainLoopTime / 1000.0) << " s         " << std::endl;

  // OpenMP loop:
  std::cout << std::endl << "OpenMP loop:" << std::endl;

  // Print OpenMP's variables.
  char *ompNumThreads, *ompSchedule;
  if ((ompNumThreads = getenv("OMP_NUM_THREADS"))) std::cout << "OMP_NUM_THREADS: " << ompNumThreads << std::endl;
  if ((ompSchedule = getenv("OMP_SCHEDULE"))) std::cout << "OMP_SCHEDULE: " << ompSchedule << std::endl;

  // Initialize the timer [1].
  auto ompLoopStart = std::chrono::steady_clock::now();

  // Loop.
#pragma omp parallel for schedule(runtime)
  for (int i = 1; i <= 1000000; i++) {
    if (i % 100 == 0) std::cout << "Iteration " << std::to_string(i) << '\r' << std::flush;
  }

  // Print the time [1, 2].
  auto ompLoopEnd = std::chrono::steady_clock::now();
  long ompLoopTime = std::chrono::duration_cast<std::chrono::milliseconds>(ompLoopEnd - ompLoopStart).count();
  std::cout << "Time: " << std::to_string(ompLoopTime / 1000.0) << " s         " << std::endl;

  return 0;
}

/*
 * Sources:
 * [1] https://stackoverflow.com/a/18685338
 * [2] https://stackoverflow.com/a/11062846
 */