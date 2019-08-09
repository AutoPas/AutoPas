/**
 * @file main.cpp
 * @date 10.02.2019
 * @author jspahl
 */

#include <chrono>
#include <iostream>
#include "autopas/autopasIncludes.h"
#include "autopas/containers/directSum/DirectSumTraversal.h"
#include "autopas/containers/linkedCells/traversals/C01CudaTraversal.h"
#include "autopas/pairwiseFunctors/LJFunctor.h"

using namespace std;
using namespace autopas;

class MyMolecule : public Particle {
 public:
  MyMolecule() : Particle(), _myvar(0) {}

  MyMolecule(std::array<double, 3> r, std::array<double, 3> v, unsigned long i, int myvar)
      : Particle(r, v, i), _myvar(myvar) {}

  void print() {
    cout << "Molecule with position: ";
    for (auto &r : getR()) {
      cout << r << ", ";
    }
    cout << "and force: ";

    for (auto &f : getF()) {
      cout << f << ", ";
    }
    cout << "ID: " << getID();
    cout << " myvar: " << _myvar << endl;
  }

  // typedef autopas::utils::SoAType<size_t, float, float, float, float, float, float>::Type SoAArraysType;
  // typedef autopas::utils::CudaSoAType<size_t, float, float, float, float, float, float>::Type CudaDeviceArraysType;

 private:
  int _myvar;
};

template <class Container>
void fillSpaceWithGrid(Container &pc, std::array<double, 3> boxMin, std::array<double, 3> boxMax, double gridsize,
                       int maxN = 10000) {
  int i = 0;

  for (double x = boxMin[0]; x < boxMax[0]; x += gridsize) {
    for (double y = boxMin[1]; y < boxMax[1]; y += gridsize) {
      for (double z = boxMin[2]; z < boxMax[2]; z += gridsize) {
        std::array<double, 3> arr({x, y, z});
        MyMolecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), i);
        pc.addParticle(m);
        if (++i >= maxN) {
          return;
        }
      }
    }
  }
}

int main(int argc, char **argv) {
  autopas::Logger::create();
  int maxIterations = 5;
  long numParticles = 10000;

  std::array<double, 3> boxMin({0., 0., 0.}), boxMax({33., 33., 33.});
  double cutoff = 3.0;
  double skin = 0.;
  double epsilon = 2.0;
  double sigma = 0.4;

  DirectSum<MyMolecule, FullParticleCell<MyMolecule>> dir(boxMin, boxMax, cutoff, skin);
  LinkedCells<MyMolecule, FullParticleCell<MyMolecule>> lc(boxMin, boxMax, cutoff, skin);

  fillSpaceWithGrid<>(dir, boxMin, boxMax, 0.8, numParticles);
  fillSpaceWithGrid<>(lc, boxMin, boxMax, 0.8, numParticles);

  typedef LJFunctor<MyMolecule, FullParticleCell<MyMolecule>> Func;

  Func func(cutoff, epsilon, sigma, 0.0);

  DirectSumTraversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::aos, false> traversalAoS(&func, cutoff);
  DirectSumTraversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::soa, false> traversalSoA(&func, cutoff);
  DirectSumTraversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::cuda, false> traversalCuda(&func, cutoff);
  DirectSumTraversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::cuda, true> traversalCudaN3(&func, cutoff);

  C08Traversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::soa, true> traversalc08N3(
      lc.getCellBlock().getCellsPerDimensionWithHalo(), &func, cutoff + skin, lc.getCellBlock().getCellLength());
  C01Traversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::cuda, false> C01Cuda(
      lc.getCellBlock().getCellsPerDimensionWithHalo(), &func, cutoff + skin, lc.getCellBlock().getCellLength());
  C01CudaTraversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::cuda, false> traversalLCcuda(
      lc.getCellBlock().getCellsPerDimensionWithHalo(), &func);
  C01CudaTraversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::cuda, true> traversalLCcudaN3(
      lc.getCellBlock().getCellsPerDimensionWithHalo(), &func);

  dir.iteratePairwise(&traversalAoS);
  lc.iteratePairwise(&traversalc08N3);

  auto start = std::chrono::high_resolution_clock::now();
  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<int64_t, micro>::rep duration;

  start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < maxIterations; ++i) {
    dir.iteratePairwise(&traversalAoS);
  }
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  cout << "DsAoS: " << maxIterations << " iterations with " << dir.getNumParticles() << " particles took: " << duration
       << " microseconds" << endl;

  start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < maxIterations; ++i) {
    dir.iteratePairwise(&traversalSoA);
  }
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  cout << "DsSoA: " << maxIterations << " iterations with " << dir.getNumParticles() << " particles took: " << duration
       << " microseconds" << endl;

  start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < maxIterations; ++i) {
    dir.iteratePairwise(&traversalCuda);
  }
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  cout << "DsCuda:" << maxIterations << " iterations with " << dir.getNumParticles() << " particles took: " << duration
       << " microseconds" << endl;

  start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < maxIterations; ++i) {
    dir.iteratePairwise(&traversalCudaN3);
  }
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  cout << "DsCudaN3:" << maxIterations << " iterations with " << dir.getNumParticles()
       << " particles took: " << duration << " microseconds" << endl;

  start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < maxIterations; ++i) {
    lc.iteratePairwise(&traversalc08N3);
  }
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  cout << "LcSoAN3:" << maxIterations << " iterations with " << lc.getNumParticles() << " particles took: " << duration
       << " microseconds" << endl;

  start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < maxIterations; ++i) {
    lc.iteratePairwise(&C01Cuda);
  }
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  cout << "LcCudaN3:" << maxIterations << " iterations with " << lc.getNumParticles() << " particles took: " << duration
       << " microseconds" << endl;

  start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < maxIterations; ++i) {
    lc.iteratePairwise(&traversalLCcuda);
  }
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  cout << "LcTraversalCudaNoN3:" << maxIterations << " iterations with " << lc.getNumParticles()
       << " particles took: " << duration << " microseconds" << endl;

  start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < maxIterations; ++i) {
    lc.iteratePairwise(&traversalLCcudaN3);
  }
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  cout << "LcTraversalCudaN3:" << maxIterations << " iterations with " << lc.getNumParticles()
       << " particles took: " << duration << " microseconds" << endl;
  return EXIT_SUCCESS;
}
