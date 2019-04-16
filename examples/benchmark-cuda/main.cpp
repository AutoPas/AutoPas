/**
 * @file main.cpp
 * @date 10.02.2019
 * @author jspahl
 */

#include <chrono>
#include <iomanip>
#include <iostream>
#include "autopas/autopasIncludes.h"
#include "autopas/containers/directSum/DirectSumTraversal.h"
#include "autopas/containers/linkedCells/traversals/C01CudaTraversal.h"
#include "autopas/containers/verletClusterLists/VerletClusterCells.h"
#include "autopas/pairwiseFunctors/LJFunctor.h"

using namespace std;
using namespace autopas;

typedef ParticleFP32 MyMolecule;
typedef MyMolecule::ParticleFloatingPointType precision;

int numSamples = 100;

template <class Container>
void fillSpaceWithGrid(Container& pc, std::array<precision, 3> boxMin, std::array<precision, 3> boxMax,
                       precision gridsize, int maxN = 10000) {
  int i = 0;

  for (precision x = boxMin[0]; x < boxMax[0]; x += gridsize) {
    for (precision y = boxMin[1]; y < boxMax[1]; y += gridsize) {
      for (precision z = boxMin[2]; z < boxMax[2]; z += gridsize) {
        std::array<precision, 3> arr({x, y, z});
        MyMolecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i));
        pc.addParticle(m);
        if (++i >= maxN) {
          return;
        }
      }
    }
  }
}

template <class Container, class Traversal, class Functor>
void run(Container& pc, Traversal* t, Functor* f) {
  auto start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < numSamples; ++i) {
    pc.iteratePairwise(f, t, false);
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  cout << pc.getContainerType() << ", " << setw(2) << t->getTraversalType() << ", " << setw(8) << duration << "ms"
       << endl;
}

template <class Container, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
void build(Container& pc, PairwiseFunctor* pairwiseFunctor, TraversalOption traversalType) {
  auto traversalSelector = pc.generateTraversalSelector();

  // auto traversal = pc.generateTraversalSelector().generateTraversal<PairwiseFunctor, DataLayout,
  // useNewton3>(traversalType, pairwiseFunctor); run(pc, traversal, pairwiseFunctor);
}

int main(int argc, char** argv) {
  autopas::Logger::create();
  long numParticles = 200000;

  if (argc == 2) {
    numSamples = stoi(argv[1]);
  }

  std::array<precision, 3> boxMin({0., 0., 0.}), boxMax({33., 33., 33.});
  precision cutoff = 1.0;
  precision epsilon = 2.0;
  precision sigma = 0.4;

  DirectSum<MyMolecule, FullParticleCell<MyMolecule>> dir(boxMin, boxMax, cutoff);
  LinkedCells<MyMolecule, FullParticleCell<MyMolecule>> lc(boxMin, boxMax, cutoff);
  VerletClusterCells<MyMolecule> vcc(boxMin, boxMax, cutoff, 0.1, 101, 256);
  VerletClusterCells<MyMolecule> vcc2(boxMin, boxMax, cutoff, 0.1, 101, 1024);

  fillSpaceWithGrid<>(dir, boxMin, boxMax, 0.3, numParticles);
  fillSpaceWithGrid<>(lc, boxMin, boxMax, 0.3, numParticles);
  fillSpaceWithGrid<>(vcc, boxMin, boxMax, 0.3, numParticles);
  fillSpaceWithGrid<>(vcc2, boxMin, boxMax, 0.3, numParticles);

  typedef LJFunctor<MyMolecule, FullParticleCell<MyMolecule>> Func;

  Func func(cutoff, epsilon, sigma, 0.0);

  DirectSumTraversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::aos, false> traversalAoS(&func);
  DirectSumTraversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::soa, false> traversalSoA(&func);
  DirectSumTraversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::cuda, false> traversalCuda(&func);
  DirectSumTraversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::cuda, true> traversalCudaN3(&func);

  C08Traversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::soa, true> traversalc08N3(
      lc.getCellBlock().getCellsPerDimensionWithHalo(), &func);
  C01Traversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::cuda, false> C01Cuda(
      lc.getCellBlock().getCellsPerDimensionWithHalo(), &func);
  C01CudaTraversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::cuda, false> traversalLCcuda(
      lc.getCellBlock().getCellsPerDimensionWithHalo(), &func);
  C01CudaTraversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::cuda, true> traversalLCcudaN3(
      lc.getCellBlock().getCellsPerDimensionWithHalo(), &func);

  VerletClusterCellsTraversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::cuda, false> traversalvcccudaNoN3(
      &func);
  VerletClusterCellsTraversal<FullParticleCell<MyMolecule>, Func, DataLayoutOption::cuda, true> traversalvcccudaN3(
      &func);

  // build<LinkedCells<MyMolecule, FullParticleCell<MyMolecule>>, Func, DataLayoutOption::soa, true>(lc, &func,
  // TraversalOption::c18);
  run(lc, &traversalc08N3, &func);
  run(lc, &traversalLCcuda, &func);
  run(lc, &traversalLCcudaN3, &func);

  run(vcc, &traversalvcccudaNoN3, &func);
  run(vcc2, &traversalvcccudaNoN3, &func);

  return EXIT_SUCCESS;
}
