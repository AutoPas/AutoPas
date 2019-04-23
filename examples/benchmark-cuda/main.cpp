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

template <class Particle, class Container>
void fillWithParticles(Container& container, std::array<size_t, 3> particlesPerDim, Particle& defaultParticle,
                       std::array<typename Particle::ParticleFloatingPointType, 3> spacing,
                       std::array<typename Particle::ParticleFloatingPointType, 3> offset) {
  size_t id = 0;
  for (unsigned int z = 0; z < particlesPerDim[2]; ++z) {
    for (unsigned int y = 0; y < particlesPerDim[1]; ++y) {
      for (unsigned int x = 0; x < particlesPerDim[0]; ++x) {
        Particle p(defaultParticle);
        p.setR({x * spacing[0] + offset[0], y * spacing[1] + offset[1], z * spacing[2] + offset[2]});
        p.setID(id++);
        container.addParticle(p);
      }
    }
  }
}

template <class Container>
void fillSpaceWithGrid(Container& pc, std::array<precision, 3> boxMin, std::array<precision, 3> boxMax,
                       precision gridsize) {
  array<size_t, 3> nParticles;
  for (int i = 0; i < 3; ++i) {
    nParticles[i] = (size_t)((boxMax[0] - boxMin[0]) / gridsize) - 1;
  }

  cout << "numParticles: " << nParticles[0] * nParticles[0] * nParticles[2] << endl;
  MyMolecule dummyParticle;
  fillWithParticles(pc, nParticles, dummyParticle, {gridsize, gridsize, gridsize},
                    {gridsize / 2, gridsize / 2, gridsize / 2});
}

template <class Container, class Traversal, class Functor>
void run(Container& pc, Traversal* t, Functor* f) {
  pc.iteratePairwise(f, t, false);

  auto start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < numSamples; ++i) {
    pc.iteratePairwise(f, t, false);
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  cout << pc.getContainerType() << ", " << setw(2) << t->getTraversalType() << ", " << setw(12) << duration << "ms"
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

  if (argc == 2) {
    numSamples = stoi(argv[1]);
  }

  std::array<precision, 3> boxMin({0., 0., 0.}), boxMax({12., 12., 12.});
  precision cutoff = 3.0;
  precision epsilon = 2.0;
  precision sigma = 0.4;
  precision gridSize = 0.3;

  DirectSum<MyMolecule, FullParticleCell<MyMolecule>> dir(boxMin, boxMax, cutoff);
  LinkedCells<MyMolecule, FullParticleCell<MyMolecule>> lc(boxMin, boxMax, cutoff);
  VerletClusterCells<MyMolecule> vcc(boxMin, boxMax, cutoff, 0.1, 1000, 64);

  fillSpaceWithGrid<>(dir, boxMin, boxMax, gridSize);
  fillSpaceWithGrid<>(lc, boxMin, boxMax, gridSize);
  fillSpaceWithGrid<>(vcc, boxMin, boxMax, gridSize);

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
  run(vcc, &traversalvcccudaN3, &func);

  return EXIT_SUCCESS;
}
