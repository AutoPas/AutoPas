/**
 * @file main.cpp
 * @date 2.01.2019
 * @author jspahl
 */

#include <iostream>
#include <vector>
#include "autopas/autopasIncludes.h"
#include "autopas/containers/directSum/DirectSumTraversal.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/pairwiseFunctors/LJFunctor.h"
#include "autopas/utils/CudaDeviceVector.h"
#include "../md-flexible/ParticleClassLibrary.h"

using namespace std;
using namespace autopas;

class MyMoleculeFP64 : public ParticleFP64 {
 public:
  MyMoleculeFP64() : ParticleFP64() {}

  MyMoleculeFP64(std::array<double, 3> r, std::array<double, 3> v, unsigned long i) : ParticleFP64(r, v, i) {}

  void print() {
    cout << "Molecule with position: ";
    for (auto &r : getR()) {
      cout << r << ", ";
    }
    cout << "and force: ";

    for (auto &f : getF()) {
      cout << f << ", ";
    }
    cout << "ID: " << getID() << endl;
  }
};

template <class ParticleCell, class Particle>
void fillSpaceWithGrid(ParticleCell &pc, std::array<typename Particle::ParticleFloatingPointType, 3> boxMin,
                       std::array<typename Particle::ParticleFloatingPointType, 3> boxMax,
                       typename Particle::ParticleFloatingPointType gridsize, int maxN = 10000) {
  int i = 0;

  for (typename Particle::ParticleFloatingPointType x = boxMin[0]; x < boxMax[0]; x += gridsize) {
    for (typename Particle::ParticleFloatingPointType y = boxMin[1]; y < boxMax[1]; y += gridsize) {
      for (typename Particle::ParticleFloatingPointType z = boxMin[2]; z < boxMax[2]; z += gridsize) {
        Particle m({x, y, z}, {0., 0., 0.}, static_cast<unsigned long>(i));
        pc.addParticle(m);
        if (++i >= maxN) {
          return;
        }
      }
    }
  }
}

template <typename Particle>
void testRun(LJFunctor<Particle, FullParticleCell<Particle>> &func, FullParticleCell<Particle> &fpc1,
             FullParticleCell<Particle> &fpc2, int num_threads = 32, bool newton3 = false) {
  cout << "NumC1: " << fpc1.numParticles() << "; NumC2: " << fpc2.numParticles() << "; threads: " << num_threads
       << "; n3: " << newton3 << "; " << endl;

  func.getCudaWrapper()->setNumThreads(num_threads);
  func.SoALoader(fpc1, fpc1._particleSoABuffer);
  func.SoALoader(fpc2, fpc2._particleSoABuffer);
  func.deviceSoALoader(fpc1._particleSoABuffer, fpc1._particleSoABufferDevice);
  func.deviceSoALoader(fpc2._particleSoABuffer, fpc2._particleSoABufferDevice);

  auto start = std::chrono::high_resolution_clock::now();

  func.CudaFunctor(fpc1._particleSoABufferDevice, newton3);
  func.deviceSoAExtractor(fpc1._particleSoABuffer, fpc1._particleSoABufferDevice);
  auto mid = std::chrono::high_resolution_clock::now();
  func.CudaFunctor(fpc1._particleSoABufferDevice, fpc2._particleSoABufferDevice, newton3);

  func.deviceSoAExtractor(fpc2._particleSoABuffer, fpc2._particleSoABufferDevice);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  auto soloT = std::chrono::duration_cast<std::chrono::microseconds>(mid - start).count();
  auto pairT = std::chrono::duration_cast<std::chrono::microseconds>(stop - mid).count();
  cout << "->" << duration << " microseconds; (" << soloT << ", " << pairT << ")" << endl;
}

template <typename Particle>
void run(int numParticles) {
  autopas::Logger::create();
  map<unsigned long, double> universalMap;
  for (int i = 0; i < numParticles; i++) {
    universalMap.emplace((unsigned long)i, 1.0);
  }
  ParticleClassLibrary PCL = ParticleClassLibrary(universalMap, universalMap, universalMap);
  std::array<typename Particle::ParticleFloatingPointType, 3> boxMin({0., 0., 0.}), boxMax({10., 10., 10.});
  typename Particle::ParticleFloatingPointType cutoff = 100.0;

  FullParticleCell<Particle> fpc1;
  FullParticleCell<Particle> fpc2;

  fillSpaceWithGrid<FullParticleCell<Particle>, Particle>(fpc1, {0, 0, 0}, {9, 9, 9}, 0.19, numParticles);
  fillSpaceWithGrid<FullParticleCell<Particle>, Particle>(fpc2, {9.2, 0, 0}, {18.2, 9, 9}, 0.19, numParticles);

  typedef LJFunctor<Particle, FullParticleCell<Particle>> Func;
  Func func(cutoff, PCL, 0.0);

  // number of threads used per block by the device
  vector<int> v = {32, 64, 128, 256, 512, 1024};
  for (auto it : v) {
    testRun(func, fpc1, fpc2, it, false);
  }
  for (auto it : v) {
    testRun(func, fpc1, fpc2, it, true);
  }
}

int main(int argc, char **argv) {
  int numParticles = 8000;
  if (argc == 2) {
    numParticles = stoi(string(argv[1]));
  }

  cout << "Test double precision" << endl;
  run<ParticleFP64>(numParticles);

  cout << "Test single precision" << endl;
  run<ParticleFP32>(numParticles);

  return EXIT_SUCCESS;
}
