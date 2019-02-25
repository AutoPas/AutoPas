/**
 * @file main.cpp
 * @date 17.01.2018
 * @author tchipev
 */

#include <iostream>
#include <vector>
#include "autopas/autopasIncludes.h"
#include "autopas/containers/directSum/DirectSumTraversal.h"
#include "autopas/pairwiseFunctors/LJFunctor.h"
#include "autopas/utils/CudaDeviceVector.h"

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

template <class ParticleCell>
void fillSpaceWithGrid(ParticleCell &pc, std::array<double, 3> boxMin, std::array<double, 3> boxMax, double gridsize,
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

void testRun(LJFunctor<MyMolecule, FullParticleCell<MyMolecule>> &func, FullParticleCell<MyMolecule> &fpc1,
             FullParticleCell<MyMolecule> &fpc2, int num_threads = 32, bool newton3 = false) {
  cout << "NumC1: " << fpc1.numParticles() << "; NumC2: " << fpc2.numParticles() << "; threads: " << num_threads
       << "; n3: " << newton3 << "; " << endl;

  func.setCudaOptions(num_threads);
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
  cout << "->" << duration << "microseconds; (" <<soloT <<", "<<pairT<<")" <<endl;
}

int main(int argc, char **argv) {
  int numParticles = 100000;
  if (argc == 2) {
    numParticles = stoi(string(argv[1]));
  }
  autopas::Logger::create();

  std::array<double, 3> boxMin({0., 0., 0.}), boxMax({10., 10., 10.});
  double cutoff = 100.0;

  DirectSum<MyMolecule, FullParticleCell<MyMolecule>> dir(boxMin, boxMax, cutoff);

  FullParticleCell<MyMolecule> fpc1;
  FullParticleCell<MyMolecule> fpc2;

  fillSpaceWithGrid(fpc1, {0, 0, 0}, {9, 9, 9}, 0.19, numParticles);
  fillSpaceWithGrid(fpc2, {9.2, 0, 0}, {18.2, 9, 9}, 0.19, numParticles);

  typedef LJFunctor<MyMolecule, FullParticleCell<MyMolecule>> Func;
  Func func(cutoff, 1.0, 1.0, 0.0);

  vector<int> v = {32, 64, 128, 256, 512, 1024};
  for (auto it : v) {
    testRun(func, fpc1, fpc2, it, false);
  }
  for (auto it : v) {
    testRun(func, fpc1, fpc2, it, true);
  }
  return EXIT_SUCCESS;
}
