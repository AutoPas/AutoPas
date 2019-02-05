/**
 * @file main.cpp
 * @date 17.01.2018
 * @author tchipev
 */

#include <iostream>
#include "autopas/autopasIncludes.h"
#include "autopas/pairwiseFunctors/LJFunctorCuda.h"
#include "autopas/pairwiseFunctors/CellFunctorCuda.h"
#include "autopas/containers/directSum/DirectSumTraversalCuda.h"
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

 private:
  int _myvar;
};


template <class ParticleCell>
void addAFewParticles(ParticleCell &pc) {
  static int i = 0;
  int iEnd = i + 4;
  for (; i < iEnd; ++i) {
    std::array<double, 3> arr({static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)});
    MyMolecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), i);
    pc.addParticle(m);
  }
}

int main() {
  autopas::Logger::create();

  std::array<double, 3> boxMin({0., 0., 0.}), boxMax({10., 10., 10.});
  double cutoff = 3.0;

  DirectSum<MyMolecule, FullParticleCell<MyMolecule>> dir(boxMin, boxMax, cutoff);
  addAFewParticles<>(dir);
  
  typedef LJFunctorCuda<MyMolecule, FullParticleCell<MyMolecule>> cudaF;
  cudaF functor(cutoff, 2.0, 3.0, 0.0);

  DirectSumTraversalCuda<FullParticleCell<MyMolecule>, cudaF, false, false> traversalAoS(&functor);
  DirectSumTraversalCuda<FullParticleCell<MyMolecule>, cudaF, true, false> traversalSoA(&functor);

  for (auto pi = dir.begin(); pi.isValid(); ++pi) {
    pi->print();
  }

  dir.iteratePairwiseAoSCuda(&functor, &traversalAoS, false);
  dir.iteratePairwiseSoACuda(&functor, &traversalSoA, false);

  for (auto pi = dir.begin(); pi.isValid(); ++pi) {
    pi->print();
  }

  autopas::CudaDeviceVector<double> cdv (32);

  cout << "Hodor" << endl;
  return EXIT_SUCCESS;
}

