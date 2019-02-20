/**
 * @file main.cpp
 * @date 17.01.2018
 * @author tchipev
 */

#include <iostream>
#include "autopas/autopasIncludes.h"
#include "autopas/containers/directSum/DirectSumTraversal.h"
#include "autopas/pairwiseFunctors/LJFunctor.h"
#include "autopas/pairwiseFunctors/LJFunctorCuda.h"
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

  typedef LJFunctor<MyMolecule, FullParticleCell<MyMolecule>> Func;
  Func func(cutoff, 1.0, 1.0, 0.0);

  DirectSumTraversal<FullParticleCell<MyMolecule>, Func, true, false, true> traversal(&func);

  for (auto pi = dir.begin(); pi.isValid(); ++pi) {
    pi->print();
  }
  cout << endl;

  dir.iteratePairwiseSoACuda(&func, &traversal, false);
  for (auto pi = dir.begin(); pi.isValid(); ++pi) {
    pi->print();
  }
  cout << endl;

  autopas::utils::CudaDeviceVector<double> cdv(32);
  std::vector<double> data = {3, 0, 0, 0, 0, 0, 1.9, 0, 0, 0, 2, 0};
  std::vector<double> res(12, -1);

  loadConstants(9., 24., 1.);
  cdv.copyHostToDevice(12, data.data());

  cdv.copyDeviceToHost(12, res.data());

  for (auto it : data) {
    std::cout << it << ", ";
  }
  cout << endl;
  for (auto it : res) {
    std::cout << it << ", ";
  }
  cout << endl;

  cout << "Hodor" << endl;
  return EXIT_SUCCESS;
}
