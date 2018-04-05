/*
 * main.cpp
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#include <iostream>
#include "autopas.h"

using namespace std;
using namespace autopas;

class MyMolecule : public Particle {
 public:
  MyMolecule() : Particle(), _myvar(0) {}

  MyMolecule(std::array<double, 3> r, std::array<double, 3> v, unsigned long i,
             int myvar)
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

// how do we recommend to use this thing?
// the user should create a wrapper around our ParticleContainer? That way they
// can compile whatever they actually need w.r.t. templates once.
// with a bit of documentation, they can

void testParticleContainerFull();

void testParticleContainerRMM();

template <class ParticleCell>
void addAFewParticles(ParticleCell &pc) {
  static int i = 0;
  int iEnd = i + 4;
  for (; i < iEnd; ++i) {
    std::array<double, 3> arr({static_cast<double>(i), static_cast<double>(i),
                               static_cast<double>(i)});
    MyMolecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), i);
    pc.addParticle(m);
  }
}

int main() {
  testParticleContainerFull();
  testParticleContainerRMM();

  std::array<double, 3> boxMin({0., 0., 0.}), boxMax({10., 10., 10.});
  double cutoff = 1.0;

  LinkedCells<MyMolecule, FullParticleCell<MyMolecule>> lc(boxMin, boxMax,
                                                           cutoff);
  //	VerletLists<MyMolecule, FullParticleCell<MyMolecule>> vl;
  DirectSum<MyMolecule, FullParticleCell<MyMolecule>> dir(boxMin, boxMax,
                                                          cutoff);

  cout << "Hodor" << endl;
  return EXIT_SUCCESS;
}

void testParticleContainerFull() {
  cout << " =========== testing ParticleContainerFull ============" << endl;
  std::array<double, 3> boxMin({0., 0., 0.}), boxMax({10., 10., 10.});
  double cutoff = 1.0;

  LinkedCells<MyMolecule, FullParticleCell<MyMolecule>> pc(boxMin, boxMax,
                                                           cutoff);
  //pc.init();  // empty full empty full empty

  // add a few particles to the second cell
  for (int i = 0; i < 4; ++i) {
    std::array<double, 3> arr({static_cast<double>(i), static_cast<double>(i),
                               static_cast<double>(i)});
    MyMolecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), i);
    pc.addParticle(m);
  }
  // add a few particles to the fourth cell
  for (int i = 4; i < 8; ++i) {
    std::array<double, 3> arr({static_cast<double>(i), static_cast<double>(i),
                               static_cast<double>(i)});
    MyMolecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), i);
    pc.addParticle(m);
  }

  for (auto pi = pc.begin(); pi.isValid(); ++pi) {
    pi->print();
  }

  cout << " =========== done testing ParticleContainerFull ============"
       << endl;
}

void testParticleContainerRMM() {
  cout << " =========== testing ParticleContainerRMM ============" << endl;
  std::array<double, 3> boxMin({0., 0., 0.}), boxMax({10., 10., 10.});
  double cutoff = 1.0;
  LinkedCells<MyMolecule, RMMParticleCell<MyMolecule>> pc(boxMin, boxMax,
                                                          cutoff);
  //pc.init();  // empty full empty full empty

  // add a few particles to the second cell
  for (int i = 0; i < 4; ++i) {
    std::array<double, 3> arr({static_cast<double>(i), static_cast<double>(i),
                               static_cast<double>(i)});
    MyMolecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), i);
    pc.addParticle(m);
  }
  // add a few particles to the fourth cell
  for (int i = 4; i < 8; ++i) {
    std::array<double, 3> arr({static_cast<double>(i), static_cast<double>(i),
                               static_cast<double>(i)});
    MyMolecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), i);
    pc.addParticle(m);
  }

  for (auto pi = pc.begin(); pi.isValid(); ++pi) {
    pi->print();
  }

  cout << " =========== done testing ParticleContainerRMM ============" << endl;
}
