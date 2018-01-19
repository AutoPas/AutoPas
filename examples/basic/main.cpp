/*
 * main.cpp
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */


#include "autopas.h"
#include <iostream>

using namespace std;
using namespace autopas;

class MyMolecule : public Particle {
public:
    MyMolecule() : Particle(), _myvar(0) {}

    MyMolecule(std::array<double, 3> r, std::array<double, 3> v, unsigned long i, int myvar) : Particle(r, v, i),
                                                                                               _myvar(myvar) {}

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
// the user should create a wrapper around our ParticleContainer? That way they can compile whatever they actually need w.r.t. templates once.
// with a bit of documentation, they can

void testParticleContainerFull();

void testParticleContainerRMM();

template<class ParticleCell>
void addAFewParticles(ParticleCell &pc) {
    static int i = 0;
    int iEnd = i + 4;
    for (; i < iEnd; ++i) {
        std::array<double, 3> arr({static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)});
        MyMolecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), i);
        pc.addParticle(m);
    }
}

int main(void) {
    testParticleContainerFull();
    testParticleContainerRMM();

    LinkedCells<MyMolecule, FullParticleCell<MyMolecule>> lc;
//	VerletLists<MyMolecule, FullParticleCell<MyMolecule>> vl;
    DirectSum<MyMolecule, FullParticleCell<MyMolecule>> dir;

    cout << "Hodor" << endl;
    return EXIT_SUCCESS;
}

void testParticleContainerFull() {
    cout << " =========== testing ParticleContainerFull ============" << endl;
    LinkedCells<MyMolecule, FullParticleCell<MyMolecule>> pc;
    pc.init(); // empty full empty full empty

    // add a few particles to the second cell
    for (int i = 0; i < 4; ++i) {
        std::array<double, 3> arr({static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)});
        MyMolecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), i);
        pc.addParticle(m, 1);
    }
    // add a few particles to the fourth cell
    for (int i = 4; i < 8; ++i) {
        std::array<double, 3> arr({static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)});
        MyMolecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), i);
        pc.addParticle(m, 3);
    }

    for (auto pi = pc.begin(); pi.isValid(); ++pi) {
        pi->print();
    }

    cout << " =========== done testing ParticleContainerFull ============" << endl;
}

void testParticleContainerRMM() {
    cout << " =========== testing ParticleContainerRMM ============" << endl;
    LinkedCells<MyMolecule, RMMParticleCell<MyMolecule>> pc;
    pc.init(); // empty full empty full empty

    // add a few particles to the second cell
    for (int i = 0; i < 4; ++i) {
        std::array<double, 3> arr({static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)});
        MyMolecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), i);
        pc.addParticle(m, 1);
    }
    // add a few particles to the fourth cell
    for (int i = 4; i < 8; ++i) {
        std::array<double, 3> arr({static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)});
        MyMolecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), i);
        pc.addParticle(m, 3);
    }

    for (auto pi = pc.begin(); pi.isValid(); ++pi) {
        pi->print();
    }

    cout << " =========== done testing ParticleContainerRMM ============" << endl;
}



