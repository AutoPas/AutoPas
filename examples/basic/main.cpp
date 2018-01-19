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

void testFullParticleCell();

void testRMMParticleCell();

void testParticleIteratorFull();

void testParticleIteratorRMM();

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
    // test FullParticleCell:
    testFullParticleCell();
    testRMMParticleCell();
    testParticleIteratorFull();
    testParticleIteratorRMM();
    testParticleContainerFull();
    testParticleContainerRMM();

    std::array<int, 3> a({1, 2, 3}), b({4, 5, 6});
    std::array<int, 3> result = arrayMath::add(a, b);
    cout << "adding {1, 2, 3} and {4, 5, 6} = {" << result[0] << ", " << result[1] << ", " << result[2] << "}" << endl;

    std::array<double, 3> ad({1.1, 2.2, 3.3}), bd({4.4, 5.5, 6.6});
    std::array<double, 3> resultd = arrayMath::add(ad, bd);
    cout << "adding {1.1, 2.2, 3.3} and {4.4, 5.5, 6.6} = {" << resultd[0] << ", " << resultd[1] << ", " << resultd[2]
         << "}" << endl;
    resultd = arrayMath::add(resultd, ad);
    cout << "adding {5.5, 7.7, 9.9} and {1.1, 2.2, 3.3} = {" << resultd[0] << ", " << resultd[1] << ", " << resultd[2]
         << "}" << endl;

    resultd = arrayMath::mul(ad, ad);
    cout << "mul {1.1, 2.2, 3.3} and {1.1, 2.2, 3.3} = {" << resultd[0] << ", " << resultd[1] << ", " << resultd[2]
         << "}" << endl;

    double rd = arrayMath::dot(ad, ad);
    cout << "dot {1.1, 2.2, 3.3} and {1.1, 2.2, 3.3} = {" << rd << "}" << endl;

    resultd = arrayMath::mulScalar(ad, 2.0);
    cout << "mulScalar {1.1, 2.2, 3.3} and 2.0 = {" << resultd[0] << ", " << resultd[1] << ", " << resultd[2] << "}"
         << endl;


    LinkedCells<MyMolecule, FullParticleCell<MyMolecule>> lc;
//	VerletLists<MyMolecule, FullParticleCell<MyMolecule>> vl;
    DirectSum<MyMolecule, FullParticleCell<MyMolecule>> dir;

    cout << "Hodor" << endl;
    return EXIT_SUCCESS;
}

void testFullParticleCell() {
    FullParticleCell<MyMolecule> fpc;
    addAFewParticles(fpc);


    SingleCellIterator<MyMolecule, FullParticleCell<MyMolecule>> it(&fpc);
    for (; it.isValid(); ++it) {
        it->print();
    }
}

void testRMMParticleCell() {
    RMMParticleCell<MyMolecule> rmmc;
    addAFewParticles(rmmc);

    cout << "===============" << endl;
    cout << "well, ids are not saved yet: " << endl;
    SingleCellIterator<MyMolecule, RMMParticleCell<MyMolecule>> it(&rmmc);
    for (; it.isValid(); ++it) {
        it->print();
    }
    cout << "===============" << endl;
}

void testParticleIteratorFull() {
    cout << " ============ testing FullCell ============ " << endl;
    std::vector<FullParticleCell<MyMolecule>> vec;
    vec.resize(1);

    cout << "one set" << endl;
    addAFewParticles(vec[0]);
    for (ParticleIterator<MyMolecule, FullParticleCell<MyMolecule>> pi(&vec);
         pi.isValid(); ++pi) {
        pi->print();
    }

    cout << "two sets" << endl;
    vec.resize(2);
    addAFewParticles(vec[1]);

    for (ParticleIterator<MyMolecule, FullParticleCell<MyMolecule>> pi(&vec);
         pi.isValid(); ++pi) {
        pi->print();
    }

    cout << "empty " << endl;
    vec.clear();
    for (ParticleIterator<MyMolecule, FullParticleCell<MyMolecule>> pi(&vec);
         pi.isValid(); ++pi) {
        pi->print();
    }

    cout << "empty full " << endl;
    vec.resize(2);
    addAFewParticles(vec[1]);
    for (ParticleIterator<MyMolecule, FullParticleCell<MyMolecule>> pi(&vec);
         pi.isValid(); ++pi) {
        pi->print();
    }

    cout << "empty full empty" << endl;
    vec.resize(3);
    for (ParticleIterator<MyMolecule, FullParticleCell<MyMolecule>> pi(&vec);
         pi.isValid(); ++pi) {
        pi->print();
    }

    cout << "empty full empty full" << endl;
    vec.resize(4);
    addAFewParticles(vec[3]);
    for (ParticleIterator<MyMolecule, FullParticleCell<MyMolecule>> pi(&vec);
         pi.isValid(); ++pi) {
        pi->print();
    }

    cout << "full empty " << endl;
    vec.clear();
    vec.resize(2);
    addAFewParticles(vec[0]);
    for (ParticleIterator<MyMolecule, FullParticleCell<MyMolecule>> pi(&vec);
         pi.isValid(); ++pi) {
        pi->print();
    }
}

void testParticleIteratorRMM() {
    cout << " ============ testing RMMCell ============ " << endl;
    std::vector<RMMParticleCell<MyMolecule>> vec;
    vec.resize(1);

    cout << "one set" << endl;
    addAFewParticles(vec[0]);
    for (ParticleIterator<MyMolecule, RMMParticleCell<MyMolecule>> pi(&vec);
         pi.isValid(); ++pi) {
        pi->print();
    }

    cout << "two sets" << endl;
    vec.resize(2);
    addAFewParticles(vec[1]);

    for (ParticleIterator<MyMolecule, RMMParticleCell<MyMolecule>> pi(&vec);
         pi.isValid(); ++pi) {
        pi->print();
    }

    cout << "empty " << endl;
    vec.clear();
    for (ParticleIterator<MyMolecule, RMMParticleCell<MyMolecule>> pi(&vec);
         pi.isValid(); ++pi) {
        pi->print();
    }

    cout << "empty full " << endl;
    vec.resize(2);
    addAFewParticles(vec[1]);
    for (ParticleIterator<MyMolecule, RMMParticleCell<MyMolecule>> pi(&vec);
         pi.isValid(); ++pi) {
        pi->print();
    }

    cout << "empty full empty" << endl;
    vec.resize(3);
    for (ParticleIterator<MyMolecule, RMMParticleCell<MyMolecule>> pi(&vec);
         pi.isValid(); ++pi) {
        pi->print();
    }

    cout << "empty full empty full" << endl;
    vec.resize(4);
    addAFewParticles(vec[3]);
    for (ParticleIterator<MyMolecule, RMMParticleCell<MyMolecule>> pi(&vec);
         pi.isValid(); ++pi) {
        pi->print();
    }

    cout << "full empty " << endl;
    vec.clear();
    vec.resize(2);
    addAFewParticles(vec[0]);
    for (ParticleIterator<MyMolecule, RMMParticleCell<MyMolecule>> pi(&vec);
         pi.isValid(); ++pi) {
        pi->print();
    }
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



