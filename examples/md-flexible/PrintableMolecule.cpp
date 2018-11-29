/**
 * @file PrintableMolecule.cpp
 * @author F. Gratl
 * @date 11/26/18
 */

#include "PrintableMolecule.h"

void PrintableMolecule::print() {
  std::cout << "Molecule with position: ";
  for (auto &r : getR()) {
    std::cout << std::setw(10) << r << ", ";
  }
  std::cout << "and force: ";

  for (auto &f : getF()) {
    std::cout << std::setw(15) << f << ", ";
  }
  std::cout << "ID: " << std::setw(5) << getID();
  std::cout << std::endl;
}
