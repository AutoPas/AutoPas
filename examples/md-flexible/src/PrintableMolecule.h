/**
 * @file PrintableMolecule.h
 * @author F. Gratl
 * @date 11/26/18
 */

#pragma once

#include <array>
#include <iomanip>
#include <iostream>

#include "autopas/molecularDynamics/MoleculeLJ.h"

/**
 * Example for a custom particle type derived from a autopas molecule type.
 */
class PrintableMolecule
    : public autopas::MoleculeLJ<> /*apparently c++17 doesnt need <> but doesnt compile without it*/ {
 public:
  /**
   * Empty Constructor.
   */
  PrintableMolecule() : autopas::MoleculeLJ<>() {}

  /**
   * Constructor.
   * @param pos Position
   * @param v Veloctiy
   * @param moleculeId Molecule ID
   * @param typeId Molecule Type ID
   */
  PrintableMolecule(std::array<double, 3> pos, std::array<double, 3> v, unsigned long moleculeId,
                    unsigned int typeId = 0)
      : autopas::MoleculeLJ<>(pos, v, moleculeId, typeId) {}

  /**
   * Print molecule properties to std out.
   */
  void print() {
    std::cout << "Molecule with position: ";
    for (auto &r : this->getR()) {
      std::cout << std::setw(10) << r << ", ";
    }
    std::cout << "and force: ";

    for (auto &f : this->getF()) {
      std::cout << std::setw(15) << f << ", ";
    }
    std::cout << "ID: " << std::setw(5) << this->getID();
    std::cout << std::endl;
  }
};
