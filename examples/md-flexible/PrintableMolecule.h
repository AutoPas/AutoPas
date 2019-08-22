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
  PrintableMolecule() : MoleculeLJ() {}

  /**
   * Constructor.
   * @param pos Position
   * @param v Veloctiy
   * @param moleculeId Molecule ID
   * @param typeId Molecule Type ID
   */
  PrintableMolecule(std::array<double, 3> pos, std::array<double, 3> v, unsigned long moleculeId,
                    unsigned int typeId = 0)
      : MoleculeLJ(pos, v, moleculeId, typeId) {}

  /**
   * Print molecule properties to std out.
   */
  void print();
};
