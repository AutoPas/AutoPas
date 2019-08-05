/**
 * @file PrintableMolecule.h
 * @author F. Gratl
 * @date 11/26/18
 */

#pragma once

#include <array>
#include <iomanip>
#include <iostream>
#include "autopas/particles/MoleculeLJ.h"

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
   * @param r Position
   * @param v Veloctiy
   * @param i Molecule ID
   * @param typeID Molecule Type ID
   */
  PrintableMolecule(std::array<double, 3> r, std::array<double, 3> v, unsigned long i, unsigned int typeID = 0)
      : MoleculeLJ(r, v, i, typeID) {}

  /**
   * Print molecule properties to std out.
   */
  void print();
};
