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
class PrintableMolecule : public autopas::MoleculeLJ {
 public:
  PrintableMolecule() : MoleculeLJ() {}

  PrintableMolecule(std::array<double, 3> r, std::array<double, 3> v, unsigned long i) : MoleculeLJ(r, v, i) {}

  /**
   * Print molecule properties to std out.
   */
  void print();
};
