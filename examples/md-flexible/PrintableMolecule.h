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
template <typename floatType>
class PrintableMoleculeBase : public autopas::MoleculeLJ {
 public:
  PrintableMoleculeBase() : autopas::MoleculeLJ() {}

  PrintableMoleculeBase(std::array<double, 3> r, std::array<double, 3> v, unsigned long i)
      : autopas::MoleculeLJ(r, v, i) {}

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

  /**
   * Override the type for the soa storage.
   */
  typedef typename autopas::utils::SoAType<unsigned long, floatType, floatType, floatType, floatType, floatType,
                                           floatType, floatType>::Type SoAArraysType;

#if defined(AUTOPAS_CUDA)
  /**
   * Override the type for the soa storage.
   */
  typedef typename autopas::utils::CudaSoAType<unsigned long, floatType, floatType, floatType, floatType, floatType,
                                               floatType, floatType>::Type CudaDeviceArraysType;
#else
  /**
   * Override the type for the soa storage.
   * empty if compiled without Cuda Support
   */
  using CudaDeviceArraysType = typename autopas::utils::CudaSoAType<>::Type;
#endif
};
