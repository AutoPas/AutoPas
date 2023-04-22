/**
* @file LJMultisiteFunctorAVXTest.h
 * @author Q. Behrami
 * @date 13/04/2023
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/LJMultisiteFunctor.h"
#include "autopas/molecularDynamics/LJMultisiteFunctorAVX.h"
#include "autopas/molecularDynamics/MultisiteMoleculeLJ.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"

/**
* Test class for the LJMultisiteFunctorAVX
*/
class LJMultisiteFunctorAVXTest : public AutoPasTestBase{
 public:
  /**
   * Constructor
   */
  LJMultisiteFunctorAVXTest() = default;

  /**
   * Generates a particle properties library for use with ::generateMolecules
   * @param PPL
   */
  void generatePPL(ParticlePropertiesLibrary<double, size_t> *PPL);
  /**
   * Generates a reproducible set of molecules of varying numbers of molecules.
   * @param molecules
   * @param offset position offset
   */
  void generateMolecules(std::vector<autopas::MultisiteMoleculeLJ> *molecules, std::array<double, 3> offset);

  /**
   * Compares the correctness of the SoACellFunctor against that of the AoSFunctor.
   * @tparam newton3 enables N3L optimizations for both functors.
   * @tparam calculateGlobals Both functors calculate global attributes such as potential energy and virial which are
   * compared in addition to force and torque.
   * @tparam applyShift applies a shift to the potential energy calculation for both functors.
   * @param molecules vector of multi-site molecules
   * @param PPL Particle Properties Library
   * @param cutoff
   */
  template <bool newton3, bool calculateGlobals, bool applyShift>
  void testSoACellAgainstAoS(std::vector<autopas::MultisiteMoleculeLJ> molecules, ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

};