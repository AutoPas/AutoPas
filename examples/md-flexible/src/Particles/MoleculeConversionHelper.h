/**
 * @file MoleculeConversionHelper.h
 * @date 15/03/2022
 * @author Samuel Newcome
*/

#pragma once

#include "MulticenteredMoleculeLJ.h"

/**
 * Helper function to convert between types of molecules.
 * Primary purpose of this function is to reduce rotational molecules to non-rotational during the initialisation of
 * md-flexible.
 * @param complexMolecule molecule with at least the components of MoleculeLJ
 * @return molecule without additional components
 */
template <class ParticleClass>
autopas::MoleculeLJ returnSimpleMolecule(ParticleClass complexMolecule) {
  autopas::MoleculeLJ simpleMolecule;
  simpleMolecule.setR(complexMolecule.getR());
  simpleMolecule.setV(complexMolecule.getV());
  simpleMolecule.setID(complexMolecule.getID());
  simpleMolecule.setTypeId(complexMolecule.getTypeId());
  return simpleMolecule;
}