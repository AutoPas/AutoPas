/**
 * @file MoleculeConversionHelper.h
 * @date 15/03/2022
 * @author Samuel Newcome
*/

#pragma once

#include "MulticenteredMoleculeLJ.h"

/**
 * Helper function to convert a rotational molecule to a simple molecule.
 * @param complexMolecule molecule with rotational components
 * @return molecule without rotational components
 */
autopas::MoleculeLJ returnSimpleMolecule(MulticenteredMoleculeLJ complexMolecule) {
  autopas::MoleculeLJ simpleMolecule;
  simpleMolecule.setR(complexMolecule.getR());
  simpleMolecule.setV(complexMolecule.getV());
  simpleMolecule.setID(complexMolecule.getID());
  simpleMolecule.setTypeId(complexMolecule.getTypeId());
  return simpleMolecule;
}