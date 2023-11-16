//
// Created by johnny on 14.11.23.
//

#include "MoleculeContainer.h"

#if defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH)
MoleculeContainer::MoleculeContainer():_molecules{} {}

void MoleculeContainer::add(mdLib::MultisiteMoleculeLJ &&multisiteMoleculeLj) {
  _molecules.push_back(multisiteMoleculeLj);
}

size_t MoleculeContainer::size() {
  return _molecules.size();
}
MoleculeType &MoleculeContainer::get(size_t moleculeID) {
  //assert moleculeID == _molecules[moleculeID].getID();
  return _molecules[moleculeID];
}
#endif
