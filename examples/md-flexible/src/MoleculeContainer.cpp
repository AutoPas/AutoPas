//
// Created by johnny on 14.11.23.
//

#include "MoleculeContainer.h"

#if defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH)
MoleculeContainer::MoleculeContainer():_molecules{} {}

void MoleculeContainer::push_back(MoleculeType &&multisiteMoleculeLj) {
  _molecules.push_back(multisiteMoleculeLj);
}

size_t MoleculeContainer::size() const {
  return _molecules.size();
}
MoleculeType &MoleculeContainer::get(size_t moleculeID) {
  //assert moleculeID == _molecules[moleculeID].getID();
  return _molecules[moleculeID];
}

const MoleculeType &MoleculeContainer::getConst(size_t moleculeID) const {
  return _molecules[moleculeID];
}
#endif
