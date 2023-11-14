//
// Created by johnny on 14.11.23.
//

#pragma once
#include "TypeDefinitions.h"

#if defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH) and (MD_FLEXIBLE_MODE == MULTISITE)
class MoleculeContainer {
 public:
  MoleculeContainer();

  void add(MoleculeType && multisiteMoleculeLj);

  size_t size();

 private:
  std::vector<MoleculeType> _molecules;
};
#endif