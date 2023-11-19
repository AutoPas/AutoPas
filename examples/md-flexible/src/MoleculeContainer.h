//
// Created by johnny on 14.11.23.
//

#pragma once
#include "TypeDefinitions.h"

#if defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH) and (MD_FLEXIBLE_MODE == MULTISITE)
class MoleculeContainer {
 public:
  MoleculeContainer();

  void push_back(MoleculeType && multisiteMoleculeLj);

  size_t size() const;

  MoleculeType& get(size_t index);

  [[nodiscard]] const MoleculeType & getConst(size_t moleculeID) const;

 private:
  std::vector<MoleculeType> _molecules;
};
#endif