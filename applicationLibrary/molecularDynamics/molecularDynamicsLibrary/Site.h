//
// Created by johnny on 14.11.23.
//

#pragma once
#include "MoleculeLJ.h"

namespace mdLib {
class Site : MoleculeLJ {
 public:
 /**
 * Constructor initializing Site completely
 * @param pos Absolute position of the site
 * @param v Velocity of the site
 * @param siteId Id uniquely identifying the site
 * @param typeId TypeId of the site.
 * @param moleculeId Identifier of the molecule this site is a part of. Sites are in a many-to-one relationship with Molecules.
   */
  Site(const std::array<double, 3> &pos, const std::array<double, 3> &v, unsigned long siteId,
       unsigned long typeId, unsigned long moleculeId);

  void setMoleculeId(unsigned long moleculeId);

  unsigned long getMoleculeId();

 private:
  /**
   * The moleculeId identifies the MultiSite Molecule this Site belongs to. It corresponds to the moleculeId of MultiSiteMoleculeLJ
   */
  unsigned long _moleculeId;
};
}