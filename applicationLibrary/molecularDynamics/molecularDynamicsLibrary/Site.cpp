//
// Created by johnny on 14.11.23.
//

#include "Site.h"

namespace mdLib {
Site::Site(const std::array<double, 3> &pos, const std::array<double, 3> &v, unsigned long siteId,
           unsigned long typeId, unsigned long moleculeId)
: MoleculeLJ(pos, v, siteId, typeId){
  _moleculeId = moleculeId;
}

void Site::setMoleculeId(unsigned long moleculeId) {
  _moleculeId = moleculeId;
}

unsigned long Site::getMoleculeId() {
  return _moleculeId;
}
}
