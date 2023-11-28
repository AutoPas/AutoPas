//
// Created by johnny on 14.11.23.
//

#include "Site.h"

namespace mdLib {
Site::Site(const std::array<double, 3> &pos,const std::array<double, 3> &v, unsigned long siteId,
           unsigned long moleculeId,size_t indexInsideMolecule)
: autopas::Particle(pos, v, siteId), _moleculeID(moleculeId), _indexInsideMolecule(indexInsideMolecule){
}

void Site::setMoleculeId(unsigned long moleculeId) {
  _moleculeID = moleculeId;
}

unsigned long Site::getMoleculeId() const{
  return _moleculeID;
}

size_t Site::getIndexInsideMolecule() const{
  return _indexInsideMolecule;
}

void Site::setIndexInsideMolecule(size_t indexInsideMolecule){
 _indexInsideMolecule = indexInsideMolecule;
}
unsigned long Site::getTypeId() const {
 return _typeId;
}
void Site::setTypeId(unsigned long typeId) {
 _typeId = typeId;
}
}
