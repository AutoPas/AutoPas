/**
* @file MoleculeLJ.cpp
* @date Originally 17/01/2018. Code pulled from header on 15/06/2023.
* @author tchipevn
*/

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

size_t Site::getTypeId() const { return _typeId; }
void Site::setTypeId(size_t typeId) { _typeId = typeId; }

std::string Site::toString() const {
 using autopas::utils::ArrayUtils::operator<<;
 std::ostringstream text;
 // clang-format off
 text << "Site"
    << "\nID                 : " << _id
    << "\nOwningMoleculeId   :" << owningMoleculeId
    << "\nIndexInsideMolecule:" << indexInsideMolecule
    << "\nPosition           : " << _r
    << "\nVelocity           : " << _v
    << "\nForce              : " << _f
    << "\nType ID            : " << _typeId
    << "\nOwnershipState     : " << _ownershipState;
 // clang-format on
 return text.str();
}
}  // namespace mdLib
