/**
* @file PositionStoringMultiSiteMolecule.h
* @date 10/10/2023
* @author Johannes Riemenschneider
 */

#include "PositionStoringMultiSiteMolecule.h"

#include "autopas/utils/Quaternion.h"

// #include "autopas/utils/Quaternion.h" //needed to update Sites automatically

namespace mdLib {
PositionStoringMultiSiteMolecule::PositionStoringMultiSiteMolecule(std::array<double, 3> r, std::array<double, 3> v,
                                                                   std::array<double, 3> angularVel, unsigned long moleculeId,
                                                                   unsigned long typeId)
    : mdLib::MoleculeLJ(r, v, moleculeId, typeId), _angularVel(angularVel), _torque({0., 0., 0.}) {}

const std::array<double, 4> &PositionStoringMultiSiteMolecule::getQuaternion() const { return _q; }

//void PositionStoringMultiSiteMolecule::setQuaternion(const std::array<double, 4> &q, const std::vector<std::array<double, 3>>& unrotated_site_positions) {
void PositionStoringMultiSiteMolecule::setQuaternion(const std::array<double, 4> &q, const ParticlePropertiesLibraryType& ppl){
  _q = q;
  //set sites based on rotation
  const std::vector<std::array<double, 3>> rotated_site_positions = autopas::utils::quaternion::rotateVectorOfPositions(_q, ppl.getSitePositions(_typeId));
  _relSitePositionsX.resize(rotated_site_positions.size());
  _relSitePositionsY.resize(rotated_site_positions.size());
  _relSitePositionsZ.resize(rotated_site_positions.size());

  //const std::vector<std::array<double, 3>> rotated_site_positions = unrotated_site_positions;
  for(int i=0; i < rotated_site_positions.size(); i++){
    _relSitePositionsX[i] = rotated_site_positions[i][0];
    _relSitePositionsY[i] = rotated_site_positions[i][1];
    _relSitePositionsZ[i] = rotated_site_positions[i][2];
  }
}

void PositionStoringMultiSiteMolecule::setOnlyQuaternion(const std::array<double, 4> &q){
  _q = q;
}

const std::vector<double> &PositionStoringMultiSiteMolecule::getRelativeSitePositionsX() const {
  return _relSitePositionsX;
}

const std::vector<double> &PositionStoringMultiSiteMolecule::getRelativeSitePositionsY() const {
  return _relSitePositionsY;
}

const std::vector<double> &PositionStoringMultiSiteMolecule::getRelativeSitePositionsZ() const {
  return _relSitePositionsZ;
}

std::array<double, 3> PositionStoringMultiSiteMolecule::getRelativeSitePosition(size_t index) const {
  return {_relSitePositionsX[index], _relSitePositionsY[index], _relSitePositionsZ[index]};
}

int PositionStoringMultiSiteMolecule::getNumberOfSites(){
  return _relSitePositionsX.size();
}

const std::array<double, 3> &PositionStoringMultiSiteMolecule::getAngularVel() const { return _angularVel; }
void PositionStoringMultiSiteMolecule::setAngularVel(const std::array<double, 3> &angularVel) { _angularVel = angularVel; }
void PositionStoringMultiSiteMolecule::addAngularVel(const std::array<double, 3> &angularVel) {
  _angularVel = autopas::utils::ArrayMath::add(_angularVel, angularVel);
}

const std::array<double, 3> &PositionStoringMultiSiteMolecule::getTorque() const { return _torque; }
void PositionStoringMultiSiteMolecule::setTorque(const std::array<double, 3> &torque) { _torque = torque; }
void PositionStoringMultiSiteMolecule::addTorque(const std::array<double, 3> &torque) {
  _torque = autopas::utils::ArrayMath::add(_torque, torque);
}
void PositionStoringMultiSiteMolecule::subTorque(const std::array<double, 3> &torque) {
  _torque = autopas::utils::ArrayMath::sub(_torque, torque);
}

std::string PositionStoringMultiSiteMolecule::toString() const {
  using autopas::utils::ArrayUtils::operator<<;
  std::ostringstream text;
  // clang-format off
     text << "PositionStoringMultiSiteMolecule"
        << "\nID                 : " << _id
        << "\nPosition           : " << _r
        << "\nVelocity           : " << _v
        << "\nForce              : " << _f
        << "\nOld Force          : " << _oldF
        //<< "\nQuaternion         : " << _q
        << "\nSite Positions X   :" << vectorToString<double>(_relSitePositionsX)
        << "\nSite Positions Y>  :" << vectorToString<double>(_relSitePositionsY)
        << "\nSite Positions  Z  :" << vectorToString<double>(_relSitePositionsZ)
        << "\nRotational Velocity: " << _angularVel
        << "\nTorque             : " << _torque
        << "\nType ID            : " << _typeId
        << "\nOwnershipState     : " << _ownershipState;
  // clang-format on
  return text.str();
}
template <typename T>
std::string PositionStoringMultiSiteMolecule::vectorToString(std::vector<T> v) const{
  using autopas::utils::ArrayUtils::operator<<;
  std::ostringstream text;
  text << "[";
  for(int i=0; i < v.size(); i++){
    text << v[i];
    if(i+1 < v.size()){
      text << ", ";
    }
  }
  return text.str();
}
}  // namespace mdLib