/**
* @file AbsoluteMultiSiteMoleculeLJ.h
* @date 10/10/2023
* @author Johannes Riemenschneider
 */

#include "AbsoluteMultiSiteMoleculeLJ.h"

namespace mdLib {
AbsoluteMultiSiteMoleculeLJ::AbsoluteMultiSiteMoleculeLJ(std::array<double, 3> r, std::array<double, 3> v,
                                        std::array<double, 3> angularVel, unsigned long moleculeId,
                                        unsigned long typeId)
   : mdLib::MoleculeLJ(r, v, moleculeId, typeId), _angularVel(angularVel), _torque({0., 0., 0.}) {}

//const std::array<double, 4> &AbsoluteMultiSiteMoleculeLJ::getQuaternion() const { return _q; }
//void AbsoluteMultiSiteMoleculeLJ::setQuaternion(const std::array<double, 4> &q) { _q = q; }

const std::vector<double> &AbsoluteMultiSiteMoleculeLJ::getAbsoluteSitePositionsX() const {
  return _absSitePositionsX;
}

const std::vector<double> &AbsoluteMultiSiteMoleculeLJ::getAbsoluteSitePositionsY() const {
  return _absSitePositionsY;
}

const std::vector<double> &AbsoluteMultiSiteMoleculeLJ::getAbsoluteSitePositionsZ() const {
  return _absSitePositionsZ;
}

const std::array<double, 3> AbsoluteMultiSiteMoleculeLJ::getAbsoluteSitePosition(size_t index) const {
  return {_absSitePositionsX[index], _absSitePositionsY[index], _absSitePositionsZ[index]};
}

int AbsoluteMultiSiteMoleculeLJ::getNumberOfSites(){
  return _absSitePositionsX.size();
}

const std::array<double, 3> &AbsoluteMultiSiteMoleculeLJ::getAngularVel() const { return _angularVel; }
void AbsoluteMultiSiteMoleculeLJ::setAngularVel(const std::array<double, 3> &angularVel) { _angularVel = angularVel; }
void AbsoluteMultiSiteMoleculeLJ::addAngularVel(const std::array<double, 3> &angularVel) {
 _angularVel = autopas::utils::ArrayMath::add(_angularVel, angularVel);
}

const std::array<double, 3> &AbsoluteMultiSiteMoleculeLJ::getTorque() const { return _torque; }
void AbsoluteMultiSiteMoleculeLJ::setTorque(const std::array<double, 3> &torque) { _torque = torque; }
void AbsoluteMultiSiteMoleculeLJ::addTorque(const std::array<double, 3> &torque) {
 _torque = autopas::utils::ArrayMath::add(_torque, torque);
}
void AbsoluteMultiSiteMoleculeLJ::subTorque(const std::array<double, 3> &torque) {
 _torque = autopas::utils::ArrayMath::sub(_torque, torque);
}

std::string AbsoluteMultiSiteMoleculeLJ::toString() const {
 using autopas::utils::ArrayUtils::operator<<;
 std::ostringstream text;
 // clang-format off
     text << "AbsoluteMultiSiteMoleculeLJ"
        << "\nID                 : " << _id
        << "\nPosition           : " << _r
        << "\nVelocity           : " << _v
        << "\nForce              : " << _f
        << "\nOld Force          : " << _oldF
        //<< "\nQuaternion         : " << _q
        << "\nSite Positions X   :" << vectorToString<double>(_absSitePositionsX)
        << "\nSite Positions Y>  :" << vectorToString<double>(_absSitePositionsY)
        << "\nSite Positions  Z  :" << vectorToString<double>(_absSitePositionsZ)
        << "\nRotational Velocity: " << _angularVel
        << "\nTorque             : " << _torque
        << "\nType ID            : " << _typeId
        << "\nOwnershipState     : " << _ownershipState;
 // clang-format on
 return text.str();
}
template <typename T>
std::string AbsoluteMultiSiteMoleculeLJ::vectorToString(std::vector<T> v) const{
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