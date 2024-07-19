/**
 * @file MoleculeLJ.cpp
 * @date Originally 17/01/2018. Code pulled from header on 15/06/2023.
 * @author tchipevn
 */

#include "MoleculeLJ_NoPPL.h"

namespace mdLib {
MoleculeLJ_NoPPL::MoleculeLJ_NoPPL(const std::array<double, 3> &pos, const std::array<double, 3> &v,
                                   unsigned long moleculeId, double squareRootEpsilon, double sigmaDiv2)
    : autopas::Particle(pos, v, moleculeId), _squareRootEpsilon{squareRootEpsilon}, _sigmaDiv2{sigmaDiv2} {}

const std::array<double, 3> &MoleculeLJ_NoPPL::getOldF() const { return _oldF; }
void MoleculeLJ_NoPPL::setOldF(const std::array<double, 3> &oldForce) { _oldF = oldForce; }

const double &MoleculeLJ_NoPPL::getSquareRootEpsilon() const { return _squareRootEpsilon; }
void MoleculeLJ_NoPPL::setSquareRootEpsilon(const double &squareRootEpsilon) { _squareRootEpsilon = squareRootEpsilon; }

const double &MoleculeLJ_NoPPL::getSigmaDiv2() const { return _sigmaDiv2; }
void MoleculeLJ_NoPPL::setSigmaDiv2(const double &sigmaDiv2) { _sigmaDiv2 = sigmaDiv2; }

double MoleculeLJ_NoPPL::getEpsilon() const { return _squareRootEpsilon * _squareRootEpsilon; }
void MoleculeLJ_NoPPL::setEpsilon(const double &epsilon) { _squareRootEpsilon = sqrt(epsilon); }

double MoleculeLJ_NoPPL::getSigma() const { return _sigmaDiv2 * 2; }
void MoleculeLJ_NoPPL::setSigma(const double &sigma) { _sigmaDiv2 = sigma / 2; }

const double &MoleculeLJ_NoPPL::getMass() const { return _mass; }
void MoleculeLJ_NoPPL::setMass(const double &mass) { _mass = mass; }

std::string MoleculeLJ_NoPPL::toString() const {
  using autopas::utils::ArrayUtils::operator<<;
  std::ostringstream text;
  // clang-format off
  text << "MoleculeLJ"
     << "\nID                 : " << _id
     << "\nPosition           : " << _r
     << "\nVelocity           : " << _v
     << "\nForce              : " << _f
     << "\nOld Force          : " << _oldF
     << "\nSqrt(Epsilon)      : " << _squareRootEpsilon
     << "\nSigma/2            : " << _sigmaDiv2
     << "\nMass               : " << _mass
     << "\nOwnershipState     : " << _ownershipState;
  // clang-format on
  return text.str();
}
}  // namespace mdLib
