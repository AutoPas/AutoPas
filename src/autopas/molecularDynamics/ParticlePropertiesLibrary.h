/**
 * @file ParticlePropertiesLibrary.h
 * @author N. Fottner
 * @date 7/4/19
 */

#pragma once

#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ExceptionHandler.h"

/**
 * This class stores the (physical) properties of molecule types, and, in the case of multi-site molecules, the location
 * of the sites relative to the the centre-of-mass.
 *
 * It also provides mixed values for (force) calculations between known types.
 */
template <typename floatType = double, typename intType = unsigned long>
class ParticlePropertiesLibrary {
 public:
  /**
   * Constructor
   * @param cutoff Cutoff for the Potential
   */
  explicit ParticlePropertiesLibrary(const double cutoff) : _cutoff(cutoff) {}

  /**
   * Copy Constructor.
   * @param particlePropertiesLibrary
   */
  ParticlePropertiesLibrary(const ParticlePropertiesLibrary &particlePropertiesLibrary) = default;

  /**
   * Copy assignment operator.
   * @param particlePropertiesLibrary
   * @return
   */
  ParticlePropertiesLibrary &operator=(const ParticlePropertiesLibrary &particlePropertiesLibrary) = default;

  /**
   * Adds the properties of a type of a LJ single-site type to the library.
   *
   * This function also precomputes all possible mixed values with already known particle types.
   * If the type id already exists the values will be overwritten.
   * @param molId
   * @param epsilon
   * @param sigma
   * @param mass
   */
  void addSimpleType(const intType molId, const floatType epsilon, const floatType sigma, const floatType mass);

  /**
   * Adds the properties of a type of a single LJ site type to the library.
   *
   * This function also precomputes all possible mixed values with already known particle types.
   * If the type id already exists the values will be overwritten.
   * @param siteId
   * @param epsilon
   * @param sigma
   * @param mass
   */
  void addSiteType(const intType siteId, const floatType epsilon, const floatType sigma, const floatType mass);

  /**
   * Adds the properties of a molecule type to the library including: position and type of all sites; and moment of
   * inertia.
   *
   * If the type id already exists the values will be overwritten.
   * @param molId
   * @param siteIds vector of IDs of sites
   * @param relPos vector of relative positions
   * @param momentOfInertia diagonalised moment of inertia as a array of 3 floats
   */
  void addMolType(const intType molId, const std::vector<intType> siteIds, const std::vector<std::array<floatType,3>> relPos, const std::array<floatType, 3> momentOfInertia);

  /**
   * Calculates the actual mixing coefficients.
   */
  void calculateMixingCoefficients();

  /**
   * Calculates the moment of inertia from the masses and positions of the site, then calculates a rotation such that
   * the moment of inertia in the rotated axes is diagonal, and applies the rotation to the sites.
   */
  void calculateMomentOfInertiaAndAdjustAxes();

  ~ParticlePropertiesLibrary() = default;

  /**
   * Returns a set of all particle types stored.
   * todo this is super ugly code
   * @return
   */
  std::set<intType> getSiteTypes() const {
    std::set<intType> typeIDs;

    for (size_t index = 0; index < _numRegisteredSiteTypes; ++index) {
      typeIDs.emplace(index);
    }

    return typeIDs;
  }

  std::set<intType> getMolTypes() const {
    std::set<intType> typeIDs;

    for (size_t index = 0; index < _numRegisteredMolTypes; ++index) {
      typeIDs.emplace(index);
    }

    return typeIDs;
  }

  std::set<intType> getTypes() const {
    return getSiteTypes();
  }


  /**
   * Getter for the site's epsilon*24.
   * @param i siteId of the site.
   * @return 24*epsilon_i
   */
  floatType get24Epsilon(intType i) const;

  /**
   * Getter for the site's squared sigma.
   * @param i siteId of the site.
   * @return sigma_i²
   */
  floatType getSigmaSquare(intType i) const;

  /**
   * Getter for the site's mass.
   * @param i siteId of molecule.
   * @return mass_i
   */
  floatType getSiteMass(intType i) const;

  /**
   * Getter for the mol's mass.
   * @param i molId of molecule.
   * @return mass_i
   */
  floatType getMolMass(intType i) const;

  /**
   * Getter for the molecule's MoI.
   * @param i molId of molecule.
   * @return moment of inertia
   */
  std::array<floatType,3> getMomentOfInertia(intType i) const;

  /**
   * Get relative site positions.
   * @param i molId of the molecule.
   * @return site positions
   */
  std::vector<std::array<floatType,3>> getSitePositions(intType i) const;

  /**
   * Get site types of the molecule.
   * @param i molId
   * @return
   */
  std::vector<intType> getSiteTypes(intType i) const;

  /**
   * Get number of sites in the molecule.
   * @param i molId
   * @return
   */
  intType getNumSites(intType i) const;

  /**
   * Returns the precomputed mixed epsilon24.
   * @param  i siteId of site one.
   * @param  j siteId of site two.
   * @return 24*epsilon_ij
   */
  inline floatType mixing24Epsilon(intType i, intType j) const {
    return _computedMixingData[i * _numRegisteredSiteTypes + j].epsilon24;
  }

  /**
   * Get complete mixing data for one pair of site types.
   * @param i siteId of site one.
   * @param j siteId of site two.
   * @return
   */
  inline auto getMixingData(intType i, intType j) const { return _computedMixingData[i * _numRegisteredSiteTypes + j]; }

  /**
   * Returns precomputed mixed squared sigma.
   * @param i siteId of site one.
   * @param j siteId of site two.
   * @return sigma_ij²
   */
  inline floatType mixingSigmaSquare(intType i, intType j) const {
    return _computedMixingData[i * _numRegisteredSiteTypes + j].sigmaSquare;
  }

  /**
   * Returns precomputed mixed shift 6.
   * @param i siteId of site one.
   * @param j siteId of site two.
   * @return shift * 6
   */
  inline floatType mixingShift6(intType i, intType j) const {
    return _computedMixingData[i * _numRegisteredSiteTypes + j].shift6;
  }

  /**
   * Calculate the shift of the lennard jones potential from given cutoff, epsilon, sigma.
   * @param epsilon24
   * @param sigmaSquare
   * @param cutoffSquared squared cutoff of the potential that should be shifted
   * @return shift multiplied by 6
   */
  static double calcShift6(double epsilon24, double sigmaSquare, double cutoffSquared);

 private:
  intType _numRegisteredSiteTypes{0};
  intType _numRegisteredMolTypes{0};
  const double _cutoff;

  std::vector<floatType> _epsilons;
  std::vector<floatType> _sigmas;
  std::vector<floatType> _siteMasses;

  // Note: this is a vector of site type Ids for the sites of a certain molecular Id
  std::vector<std::vector<intType>> _siteIds;
  // This is a vector (indexed by mol ID) of vectors of site positions (which are 3D arrays)
  std::vector<std::vector<std::array<floatType,3>>> _relativeSitePositions;
  std::vector<floatType> _molMasses;
  std::vector<std::array<floatType,3>> _momentOfInertias;
  std::vector<size_t> _numSites;

  struct PackedMixingData {
    floatType epsilon24;
    floatType sigmaSquare;
    floatType shift6;
  };

  std::vector<PackedMixingData, autopas::AlignedAllocator<PackedMixingData>> _computedMixingData;
};

template <typename floatType, typename intType>
void ParticlePropertiesLibrary<floatType, intType>::addSimpleType(const intType molId, const floatType epsilon, const floatType sigma, const floatType mass) {
  if (_numRegisteredSiteTypes != molId) {
    autopas::utils::ExceptionHandler::exception(
        "ParticlePropertiesLibrary::addSimpleType(): trying to register a type with id {}. Please register types "
        "consecutively, starting at id 0. Currently there are {} registered types.",
        molId, _numRegisteredSiteTypes);
  }
  ++_numRegisteredSiteTypes;
  _epsilons.emplace_back(epsilon);
  _sigmas.emplace_back(sigma);
  _siteMasses.emplace_back(mass);
}

template <typename floatType, typename intType>
void ParticlePropertiesLibrary<floatType, intType>::addSiteType(intType siteID, floatType epsilon, floatType sigma, floatType mass) {
  if (_numRegisteredSiteTypes != siteID) {
    autopas::utils::ExceptionHandler::exception(
        "ParticlePropertiesLibrary::addSiteType(): trying to register a site type with id {}. Please register types "
        "consecutively, starting at id 0. Currently there are {} registered types.",
        siteID, _numRegisteredSiteTypes);
  }
  ++_numRegisteredSiteTypes;
  _epsilons.emplace_back(epsilon);
  _sigmas.emplace_back(sigma);
  _siteMasses.emplace_back(mass);
}

template <typename floatType, typename intType>
void ParticlePropertiesLibrary<floatType, intType>::addMolType(const intType molId, const std::vector<intType> siteIds,
                                                               const std::vector<std::array<floatType,3>> relPos,
                                                               const std::array<floatType, 3> momentOfInertia) {
  if (_numRegisteredMolTypes != molId) {
    autopas::utils::ExceptionHandler::exception(
        "ParticlePropertiesLibrary::addMolType(): trying to register a molecule type with id {}. Please register types "
        "consecutively, starting at id 0. Currently there are {} registered types.",
        molId, _numRegisteredSiteTypes);
  }
  ++_numRegisteredMolTypes;
  _siteIds.emplace_back(siteIds);
  _relativeSitePositions.emplace_back(relPos);
  _numSites.emplace_back(siteIds.size());
  floatType molMass{0.};
  for (intType site = 0; site < siteIds.size(); ++site) {
    molMass += _siteMasses[siteIds[site]];
  }
  _molMasses.emplace_back(molMass);
  _momentOfInertias.emplace_back(momentOfInertia);
}

template <typename floatType, typename intType>
void ParticlePropertiesLibrary<floatType, intType>::calculateMixingCoefficients() {
  _computedMixingData.resize(_numRegisteredSiteTypes * _numRegisteredSiteTypes);

  auto cutoffSquare = _cutoff * _cutoff;

  for (size_t firstIndex = 0ul; firstIndex < _numRegisteredSiteTypes; ++firstIndex) {
    for (size_t secondIndex = 0ul; secondIndex < _numRegisteredSiteTypes; ++secondIndex) {
      auto globalIndex = _numRegisteredSiteTypes * firstIndex + secondIndex;

      // epsilon
      floatType epsilon24 = 24 * sqrt(_epsilons[firstIndex] * _epsilons[secondIndex]);
      _computedMixingData[globalIndex].epsilon24 = epsilon24;

      // sigma
      floatType sigma = (_sigmas[firstIndex] + _sigmas[secondIndex]) / 2.0;
      floatType sigmaSquare = sigma * sigma;
      _computedMixingData[globalIndex].sigmaSquare = sigmaSquare;

      // shift6
      floatType shift6 = calcShift6(epsilon24, sigmaSquare, cutoffSquare);
      _computedMixingData[globalIndex].shift6 = shift6;
    }
  }
}

template <typename floatType, typename intType>
void ParticlePropertiesLibrary<floatType, intType>::calculateMomentOfInertiaAndAdjustAxes() {
  // todo fix this
  _momentOfInertias.reserve(_numRegisteredMolTypes);

  for (intType molId = 0; molId < _numRegisteredMolTypes; ++molId) {
    // 1. Calculate Moment of Inertia as 3x3 array
    Eigen::MatrixXd MoI =  Eigen::MatrixXd::Zero(3);
    for (intType site = 0; site < _numSites[molId]; ++site) {
      MoI(0,0) += (_relativeSitePositions[molId][site][1]*_relativeSitePositions[molId][site][1] +
                    _relativeSitePositions[molId][site][2]*_relativeSitePositions[molId][site][2])
                   * _siteMasses[_siteIds[molId][site]];
      MoI(1,1) += (_relativeSitePositions[molId][site][0]*_relativeSitePositions[molId][site][0] +
                    _relativeSitePositions[molId][site][2]*_relativeSitePositions[molId][site][2])
                   * _siteMasses[_siteIds[molId][site]];
      MoI(2,2) += (_relativeSitePositions[molId][site][0]*_relativeSitePositions[molId][site][0] +
                    _relativeSitePositions[molId][site][1]*_relativeSitePositions[molId][site][1])
                   * _siteMasses[_siteIds[molId][site]];
      MoI(0,1) -= _relativeSitePositions[molId][site][0] * _relativeSitePositions[molId][site][1] *
                   _siteMasses[_siteIds[molId][site]];
      MoI(0,2) -= _relativeSitePositions[molId][site][0] * _relativeSitePositions[molId][site][2] *
                   _siteMasses[_siteIds[molId][site]];
      MoI(1,2) -= _relativeSitePositions[molId][site][1] * _relativeSitePositions[molId][site][2] *
                   _siteMasses[_siteIds[molId][site]];
    }
    MoI(1,0) = MoI(0,1);
    MoI(2,0) = MoI(0,2);
    MoI(2,1) = MoI(1,2);
    // 2. Get eigenvalues & -vectors of MoI
    Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(MoI);
    const auto eigenvalues = eigenSolver.eigenvalues();
    const auto eigenvectors = eigenSolver.eigenvectors();

    // 3. Construct diagonal MoI matrix as array of length 3
    const std::array<floatType,3> MoI_diag = {eigenvalues.template cast<floatType>()[0],
        eigenvalues.template cast<floatType>()[1], eigenvalues.template cast<floatType>()[2]};

    _momentOfInertias.push_back(MoI_diag);

    // 4. Rotate sites of mol using rotation matrix formed from eigenvectors to fit new orientation
    for (intType site = 0; site < _numSites[molId]; ++site) {
      const Eigen::Vector3d oldSitePos(_relativeSitePositions[molId][site].data());
      const Eigen::Vector3d newSitePos = eigenvectors * oldSitePos; // todo check
      _relativeSitePositions[molId][site][0] = newSitePos(0);
      _relativeSitePositions[molId][site][1] = newSitePos(1);
      _relativeSitePositions[molId][site][2] = newSitePos(2);
    }
  }
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::getSiteMass(intType i) const {
  return _siteMasses[i];
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::getMolMass(intType i) const {
  return _molMasses[i];
}

template <typename floatType, typename intType>
std::array<floatType,3> ParticlePropertiesLibrary<floatType, intType>::getMomentOfInertia(intType i) const {
  return _momentOfInertias[i];
}

template <typename floatType, typename intType>
std::vector<std::array<floatType,3>> ParticlePropertiesLibrary<floatType, intType>::getSitePositions(intType i) const {
  return _relativeSitePositions[i];
}

template <typename floatType, typename intType>
std::vector<intType> ParticlePropertiesLibrary<floatType, intType>::getSiteTypes(intType i) const {
  return _siteIds[i];
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::get24Epsilon(intType i) const {
  return _computedMixingData[i * _numRegisteredSiteTypes + i].epsilon24;
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::getSigmaSquare(intType i) const {
  return _computedMixingData[i * _numRegisteredSiteTypes + i].sigmaSquare;
}

template<typename floatType, typename intType>
intType ParticlePropertiesLibrary<floatType, intType>::getNumSites(intType i) const {
  return _numSites[i];
}

template <typename floatType, typename intType>
double ParticlePropertiesLibrary<floatType, intType>::calcShift6(double epsilon24, double sigmaSquare,
                                                                 double cutoffSquare) {
  auto sigmaPow2DivCutoff = sigmaSquare / (cutoffSquare);
  auto sigmaPow6DivCutoff = sigmaPow2DivCutoff * sigmaPow2DivCutoff * sigmaPow2DivCutoff;
  auto shift6 = epsilon24 * (sigmaPow6DivCutoff - sigmaPow6DivCutoff * sigmaPow6DivCutoff);
  return shift6;
}
