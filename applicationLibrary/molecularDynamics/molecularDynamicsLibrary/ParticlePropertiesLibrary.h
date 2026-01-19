/**
 * @file ParticlePropertiesLibrary.h
 * @author N. Fottner
 * @date 7/4/19
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <vector>

#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ExceptionHandler.h"

/**
 * This class stores the (physical) properties of molecule types, and, in the case of multi-site molecules, the location
 * of the sites relative to the the center-of-mass.
 *
 * It also provides mixed values for (force) calculations between known types.
 *
 * ToDo: Add a function that computes the diagonalized Moment of Inertia and adjusts the rotation of the positions of
 * the sites such that they correctly match the new axes.
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
   * Registers a new single site type to the library with a given mass.
   * @note New sites must be registered with consecutive siteIds.
   * @note This only registers the site. Potential specific parameters must be added afterwards by calling e.g.
   * `addLJParametersToSite()` for a Lennard-Jones Site.
   *
   * @param siteId
   * @param mass
   */
  void addSiteType(const intType siteId, const floatType mass);

  /**
   * Adds the LJ properties of a single site type to the library.
   *
   * Checks if a site with given siteId was already registered.
   * Old values will be overwritten.
   * @param siteId
   * @param epsilon
   * @param sigma
   */
  void addLJParametersToSite(const intType siteId, const floatType epsilon, const floatType sigma);

  /**
   * Adds the AT properties of a single site type to the library.
   *
   * Checks if a site with given siteId was already registered.
   * Old values will be overwritten.
   * @param siteId
   * @param nu
   */
  void addATParametersToSite(const intType siteId, const floatType nu);

  /**
   * Adds the properties of a molecule type to the library including: position and type of all sites, as well as the
   * diagonalized moment of inertia.
   *
   * If md-flexible has been compiled for single-site molecules, calls to this function result in an error.
   * If the type id already exists the values will be overwritten.
   *
   * @param molId
   * @param siteIds vector of IDs of sites
   * @param relPos vector of relative positions
   * @param momentOfInertia diagonalized moment of inertia as a array of 3 floats
   */
  void addMolType(const intType molId, const std::vector<intType> siteIds,
                  const std::vector<std::array<floatType, 3>> relPos, const std::array<floatType, 3> momentOfInertia);

  /**
   * Calculates the actual mixing coefficients.
   */
  void calculateMixingCoefficients();

  ~ParticlePropertiesLibrary() = default;

  /**
   * Returns the number of registered site / single-site molecule types.
   * @return Number of registered site types.
   */
  [[nodiscard]] int getNumberRegisteredSiteTypes() const { return _numRegisteredSiteTypes; }

  /**
   * Returns the number of registered multi-site molecule types.
   *
   * Throws an error if support for multi-site molecules has not been compiled.
   *
   * @return Number of registered multi-site molecule types.
   */
  [[nodiscard]] int getNumberRegisteredMolTypes() const {
#if MD_FLEXIBLE_MODE == SINGLESITE
    AutoPasLog(
        WARN,
        "ParticlePropertiesLibrary::getNumberRegisteredMolTypes(): trying to get the number of registered multi-site"
        "molecule types when md-flexible has been compiled without support for multi-site molecules. Please compile "
        "with the CMake argument '-DMD_FLEXIBLE_MODE=MULTISITE'.");
#endif
    return _numRegisteredMolTypes;
  }

  /**
   * Getter for the site's epsilon.
   * @param i Type Id of the site or single-site molecule.
   * @return epsilon_i
   */
  floatType getEpsilon(intType i) const;

  /**
   * Getter for the site's sigma.
   * @param i Type Id of the site or single-site molecule.
   * @return sigma_i
   */
  floatType getSigma(intType i) const;

  /**
   * Getter for the site's nu.
   * @param i Type Id of the site or single-site molecule.
   * @return nu_i
   */
  floatType getNu(intType i) const;

  /**
   * Getter for the site's mass.
   * @param i Type Id of the site or single-site molecule.
   * @return mass_i
   */
  floatType getSiteMass(intType i) const;

  /**
   * Getter for a molecules' mass.
   *
   * For single site molecules, this automatically gets the mass of the site with the given ID.
   * For multi site molecules, this gets the total mass of the molecule.
   *
   * @param i Type Id of a multi-site molecule.
   * @return mass_i
   */
  floatType getMolMass(intType i) const;

  /**
   * Getter for the multi-site molecule's diagonalized Moment of Inertia.
   *
   * Throws an error if support for multi-site molecules has not been compiled.
   *
   * @param i Type Id of a multi-site molecule.
   * @return diagonalized moment of inertia, with each element of the returned array corresponding to a diagonal element
   * of the MoI matrix/tensor.
   */
  std::array<floatType, 3> getMomentOfInertia(intType i) const;

  /**
   * Get relative site positions to a multi-site molecule's center-of-mass. These site positions must be appropriately
   * translated and rotated to be used.
   *
   * Throws an error if support for multi-site molecules has not been compiled.
   *
   * @param i Type Id of a multi-site molecule.
   * @return untranslated, non-rotated site positions
   */
  std::vector<std::array<floatType, 3>> getSitePositions(intType i) const;

  /**
   * Get site types of a multi-site molecule.
   *
   * Throws an error if support for multi-site molecules has not been compiled.
   *
   * @param i Type Id of a multi-site molecule.
   * @return site type Ids.
   */
  std::vector<intType> getSiteTypes(intType i) const;

  /**
   * Get number of sites of a multi-site molecule.
   *
   * Throws an error if support for multi-site molecules has not been compiled.
   *
   * @param i Type Id of a multi-site molecule.
   * @return number of sites in molecule's of type Id
   */
  intType getNumSites(intType i) const;

  /**
   * Get the largest sigma of any site of a multi-site molecule.
   *
   * Throws an error if support for multi-site molecules has not been compiled.
   *
   * @param i Type Id of a multi-site molecule.
   * @return largest sigma of that molecule's sites
   */
  floatType getMoleculesLargestSigma(intType i) const;

  /**
   * Returns the precomputed mixed epsilon * 24.
   * @param  i Id of site one.
   * @param  j Id of site two.
   * @return 24*epsilon_ij
   */
  floatType getMixing24Epsilon(intType i, intType j) const {
    return _computedLJMixingData[i * _numRegisteredSiteTypes + j].epsilon24;
  }

  /**
   * Get complete mixing data for one pair of LJ site types.
   * @param i Id of site one.
   * @param j Id of site two.
   * @return
   */
  auto getLJMixingData(intType i, intType j) const {
    return _computedLJMixingData[i * _numRegisteredSiteTypes + j];
  }

  /**
   * Get a pointer to Mixing Data for one pair of LJ site types.
   * @param i Id of site one.
   * @param j Id of site two.
   * @return
   */
  const double *getLJMixingDataPtr(intType i, intType j) {
    return reinterpret_cast<const double *>(&_computedLJMixingData[i * _numRegisteredSiteTypes + j]);
  }

  /**
   * Returns precomputed mixed squared sigma for one pair of site types.
   * @param i Id of site one.
   * @param j Id of site two.
   * @return sigma_ij²
   */
  floatType getMixingSigmaSquared(intType i, intType j) const {
    return _computedLJMixingData[i * _numRegisteredSiteTypes + j].sigmaSquared;
  }

  /**
   * Returns precomputed mixed shift * 6 for one pair of site types.
   * @param i siteId of site one.
   * @param j siteId of site two.
   * @return shift * 6
   */
  floatType getMixingShift6(intType i, intType j) const {
    return _computedLJMixingData[i * _numRegisteredSiteTypes + j].shift6;
  }

  /**
   * Calculate the shift multiplied 6 of the lennard jones potential from given cutoff, epsilon, sigma.
   * The shift * 6 is then added to the total potential energy for every pairwise interaction within the cutoff.
   * @param epsilon24 epsilon * 24
   * @param sigmaSquared sigma squared
   * @param cutoffSquared squared cutoff of the lennard-jones potential
   * @return shift multiplied by 6
   */
  static double calcShift6(double epsilon24, double sigmaSquared, double cutoffSquared);

  /**
   * Returns the precomputed mixed epsilon * 24.
   * @param  i Id of site one.
   * @param  j Id of site two.
   * @param  k Id of site three.
   * @return nu_ijk
   */
  floatType getMixingNu(intType i, intType j, intType k) const {
    return
    _computedATMixingData[i * _numRegisteredSiteTypes * _numRegisteredSiteTypes + j * _numRegisteredSiteTypes + k].nu;
  }

  /**
   * Get complete mixing data for one triplet of AT site types.
   * @param i Id of site one.
   * @param j Id of site two.
   * @param k Id of site three.
   * @return
   */
  auto getATMixingData(intType i, intType j, intType k) const {
    return
    _computedATMixingData[i * _numRegisteredSiteTypes * _numRegisteredSiteTypes + j * _numRegisteredSiteTypes + k];
  }

  struct PackedLJMixingData {
    floatType epsilon24;
    floatType sigmaSquared;
    floatType shift6;
  };

 [[nodiscard]] std::vector<PackedLJMixingData, autopas::AlignedAllocator<PackedLJMixingData>>
  getComputedLJMixingData() const {
    return _computedLJMixingData;
 }

  [[nodiscard]] intType getNumRegisteredSiteTypes() const { return _numRegisteredSiteTypes; }

 private:
  intType _numRegisteredSiteTypes{0};
  intType _numRegisteredMolTypes{0};
  const double _cutoff;

  std::vector<floatType> _epsilons;
  std::vector<floatType> _sigmas;
  std::vector<floatType> _siteMasses;
  std::vector<floatType> _nus;  // Factor for AxilrodTeller potential

  // Note: this is a vector of site type Ids for the sites of a certain molecular Id
  std::vector<std::vector<intType>> _siteIds;
  // This is a vector (indexed by mol ID) of vectors of site positions (which are 3D arrays)
  std::vector<std::vector<std::array<floatType, 3>>> _relativeSitePositions;
  std::vector<floatType> _molMasses;
  std::vector<std::array<floatType, 3>> _momentOfInertias;
  std::vector<size_t> _numSites;
  std::vector<floatType> _moleculesLargestSigma;

  // Allocate memory for the respective parameters
  bool _storeLJData{false};
  bool _storeATData{false};

  struct PackedLJMixingData {
    floatType epsilon24;
    floatType sigmaSquared;
    floatType shift6;
  };

  struct PackedATMixingData {
    floatType nu;
  };

  std::vector<PackedLJMixingData, autopas::AlignedAllocator<PackedLJMixingData>> _computedLJMixingData;
  std::vector<PackedATMixingData, autopas::AlignedAllocator<PackedATMixingData>> _computedATMixingData;
};

template <typename floatType, typename intType>
void ParticlePropertiesLibrary<floatType, intType>::addSiteType(intType siteID, floatType mass) {
  if (_numRegisteredSiteTypes != siteID) {
    autopas::utils::ExceptionHandler::exception(
        "ParticlePropertiesLibrary::addSiteType(): trying to register a site type with id {}. Please "
        "register types "
        "consecutively, starting at id 0. Currently there are {} registered types.",
        siteID, _numRegisteredSiteTypes);
  }
  ++_numRegisteredSiteTypes;
  _siteMasses.emplace_back(mass);

  // Allocate memory for all parameters of used models
  if (_storeLJData) {
    _sigmas.emplace_back(0.0);
    _epsilons.emplace_back(0.0);
  }
  if (_storeATData) {
    _nus.emplace_back(0.0);
  }
}

template <typename floatType, typename intType>
void ParticlePropertiesLibrary<floatType, intType>::addLJParametersToSite(intType siteID, floatType epsilon,
                                                                          floatType sigma) {
  if (siteID >= _numRegisteredSiteTypes) {
    autopas::utils::ExceptionHandler::exception(
        "ParticlePropertiesLibrary::addLJParametersToSite(): Trying to set lennard-jones parameters for a site type "
        "with id {},"
        " which has not been registered yet. Currently there are {} registered types.",
        siteID, _numRegisteredSiteTypes);
  }
  _storeLJData = true;
  if (_epsilons.size() != _numRegisteredSiteTypes) {
    _epsilons.resize(_numRegisteredSiteTypes);
    _sigmas.resize(_numRegisteredSiteTypes);
  }
  _epsilons[siteID] = epsilon;
  _sigmas[siteID] = sigma;
}

template <typename floatType, typename intType>
void ParticlePropertiesLibrary<floatType, intType>::addATParametersToSite(intType siteID, floatType nu) {
  if (siteID >= _numRegisteredSiteTypes) {
    autopas::utils::ExceptionHandler::exception(
        "ParticlePropertiesLibrary::addATParametersToSite(): Trying to set the axilrod-teller parameter for a site "
        "type with id {},"
        " which has not been registered yet. Currently there are {} registered types.",
        siteID, _numRegisteredSiteTypes);
  }
  _storeATData = true;
  if (_nus.size() != _numRegisteredSiteTypes) {
    _nus.resize(_numRegisteredSiteTypes);
  }
  _nus[siteID] = nu;
}

template <typename floatType, typename intType>
void ParticlePropertiesLibrary<floatType, intType>::addMolType(const intType molId, const std::vector<intType> siteIds,
                                                               const std::vector<std::array<floatType, 3>> relPos,
                                                               const std::array<floatType, 3> momentOfInertia) {
  // Error handling
#if MD_FLEXIBLE_MODE == SINGLESITE
  AutoPasLog(WARN,
             "ParticlePropertiesLibrary::addMolType(): trying to register a multi-site molecule type when md-flexible "
             "has been compiled without support for multi-site molecules. Please compile with the CMake argument "
             "'-DMD_FLEXIBLE_MODE=MULTISITE'.");
#endif
  if (_numRegisteredMolTypes != molId) {
    autopas::utils::ExceptionHandler::exception(
        "ParticlePropertiesLibrary::addMolType(): trying to register a molecule type with id {}. Please register types "
        "consecutively, starting at id 0. Currently there are {} registered types.",
        molId, _numRegisteredSiteTypes);
  }
  if (std::any_of(siteIds.cbegin(), siteIds.cend(), [this](intType i) { return i >= this->_numRegisteredSiteTypes; })) {
    autopas::utils::ExceptionHandler::exception(
        "ParticlePropertiesLibrary::addMolType(): trying to register a molecule type with an unregistered site type "
        "Id.");
  }
  if (siteIds.size() != relPos.size()) {
    autopas::utils::ExceptionHandler::exception(
        "ParticlePropertiesLibrary::addMolType(): trying to register a molecule type with vectors of site IDs and site"
        "positions that do not match in size.");
  }

  // Add molecule type if there are no errors.
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

  floatType molLargestSigma{0.};
  for (size_t site = 0; site < siteIds.size(); site++) {
    molLargestSigma = std::max(molLargestSigma, _sigmas[siteIds[site]]);
  }
  _moleculesLargestSigma.emplace_back(molLargestSigma);
}

template <typename floatType, typename intType>
void ParticlePropertiesLibrary<floatType, intType>::calculateMixingCoefficients() {
  if (_numRegisteredSiteTypes == 0) {
    autopas::utils::ExceptionHandler::AutoPasException(
        "ParticlePropertiesLibrary::calculateMixingCoefficients was called without any site types being registered!");
  }

  // There are Lennard-Jones Sites
  if (_storeLJData) {
    const auto cutoffSquared = _cutoff * _cutoff;
    _computedLJMixingData.resize(_numRegisteredSiteTypes * _numRegisteredSiteTypes);

    for (size_t firstIndex = 0ul; firstIndex < _numRegisteredSiteTypes; ++firstIndex) {
      for (size_t secondIndex = 0ul; secondIndex < _numRegisteredSiteTypes; ++secondIndex) {
        auto globalIndex = _numRegisteredSiteTypes * firstIndex + secondIndex;

        // epsilon
        const floatType epsilon24 = 24 * sqrt(_epsilons[firstIndex] * _epsilons[secondIndex]);
        _computedLJMixingData[globalIndex].epsilon24 = epsilon24;

        // sigma
        const floatType sigma = (_sigmas[firstIndex] + _sigmas[secondIndex]) / 2.0;
        const floatType sigmaSquared = sigma * sigma;
        _computedLJMixingData[globalIndex].sigmaSquared = sigmaSquared;

        // shift6
        const floatType shift6 = calcShift6(epsilon24, sigmaSquared, cutoffSquared);
        _computedLJMixingData[globalIndex].shift6 = shift6;
      }
    }
  }

  if (_storeATData) {
    _computedATMixingData.resize(_numRegisteredSiteTypes * _numRegisteredSiteTypes * _numRegisteredSiteTypes);
    for (size_t firstIndex = 0ul; firstIndex < _numRegisteredSiteTypes; ++firstIndex) {
      for (size_t secondIndex = 0ul; secondIndex < _numRegisteredSiteTypes; ++secondIndex) {
        for (size_t thirdIndex = 0ul; thirdIndex < _numRegisteredSiteTypes; ++thirdIndex) {
          const auto globalIndex3B = _numRegisteredSiteTypes * _numRegisteredSiteTypes * firstIndex +
                                     _numRegisteredSiteTypes * secondIndex + thirdIndex;
          // geometric mixing as used in e.g. https://doi.org/10.1063/1.3567308
          const floatType mixedNu = cbrt(_nus[firstIndex] * _nus[secondIndex] * _nus[thirdIndex]);
          _computedATMixingData[globalIndex3B].nu = mixedNu;
        }
      }
    }
  }
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::getSiteMass(intType i) const {
  return _siteMasses[i];
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::getMolMass(intType i) const {
#if MD_FLEXIBLE_MODE == MULTISITE
  return _molMasses[i];
#else
  return _siteMasses[i];
#endif
}

template <typename floatType, typename intType>
std::array<floatType, 3> ParticlePropertiesLibrary<floatType, intType>::getMomentOfInertia(intType i) const {
#if MD_FLEXIBLE_MODE == SINGLESITE
  AutoPasLog(WARN,
             "ParticlePropertiesLibrary::getMomentOfInertia(): trying to get the Moment of Inertia of a multi-site "
             "molecule type when md-flexible has been compiled without support for multi-site molecules. Please "
             "compile with the CMake argument '-DMD_FLEXIBLE_MODE=MULTISITE'.");
#endif
  return _momentOfInertias[i];
}

template <typename floatType, typename intType>
std::vector<std::array<floatType, 3>> ParticlePropertiesLibrary<floatType, intType>::getSitePositions(intType i) const {
#if MD_FLEXIBLE_MODE == SINGLESITE
  AutoPasLog(WARN,
             "ParticlePropertiesLibrary::getSitePositions(): trying to get the site positions of a multi-site molecule "
             "type when md-flexible has been compiled without support for multi-site molecules. Please compile with "
             "the CMake argument '-DMD_FLEXIBLE_MODE=MULTISITE'.");
#endif
  return _relativeSitePositions[i];
}

template <typename floatType, typename intType>
std::vector<intType> ParticlePropertiesLibrary<floatType, intType>::getSiteTypes(intType i) const {
#if MD_FLEXIBLE_MODE == SINGLESITE
  AutoPasLog(WARN,
             "ParticlePropertiesLibrary::getSiteTypes(): trying to get the site types of a multi-site molecule type "
             "when md-flexible has been compiled without support for multi-site molecules. Please compile with the "
             "CMake argument '-DMD_FLEXIBLE_MODE=MULTISITE'.");
#endif
  return _siteIds[i];
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::getEpsilon(intType i) const {
  return _epsilons[i];
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::getSigma(intType i) const {
  return _sigmas[i];
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::getNu(intType i) const {
  return _nus[i];
}

template <typename floatType, typename intType>
intType ParticlePropertiesLibrary<floatType, intType>::getNumSites(intType i) const {
#if MD_FLEXIBLE_MODE == SINGLESITE
  AutoPasLog(WARN,
             "ParticlePropertiesLibrary::getNumSites(): trying to get the number of sites of a multi-site molecule "
             "type when md-flexible has been compiled without support for multi-site molecules. Please compile with "
             "the CMake argument '-DMD_FLEXIBLE_MODE=MULTISITE'.");
#endif
  return _numSites[i];
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::getMoleculesLargestSigma(intType i) const {
#if MD_FLEXIBLE_MODE == SINGLESITE
  AutoPasLog(WARN,
             "ParticlePropertiesLibrary::getNumSites(): trying to get the number of sites of a multi-site molecule "
             "type when md-flexible has been compiled without support for multi-site molecules. Please compile with "
             "the CMake argument '-DMD_FLEXIBLE_MODE=MULTISITE'.");
#endif
  return _moleculesLargestSigma[i];
}

template <typename floatType, typename intType>
double ParticlePropertiesLibrary<floatType, intType>::calcShift6(double epsilon24, double sigmaSquared,
                                                                 double cutoffSquared) {
  const auto sigmaDivCutoffPow2 = sigmaSquared / cutoffSquared;
  const auto sigmaDivCutoffPow6 = sigmaDivCutoffPow2 * sigmaDivCutoffPow2 * sigmaDivCutoffPow2;
  const auto shift6 = epsilon24 * (sigmaDivCutoffPow6 - sigmaDivCutoffPow6 * sigmaDivCutoffPow6);
  return shift6;
}
