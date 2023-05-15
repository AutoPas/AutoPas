/**
 * @file RegularGridDecomposition.cpp
 * @author J. Körner
 * @date 19.04.2021
 */
#include "RegularGridDecomposition.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <numeric>

#include "DomainTools.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Quaternion.h"
#include "autopas/utils/Math.h"
#include "src/ParticleCommunicator.h"
#include "src/TypeDefinitions.h"

RegularGridDecomposition::RegularGridDecomposition(const MDFlexConfig &configuration)
    : _loadBalancerOption(configuration.loadBalancer.value),
      _cutoffWidth(configuration.cutoff.value),
      _skinWidthPerTimestep(configuration.verletSkinRadiusPerTimestep.value),
      _rebuildFrequency(configuration.verletRebuildFrequency.value),
      _globalBoxMin(configuration.boxMin.value),
      _globalBoxMax(configuration.boxMax.value),
      _boundaryType(configuration.boundaryOption.value),
      _mpiCommunicationNeeded(
#if defined(AUTOPAS_INCLUDE_MPI)
          true
#else
          false
#endif
      ) {
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &_subdomainCount);
  // if there is only one rank, but we are running with MPI we actually do not need it.
  if (_subdomainCount == 1) {
    _mpiCommunicationNeeded = false;
  }

  _decomposition = DomainTools::generateDecomposition(_subdomainCount, configuration.subdivideDimension.value);

  // For reflection, we interact particles with their 'mirror' only if it would cause a repulsive force. This occurs
  // only for particles within sixthRootOfTwo * sigma / 2.0 of the boundary. For minimal redundant iterating over
  // particles, we iterate once over particles within the largest possible range where a particle might experience a
  // repulsion, i.e. sixthRootOfTwo * maxSigma / 2.0.
  double maxSigma{0};
  for (const auto &[_, sigma] : configuration.sigmaMap.value) {
    maxSigma = std::max(maxSigma, sigma);
  }
  _maxReflectiveSkin = sixthRootOfTwo * maxSigma / 2.;

  // initialize _communicator and _domainIndex
  initializeMPICommunicator();
  // initialize _planarCommunicators and domainID
  initializeLocalDomain();
  // initialize _localBoxMin/Max
  initializeLocalBox();
  // initialize _neighborDomainIndices
  initializeNeighborIndices();

#if defined(MD_FLEXIBLE_ENABLE_ALLLBL)
  if (_loadBalancerOption == LoadBalancerOption::all) {
    _allLoadBalancer = std::make_unique<ALL::ALL<double, double>>(ALL::TENSOR, _dimensionCount, 0);
    _allLoadBalancer->setCommunicator(_communicator);

    const double minDomainSize = 2 * (_cutoffWidth + _skinWidthPerTimestep * _rebuildFrequency);
    _allLoadBalancer->setMinDomainSize({minDomainSize, minDomainSize, minDomainSize});
    _allLoadBalancer->setup();
  }
#else
  if (_domainIndex == 0 and _loadBalancerOption == LoadBalancerOption::all) {
    std::cout << "ALL load balancer has been disabled during compile time. Load balancing will be turned off."
              << std::endl;
  }
#endif
}

RegularGridDecomposition::~RegularGridDecomposition() = default;

void RegularGridDecomposition::update(const double &work) {
  if (_mpiCommunicationNeeded) {
    switch (_loadBalancerOption) {
      case LoadBalancerOption::invertedPressure: {
        balanceWithInvertedPressureLoadBalancer(work);
        break;
      }
#if defined(MD_FLEXIBLE_ENABLE_ALLLBL)
      case LoadBalancerOption::all: {
        balanceWithAllLoadBalancer(work);
        break;
      }
#endif
      default: {
        // do nothing
      }
    }
  }
}

int RegularGridDecomposition::getNumberOfSubdomains() const {
  return std::accumulate(_decomposition.begin(), _decomposition.end(), 1, std::multiplies<>());
}

void RegularGridDecomposition::initializeMPICommunicator() {
  // have to cast away const because MPI requires it
  autopas::AutoPas_MPI_Cart_create(AUTOPAS_MPI_COMM_WORLD, _dimensionCount, _decomposition.data(),
                                   const_cast<int *>(_periods.data()), true, &_communicator);
  autopas::AutoPas_MPI_Comm_rank(_communicator, &_domainIndex);
}

void RegularGridDecomposition::initializeLocalDomain() {
  // have to cast away const because MPI requires it
  autopas::AutoPas_MPI_Cart_get(_communicator, _dimensionCount, _decomposition.data(),
                                const_cast<int *>(_periods.data()), _domainId.data());

  // Create planar communicators used for diffuse load balancing.
  for (int i = 0; i < _dimensionCount; ++i) {
    if (_mpiCommunicationNeeded) {
      const int key = _decomposition[(i + 1) % _dimensionCount] * _domainId[(i + 2) % _dimensionCount] +
                      _domainId[(i + 1) % _dimensionCount];
      autopas::AutoPas_MPI_Comm_split(_communicator, _domainId[i], key, &_planarCommunicators[i]);
    } else {
      _planarCommunicators[i] = AUTOPAS_MPI_COMM_WORLD;
    }
  }
}

void RegularGridDecomposition::initializeLocalBox() {
  for (int i = 0; i < 3; ++i) {
    double localBoxWidth = (_globalBoxMax[i] - _globalBoxMin[i]) / static_cast<double>(_decomposition[i]);

    _localBoxMin[i] = _domainId[i] * localBoxWidth + _globalBoxMin[i];
    _localBoxMax[i] = (_domainId[i] + 1) * localBoxWidth + _globalBoxMin[i];
  }
}

void RegularGridDecomposition::initializeNeighborIndices() {
  for (int i = 0; i < 3; ++i) {
    auto neighborIndex = i * 2;
    auto precedingNeighborId = _domainId;
    precedingNeighborId[i] = (--precedingNeighborId[i] + _decomposition[i]) % _decomposition[i];
    _neighborDomainIndices[neighborIndex] = DomainTools::convertIdToIndex(precedingNeighborId, _decomposition);

    ++neighborIndex;
    auto succeedingNeighborId = _domainId;
    succeedingNeighborId[i] = (++succeedingNeighborId[i] + _decomposition[i]) % _decomposition[i];
    _neighborDomainIndices[neighborIndex] = DomainTools::convertIdToIndex(succeedingNeighborId, _decomposition);
  }
}

bool RegularGridDecomposition::isInsideLocalDomain(const std::array<double, 3> &coordinates) const {
  return DomainTools::isInsideDomain(coordinates, _localBoxMin, _localBoxMax);
}

std::array<int, 6> RegularGridDecomposition::getExtentOfSubdomain(const int subdomainIndex) const {
  return DomainTools::getExtentOfSubdomain(subdomainIndex, _decomposition);
}

void RegularGridDecomposition::exchangeHaloParticles(AutoPasType &autoPasContainer) {
  using autopas::utils::Math::isNear;
  std::vector<ParticleType> haloParticles{};

  for (int dimensionIndex = 0; dimensionIndex < _dimensionCount; ++dimensionIndex) {
    // completely bypass Halo particle exchange in this dimension if boundaries in this direction are not periodic
    // *and* if both local boundaries are the global boundaries in this dimension
    if (_boundaryType[dimensionIndex] != options::BoundaryTypeOption::periodic and
        isNear(_localBoxMin[dimensionIndex], _globalBoxMin[dimensionIndex]) and
        isNear(_localBoxMax[dimensionIndex], _globalBoxMax[dimensionIndex])) {
      continue;
    }

    auto particlesForLeftNeighbor = collectHaloParticlesForLeftNeighbor(autoPasContainer, dimensionIndex);
    auto particlesForRightNeighbor = collectHaloParticlesForRightNeighbor(autoPasContainer, dimensionIndex);

    const double leftHaloMin = _localBoxMin[dimensionIndex] - _skinWidthPerTimestep * _rebuildFrequency;
    const double leftHaloMax = _localBoxMin[dimensionIndex] + _cutoffWidth + _skinWidthPerTimestep * _rebuildFrequency;
    const double rightHaloMin = _localBoxMax[dimensionIndex] - _cutoffWidth - _skinWidthPerTimestep * _rebuildFrequency;
    const double rightHaloMax = _localBoxMax[dimensionIndex] + _skinWidthPerTimestep * _rebuildFrequency;

    for (const auto &particle : haloParticles) {
      std::array<double, _dimensionCount> position = particle.getR();
      if (position[dimensionIndex] >= leftHaloMin and position[dimensionIndex] < leftHaloMax) {
        particlesForLeftNeighbor.push_back(particle);

        // Apply boundary condition
        if (isNear(_localBoxMin[dimensionIndex], _globalBoxMin[dimensionIndex])) {
          position[dimensionIndex] =
              position[dimensionIndex] + (_globalBoxMax[dimensionIndex] - _globalBoxMin[dimensionIndex]);
          particlesForLeftNeighbor.back().setR(position);
        }
      } else if (position[dimensionIndex] >= rightHaloMin and position[dimensionIndex] < rightHaloMax) {
        particlesForRightNeighbor.push_back(particle);

        // Apply boundary condition
        if (isNear(_localBoxMax[dimensionIndex], _globalBoxMax[dimensionIndex])) {
          position[dimensionIndex] =
              position[dimensionIndex] - (_globalBoxMax[dimensionIndex] - _globalBoxMin[dimensionIndex]);
          particlesForRightNeighbor.back().setR(position);
        }
      }
    }
    // See documentation for _neighborDomainIndices to explain the indexing
    const int leftNeighbor = _neighborDomainIndices[(dimensionIndex * 2) % _neighborCount];
    const int rightNeighbor = _neighborDomainIndices[(dimensionIndex * 2 + 1) % _neighborCount];
    const auto receivedHaloParticles = sendAndReceiveParticlesLeftAndRight(
        particlesForLeftNeighbor, particlesForRightNeighbor, leftNeighbor, rightNeighbor);
    haloParticles.insert(haloParticles.end(), receivedHaloParticles.begin(), receivedHaloParticles.end());
  }
  autoPasContainer.addHaloParticles(haloParticles);
}

void RegularGridDecomposition::exchangeMigratingParticles(AutoPasType &autoPasContainer,
                                                          std::vector<ParticleType> &emigrants) {
  using autopas::utils::Math::isNear;
  for (int dimensionIndex = 0; dimensionIndex < _dimensionCount; ++dimensionIndex) {
    // if this rank spans the whole dimension but is not periodic -> skip.
    if (_boundaryType[dimensionIndex] != options::BoundaryTypeOption::periodic and
        isNear(_localBoxMin[dimensionIndex], _globalBoxMin[dimensionIndex]) and
        isNear(_localBoxMax[dimensionIndex], _globalBoxMax[dimensionIndex]))
      continue;

    // If the ALL load balancer is used, it may happen that particles migrate to a non-adjacent domain.
    // Therefore, we need to migrate particles as many times as there are grid cells along the dimension.
    const int maximumSendSteps = _loadBalancerOption == LoadBalancerOption::all ? _decomposition[dimensionIndex] : 1;

    for (int gridIndex = 0; gridIndex < maximumSendSteps; ++gridIndex) {
      // See documentation for _neighborDomainIndices to explain the indexing
      const int leftNeighbor = _neighborDomainIndices[(dimensionIndex * 2) % _neighborCount];
      const int rightNeighbor = _neighborDomainIndices[(dimensionIndex * 2 + 1) % _neighborCount];

      const auto &[particlesForLeftNeighbor, particlesForRightNeighbor, remainingEmigrants] =
          categorizeParticlesIntoLeftAndRightNeighbor(emigrants, dimensionIndex);
      emigrants = remainingEmigrants;

      const auto immigrants = sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbor, particlesForRightNeighbor,
                                                                  leftNeighbor, rightNeighbor);
#ifdef AUTOPAS_OPENMP
#pragma omp declare reduction(vecMergeParticle : std::remove_reference_t<decltype(emigrants)> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
// make sure each buffer gets filled equally while not inducing scheduling overhead
#pragma omp parallel for reduction(vecMergeParticle \
                                   : emigrants),    \
    schedule(static, std::max(1ul, immigrants.size() / omp_get_max_threads()))
#endif
      // we can't use range based for loops here because clang accepts this only starting with version 11
      for (size_t i = 0; i < immigrants.size(); ++i) {
        const auto &particle = immigrants[i];
        if (isInsideLocalDomain(particle.getR())) {
          autoPasContainer.addParticle(particle);
        } else {
          emigrants.push_back(particle);
        }
      }
    }
  }
  // sanity check: if the simulation is a closed system there should be no emigrants left at this point.
  if (std::all_of(_boundaryType.begin(), _boundaryType.end(),
                  [](const auto &boundary) {
                    return boundary == options::BoundaryTypeOption::periodic or
                           boundary == options::BoundaryTypeOption::reflective;
                  }) and
      not emigrants.empty()) {
    using autopas::utils::ArrayUtils::operator<<;
    std::stringstream ss;
    ss << "Rank " << _domainIndex << ": All boundaries are periodic or reflective but " << emigrants.size()
       << " migrants could not be re-inserted:\n"
       << autopas::utils::ArrayUtils::to_string(emigrants, "\n", {"", ""}) << "\n\n"
       << "Local box min: " << _localBoxMin << "\n"
       << "Local box max: " << _localBoxMax;
    throw std::runtime_error(ss.str());
  }
}

void RegularGridDecomposition::reflectParticlesAtBoundaries(AutoPasType &autoPasContainer,
                                                            ParticlePropertiesLibraryType &particlePropertiesLib) {
  using autopas::utils::Math::isNear;
  std::array<double, _dimensionCount> reflSkinMin{}, reflSkinMax{};

  for (int dimensionIndex = 0; dimensionIndex < _dimensionCount; ++dimensionIndex) {
    // skip if boundary is not reflective
    if (_boundaryType[dimensionIndex] != options::BoundaryTypeOption::reflective) continue;

    auto reflect = [&](bool isUpper) {
      const auto boundaryPosition = isUpper ? reflSkinMax[dimensionIndex] : reflSkinMin[dimensionIndex];

      for (auto p = autoPasContainer.getRegionIterator(reflSkinMin, reflSkinMax, autopas::IteratorBehavior::owned);
           p.isValid(); ++p) {
        // Check that particle is within 6th root of 2 * sigma
        const auto position = p->getR();
        const auto distanceToBoundary = std::abs(position[dimensionIndex] - boundaryPosition);

        // For single-site molecules, we discard molecules further than sixthRootOfTwo * sigma, and are left only with
        // molecules who will experience repulsion from the boundary.

        // For multi-site molecules, we discard molecules with center-of-mass further than sixthRootOfTwo * the largest sigma of any site of
        // that molecule. Some molecules may experience attraction, and this is only stopped after calculation of force
        // with mirror particle.
        //
        // Note, there is a scenario where a molecule has center-of-mass further than sixthRootOfTwo * sigma, but has a
        // site closer than this distance, with a large enough epsilon, that repulsion would occur. For computational cost
        // reasons, this scenario is neglected - no repulsion occurs. This *should*, in theory, with an appropriate step-size
        // and molecular model, not cause any problems.

        // To produce a "mirrored" multi-site molecule that could be used with the multi-site molecule functor would be very nasty.
        // As such, we reimplement the kernel of the lennard-jones force here.
        // We also use this for single-site molecules, primarily for consistency but it is also suspected to be cheaper than creating
        // a mirror particle for use with the actual functor.


        // Calculates force acting on site from another site
        const auto LJKernel = [](const std::array<double, 3> sitePosition, const std::array<double, 3> mirrorSitePosition, const double sigmaSquared, const double epsilon24) {
          const auto displacement = autopas::utils::ArrayMath::sub(sitePosition, mirrorSitePosition);
          const auto distanceSquared = autopas::utils::ArrayMath::dot(displacement, displacement);

          const auto inverseDistanceSquared = 1. / distanceSquared;
          const auto lj2 = sigmaSquared * inverseDistanceSquared;
          const auto lj6 = lj2 * lj2 * lj2;
          const auto lj12 = lj6 * lj6;
          const auto lj12m6 = lj12 - lj6;
          const auto scalarMultiple = epsilon24 * (lj12 + lj12m6) * inverseDistanceSquared;

          return autopas::utils::ArrayMath::mulScalar(displacement, scalarMultiple);
        };

        const bool reflectMoleculeFlag =
#if MD_FLEXIBLE_MODE==MULTISITE
            distanceToBoundary < sixthRootOfTwo * particlePropertiesLib.getMoleculesLargestSigma(p->getTypeId()) / 2.;
#else
            distanceToBoundary < sixthRootOfTwo * particlePropertiesLib.getSigma(p->getTypeId()) / 2.;
#endif
        if (reflectMoleculeFlag) {
#if MD_FLEXIBLE_MODE==MULTISITE
          // Keep track of current force and torque to see if molecule is repulsed, and, if not, reset the force.
          const auto currentForce = p->getF();
          const auto currentTorque = p->getTorque();

          // load site positions and types
          const auto unrotatedSitePositions = particlePropertiesLib.getSitePositions(p->getTypeId());
          const auto rotatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(p->getQ(), unrotatedSitePositions);
          const auto exactSitePositions = [rotatedSitePositions, position]() {
            std::vector<std::array<double, 3>> returnedPositions{};
            returnedPositions.reserve(rotatedSitePositions.size());
            for (const auto & rotatedSitePosition : rotatedSitePositions) {
              returnedPositions.push_back(autopas::utils::ArrayMath::add(rotatedSitePosition, position));
            }
            return returnedPositions;
          } ();
          const auto siteTypes = particlePropertiesLib.getSiteTypes(p->getTypeId());

          // get positions of opposing mirror sites
          const auto exactMirrorSitePositions = [exactSitePositions, boundaryPosition, dimensionIndex]() {
            std::vector<std::array<double, 3>> returnedPositions{};
            returnedPositions.reserve(exactSitePositions.size());
            for (const auto & exactSitePosition : exactSitePositions) {
              auto mirrorPosition = exactSitePosition;
              const auto displacementToBoundary = boundaryPosition - exactSitePosition[dimensionIndex];
              mirrorPosition[dimensionIndex] += 2 * displacementToBoundary;
              returnedPositions.push_back(mirrorPosition);
            }
            return returnedPositions;
          } ();

          // Add forces + torques for molecule-to-molecule interaction
          for (int site = 0; site < particlePropertiesLib.getNumSites(p->getTypeId()); site++) {
            for (int mirrorSite = 0; mirrorSite < particlePropertiesLib.getNumSites(p->getTypeId()); mirrorSite++) {
              const auto sigmaSquared = particlePropertiesLib.getMixingSigmaSquared(siteTypes[site], siteTypes[site]);
              const auto epsilon24 = particlePropertiesLib.getMixing24Epsilon(siteTypes[site], siteTypes[site]);
              const auto force = LJKernel(exactSitePositions[site], exactMirrorSitePositions[mirrorSite], sigmaSquared, epsilon24);
              p->addF(force);
              p->addTorque((autopas::utils::ArrayMath::cross(rotatedSitePositions[site], force)));
            }
          }
#else
          const auto siteType = p->getTypeId();
          const auto mirrorPosition = [position, boundaryPosition, dimensionIndex]() {
            const auto displacementToBoundary = boundaryPosition - position[dimensionIndex];
            auto returnedPosition = position;
            returnedPosition[dimensionIndex] += 2 * displacementToBoundary;
            return returnedPosition;
          } ();
          const auto sigmaSquared = particlePropertiesLib.getMixingSigmaSquared(siteType, siteType);
          const auto epsilon24 = particlePropertiesLib.getMixing24Epsilon(siteType, siteType);
          const auto force = LJKernel(position, mirrorPosition, sigmaSquared, epsilon24);
          p->addF(force);
#endif

#if MD_FLEXIBLE_MODE==MULTISITE
          // test if attraction has occurred
          const bool reflectionIsAttractive = isUpper ? p->getF()[dimensionIndex] - currentForce[dimensionIndex] > 0 :
                                                      p->getF()[dimensionIndex] - currentForce[dimensionIndex] < 0;
          // reset force if no attraction has occurred
          if (reflectionIsAttractive) {
            p->setF(currentForce);
            p->setTorque(currentTorque);
          }
#endif
        }
      }
    };

    // apply if we are at a global boundary on lower end of the dimension
    if (isNear(_localBoxMin[dimensionIndex], _globalBoxMin[dimensionIndex])) {
      reflSkinMin = _globalBoxMin;
      reflSkinMax = _globalBoxMax;
      reflSkinMax[dimensionIndex] = _globalBoxMin[dimensionIndex] + _maxReflectiveSkin;

      reflect(false);
    }
    // apply if we are at a global boundary on upper end of the dimension
    if (isNear(_localBoxMax[dimensionIndex], _globalBoxMax[dimensionIndex])) {
      reflSkinMin = _globalBoxMin;
      reflSkinMax = _globalBoxMax;
      reflSkinMin[dimensionIndex] = _globalBoxMax[dimensionIndex] - _maxReflectiveSkin;

      reflect(true);
    }
  }
}

std::vector<ParticleType> RegularGridDecomposition::sendAndReceiveParticlesLeftAndRight(
    const std::vector<ParticleType> &particlesToLeft, const std::vector<ParticleType> &particlesToRight,
    int leftNeighbor, int rightNeighbor) {
  std::vector<ParticleType> receivedParticles{};
  // only actually send / receive if we are not sending / receiving to ourselves
  if (_mpiCommunicationNeeded and leftNeighbor != _domainIndex) {
    ParticleCommunicator particleCommunicator(_communicator);

    particleCommunicator.sendParticles(particlesToLeft, leftNeighbor);
    particleCommunicator.sendParticles(particlesToRight, rightNeighbor);

    particleCommunicator.receiveParticles(receivedParticles, leftNeighbor);
    particleCommunicator.receiveParticles(receivedParticles, rightNeighbor);

    particleCommunicator.waitForSendRequests();
  } else {
    receivedParticles.insert(receivedParticles.end(), particlesToLeft.begin(), particlesToLeft.end());
    receivedParticles.insert(receivedParticles.end(), particlesToRight.begin(), particlesToRight.end());
  }
  return receivedParticles;
}

std::vector<ParticleType> RegularGridDecomposition::collectHaloParticlesForLeftNeighbor(AutoPasType &autoPasContainer,
                                                                                        size_t direction) {
  using autopas::utils::Math::isNear;
  using namespace autopas::utils::ArrayMath::literals;

  std::vector<ParticleType> haloParticles{};
  // Calculate halo box for left neighbor
  const auto skinWidth = _skinWidthPerTimestep * _rebuildFrequency;
  const std::array<double, _dimensionCount> boxMin = _localBoxMin - skinWidth;
  const std::array<double, _dimensionCount> boxMax = [&]() {
    auto boxMax = _localBoxMax + skinWidth;
    boxMax[direction] = _localBoxMin[direction] + _cutoffWidth + skinWidth;
    return boxMax;
  }();

  // Collect the halo particles for the left neighbor
  for (auto particle = autoPasContainer.getRegionIterator(boxMin, boxMax, autopas::IteratorBehavior::owned);
       particle.isValid(); ++particle) {
    haloParticles.push_back(*particle);

    // if the particle is outside the global box move it to the other side (periodic boundary)
    if (isNear(_localBoxMin[direction], _globalBoxMin[direction])) {
      auto position = particle->getR();
      position[direction] = position[direction] + (_globalBoxMax[direction] - _globalBoxMin[direction]);
      haloParticles.back().setR(position);
    }
  }
  return haloParticles;
}

std::vector<ParticleType> RegularGridDecomposition::collectHaloParticlesForRightNeighbor(AutoPasType &autoPasContainer,
                                                                                         size_t direction) {
  using autopas::utils::Math::isNear;
  using namespace autopas::utils::ArrayMath::literals;

  std::vector<ParticleType> haloParticles;
  // Calculate left halo box of right neighbor
  const auto skinWidth = _skinWidthPerTimestep * _rebuildFrequency;
  const std::array<double, _dimensionCount> boxMax = _localBoxMax + skinWidth;
  const std::array<double, _dimensionCount> boxMin = [&]() {
    auto boxMin = _localBoxMin - skinWidth;
    boxMin[direction] = _localBoxMax[direction] - _cutoffWidth - skinWidth;
    return boxMin;
  }();

  // Collect the halo particles for the right neighbor
  for (auto particle = autoPasContainer.getRegionIterator(boxMin, boxMax, autopas::IteratorBehavior::owned);
       particle.isValid(); ++particle) {
    haloParticles.push_back(*particle);

    // if the particle is outside the global box move it to the other side (periodic boundary)
    if (isNear(_localBoxMax[direction], _globalBoxMax[direction])) {
      auto position = particle->getR();
      position[direction] = position[direction] - (_globalBoxMax[direction] - _globalBoxMin[direction]);
      haloParticles.back().setR(position);
    }
  }
  return haloParticles;
}

std::tuple<std::vector<ParticleType>, std::vector<ParticleType>, std::vector<ParticleType>>
RegularGridDecomposition::categorizeParticlesIntoLeftAndRightNeighbor(const std::vector<ParticleType> &particles,
                                                                      size_t direction) {
  using autopas::utils::Math::isNear;
  using namespace autopas::utils::ArrayMath::literals;
  const std::array<double, _dimensionCount> globalBoxLength = _globalBoxMax - _globalBoxMin;

  // The chosen size is the best guess based on the particles vector being distributed into three other vectors.
  const auto sizeEstimate = particles.size() / 3;
  std::vector<ParticleType> leftNeighborParticles;
  std::vector<ParticleType> rightNeighborParticles;
  std::vector<ParticleType> uncategorizedParticles;
  leftNeighborParticles.reserve(sizeEstimate);
  rightNeighborParticles.reserve(sizeEstimate);
  uncategorizedParticles.reserve(sizeEstimate);

  for (const auto &particle : particles) {
    auto position = particle.getR();
    // if the particle is left of the box
    if (position[direction] < _localBoxMin[direction]) {
      leftNeighborParticles.push_back(particle);

      // if the particle is outside the global box move it to the other side (periodic boundary)
      if (isNear(_localBoxMin[direction], _globalBoxMin[direction])) {
        // TODO: check if this failsafe is really reasonable.
        //  It should only trigger if a particle's position was already inside the box?
        const auto periodicPosition = position[direction] + globalBoxLength[direction];
        const auto justInsideOfBox = std::nextafter(_globalBoxMax[direction], _globalBoxMin[direction]);
        position[direction] = std::min(justInsideOfBox, periodicPosition);
        leftNeighborParticles.back().setR(position);
      }
      // if the particle is right of the box
    } else if (position[direction] >= _localBoxMax[direction]) {
      rightNeighborParticles.push_back(particle);

      // if the particle is outside the global box move it to the other side (periodic boundary)
      if (isNear(_localBoxMax[direction], _globalBoxMax[direction])) {
        // TODO: check if this failsafe is really reasonable.
        //  It should only trigger if a particle's position was already inside the box?
        const auto periodicPosition = position[direction] - globalBoxLength[direction];
        position[direction] = std::max(_globalBoxMin[direction], periodicPosition);
        rightNeighborParticles.back().setR(position);
      }
    } else {
      uncategorizedParticles.push_back(particle);
    }
  }
  return {leftNeighborParticles, rightNeighborParticles, uncategorizedParticles};
}

void RegularGridDecomposition::balanceWithInvertedPressureLoadBalancer(double work) {
  using autopas::utils::Math::isNear;
  // This is a dummy variable which is not being used. It is required by the non-blocking MPI_Send calls.
  autopas::AutoPas_MPI_Request dummyRequest{};

  auto oldLocalBoxMin = _localBoxMin;
  auto oldLocalBoxMax = _localBoxMax;

  std::array<double, 3> averageWorkInPlane{};

  for (size_t i = 0; i < _dimensionCount; ++i) {
    const int domainCountInPlane =
        _decomposition[(i + 1) % _dimensionCount] * _decomposition[(i + 2) % _dimensionCount];

    // Calculate the average work in the process grid plane
    averageWorkInPlane[i] = work;
    if (domainCountInPlane > 1) {
      autopas::AutoPas_MPI_Allreduce(AUTOPAS_MPI_IN_PLACE, &averageWorkInPlane[i], 1, AUTOPAS_MPI_DOUBLE,
                                     AUTOPAS_MPI_SUM, _planarCommunicators[i]);
      averageWorkInPlane[i] = averageWorkInPlane[i] / domainCountInPlane;
    }

    // Get neighbour indices
    const int leftNeighbor = _neighborDomainIndices[i * 2];
    const int rightNeighbor = _neighborDomainIndices[i * 2 + 1];

    // Send average work in plane to neighbours
    if (not isNear(_localBoxMin[i], _globalBoxMin[i])) {
      autopas::AutoPas_MPI_Isend(&averageWorkInPlane[i], 1, AUTOPAS_MPI_DOUBLE, leftNeighbor, 0, _communicator,
                                 &dummyRequest);
      autopas::AutoPas_MPI_Isend(&oldLocalBoxMax[i], 1, AUTOPAS_MPI_DOUBLE, leftNeighbor, 0, _communicator,
                                 &dummyRequest);
    }

    if (not isNear(_localBoxMax[i], _globalBoxMax[i])) {
      autopas::AutoPas_MPI_Isend(&averageWorkInPlane[i], 1, AUTOPAS_MPI_DOUBLE, rightNeighbor, 0, _communicator,
                                 &dummyRequest);
      autopas::AutoPas_MPI_Isend(&oldLocalBoxMin[i], 1, AUTOPAS_MPI_DOUBLE, rightNeighbor, 0, _communicator,
                                 &dummyRequest);
    }
  }

  for (size_t i = 0; i < _dimensionCount; ++i) {
    // Get neighbour indices
    const int leftNeighbor = _neighborDomainIndices[i * 2];
    const int rightNeighbor = _neighborDomainIndices[i * 2 + 1];

    double neighborPlaneWork{}, neighborBoundary{}, balancedPosition{};
    if (not isNear(_localBoxMin[i], _globalBoxMin[i])) {
      // Receive average work from neighbour planes.
      autopas::AutoPas_MPI_Recv(&neighborPlaneWork, 1, AUTOPAS_MPI_DOUBLE, leftNeighbor, 0, _communicator,
                                AUTOPAS_MPI_STATUS_IGNORE);
      autopas::AutoPas_MPI_Recv(&neighborBoundary, 1, AUTOPAS_MPI_DOUBLE, leftNeighbor, 0, _communicator,
                                AUTOPAS_MPI_STATUS_IGNORE);

      // Calculate balanced positions and only shift by half the resulting distance to prevent localBox min to be larger
      // than localBoxMax.
      balancedPosition = DomainTools::balanceAdjacentDomains(
          neighborPlaneWork, averageWorkInPlane[i], neighborBoundary, oldLocalBoxMax[i],
          2 * (_cutoffWidth + _skinWidthPerTimestep * _rebuildFrequency));
      _localBoxMin[i] += (balancedPosition - _localBoxMin[i]) / 2;
    }

    if (not isNear(_localBoxMax[i], _globalBoxMax[i])) {
      // Receive average work from neighbour planes.
      autopas::AutoPas_MPI_Recv(&neighborPlaneWork, 1, AUTOPAS_MPI_DOUBLE, rightNeighbor, 0, _communicator,
                                AUTOPAS_MPI_STATUS_IGNORE);
      autopas::AutoPas_MPI_Recv(&neighborBoundary, 1, AUTOPAS_MPI_DOUBLE, rightNeighbor, 0, _communicator,
                                AUTOPAS_MPI_STATUS_IGNORE);

      // Calculate balanced positions and only shift by half the resulting distance to prevent localBox min to be larger
      // than localBoxMax.
      balancedPosition = DomainTools::balanceAdjacentDomains(
          averageWorkInPlane[i], neighborPlaneWork, oldLocalBoxMin[i], neighborBoundary,
          2 * (_cutoffWidth + _skinWidthPerTimestep * _rebuildFrequency));
      _localBoxMax[i] += (balancedPosition - _localBoxMax[i]) / 2;
    }
  }
}

#if defined(MD_FLEXIBLE_ENABLE_ALLLBL)
void RegularGridDecomposition::balanceWithAllLoadBalancer(const double &work) {
  std::vector<ALL::Point<double>> domain(2, ALL::Point<double>(3));

  for (int i = 0; i < 3; ++i) {
    domain[0][i] = _localBoxMin[i];
    domain[1][i] = _localBoxMax[i];
  }

  _allLoadBalancer->setVertices(domain);
  _allLoadBalancer->setWork(work);
  _allLoadBalancer->balance();

  std::vector<ALL::Point<double>> updatedVertices = _allLoadBalancer->getVertices();

  for (int i = 0; i < 3; ++i) {
    _localBoxMin[i] = updatedVertices[0][i];
    _localBoxMax[i] = updatedVertices[1][i];
  }
}
#endif
