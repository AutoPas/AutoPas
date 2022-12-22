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
#include "src/ParticleCommunicator.h"

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

  // initialize _communicator and _domainIndex
  initializeMPICommunicator();
  // initialize _planarCommunicators and domainID
  initializeLocalDomain();
  // initialize _localBoxMin/Max
  initializeLocalBox();
  // initialize _neighborDomainIndices
  initializeNeighborIndices();

#if defined(AUTOPAS_ENABLE_ALLLBL)
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
#if defined(AUTOPAS_ENABLE_ALLLBL)
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
  std::vector<ParticleType> haloParticles{};

  for (int dimensionIndex = 0; dimensionIndex < _dimensionCount; ++dimensionIndex) {
    // completely bypass Halo particle exchange in this dimension if boundaries in this direction are not periodic
    // *and* if both local boundaries are the global boundaries in this dimension
    if (_boundaryType[dimensionIndex] != options::BoundaryTypeOption::periodic and
        _localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex] and
        _localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
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
        if (_localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex]) {
          position[dimensionIndex] =
              position[dimensionIndex] + (_globalBoxMax[dimensionIndex] - _globalBoxMin[dimensionIndex]);
          particlesForLeftNeighbor.back().setR(position);
        }
      } else if (position[dimensionIndex] >= rightHaloMin and position[dimensionIndex] < rightHaloMax) {
        particlesForRightNeighbor.push_back(particle);

        // Apply boundary condition
        if (_localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
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
  for (const auto &particle : haloParticles) {
    autoPasContainer.addHaloParticle(particle);
  }
}

void RegularGridDecomposition::exchangeMigratingParticles(AutoPasType &autoPasContainer,
                                                          std::vector<ParticleType> &emigrants) {
  for (int dimensionIndex = 0; dimensionIndex < _dimensionCount; ++dimensionIndex) {
    // if this rank spans the whole dimension but is not periodic -> skip.
    if (_boundaryType[dimensionIndex] != options::BoundaryTypeOption::periodic and
        _localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex] and
        _localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex])
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

      for (const auto &particle : immigrants) {
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

void RegularGridDecomposition::reflectParticlesAtBoundaries(AutoPasType &autoPasContainer) {
  std::array<double, _dimensionCount> reflSkinMin{}, reflSkinMax{};

  for (int dimensionIndex = 0; dimensionIndex < _dimensionCount; ++dimensionIndex) {
    // skip if boundary is not reflective
    if (_boundaryType[dimensionIndex] != options::BoundaryTypeOption::reflective) continue;

    auto reflect = [&](bool isUpper) {
      for (auto p = autoPasContainer.getRegionIterator(reflSkinMin, reflSkinMax, autopas::IteratorBehavior::owned);
           p.isValid(); ++p) {
        auto vel = p->getV();
        // reverse velocity in dimension if towards boundary
        if ((vel[dimensionIndex] < 0) xor isUpper) {
          vel[dimensionIndex] *= -1;
        }
        p->setV(vel);
      }
    };

    // apply if we are at a global boundary on lower end of the dimension
    if (_localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex]) {
      reflSkinMin = _globalBoxMin;
      reflSkinMax = _globalBoxMax;
      reflSkinMax[dimensionIndex] = _globalBoxMin[dimensionIndex] + autoPasContainer.getVerletSkinPerTimestep() / 2;

      reflect(false);
    }
    // apply if we are at a global boundary on upper end of the dimension
    if (_localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
      reflSkinMin = _globalBoxMin;
      reflSkinMax = _globalBoxMax;
      reflSkinMin[dimensionIndex] = _globalBoxMax[dimensionIndex] - autoPasContainer.getVerletSkinPerTimestep() / 2;

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
  std::vector<ParticleType> haloParticles{};
  // Calculate halo box for left neighbor
  const auto skinWidth = _skinWidthPerTimestep * _rebuildFrequency;
  const std::array<double, _dimensionCount> boxMin = autopas::utils::ArrayMath::subScalar(_localBoxMin, skinWidth);
  const std::array<double, _dimensionCount> boxMax = [&]() {
    auto boxMax = autopas::utils::ArrayMath::addScalar(_localBoxMax, skinWidth);
    boxMax[direction] = _localBoxMin[direction] + _cutoffWidth + skinWidth;
    return boxMax;
  }();

  // Collect the halo particles for the left neighbor
  for (auto particle = autoPasContainer.getRegionIterator(boxMin, boxMax, autopas::IteratorBehavior::owned);
       particle.isValid(); ++particle) {
    haloParticles.push_back(*particle);

    // if the particle is outside the global box move it to the other side (periodic boundary)
    if (_localBoxMin[direction] == _globalBoxMin[direction]) {
      auto position = particle->getR();
      position[direction] = position[direction] + (_globalBoxMax[direction] - _globalBoxMin[direction]);
      haloParticles.back().setR(position);
    }
  }
  return haloParticles;
}

std::vector<ParticleType> RegularGridDecomposition::collectHaloParticlesForRightNeighbor(AutoPasType &autoPasContainer,
                                                                                         size_t direction) {
  std::vector<ParticleType> haloParticles;
  // Calculate left halo box of right neighbor
  const auto skinWidth = _skinWidthPerTimestep * _rebuildFrequency;
  const std::array<double, _dimensionCount> boxMax = autopas::utils::ArrayMath::addScalar(_localBoxMax, skinWidth);
  const std::array<double, _dimensionCount> boxMin = [&]() {
    auto boxMin = autopas::utils::ArrayMath::subScalar(_localBoxMin, skinWidth);
    boxMin[direction] = _localBoxMax[direction] - _cutoffWidth - skinWidth;
    return boxMin;
  }();

  // Collect the halo particles for the right neighbor
  for (auto particle = autoPasContainer.getRegionIterator(boxMin, boxMax, autopas::IteratorBehavior::owned);
       particle.isValid(); ++particle) {
    haloParticles.push_back(*particle);

    // if the particle is outside the global box move it to the other side (periodic boundary)
    if (_localBoxMax[direction] == _globalBoxMax[direction]) {
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
  const std::array<double, _dimensionCount> globalBoxLength =
      autopas::utils::ArrayMath::sub(_globalBoxMax, _globalBoxMin);

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
      if (_localBoxMin[direction] == _globalBoxMin[direction]) {
        // TODO: check if this failsafe is really reasonable.
        //  It should only trigger if a particle's position was already inside the box?
        const auto periodicPosition = position[direction] + globalBoxLength[direction];
        const auto justInsideOfBox = std::nextafter(_globalBoxMax[direction], _globalBoxMin[direction]);
        position[direction] = std::min(justInsideOfBox, periodicPosition);
        leftNeighborParticles.back().setR(position);
      }
    } else  // if the particle is right of the box
      if (position[direction] >= _localBoxMax[direction]) {
        rightNeighborParticles.push_back(particle);

        // if the particle is outside the global box move it to the other side (periodic boundary)
        if (_localBoxMax[direction] == _globalBoxMax[direction]) {
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
    if (_localBoxMin[i] != _globalBoxMin[i]) {
      autopas::AutoPas_MPI_Isend(&averageWorkInPlane[i], 1, AUTOPAS_MPI_DOUBLE, leftNeighbor, 0, _communicator,
                                 &dummyRequest);
      autopas::AutoPas_MPI_Isend(&oldLocalBoxMax[i], 1, AUTOPAS_MPI_DOUBLE, leftNeighbor, 0, _communicator,
                                 &dummyRequest);
    }

    if (_localBoxMax[i] != _globalBoxMax[i]) {
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
    if (_localBoxMin[i] != _globalBoxMin[i]) {
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

    if (_localBoxMax[i] != _globalBoxMax[i]) {
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

#if defined(AUTOPAS_ENABLE_ALLLBL)
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
