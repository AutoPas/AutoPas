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
      _skinWidth(configuration.verletSkinRadius.value),
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
  // if there is only one rank but we are running with MPI we actually do not need it.
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

    const double minDomainSize = 2 * (_cutoffWidth + _skinWidth);
    _allLoadBalancer->setMinDomainSize({minDomainSize, minDomainSize, minDomainSize});
    _allLoadBalancer->setup();
  }
#else
  if (_domainIndex == 0 and _loadBalancerOption == LoadBalancerOption::all) {
    std::cout << "ALL loadbalancer has been disabled during compile time. Load balancing will be turned off."
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
  return std::accumulate(_decomposition.begin(), _decomposition.end(), 1, std::multiplies<int>());
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
    auto preceedingNeighborId = _domainId;
    preceedingNeighborId[i] = (--preceedingNeighborId[i] + _decomposition[i]) % _decomposition[i];
    _neighborDomainIndices[neighborIndex] = DomainTools::convertIdToIndex(preceedingNeighborId, _decomposition);

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
        _localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex])
      continue;

    std::vector<ParticleType> particlesForLeftNeighbor{};
    std::vector<ParticleType> particlesForRightNeighbor{};

    collectHaloParticlesForLeftNeighbor(autoPasContainer, dimensionIndex, particlesForLeftNeighbor);
    collectHaloParticlesForRightNeighbor(autoPasContainer, dimensionIndex, particlesForRightNeighbor);

    double leftHaloMin = _localBoxMin[dimensionIndex] - _skinWidth;
    double leftHaloMax = _localBoxMin[dimensionIndex] + _cutoffWidth + _skinWidth;
    double rightHaloMin = _localBoxMax[dimensionIndex] - _cutoffWidth - _skinWidth;
    double rightHaloMax = _localBoxMax[dimensionIndex] + _skinWidth;

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
    int leftNeighbor = _neighborDomainIndices[(dimensionIndex * 2) % _neighborCount];
    int rightNeighbor = _neighborDomainIndices[(dimensionIndex * 2 + 1) % _neighborCount];
    sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbor, particlesForRightNeighbor, leftNeighbor,
                                        rightNeighbor, haloParticles);
  }
  for (const auto &particle : haloParticles) {
    autoPasContainer.addHaloParticle(particle);
  }
}

void RegularGridDecomposition::exchangeMigratingParticles(AutoPasType &autoPasContainer,
                                                          std::vector<ParticleType> &emigrants) {
  for (int dimensionIndex = 0; dimensionIndex < _dimensionCount; ++dimensionIndex) {
    // if this rank spans the whole dimension but it is not periodic -> skip.
    if (_boundaryType[dimensionIndex] != options::BoundaryTypeOption::periodic and
        _localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex] and
        _localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex])
      continue;

    // If the ALL load balancer is used, it may happen that particles migrate to a non adjacent domain.
    // Therefore we need to migrate particles as many times as there are grid cells along the dimension.
    int maximumSendSteps = _loadBalancerOption == LoadBalancerOption::all ? _decomposition[dimensionIndex] : 1;

    for (int gridIndex = 0; gridIndex < maximumSendSteps; ++gridIndex) {
      std::vector<ParticleType> immigrants, remainingEmigrants;
      std::vector<ParticleType> particlesForLeftNeighbor;
      std::vector<ParticleType> particlesForRightNeighbor;

      // See documentation for _neighborDomainIndices to explain the indexing
      int leftNeighbor = _neighborDomainIndices[(dimensionIndex * 2) % _neighborCount];
      int rightNeighbor = _neighborDomainIndices[(dimensionIndex * 2 + 1) % _neighborCount];

      categorizeParticlesIntoLeftAndRightNeighbor(emigrants, dimensionIndex, particlesForLeftNeighbor,
                                                  particlesForRightNeighbor, remainingEmigrants);
      emigrants = remainingEmigrants;

      sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbor, particlesForRightNeighbor, leftNeighbor,
                                          rightNeighbor, immigrants);

      for (const auto &particle : immigrants) {
        if (isInsideLocalDomain(particle.getR())) {
          autoPasContainer.addParticle(particle);
        } else {
          emigrants.push_back(particle);
        }
      }
    }
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
      reflSkinMax[dimensionIndex] = _globalBoxMin[dimensionIndex] + autoPasContainer.getVerletSkin() / 2;

      reflect(false);
    }
    // apply if we are at a global boundary on upper end of the dimension
    if (_localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
      reflSkinMin = _globalBoxMin;
      reflSkinMax = _globalBoxMax;
      reflSkinMin[dimensionIndex] = _globalBoxMax[dimensionIndex] - autoPasContainer.getVerletSkin() / 2;

      reflect(true);
    }
  }
}

void RegularGridDecomposition::sendAndReceiveParticlesLeftAndRight(std::vector<ParticleType> &particlesToLeft,
                                                                   std::vector<ParticleType> &particlesToRight,
                                                                   const int &leftNeighbor, const int &rightNeighbor,
                                                                   std::vector<ParticleType> &receivedParticles) {
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
}

void RegularGridDecomposition::collectHaloParticlesForLeftNeighbor(AutoPasType &autoPasContainer,
                                                                   const size_t &direction,
                                                                   std::vector<ParticleType> &haloParticles) {
  // Calculate halo box for left neighbor
  const std::array<double, _dimensionCount> boxMin = autopas::utils::ArrayMath::subScalar(_localBoxMin, _skinWidth);
  const std::array<double, _dimensionCount> boxMax = [&]() {
    auto boxMax = autopas::utils::ArrayMath::addScalar(_localBoxMax, _skinWidth);
    boxMax[direction] = _localBoxMin[direction] + _cutoffWidth + _skinWidth;
    return boxMax;
  }();

  // Collect the halo particles for the left neighbor
  for (auto particle = autoPasContainer.getRegionIterator(boxMin, boxMax, autopas::IteratorBehavior::owned);
       particle.isValid(); ++particle) {
    std::array<double, _dimensionCount> position = particle->getR();
    haloParticles.push_back(*particle);

    // Apply boundary condition
    if (_localBoxMin[direction] == _globalBoxMin[direction]) {
      position[direction] = position[direction] + (_globalBoxMax[direction] - _globalBoxMin[direction]);
      haloParticles.back().setR(position);
    }
  }
}

void RegularGridDecomposition::collectHaloParticlesForRightNeighbor(AutoPasType &autoPasContainer,
                                                                    const size_t &direction,
                                                                    std::vector<ParticleType> &haloParticles) {
  // Calculate left halo box of right neighbor
  const std::array<double, _dimensionCount> boxMax = autopas::utils::ArrayMath::addScalar(_localBoxMax, _skinWidth);
  const std::array<double, _dimensionCount> boxMin = [&]() {
    auto boxMin = autopas::utils::ArrayMath::subScalar(_localBoxMin, _skinWidth);
    boxMin[direction] = _localBoxMax[direction] - _cutoffWidth - _skinWidth;
    return boxMax;
  }();

  // Collect the halo particles for the right neighbor
  for (auto particle = autoPasContainer.getRegionIterator(boxMin, boxMax, autopas::IteratorBehavior::owned);
       particle.isValid(); ++particle) {
    std::array<double, _dimensionCount> position = particle->getR();
    haloParticles.push_back(*particle);

    // Apply boundary condition
    if (_localBoxMax[direction] == _globalBoxMax[direction]) {
      position[direction] = position[direction] - (_globalBoxMax[direction] - _globalBoxMin[direction]);
      haloParticles.back().setR(position);
    }
  }
}

void RegularGridDecomposition::categorizeParticlesIntoLeftAndRightNeighbor(
    const std::vector<ParticleType> &particles, const size_t &direction,
    std::vector<ParticleType> &leftNeighborParticles, std::vector<ParticleType> &rightNeighborParticles,
    std::vector<ParticleType> &uncategorizedParticles) {
  const std::array<double, _dimensionCount> globalBoxLength =
      autopas::utils::ArrayMath::sub(_globalBoxMax, _globalBoxMin);

  /**
   * The chosen size is the best guess based on the particles vector being distributed into three other vectors.
   */
  leftNeighborParticles.reserve(particles.size() / 3);
  rightNeighborParticles.reserve(particles.size() / 3);
  uncategorizedParticles.reserve(particles.size() / 3);

  for (const auto &particle : particles) {
    std::array<double, _dimensionCount> position = particle.getR();
    if (position[direction] < _localBoxMin[direction]) {
      leftNeighborParticles.push_back(particle);

      // Apply boundary condition
      if (_localBoxMin[direction] == _globalBoxMin[direction]) {
        position[direction] = std::min(std::nextafter(_globalBoxMax[direction], _globalBoxMin[direction]),
                                       position[direction] + globalBoxLength[direction]);
        leftNeighborParticles.back().setR(position);
      }
    } else if (position[direction] >= _localBoxMax[direction]) {
      rightNeighborParticles.push_back(particle);

      // Apply boundary condition
      if (_localBoxMax[direction] == _globalBoxMax[direction]) {
        position[direction] = std::max(_globalBoxMin[direction], position[direction] - globalBoxLength[direction]);
        rightNeighborParticles.back().setR(position);
      }
    } else {
      uncategorizedParticles.push_back(particle);
    }
  }
}

void RegularGridDecomposition::balanceWithInvertedPressureLoadBalancer(const double &work) {
  // This is a dummy variable which is not being used. It is required by the non-blocking MPI_Send calls.
  autopas::AutoPas_MPI_Request dummyRequest;

  auto oldLocalBoxMin = _localBoxMin;
  auto oldLocalBoxMax = _localBoxMax;

  std::array<double, 3> averageWorkInPlane{};

  for (int i = 0; i < _dimensionCount; ++i) {
    const int domainCountInPlane =
        _decomposition[(i + 1) % _dimensionCount] * _decomposition[(i + 2) % _dimensionCount];

    // Calculate the average work in the process grid plane
    averageWorkInPlane[i] = work;
    if (domainCountInPlane > 1) {
      autopas::AutoPas_MPI_Allreduce(&work, &averageWorkInPlane[i], 1, AUTOPAS_MPI_DOUBLE, AUTOPAS_MPI_SUM,
                                     _planarCommunicators[i]);
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

  for (int i = 0; i < _dimensionCount; ++i) {
    // Get neighbour indices
    const int leftNeighbor = _neighborDomainIndices[i * 2];
    const int rightNeighbor = _neighborDomainIndices[i * 2 + 1];

    double neighborPlaneWork, neighborBoundary, balancedPosition;
    if (_localBoxMin[i] != _globalBoxMin[i]) {
      // Receive average work from neighbour planes.
      autopas::AutoPas_MPI_Recv(&neighborPlaneWork, 1, AUTOPAS_MPI_DOUBLE, leftNeighbor, 0, _communicator,
                                AUTOPAS_MPI_STATUS_IGNORE);
      autopas::AutoPas_MPI_Recv(&neighborBoundary, 1, AUTOPAS_MPI_DOUBLE, leftNeighbor, 0, _communicator,
                                AUTOPAS_MPI_STATUS_IGNORE);

      // Calculate balanced positions and only shift by half the resulting distance to prevent localBox min to be larger
      // than localBoxMax.
      balancedPosition = DomainTools::balanceAdjacentDomains(neighborPlaneWork, averageWorkInPlane[i], neighborBoundary,
                                                             oldLocalBoxMax[i], 2 * (_cutoffWidth + _skinWidth));
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
      balancedPosition =
          DomainTools::balanceAdjacentDomains(averageWorkInPlane[i], neighborPlaneWork, oldLocalBoxMin[i],
                                              neighborBoundary, 2 * (_cutoffWidth + _skinWidth));
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
