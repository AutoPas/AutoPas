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
#include "src/ParticleSerializationTools.h"

RegularGridDecomposition::RegularGridDecomposition(const MDFlexConfig &configuration)
    : _cutoffWidth(configuration.cutoff.value), _skinWidth(configuration.verletSkinRadius.value) {
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &_subdomainCount);

#if defined(AUTOPAS_INCLUDE_MPI)
  _mpiCommunicationNeeded = true;
#else
  _mpiCommunicationNeeded = false;
#endif

  if (_subdomainCount == 1) {
    _mpiCommunicationNeeded = false;
  }

  DomainTools::generateDecomposition(_subdomainCount, configuration.subdivideDimension.value, _decomposition);

  initializeMPICommunicator();

  initializeLocalDomain();

  initializeGlobalBox(configuration.boxMin.value, configuration.boxMax.value);

  initializeLocalBox();

  initializeNeighborIds();

  _loadBalancer = configuration.loadBalancer.value;

  if (_loadBalancer == LoadBalancerOption::all) {
    _allLoadBalancer = std::make_unique<ALL::ALL<double, double>>(ALL::TENSOR, _dimensionCount, 0);
    _allLoadBalancer->setCommunicator(_communicator);

    const double minDomainSize = 2 * (_cutoffWidth + _skinWidth);
    _allLoadBalancer->setMinDomainSize({minDomainSize, minDomainSize, minDomainSize});
    _allLoadBalancer->setup();
  }
}

RegularGridDecomposition::~RegularGridDecomposition() {}

void RegularGridDecomposition::update(const double &work) {
  if (_mpiCommunicationNeeded) {
    switch (_loadBalancer) {
      case LoadBalancerOption::invertedPressure: {
        balanceWithInvertedPressureLoadBalancer(work);
        break;
      }
      case LoadBalancerOption::all: {
        balanceWithAllLoadBalancer(work);
        break;
      }
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
  std::vector<int> periods(3, 1);
  autopas::AutoPas_MPI_Cart_create(AUTOPAS_MPI_COMM_WORLD, 3, _decomposition.data(), periods.data(), true,
                                   &_communicator);
  autopas::AutoPas_MPI_Comm_rank(_communicator, &_domainIndex);
}

void RegularGridDecomposition::initializeLocalDomain() {
  _domainId = {0, 0, 0};
  autopas::AutoPas_MPI_Comm_rank(_communicator, &_domainIndex);

  std::vector<int> periods(3, 1);
  autopas::AutoPas_MPI_Cart_get(_communicator, 3, _decomposition.data(), periods.data(), _domainId.data());

  // Create planar communicators used for diffuse load balancing.
  for (int i = 0; i < _dimensionCount; ++i) {
    if (_mpiCommunicationNeeded) {
      const size_t key = _decomposition[(i + 1) % _dimensionCount] * _domainId[(i + 2) % _dimensionCount] +
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

    if (_domainId[i] == 0) {
      _localBoxMin[i] = _globalBoxMin[i];
    } else if (_domainId[i] == _decomposition[i] - 1) {
      _localBoxMax[i] = _globalBoxMax[i];
    }
  }
}

void RegularGridDecomposition::initializeNeighborIds() {
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

void RegularGridDecomposition::initializeGlobalBox(const std::array<double, 3> &globalBoxMin,
                                                   const std::array<double, 3> &globalBoxMax) {
  for (int i = 0; i < 3; ++i) {
    _globalBoxMin[i] = globalBoxMin[i];
    _globalBoxMax[i] = globalBoxMax[i];
  }
}

bool RegularGridDecomposition::isInsideLocalDomain(const std::array<double, 3> &coordinates) const {
  return DomainTools::isInsideDomain(coordinates, _localBoxMin, _localBoxMax);
}

std::array<int, 6> RegularGridDecomposition::getExtentOfSubdomain(const int subdomainIndex) const {
  return DomainTools::getExtentOfSubdomain(subdomainIndex, _decomposition);
}

void RegularGridDecomposition::exchangeHaloParticles(SharedAutoPasContainer &autoPasContainer) {
  std::vector<ParticleType> haloParticles{};

  for (int dimensionIndex = 0; dimensionIndex < _dimensionCount; ++dimensionIndex) {
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

      // See documentation for _neighborDomainIndices to explain the indexing
      int leftNeighbor = _neighborDomainIndices[(dimensionIndex * 2) % _neighborCount];
      int rightNeighbor = _neighborDomainIndices[(dimensionIndex * 2 + 1) % _neighborCount];
      sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbor, particlesForRightNeighbor, leftNeighbor,
                                          rightNeighbor, haloParticles);
    }
  }

  for (const auto &particle : haloParticles) {
    autoPasContainer->addHaloParticle(particle);
  }
}

void RegularGridDecomposition::exchangeMigratingParticles(SharedAutoPasContainer &autoPasContainer,
                                                          std::vector<ParticleType> &emigrants) {
  for (int dimensionIndex = 0; dimensionIndex < _dimensionCount; ++dimensionIndex) {
    // If the ALL load balancer is used, it may happen that particles migrate to a non adjacent domain.
    // Therefore we need to migrate particles as many times as there are grid cells along the dimension.
    int maximumSendSteps = _loadBalancer == LoadBalancerOption::all ? _decomposition[dimensionIndex] : 1;

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
          autoPasContainer->addParticle(particle);
        } else {
          emigrants.push_back(particle);
        }
      }

      immigrants.clear();
    }
  }

  // Add remaining emigrants to current container
  for (const auto &particle : emigrants) {
    autoPasContainer->addParticle(particle);
  }
}

void RegularGridDecomposition::sendParticles(const std::vector<ParticleType> &particles, const int &receiver) {
  std::vector<char> buffer;

  for (const auto &particle : particles) {
    ParticleSerializationTools::serializeParticle(particle, buffer);
  }

  sendDataToNeighbor(buffer, receiver);
}

void RegularGridDecomposition::receiveParticles(std::vector<ParticleType> &receivedParticles, const int &source) {
  std::vector<char> receiveBuffer;

  receiveDataFromNeighbor(source, receiveBuffer);

  if (!receiveBuffer.empty()) {
    ParticleSerializationTools::deserializeParticles(receiveBuffer, receivedParticles);
  }
}

void RegularGridDecomposition::sendDataToNeighbor(std::vector<char> sendBuffer, const int &neighbor) {
  _sendBuffers.push_back(sendBuffer);

  autopas::AutoPas_MPI_Request sendRequest;
  _sendRequests.push_back(sendRequest);

  autopas::AutoPas_MPI_Isend(_sendBuffers.back().data(), _sendBuffers.back().size(), AUTOPAS_MPI_CHAR, neighbor, 0,
                             _communicator, &_sendRequests.back());
}

void RegularGridDecomposition::receiveDataFromNeighbor(const int &neighbor, std::vector<char> &receiveBuffer) {
  autopas::AutoPas_MPI_Status status;
  autopas::AutoPas_MPI_Probe(neighbor, 0, _communicator, &status);

  int receiveBufferSize;
  autopas::AutoPas_MPI_Get_count(&status, AUTOPAS_MPI_CHAR, &receiveBufferSize);
  receiveBuffer.resize(receiveBufferSize);

  autopas::AutoPas_MPI_Recv(receiveBuffer.data(), receiveBufferSize, AUTOPAS_MPI_CHAR, neighbor, 0, _communicator,
                            AUTOPAS_MPI_STATUS_IGNORE);
}

void RegularGridDecomposition::sendAndReceiveParticlesLeftAndRight(std::vector<ParticleType> &particlesToLeft,
                                                                   std::vector<ParticleType> &particlesToRight,
                                                                   const int &leftNeighbor, const int &rightNeighbor,
                                                                   std::vector<ParticleType> &receivedParticles) {
  if (_mpiCommunicationNeeded and leftNeighbor != _domainIndex) {
    sendParticles(particlesToLeft, leftNeighbor);
    sendParticles(particlesToRight, rightNeighbor);

    receiveParticles(receivedParticles, leftNeighbor);
    receiveParticles(receivedParticles, rightNeighbor);

    waitForSendRequests();
  } else {
    receivedParticles.insert(receivedParticles.end(), particlesToLeft.begin(), particlesToLeft.end());
    receivedParticles.insert(receivedParticles.end(), particlesToRight.begin(), particlesToRight.end());
  }
}

void RegularGridDecomposition::waitForSendRequests() {
  std::vector<autopas::AutoPas_MPI_Status> sendStates;
  sendStates.resize(_sendRequests.size());
  autopas::AutoPas_MPI_Waitall(_sendRequests.size(), _sendRequests.data(), sendStates.data());
  _sendRequests.clear();
  _sendBuffers.clear();
}

void RegularGridDecomposition::collectHaloParticlesForLeftNeighbor(SharedAutoPasContainer &autoPasContainer,
                                                                   const size_t &direction,
                                                                   std::vector<ParticleType> &haloParticles) {
  std::array<double, _dimensionCount> boxMin, boxMax;

  // Calculate halo box for left neighbor
  for (int i = 0; i < _dimensionCount; ++i) {
    boxMin[i] = _localBoxMin[i] - _skinWidth;
    boxMax[i] = _skinWidth + _localBoxMax[i];
  }
  boxMax[direction] = _localBoxMin[direction] + _cutoffWidth + _skinWidth;

  // Collect the halo particles for the left neighbor
  for (auto particle = autoPasContainer->getRegionIterator(boxMin, boxMax, autopas::IteratorBehavior::owned);
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

void RegularGridDecomposition::collectHaloParticlesForRightNeighbor(SharedAutoPasContainer &autoPasContainer,
                                                                    const size_t &direction,
                                                                    std::vector<ParticleType> &haloParticles) {
  std::array<double, _dimensionCount> boxMin, boxMax;

  // Calculate left halo box of right neighbor
  for (int i = 0; i < _dimensionCount; ++i) {
    boxMin[i] = _localBoxMin[i] - _skinWidth;
    boxMax[i] = _localBoxMax[i] + _skinWidth;
  }
  boxMin[direction] = _localBoxMax[direction] - _cutoffWidth - _skinWidth;

  // Collect the halo particles for the right neighbor
  for (auto particle = autoPasContainer->getRegionIterator(boxMin, boxMax, autopas::IteratorBehavior::owned);
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
  // This is a dummy variable which is not being used.
  autopas::AutoPas_MPI_Request dummyRequest;

  auto oldLocalBoxMin = _localBoxMin;
  auto oldLocalBoxMax = _localBoxMax;

  std::array<double, 3> distributedWorkInPlane{};

  for (int i = 0; i < _dimensionCount; ++i) {
    const int domainCountInPlane =
        _decomposition[(i + 1) % _dimensionCount] * _decomposition[(i + 2) % _dimensionCount];

    distributedWorkInPlane[i] = work;
    if (domainCountInPlane > 1) {
      autopas::AutoPas_MPI_Allreduce(&work, &distributedWorkInPlane[i], 1, AUTOPAS_MPI_DOUBLE, AUTOPAS_MPI_SUM,
                                     _planarCommunicators[i]);
      distributedWorkInPlane[i] = distributedWorkInPlane[i] / domainCountInPlane;
    }

    const int leftNeighbor = _neighborDomainIndices[i * 2];
    const int rightNeighbor = _neighborDomainIndices[i * 2 + 1];

    if (_localBoxMin[i] != _globalBoxMin[i]) {
      autopas::AutoPas_MPI_Isend(&distributedWorkInPlane[i], 1, AUTOPAS_MPI_DOUBLE, leftNeighbor, 0, _communicator,
                                 &dummyRequest);
      autopas::AutoPas_MPI_Isend(&oldLocalBoxMax[i], 1, AUTOPAS_MPI_DOUBLE, leftNeighbor, 0, _communicator,
                                 &dummyRequest);
    }

    if (_localBoxMax[i] != _globalBoxMax[i]) {
      autopas::AutoPas_MPI_Isend(&distributedWorkInPlane[i], 1, AUTOPAS_MPI_DOUBLE, rightNeighbor, 0, _communicator,
                                 &dummyRequest);
      autopas::AutoPas_MPI_Isend(&oldLocalBoxMin[i], 1, AUTOPAS_MPI_DOUBLE, rightNeighbor, 0, _communicator,
                                 &dummyRequest);
    }
  }

  for (int i = 0; i < _dimensionCount; ++i) {
    const int leftNeighbor = _neighborDomainIndices[i * 2];
    const int rightNeighbor = _neighborDomainIndices[i * 2 + 1];

    double neighborPlaneWork, neighborBoundary, balancedPosition;
    if (_localBoxMin[i] != _globalBoxMin[i]) {
      autopas::AutoPas_MPI_Recv(&neighborPlaneWork, 1, AUTOPAS_MPI_DOUBLE, leftNeighbor, 0, _communicator,
                                AUTOPAS_MPI_STATUS_IGNORE);
      autopas::AutoPas_MPI_Recv(&neighborBoundary, 1, AUTOPAS_MPI_DOUBLE, leftNeighbor, 0, _communicator,
                                AUTOPAS_MPI_STATUS_IGNORE);

      balancedPosition =
          DomainTools::balanceAdjacentDomains(neighborPlaneWork, distributedWorkInPlane[i], neighborBoundary,
                                              oldLocalBoxMax[i], 2 * (_cutoffWidth + _skinWidth));
      _localBoxMin[i] += (balancedPosition - _localBoxMin[i]) / 2;
    }

    if (_localBoxMax[i] != _globalBoxMax[i]) {
      double neighborPlaneWork, neighborBoundary;
      autopas::AutoPas_MPI_Recv(&neighborPlaneWork, 1, AUTOPAS_MPI_DOUBLE, rightNeighbor, 0, _communicator,
                                AUTOPAS_MPI_STATUS_IGNORE);
      autopas::AutoPas_MPI_Recv(&neighborBoundary, 1, AUTOPAS_MPI_DOUBLE, rightNeighbor, 0, _communicator,
                                AUTOPAS_MPI_STATUS_IGNORE);

      balancedPosition =
          DomainTools::balanceAdjacentDomains(distributedWorkInPlane[i], neighborPlaneWork, oldLocalBoxMin[i],
                                              neighborBoundary, 2 * (_cutoffWidth + _skinWidth));
      _localBoxMax[i] += (balancedPosition - _localBoxMax[i]) / 2;
    }
  }
}

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
