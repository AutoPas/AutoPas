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

RegularGridDecomposition::RegularGridDecomposition(const std::array<double, 3> &globalBoxMin,
                                                   const std::array<double, 3> &globalBoxMax,
                                                   const std::array<bool, 3> subdivideDimension,
                                                   const double &cutoffWidth, const double &skinWidth)
    : _cutoffWidth(cutoffWidth), _skinWidth(skinWidth) {
#if defined(AUTOPAS_INCLUDE_MPI)
  _mpiIsEnabled = true;
#else
  _mpiIsEnabled = false;
#endif

  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &_subdomainCount);

  if (_subdomainCount == 1) {
    _mpiIsEnabled = false;
  }

  if (_mpiIsEnabled) {
    std::cout << "MPI will be used." << std::endl;
  } else {
    std::cout << "MPI will not be used." << std::endl;
  }

  DomainTools::generateDecomposition(_subdomainCount, subdivideDimension, _decomposition);

  std::cout << "Decomposition: " << autopas::utils::ArrayUtils::to_string(_decomposition) << std::endl;

  initializeMPICommunicator();

  initializeLocalDomain();

  initializeGlobalBox(globalBoxMin, globalBoxMax);

  initializeLocalBox();

  initializeNeighbourIds();
}

RegularGridDecomposition::~RegularGridDecomposition() {}

void RegularGridDecomposition::update(const double &work) {
  if (_mpiIsEnabled) {
    // This corresponds to the work for each direction (dimension * 2).
    const double distributedWork = calculateDistributedWork(work);

    // This is a dummy variable which is not being used.
    autopas::AutoPas_MPI_Request dummyRequest;

    for (int i = 0; i < _dimensionCount; ++i) {
      // Calculate average distributed work in planar communicator.
      const int domainCountAlongDimension = _decomposition[(i + 1) % _dimensionCount];

      double distributedWorkInPlane = distributedWork;
      if (domainCountAlongDimension > 1) {
        autopas::AutoPas_MPI_Allreduce(&distributedWork, &distributedWorkInPlane, 1, AUTOPAS_MPI_DOUBLE,
                                       AUTOPAS_MPI_SUM, _planarCommunicators[i]);
        distributedWorkInPlane = distributedWorkInPlane / domainCountAlongDimension;
      }

      const int leftNeighbour = _neighbourDomainIndices[i * 2];
      const int rightNeighbour = _neighbourDomainIndices[i * 2 + 1];

      if (_localBoxMin[i] != _globalBoxMin[i]) {
        autopas::AutoPas_MPI_Isend(&distributedWorkInPlane, 1, AUTOPAS_MPI_DOUBLE, leftNeighbour, 0, _communicator,
                                   &dummyRequest);
        autopas::AutoPas_MPI_Isend(&_localBoxMax[i], 1, AUTOPAS_MPI_DOUBLE, leftNeighbour, 0, _communicator,
                                   &dummyRequest);
      }

      if (_localBoxMax[i] != _globalBoxMax[i]) {
        autopas::AutoPas_MPI_Isend(&distributedWorkInPlane, 1, AUTOPAS_MPI_DOUBLE, rightNeighbour, 0, _communicator,
                                   &dummyRequest);
        autopas::AutoPas_MPI_Isend(&_localBoxMin[i], 1, AUTOPAS_MPI_DOUBLE, rightNeighbour, 0, _communicator,
                                   &dummyRequest);
      }

      double balancedPosition;
      if (_localBoxMin[i] != _globalBoxMin[i]) {
        double neighbourPlaneWork, neighbourBoundary;
        autopas::AutoPas_MPI_Recv(&neighbourPlaneWork, 1, AUTOPAS_MPI_DOUBLE, leftNeighbour, 0, _communicator,
                                  AUTOPAS_MPI_STATUS_IGNORE);
        autopas::AutoPas_MPI_Recv(&neighbourBoundary, 1, AUTOPAS_MPI_DOUBLE, leftNeighbour, 0, _communicator,
                                  AUTOPAS_MPI_STATUS_IGNORE);

        balancedPosition =
            DomainTools::balanceAdjacentDomains(neighbourPlaneWork, distributedWorkInPlane, neighbourBoundary,
                                                _localBoxMax[i], 2 * (_cutoffWidth + _skinWidth));
        _localBoxMin[i] += (balancedPosition - _localBoxMin[i]) / 2;
      }

      if (_localBoxMax[i] != _globalBoxMax[i]) {
        double neighbourPlaneWork, neighbourBoundary;
        autopas::AutoPas_MPI_Recv(&neighbourPlaneWork, 1, AUTOPAS_MPI_DOUBLE, rightNeighbour, 0, _communicator,
                                  AUTOPAS_MPI_STATUS_IGNORE);
        autopas::AutoPas_MPI_Recv(&neighbourBoundary, 1, AUTOPAS_MPI_DOUBLE, rightNeighbour, 0, _communicator,
                                  AUTOPAS_MPI_STATUS_IGNORE);

        balancedPosition =
            DomainTools::balanceAdjacentDomains(distributedWorkInPlane, neighbourPlaneWork, _localBoxMin[i],
                                                neighbourBoundary, 2 * (_cutoffWidth + _skinWidth));
        _localBoxMax[i] += (balancedPosition - _localBoxMax[i]) / 2;
      }
    }
  }
}

void RegularGridDecomposition::initializeMPICommunicator() {
  std::vector<int> periods(3, 1);
  autopas::AutoPas_MPI_Cart_create(AUTOPAS_MPI_COMM_WORLD, 3, _decomposition.data(), periods.data(), true,
                                   &_communicator);
  autopas::AutoPas_MPI_Comm_rank(_communicator, &_domainIndex);

  // Create planar communicators used for diffuse load balancing.
  for (int i = 0; i < _dimensionCount; ++i) {
    if (_mpiIsEnabled) {
      autopas::AutoPas_MPI_Comm_split(_communicator, _domainId[i], _domainId[i], &_planarCommunicators[i]);
    } else {
      _planarCommunicators[i] = AUTOPAS_MPI_COMM_WORLD;
    }
  }
}

void RegularGridDecomposition::initializeLocalDomain() {
  _domainId = {0, 0, 0};
  autopas::AutoPas_MPI_Comm_rank(_communicator, &_domainIndex);

  std::vector<int> periods(3, 1);
  autopas::AutoPas_MPI_Cart_get(_communicator, 3, _decomposition.data(), periods.data(), _domainId.data());
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

void RegularGridDecomposition::initializeNeighbourIds() {
  for (int i = 0; i < 3; ++i) {
    auto neighbourIndex = i * 2;
    auto preceedingNeighbourId = _domainId;
    preceedingNeighbourId[i] = (--preceedingNeighbourId[i] + _decomposition[i]) % _decomposition[i];
    _neighbourDomainIndices[neighbourIndex] = DomainTools::convertIdToIndex(preceedingNeighbourId, _decomposition);

    ++neighbourIndex;
    auto succeedingNeighbourId = _domainId;
    succeedingNeighbourId[i] = (++succeedingNeighbourId[i] + _decomposition[i]) % _decomposition[i];
    _neighbourDomainIndices[neighbourIndex] = DomainTools::convertIdToIndex(succeedingNeighbourId, _decomposition);
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
  for (int i = 0; i < _dimensionCount; ++i) {
    std::vector<ParticleType> particlesForLeftNeighbour, particlesForRightNeighbour, haloParticles;

    std::array<double, _dimensionCount> haloBoxMin, haloBoxMax;
    haloBoxMin = {_localBoxMin[0] - _skinWidth, _localBoxMin[1] - _skinWidth, _localBoxMin[2] - _skinWidth};
    for (int i = 0; i < _dimensionCount; ++i) {
      haloBoxMax[i] = _localBoxMin[i] + _skinWidth + (_localBoxMax[i] - _localBoxMin[i]);
    }
    haloBoxMax[i] = _localBoxMin[i] + _cutoffWidth + _skinWidth;

    for (auto particle = autoPasContainer->getRegionIterator(haloBoxMin, haloBoxMax, autopas::IteratorBehavior::owned);
         particle.isValid(); ++particle) {
      std::array<double, _dimensionCount> position = particle->getR();
      particlesForLeftNeighbour.push_back(*particle);

      // Apply boundary condition
      if (_localBoxMin[i] == _globalBoxMin[i]) {
        position[i] = position[i] + (_globalBoxMax[i] - _globalBoxMin[i]);
        particlesForLeftNeighbour.back().setR(position);
      }
    }

    for (int i = 0; i < _dimensionCount; ++i) {
      haloBoxMin[i] = _localBoxMax[i] - _skinWidth - (_localBoxMax[i] - _localBoxMin[i]);
    }
    haloBoxMin[i] = _localBoxMax[i] - _cutoffWidth - _skinWidth;
    haloBoxMax = {_localBoxMax[0] + _skinWidth, _localBoxMax[1] + _skinWidth, _localBoxMax[2] + _skinWidth};

    for (auto particle = autoPasContainer->getRegionIterator(haloBoxMin, haloBoxMax, autopas::IteratorBehavior::owned);
         particle.isValid(); ++particle) {
      std::array<double, _dimensionCount> position = particle->getR();
      particlesForRightNeighbour.push_back(*particle);

      // Apply boundary condition
      if (_localBoxMax[i] == _globalBoxMax[i]) {
        position[i] = position[i] - (_globalBoxMax[i] - _globalBoxMin[i]);
        particlesForRightNeighbour.back().setR(position);
      }
    }

    // See documentation for _neighbourDomainIndices to explain the indexing
    int leftNeighbour = _neighbourDomainIndices[(i * 2) % _neighbourCount];
    int rightNeighbour = _neighbourDomainIndices[(i * 2 + 1) % _neighbourCount];
    sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbour, particlesForRightNeighbour, leftNeighbour,
                                        rightNeighbour, haloParticles);

    for (auto &particle : haloParticles) {
      autoPasContainer->addOrUpdateHaloParticle(particle);
    }

    particlesForLeftNeighbour.clear();
    particlesForRightNeighbour.clear();

    // index of next dimension
    int j = (i + 1) % _dimensionCount;
    double leftHaloMin = _localBoxMin[j] - _skinWidth;
    double leftHaloMax = _localBoxMin[j] + _cutoffWidth + _skinWidth;
    double rightHaloMin = _localBoxMax[j] - _cutoffWidth - _skinWidth;
    double rightHaloMax = _localBoxMax[j] + _skinWidth;

    for (const auto &particle : haloParticles) {
      std::array<double, _dimensionCount> position = particle.getR();
      if (position[j] >= leftHaloMin && position[j] < leftHaloMax) {
        particlesForLeftNeighbour.push_back(particle);

        // Apply boundary condition
        if (_localBoxMin[j] == _globalBoxMin[j]) {
          position[j] = position[j] + (_globalBoxMax[j] - _globalBoxMin[j]);
          particlesForLeftNeighbour.back().setR(position);
        }
      } else if (position[j] >= rightHaloMin && position[j] < rightHaloMax) {
        particlesForRightNeighbour.push_back(particle);

        // Apply boundary condition
        if (_localBoxMax[j] == _globalBoxMax[j]) {
          position[j] = position[j] - (_globalBoxMax[j] - _globalBoxMin[j]);
          particlesForRightNeighbour.back().setR(position);
        }
      }
    }

    haloParticles.clear();

    // See documentation for _neighbourDomainIndices to explain the indexing
    leftNeighbour = _neighbourDomainIndices[(j * 2) % _neighbourCount];
    rightNeighbour = _neighbourDomainIndices[(j * 2 + 1) % _neighbourCount];
    sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbour, particlesForRightNeighbour, leftNeighbour,
                                        rightNeighbour, haloParticles);

    for (auto &particle : haloParticles) {
      autoPasContainer->addOrUpdateHaloParticle(particle);
    }
  }
}

void RegularGridDecomposition::exchangeMigratingParticles(SharedAutoPasContainer &autoPasContainer,
                                                          std::vector<ParticleType> &emigrants) {
  const std::array<double, _dimensionCount> globalBoxMin = {_globalBoxMin[0], _globalBoxMin[1], _globalBoxMin[2]};
  const std::array<double, _dimensionCount> globalBoxMax = {_globalBoxMax[0], _globalBoxMax[1], _globalBoxMax[2]};
  const std::array<double, _dimensionCount> globalBoxLength =
      autopas::utils::ArrayMath::sub(globalBoxMax, globalBoxMin);

  for (int i = 0; i < _dimensionCount; ++i) {
    std::vector<ParticleType> immigrants, migrants, remainingEmigrants;

    std::vector<ParticleType> particlesForLeftNeighbour;
    std::vector<ParticleType> particlesForRightNeighbour;

    // See documentation for _neighbourDomainIndices to explain the indexing
    int leftNeighbour = _neighbourDomainIndices[(i * 2) % _neighbourCount];
    int rightNeighbour = _neighbourDomainIndices[(i * 2 + 1) % _neighbourCount];

    std::array<double, _dimensionCount> position;
    for (const auto &particle : emigrants) {
      position = particle.getR();
      if (position[i] < _localBoxMin[i]) {
        particlesForLeftNeighbour.push_back(particle);

        // Apply boundary condition
        if (_localBoxMin[i] == _globalBoxMin[i]) {
          position[i] = std::min(std::nextafter(_globalBoxMax[i], _globalBoxMin[i]), position[i] + globalBoxLength[i]);
          particlesForLeftNeighbour.back().setR(position);
        }
      } else if (position[i] >= _localBoxMax[i]) {
        particlesForRightNeighbour.push_back(particle);

        // Apply boundary condition
        if (_localBoxMax[i] == _globalBoxMax[i]) {
          position[i] = std::max(_globalBoxMin[i], position[i] - globalBoxLength[i]);
          particlesForRightNeighbour.back().setR(position);
        }
      } else {
        remainingEmigrants.push_back(particle);
      }
    }
    emigrants = remainingEmigrants;

    sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbour, particlesForRightNeighbour, leftNeighbour,
                                        rightNeighbour, immigrants);

    particlesForLeftNeighbour.clear();
    particlesForRightNeighbour.clear();

    for (const auto &particle : immigrants) {
      if (isInsideLocalDomain(particle.getR())) {
        autoPasContainer->addParticle(particle);
      } else {
        migrants.push_back(particle);
      }
    }

    immigrants.clear();

    // index of next dimension
    int j = (i + 1) % _dimensionCount;

    // See documentation for _neighbourDomainIndices to explain the indexing
    leftNeighbour = _neighbourDomainIndices[(j * 2) % _neighbourCount];
    rightNeighbour = _neighbourDomainIndices[(j * 2 + 1) % _neighbourCount];

    for (const auto &particle : migrants) {
      std::array<double, _dimensionCount> position = particle.getR();
      if (position[j] < _localBoxMin[j]) {
        particlesForLeftNeighbour.push_back(particle);

        // Apply boundary condition
        if (_localBoxMin[j] == _globalBoxMin[j]) {
          position[j] = std::min(std::nextafter(_globalBoxMax[j], _globalBoxMin[j]), position[j] + globalBoxLength[j]);
          particlesForLeftNeighbour.back().setR(position);
        }
      } else if (position[j] >= _localBoxMax[j]) {
        particlesForRightNeighbour.push_back(particle);

        // Apply boundary condition
        if (_localBoxMax[j] == _globalBoxMax[j]) {
          position[j] = std::max(_globalBoxMin[j], position[j] - globalBoxLength[j]);
          particlesForRightNeighbour.back().setR(position);
        }
      }
    }

    sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbour, particlesForRightNeighbour, leftNeighbour,
                                        rightNeighbour, immigrants);

    for (auto &particle : immigrants) {
      autoPasContainer->addParticle(particle);
    }

    waitForSendRequests();
  }

  // Add remaining emigrants to current container
  for (auto &particle : emigrants) {
    autoPasContainer->addParticle(particle);
  }
}

void RegularGridDecomposition::sendParticles(const std::vector<ParticleType> &particles, const int &receiver) {
  std::vector<char> buffer;

  for (auto &particle : particles) {
    ParticleSerializationTools::serializeParticle(particle, buffer);
  }

  size_t sizeOfParticleAttributes = sizeof(ParticleAttributes);

  sendDataToNeighbour(buffer, receiver);
}

void RegularGridDecomposition::receiveParticles(std::vector<ParticleType> &receivedParticles, const int &source) {
  std::vector<char> receiveBuffer;

  receiveDataFromNeighbour(source, receiveBuffer);

  if (!receiveBuffer.empty()) {
    ParticleSerializationTools::deserializeParticles(receiveBuffer, receivedParticles);
  }
}

void RegularGridDecomposition::sendDataToNeighbour(std::vector<char> sendBuffer, const int &neighbour) {
  _sendBuffers.push_back(sendBuffer);

  autopas::AutoPas_MPI_Request sendRequest;
  _sendRequests.push_back(sendRequest);

  autopas::AutoPas_MPI_Isend(_sendBuffers.back().data(), _sendBuffers.back().size(), AUTOPAS_MPI_CHAR, neighbour, 0,
                             _communicator, &_sendRequests.back());
}

void RegularGridDecomposition::receiveDataFromNeighbour(const int &neighbour, std::vector<char> &receiveBuffer) {
  autopas::AutoPas_MPI_Status status;
  autopas::AutoPas_MPI_Probe(neighbour, 0, _communicator, &status);

  int receiveBufferSize;
  autopas::AutoPas_MPI_Get_count(&status, AUTOPAS_MPI_CHAR, &receiveBufferSize);
  receiveBuffer.resize(receiveBufferSize);

  autopas::AutoPas_MPI_Recv(receiveBuffer.data(), receiveBufferSize, AUTOPAS_MPI_CHAR, neighbour, 0, _communicator,
                            AUTOPAS_MPI_STATUS_IGNORE);
}

void RegularGridDecomposition::sendAndReceiveParticlesLeftAndRight(const std::vector<ParticleType> &particlesToLeft,
                                                                   const std::vector<ParticleType> &particlesToRight,
                                                                   const int &leftNeighbour, const int &rightNeighbour,
                                                                   std::vector<ParticleType> &receivedParticles) {
  if (_mpiIsEnabled && leftNeighbour != _domainIndex) {
    sendParticles(particlesToLeft, leftNeighbour);
    sendParticles(particlesToRight, rightNeighbour);

    receiveParticles(receivedParticles, leftNeighbour);
    receiveParticles(receivedParticles, rightNeighbour);

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

int RegularGridDecomposition::convertIdToIndex(const std::array<int, 3> &domainId) {
  int neighbourDomainIndex = 0;

  for (int i = 0; i < 3; ++i) {
    int accumulatedTail = 1;

    if (i < _decomposition.size() - 1) {
      accumulatedTail =
          std::accumulate(_decomposition.begin() + i + 1, _decomposition.end(), 1, std::multiplies<int>());
    }

    neighbourDomainIndex += accumulatedTail * domainId[i];
  }

  return neighbourDomainIndex;
}

double RegularGridDecomposition::calculateDistributedWork(const double work) {
  size_t numberOfShiftableBoundaries = 0;
  for (int i = 0; i < _dimensionCount; ++i) {
    if (_localBoxMin[i] != _globalBoxMin[i]) {
      numberOfShiftableBoundaries++;
    }
    if (_localBoxMax[i] != _globalBoxMax[i]) {
      numberOfShiftableBoundaries++;
    }
  }

  return work / (numberOfShiftableBoundaries != 0 ? numberOfShiftableBoundaries : 1);
}
