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
                                                   const std::array<bool, 3> &subdivideDimension, double cutoffWidth,
                                                   double skinWidthPerTimestep, double rebuildFrequency,
                                                   const std::array<options::BoundaryTypeOption, 3> &boundaryConditions)
    : _cutoffWidth(cutoffWidth),
      _skinWidthPerTimestep(skinWidthPerTimestep),
      _skinWidth(_skinWidthPerTimestep*rebuildFrequency), 
      _globalBoxMin(globalBoxMin),
      _globalBoxMax(globalBoxMax),
      _boundaryType(boundaryConditions) {
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &_subdomainCount);

  int rank;
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);

#if defined(AUTOPAS_INCLUDE_MPI)
  _mpiCommunicationNeeded = true;
#else
  _mpiCommunicationNeeded = false;
#endif

  if (_subdomainCount == 1) {
    _mpiCommunicationNeeded = false;
  }

  DomainTools::generateDecomposition(_subdomainCount, subdivideDimension, _decomposition);

  initializeMPICommunicator();

  initializeLocalDomain();

  initializeLocalBox();

  initializeNeighbourIds();
}

RegularGridDecomposition::~RegularGridDecomposition() {}

int RegularGridDecomposition::getNumberOfSubdomains() const {
  return std::accumulate(_decomposition.begin(), _decomposition.end(), 1, std::multiplies<int>());
}

void RegularGridDecomposition::update() { updateLocalBox(); }

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
}

void RegularGridDecomposition::initializeLocalBox() { updateLocalBox(); }

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

void RegularGridDecomposition::updateLocalBox() {
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

bool RegularGridDecomposition::isInsideLocalDomain(const std::array<double, 3> &coordinates) const {
  return DomainTools::isInsideDomain(coordinates, _localBoxMin, _localBoxMax);
}

std::array<int, 6> RegularGridDecomposition::getExtentOfSubdomain(const int subdomainIndex) const {
  return DomainTools::getExtentOfSubdomain(subdomainIndex, _decomposition);
}

void RegularGridDecomposition::exchangeHaloParticles(SharedAutoPasContainer &autoPasContainer) {
  for (int i = 0; i < _dimensionCount; ++i) {
    std::vector<ParticleType> haloParticles{};
    for (int j = i; j < _dimensionCount; ++j) {
      const size_t dimensionIndex = j % _dimensionCount;

      // completely bypass Halo particle exchange in this dimension if boundaries in this direction are not periodic
      // *and* if both local boundaries are the global boundaries in this dimension
      if (_boundaryType[dimensionIndex] != options::BoundaryTypeOption::periodic and
          _localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex] and
          _localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex])
        continue;

      std::vector<ParticleType> particlesForLeftNeighbour{};
      std::vector<ParticleType> particlesForRightNeighbour{};

      if (j == i) {
        collectHaloParticlesForLeftNeighbour(autoPasContainer, dimensionIndex, particlesForLeftNeighbour);
        collectHaloParticlesForRightNeighbour(autoPasContainer, dimensionIndex, particlesForRightNeighbour);
      }

      double leftHaloMin = _localBoxMin[dimensionIndex] - _skinWidth;
      double leftHaloMax = _localBoxMin[dimensionIndex] + _cutoffWidth + _skinWidth;
      double rightHaloMin = _localBoxMax[dimensionIndex] - _cutoffWidth - _skinWidth;
      double rightHaloMax = _localBoxMax[dimensionIndex] + _skinWidth;

      for (const auto &particle : haloParticles) {
        auto position = particle.getR();

        // check left boundary is a periodic (global) boundary
        if (_boundaryType[dimensionIndex] == options::BoundaryTypeOption::periodic and
            _localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex]) {
          if (position[dimensionIndex] >= leftHaloMin and position[dimensionIndex] < leftHaloMax) {
            particlesForLeftNeighbour.push_back(particle);

            // Apply boundary condition
            if (_localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex]) {
              position[dimensionIndex] =
                  position[dimensionIndex] + (_globalBoxMax[dimensionIndex] - _globalBoxMin[dimensionIndex]);
              particlesForLeftNeighbour.back().setR(position);
            }
          }
        }

        // check right boundary is a periodic (global) boundary
        if (_boundaryType[dimensionIndex] == options::BoundaryTypeOption::periodic and
            _localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
          if (position[dimensionIndex] >= rightHaloMin and position[dimensionIndex] < rightHaloMax) {
            particlesForRightNeighbour.push_back(particle);

            // Apply boundary condition
            if (_localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
              position[dimensionIndex] =
                  position[dimensionIndex] - (_globalBoxMax[dimensionIndex] - _globalBoxMin[dimensionIndex]);
              particlesForRightNeighbour.back().setR(position);
            }
          }
        }
      }

      // See documentation for _neighbourDomainIndices to explain the indexing
      int leftNeighbour = _neighbourDomainIndices[(dimensionIndex * 2) % _neighbourCount];
      int rightNeighbour = _neighbourDomainIndices[(dimensionIndex * 2 + 1) % _neighbourCount];
      sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbour, particlesForRightNeighbour, leftNeighbour,
                                          rightNeighbour, haloParticles);
    }

    for (const auto &particle : haloParticles) {
      autoPasContainer->addHaloParticle(particle);
    }
  }
}

void RegularGridDecomposition::exchangeMigratingParticles(SharedAutoPasContainer &autoPasContainer) {
  auto emigrants = autoPasContainer->updateContainer();

  const std::array<double, _dimensionCount> globalBoxLength =
      autopas::utils::ArrayMath::sub(_globalBoxMax, _globalBoxMin);

  for (int i = 0; i < _dimensionCount; ++i) {
    for (int j = i; j < _dimensionCount; ++j) {
      const size_t dimensionIndex = j % _dimensionCount;
      // completely bypass particle exchange in this dimension if boundaries in this direction are not periodic *and*
      // if both local boundaries are the global boundaries in this dimension
      if (_boundaryType[dimensionIndex] != options::BoundaryTypeOption::periodic and
          _localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex] and
          _localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex])
        continue;

      std::vector<ParticleType> immigrants, remainingEmigrants;
      std::vector<ParticleType> particlesForLeftNeighbour;
      std::vector<ParticleType> particlesForRightNeighbour;

      // See documentation for _neighbourDomainIndices to explain the indexing
      int leftNeighbour = _neighbourDomainIndices[(dimensionIndex * 2) % _neighbourCount];
      int rightNeighbour = _neighbourDomainIndices[(dimensionIndex * 2 + 1) % _neighbourCount];

      for (const auto &particle : emigrants) {
        auto position = particle.getR();
        // check left boundary is a periodic (global) boundary
        if (_boundaryType[dimensionIndex] == options::BoundaryTypeOption::periodic and
            _localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex]) {
          if (position[dimensionIndex] < _localBoxMin[dimensionIndex]) {
            particlesForLeftNeighbour.push_back(particle);

            // Apply boundary condition
            if (_localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex]) {
              position[dimensionIndex] =
                  std::min(std::nextafter(_globalBoxMax[dimensionIndex], _globalBoxMin[dimensionIndex]),
                           position[dimensionIndex] + globalBoxLength[dimensionIndex]);
              particlesForLeftNeighbour.back().setR(position);
            }
          }
        }
        // check right boundary is a periodic (global) boundary
        if (_boundaryType[dimensionIndex] == options::BoundaryTypeOption::periodic and
            _localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
          if (position[dimensionIndex] >= _localBoxMax[dimensionIndex]) {
            particlesForRightNeighbour.push_back(particle);

            // Apply boundary condition
            if (_localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
              position[dimensionIndex] =
                  std::max(_globalBoxMin[dimensionIndex], position[dimensionIndex] - globalBoxLength[dimensionIndex]);
              particlesForRightNeighbour.back().setR(position);
            }
          } else {
            remainingEmigrants.push_back(particle);
          }
        }
      }
      emigrants = remainingEmigrants;

      sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbour, particlesForRightNeighbour, leftNeighbour,
                                          rightNeighbour, immigrants);

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
}

void RegularGridDecomposition::reflectParticlesAtBoundaries(SharedAutoPasContainer &autoPasContainer) {
  std::array<double, _dimensionCount> reflSkinMin{}, reflSkinMax{};

  for (int dimensionIndex = 0; dimensionIndex < _dimensionCount; ++dimensionIndex) {
    // skip if boundary is not reflective
    if (_boundaryType[dimensionIndex] != options::BoundaryTypeOption::reflective) continue;

    auto reflect = [&](bool isUpper) {
      for (auto p = autoPasContainer->getRegionIterator(reflSkinMin, reflSkinMax, autopas::IteratorBehavior::owned);
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
      reflSkinMax[dimensionIndex] = _globalBoxMin[dimensionIndex] + autoPasContainer->verletSkin() / 2;

      reflect(false);
    }
    // apply if we are at a global boundary on upper end of the dimension
    if (_localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
      reflSkinMin = _globalBoxMin;
      reflSkinMax = _globalBoxMax;
      reflSkinMin[dimensionIndex] = _globalBoxMax[dimensionIndex] - autoPasContainer->verletSkin() / 2;

      reflect(true);
    }
  }
}

void RegularGridDecomposition::sendParticles(const std::vector<ParticleType> &particles, const int &receiver) {
  std::vector<char> buffer;

  for (const auto &particle : particles) {
    ParticleSerializationTools::serializeParticle(particle, buffer);
  }

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

void RegularGridDecomposition::sendAndReceiveParticlesLeftAndRight(std::vector<ParticleType> &particlesToLeft,
                                                                   std::vector<ParticleType> &particlesToRight,
                                                                   const int &leftNeighbour, const int &rightNeighbour,
                                                                   std::vector<ParticleType> &receivedParticles) {
  if (_mpiCommunicationNeeded and leftNeighbour != _domainIndex) {
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

void RegularGridDecomposition::collectHaloParticlesForLeftNeighbour(SharedAutoPasContainer &autoPasContainer,
                                                                    const size_t &direction,
                                                                    std::vector<ParticleType> &haloParticles) {
  std::array<double, _dimensionCount> boxMin, boxMax;

  // Calculate halo box for left neighbour
  for (int i = 0; i < _dimensionCount; ++i) {
    boxMin[i] = _localBoxMin[i] - _skinWidth;
    boxMax[i] = _skinWidth + _localBoxMax[i];
  }
  boxMax[direction] = _localBoxMin[direction] + _cutoffWidth + _skinWidth;

  // Collect the halo particles for the left neighbour
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

void RegularGridDecomposition::collectHaloParticlesForRightNeighbour(SharedAutoPasContainer &autoPasContainer,
                                                                     const size_t &direction,
                                                                     std::vector<ParticleType> &haloParticles) {
  std::array<double, _dimensionCount> boxMin, boxMax;

  // Calculate left halo box of right neighbour
  for (int i = 0; i < _dimensionCount; ++i) {
    boxMin[i] = _localBoxMin[i] - _skinWidth;
    boxMax[i] = _localBoxMax[i] + _skinWidth;
  }
  boxMin[direction] = _localBoxMax[direction] - _cutoffWidth - _skinWidth;

  // Collect the halo particles for the right neighbour
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

void RegularGridDecomposition::categorizeParticlesIntoLeftAndRightNeighbour(
    const std::vector<ParticleType> &particles, const size_t &direction,
    std::vector<ParticleType> &leftNeighbourParticles, std::vector<ParticleType> &rightNeighbourParticles,
    std::vector<ParticleType> &uncategorizedParticles) {
  for (const auto &particle : particles) {
    std::array<double, _dimensionCount> position = particle.getR();
    if (position[direction] < _localBoxMin[direction]) {
      leftNeighbourParticles.push_back(particle);

      // Apply boundary condition
      if (_localBoxMin[direction] == _globalBoxMin[direction]) {
        position[direction] = std::min(std::nextafter(_globalBoxMax[direction], _globalBoxMin[direction]),
                                       position[direction] + _globalBoxMax[direction] - _globalBoxMin[direction]);
        leftNeighbourParticles.back().setR(position);
      }
    } else if (position[direction] >= _localBoxMax[direction]) {
      rightNeighbourParticles.push_back(particle);

      // Apply boundary condition
      if (_localBoxMax[direction] == _globalBoxMax[direction]) {
        position[direction] = std::max(_globalBoxMin[direction],
                                       position[direction] - (_globalBoxMax[direction] - _globalBoxMin[direction]));
        rightNeighbourParticles.back().setR(position);
      }
    } else {
      uncategorizedParticles.push_back(particle);
    }
  }
}
