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

RegularGridDecomposition::RegularGridDecomposition(const std::array<double, 3> &globalBoxMin,
                                                   const std::array<double, 3> &globalBoxMax,
                                                   const std::array<bool, 3> &subdivideDimension,
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
      // Create planar communicator along shift axis. The shift axis is vertical to the dimension's direction.
      autopas::AutoPas_MPI_Comm planarCommunicator;
      autopas::AutoPas_MPI_Comm_split(_communicator, _domainId[i], _domainId[i], &planarCommunicator);

      // Calculate average distributed work in planar communicator.
      double distributedWorkInPlane;
      autopas::AutoPas_MPI_Allreduce(&distributedWork, &distributedWorkInPlane, 1, AUTOPAS_MPI_DOUBLE, AUTOPAS_MPI_SUM,
                                     planarCommunicator);
      distributedWorkInPlane = distributedWorkInPlane / _decomposition[i + 1];

      const int leftNeighbour = _neighbourDomainIndices[(i * 2) % _neighbourCount];
      const int rightNeighbour = _neighbourDomainIndices[(i * 2 + 1) % _neighbourCount];
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

      if (_localBoxMin[i] != _globalBoxMin[i]) {
        double neighbourPlaneWork, neighbourBoundary;
        autopas::AutoPas_MPI_Recv(&neighbourPlaneWork, 1, AUTOPAS_MPI_DOUBLE, leftNeighbour, 0, _communicator,
                                  AUTOPAS_MPI_STATUS_IGNORE);
        autopas::AutoPas_MPI_Recv(&neighbourBoundary, 1, AUTOPAS_MPI_DOUBLE, leftNeighbour, 0, _communicator,
                                  AUTOPAS_MPI_STATUS_IGNORE);

        _localBoxMin[i] = DomainTools::balanceAdjacentDomains(neighbourPlaneWork, distributedWorkInPlane,
                                                              neighbourBoundary, _localBoxMax[i]);
      }
      if (_localBoxMax[i] != _globalBoxMax[i]) {
        double neighbourPlaneWork, neighbourBoundary;
        autopas::AutoPas_MPI_Recv(&neighbourPlaneWork, 1, AUTOPAS_MPI_DOUBLE, rightNeighbour, 0, _communicator,
                                  AUTOPAS_MPI_STATUS_IGNORE);
        autopas::AutoPas_MPI_Recv(&neighbourBoundary, 1, AUTOPAS_MPI_DOUBLE, rightNeighbour, 0, _communicator,
                                  AUTOPAS_MPI_STATUS_IGNORE);

        _localBoxMax[i] = DomainTools::balanceAdjacentDomains(distributedWorkInPlane, neighbourPlaneWork,
                                                              _localBoxMin[i], neighbourBoundary);
      }
    }
  }
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
    _neighbourDomainIndices[neighbourIndex] = convertIdToIndex(preceedingNeighbourId);

    ++neighbourIndex;
    auto succeedingNeighbourId = _domainId;
    succeedingNeighbourId[i] = (++succeedingNeighbourId[i] + _decomposition[i]) % _decomposition[i];
    _neighbourDomainIndices[neighbourIndex] = convertIdToIndex(succeedingNeighbourId);
  }
}

void RegularGridDecomposition::initializeGlobalBox(const std::array<double, 3> &globalBoxMin,
                                                   const std::array<double, 3> &globalBoxMax) {
  for (int i = 0; i < 3; ++i) {
    _globalBoxMin[i] = globalBoxMin[i];
    _globalBoxMax[i] = globalBoxMax[i];
  }
}

bool RegularGridDecomposition::isInsideLocalDomain(const std::array<double, 3> &coordinates) {
  return DomainTools::isInsideDomain(coordinates, _localBoxMin, _localBoxMax);
}

void RegularGridDecomposition::exchangeHaloParticles(SharedAutoPasContainer &autoPasContainer) {
  for (int i = 0; i < _dimensionCount; ++i) {
    std::vector<ParticleType> haloParticles{};
    for (int j = i; j < _dimensionCount; ++j) {
      const size_t dimensionIndex = j % _dimensionCount;

      std::vector<ParticleType> particlesForLeftNeighbour, particlesForRightNeighbour;
      std::array<double, _dimensionCount> haloBoxMin, haloBoxMax;
      if (j == i) {
        for (int k = 0; k < _dimensionCount; ++k) {
          haloBoxMin[k] = _localBoxMin[k] - _skinWidth;
          haloBoxMax[k] = _skinWidth + _localBoxMax[k];
        }
        haloBoxMax[dimensionIndex] = _localBoxMin[dimensionIndex] + _cutoffWidth + _skinWidth;

        for (auto particleIterator =
                 autoPasContainer->getRegionIterator(haloBoxMin, haloBoxMax, autopas::IteratorBehavior::owned);
             particleIterator.isValid(); ++particleIterator) {
          std::array<double, _dimensionCount> position = particleIterator->getR();
          particlesForLeftNeighbour.push_back(*particleIterator);

          // Apply boundary condition
          if (_localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex]) {
            position[dimensionIndex] =
                position[dimensionIndex] + (_globalBoxMax[dimensionIndex] - _globalBoxMin[dimensionIndex]);
            particlesForLeftNeighbour.back().setR(position);
          }
        }

        for (int k = 0; k < _dimensionCount; ++k) {
          haloBoxMin[k] = _localBoxMax[k] - _skinWidth - (_localBoxMax[k] - _localBoxMin[k]);
          haloBoxMax[k] = _localBoxMax[k] + _skinWidth;
        }
        haloBoxMin[dimensionIndex] = _localBoxMax[dimensionIndex] - _cutoffWidth - _skinWidth;

        for (auto particleIterator =
                 autoPasContainer->getRegionIterator(haloBoxMin, haloBoxMax, autopas::IteratorBehavior::owned);
             particleIterator.isValid(); ++particleIterator) {
          std::array<double, _dimensionCount> position = particleIterator->getR();
          particlesForRightNeighbour.push_back(*particleIterator);

          // Apply boundary condition
          if (_localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
            position[dimensionIndex] =
                position[dimensionIndex] - (_globalBoxMax[dimensionIndex] - _globalBoxMin[dimensionIndex]);
            particlesForRightNeighbour.back().setR(position);
          }
        }
      }

      double leftHaloMin = _localBoxMin[dimensionIndex] - _skinWidth;
      double leftHaloMax = _localBoxMin[dimensionIndex] + _cutoffWidth + _skinWidth;
      double rightHaloMin = _localBoxMax[dimensionIndex] - _cutoffWidth - _skinWidth;
      double rightHaloMax = _localBoxMax[dimensionIndex] + _skinWidth;

      for (const auto &particle : haloParticles) {
        std::array<double, _dimensionCount> position = particle.getR();
        if (position[dimensionIndex] >= leftHaloMin and position[dimensionIndex] < leftHaloMax) {
          particlesForLeftNeighbour.push_back(particle);

          // Apply boundary condition
          if (_localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex]) {
            position[dimensionIndex] =
                position[dimensionIndex] + (_globalBoxMax[dimensionIndex] - _globalBoxMin[dimensionIndex]);
            particlesForLeftNeighbour.back().setR(position);
          }
        } else if (position[dimensionIndex] >= rightHaloMin and position[dimensionIndex] < rightHaloMax) {
          particlesForRightNeighbour.push_back(particle);

          // Apply boundary condition
          if (_localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
            position[dimensionIndex] =
                position[dimensionIndex] - (_globalBoxMax[dimensionIndex] - _globalBoxMin[dimensionIndex]);
            particlesForRightNeighbour.back().setR(position);
          }
        }
      }

      // See documentation for _neighbourDomainIndices to explain the indexing
      int leftNeighbour = _neighbourDomainIndices[(dimensionIndex * 2) % _neighbourCount];
      int rightNeighbour = _neighbourDomainIndices[(dimensionIndex * 2 + 1) % _neighbourCount];
      sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbour, particlesForRightNeighbour, leftNeighbour,
                                          rightNeighbour, haloParticles);
    }

    for (auto &particle : haloParticles) {
      autoPasContainer->addOrUpdateHaloParticle(particle);
    }
  }
}

void RegularGridDecomposition::exchangeMigratingParticles(SharedAutoPasContainer &autoPasContainer,
                                                          std::vector<ParticleType> &emigrants, const bool &updated) {
  if (updated) {
    const std::array<double, _dimensionCount> globalBoxMin = {_globalBoxMin[0], _globalBoxMin[1], _globalBoxMin[2]};
    const std::array<double, _dimensionCount> globalBoxMax = {_globalBoxMax[0], _globalBoxMax[1], _globalBoxMax[2]};
    const std::array<double, _dimensionCount> globalBoxLength =
        autopas::utils::ArrayMath::sub(globalBoxMax, globalBoxMin);

    for (int i = 0; i < _dimensionCount; ++i) {
      for (int j = i; j < _dimensionCount; ++j) {
        const size_t dimensionIndex = j % _dimensionCount;

        std::vector<ParticleType> immigrants, remainingEmigrants;
        std::vector<ParticleType> particlesForLeftNeighbour;
        std::vector<ParticleType> particlesForRightNeighbour;

        // See documentation for _neighbourDomainIndices to explain the indexing
        int leftNeighbour = _neighbourDomainIndices[(dimensionIndex * 2) % _neighbourCount];
        int rightNeighbour = _neighbourDomainIndices[(dimensionIndex * 2 + 1) % _neighbourCount];

        std::array<double, _dimensionCount> position;
        for (const auto &particle : emigrants) {
          position = particle.getR();
          if (position[dimensionIndex] < _localBoxMin[dimensionIndex]) {
            particlesForLeftNeighbour.push_back(particle);

            // Apply boundary condition
            if (_localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex]) {
              position[dimensionIndex] =
                  std::min(std::nextafter(_globalBoxMax[dimensionIndex], _globalBoxMin[dimensionIndex]),
                           position[dimensionIndex] + globalBoxLength[dimensionIndex]);
              particlesForLeftNeighbour.back().setR(position);
            }
          } else if (position[dimensionIndex] >= _localBoxMax[dimensionIndex]) {
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
}

void RegularGridDecomposition::sendAndReceiveParticlesLeftAndRight(const std::vector<ParticleType> &particlesToLeft,
                                                                   const std::vector<ParticleType> &particlesToRight,
                                                                   const int &leftNeighbour, const int &rightNeighbour,
                                                                   std::vector<ParticleType> &receivedParticles) {
  if (_mpiIsEnabled && leftNeighbour != _domainIndex) {
    ParticleCommunicator particleCommunicator(_communicator);
    particleCommunicator.sendParticles(particlesToLeft, leftNeighbour);
    particleCommunicator.sendParticles(particlesToRight, rightNeighbour);

    particleCommunicator.receiveParticles(receivedParticles, leftNeighbour);
    particleCommunicator.receiveParticles(receivedParticles, rightNeighbour);

    particleCommunicator.waitForSendRequests();
  } else {
    receivedParticles.insert(receivedParticles.end(), particlesToLeft.begin(), particlesToLeft.end());
    receivedParticles.insert(receivedParticles.end(), particlesToRight.begin(), particlesToRight.end());
  }
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
