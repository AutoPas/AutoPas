/**
 * @file RegularGrid.cpp
 * @author J. Körner
 * @date 19.04.2021
 */
#include "RegularGrid.h"

#include <algorithm>
#include <functional>
#include <math.h>
#include <numeric>


#include "autopas/utils/ArrayUtils.h"
#include "DomainTools.h"
#include "src/ParticleSerializationTools.h"

namespace {
/**
 * Calculates the prime factorization of a number.
 * @param number The number for which the factirization will be calculated.
 * @oPrimeFactors The container, where the prime factos will be appended.
 */
void calculatePrimeFactors(unsigned int number, std::list<unsigned int> &oPrimeFactors) {
  while (number % 2 == 0) {
    oPrimeFactors.push_back(2);
    number = number / 2;
  }

  for (unsigned int i = 3; i <= number; i = i + 2) {
    while (number % i == 0) {
      oPrimeFactors.push_back(i);
      number = number / i;
    }
  }
}
}

RegularGrid::RegularGrid(const int &dimensionCount, const std::vector<double> &globalBoxMin,
                         const std::vector<double> &globalBoxMax) {
  _dimensionCount = dimensionCount;

  MPI_Comm_size(MPI_COMM_WORLD, &_subdomainCount);

  initializeDecomposition();

  initializeMPICommunicator();

  initializeLocalDomain();

  initializeGlobalBox(globalBoxMin, globalBoxMax);

  initializeLocalBox();

  initializeNeighbourIds();

  std::cout << "Rank: " << _domainIndex << ", " << autopas::utils::ArrayUtils::to_string(_neighbourDomainIndices)
            << std::endl;
}

RegularGrid::~RegularGrid() { }

void RegularGrid::update() { updateLocalBox(); }

void RegularGrid::setHaloWidth(double width){
  _haloWidth = width; 
}

void RegularGrid::initializeDecomposition() {
  std::list<unsigned int> primeFactors;
  calculatePrimeFactors(_subdomainCount, primeFactors);

  while (primeFactors.size() > _dimensionCount) {
    primeFactors.sort();
    auto firstElement = primeFactors.front();
    primeFactors.pop_front();
    primeFactors.front() *= firstElement;
  }

  _decomposition.resize(_dimensionCount);

  for (auto &dimensionSize : _decomposition) {
    if (primeFactors.size() > 0) {
      dimensionSize = primeFactors.front();
      primeFactors.pop_front();
    } else {
      dimensionSize = 1;
    }
  }
}

void RegularGrid::initializeMPICommunicator() {
  std::vector<int> periods(_dimensionCount, 1);
  MPI_Cart_create(MPI_COMM_WORLD, _dimensionCount, _decomposition.data(), periods.data(), true, &_communicator);
  MPI_Comm_rank(_communicator, &_domainIndex);
}

void RegularGrid::initializeLocalDomain() {
  _domainId.resize(_dimensionCount);
  MPI_Comm_rank(_communicator, &_domainIndex);

  std::vector<int> periods(_dimensionCount, 1);
  MPI_Cart_get(_communicator, _dimensionCount, _decomposition.data(), periods.data(), _domainId.data());
}

void RegularGrid::initializeLocalBox() {
  _localBoxMin.resize(_dimensionCount);
  _localBoxMax.resize(_dimensionCount);
  updateLocalBox();

  for (int i = 0; i < _localBoxMin.size(); ++i){
    _haloWidth += (_localBoxMax[i] - _localBoxMin[i]) / 20.0;
  }

  _haloWidth = _haloWidth / _localBoxMin.size();
}

void RegularGrid::initializeNeighbourIds() {
  _neighbourDomainIndices.resize(_dimensionCount * 2);

  for (int i = 0; i < _dimensionCount; ++i) {
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

int RegularGrid::convertIdToIndex(const std::vector<int> &domainId) {
  int neighbourDomainIndex = 0;

  for (int i = 0; i < _dimensionCount; ++i) {
    int accumulatedTail = 1;

    if (i < _decomposition.size() - 1) {
      accumulatedTail =
          std::accumulate(_decomposition.begin() + i + 1, _decomposition.end(), 1, std::multiplies<int>());
    }

    neighbourDomainIndex += accumulatedTail * domainId[i];
  }

  return neighbourDomainIndex;
}

void RegularGrid::updateLocalBox() {
  for (int i = 0; i < _dimensionCount; ++i) {
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

void RegularGrid::initializeGlobalBox(const std::vector<double> &globalBoxMin,
                                      const std::vector<double> &globalBoxMax) {
  _globalBoxMin.resize(_dimensionCount);
  _globalBoxMax.resize(_dimensionCount);
  for (int i = 0; i < _dimensionCount; ++i) {
    _globalBoxMin[i] = globalBoxMin[i];
    _globalBoxMax[i] = globalBoxMax[i];
  }
}

void RegularGrid::exchangeHaloParticles(SharedAutoPasContainer &autoPasContainer) {
  int dimensionCount = _localBoxMin.size();
  int neighbourCount = dimensionCount * 2;

  std::vector<ParticleType> particlesForLeftNeighbour;
  std::vector<ParticleType> particlesForRightNeighbour;
  std::vector<ParticleType> haloParticles;

  double particlePosition;

  // @todo: these are named in a misinforming manner. They do not refer to the current domain's halos,
  // but to the halo boundaries of the respective neighbour
  double leftNeighbourHaloEnd, rightNeighbourHaloStart;
  double bottomNeighbourHaloEnd, topNeighbourHaloStart;

  int leftNeighbour, rightNeighbour;

  int nextDimensionIndex;

  for (int i = 0; i < dimensionCount; ++i) {
    leftNeighbour = _neighbourDomainIndices[(i * 2) % neighbourCount];
    rightNeighbour = _neighbourDomainIndices[(i * 2 + 1) % neighbourCount];

    leftNeighbourHaloEnd = _localBoxMin[i] + _haloWidth;
    rightNeighbourHaloStart = _localBoxMax[i] - _haloWidth;

    nextDimensionIndex = (i + 1) % dimensionCount;

    if (leftNeighbour != _domainIndex && rightNeighbour != _domainIndex) {
      for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
        particlePosition = particle->getR()[i];
        if (particlePosition < leftNeighbourHaloEnd) {
          particlesForLeftNeighbour.push_back(*particle);
        } else if (particlePosition > rightNeighbourHaloStart) {
          particlesForRightNeighbour.push_back(*particle);
        }
      }

      sendParticles(particlesForLeftNeighbour, leftNeighbour);
      sendParticles(particlesForRightNeighbour, rightNeighbour);

      receiveParticles(haloParticles, leftNeighbour);
      receiveParticles(haloParticles, rightNeighbour);

      waitForSendRequests();
    }

    particlesForLeftNeighbour.clear();
    particlesForRightNeighbour.clear();

    leftNeighbourHaloEnd = _localBoxMin[nextDimensionIndex] + _haloWidth;
    rightNeighbourHaloStart = _localBoxMax[nextDimensionIndex] - _haloWidth;

    leftNeighbour = _neighbourDomainIndices[(leftNeighbour + 2) % neighbourCount];
    rightNeighbour = _neighbourDomainIndices[(rightNeighbour + 2) % neighbourCount];

    if (leftNeighbour != _domainIndex && rightNeighbour != _domainIndex) {
      for (auto &particle : haloParticles) {
        particlePosition = particle.getR()[nextDimensionIndex];
        if (particlePosition < leftNeighbourHaloEnd) {
          particlesForLeftNeighbour.push_back(particle);
        } else if (particlePosition > rightNeighbourHaloStart) {
          particlesForRightNeighbour.push_back(particle);
        }
      }

      sendParticles(particlesForLeftNeighbour, leftNeighbour);
      sendParticles(particlesForRightNeighbour, rightNeighbour);

      receiveParticles(haloParticles, leftNeighbour);
      receiveParticles(haloParticles, rightNeighbour);

      waitForSendRequests();
    }

    particlesForLeftNeighbour.clear();
    particlesForRightNeighbour.clear();
  }

  std::vector<double> haloParticlePosition = {0.0, 0.0, 0.0};
  for (auto &particle : haloParticles) {
    haloParticlePosition[0] = particle.getR()[0];
    haloParticlePosition[1] = particle.getR()[1];
    haloParticlePosition[2] = particle.getR()[2];

    if (!isInsideLocalDomain(haloParticlePosition)) {
      autoPasContainer->addOrUpdateHaloParticle(particle);
    }
  }
}

void RegularGrid::exchangeMigratingParticles(SharedAutoPasContainer &autoPasContainer) {
  auto [emigrants, updated] = autoPasContainer->updateContainer();
  int dimensionCount = _localBoxMin.size();

  std::vector<ParticleType> immigrants;
  std::vector<ParticleType> migrants;

  std::vector<ParticleType> particlesForLeftNeighbour;
  std::vector<ParticleType> particlesForRightNeighbour;

  int nextDimensionIndex;
  int leftNeighbour, rightNeighbour;
  int neighbourCount = dimensionCount * 2;
  double particlePosition;

  for (int i = 0; i < dimensionCount; ++i) {
    leftNeighbour = _neighbourDomainIndices[(i * 2) % neighbourCount];
    rightNeighbour = _neighbourDomainIndices[(i * 2 + 1) % neighbourCount];

    if (leftNeighbour != _domainIndex && rightNeighbour != _domainIndex) {
      for (auto &particle : emigrants) {
        particlePosition = particle.getR()[i];
        if (particlePosition - _localBoxMin[i] < 0) {
          particlesForLeftNeighbour.push_back(particle);
        } else if (particlePosition - _localBoxMax[i] > 0) {
          particlesForRightNeighbour.push_back(particle);
        }
      }

      sendParticles(particlesForLeftNeighbour, leftNeighbour);
      sendParticles(particlesForRightNeighbour, rightNeighbour);

      receiveParticles(migrants, leftNeighbour);
      receiveParticles(migrants, rightNeighbour);

      waitForSendRequests();
    }

    particlesForLeftNeighbour.clear();
    particlesForRightNeighbour.clear();

    leftNeighbour = _neighbourDomainIndices[(leftNeighbour + 2) % neighbourCount];
    rightNeighbour = _neighbourDomainIndices[(rightNeighbour + 2) % neighbourCount];

    if (leftNeighbour != _domainIndex && rightNeighbour != _domainIndex) {
      nextDimensionIndex = (i + 1) % dimensionCount;

      for (auto &particle : migrants) {
        particlePosition = particle.getR()[nextDimensionIndex];
        if (particlePosition - _localBoxMin[nextDimensionIndex] < 0) {
          particlesForLeftNeighbour.push_back(particle);
        } else if (particlePosition - _localBoxMax[nextDimensionIndex] > 0) {
          particlesForRightNeighbour.push_back(particle);
        } else {
          immigrants.push_back(particle);
        }
      }

      migrants.clear();

      sendParticles(particlesForLeftNeighbour, leftNeighbour);
      sendParticles(particlesForRightNeighbour, rightNeighbour);

      receiveParticles(immigrants, leftNeighbour);
      receiveParticles(immigrants, rightNeighbour);

      waitForSendRequests();
    }

    particlesForLeftNeighbour.clear();
    particlesForRightNeighbour.clear();
  }

  for (auto &particle : immigrants) {
    autoPasContainer->addParticle(particle);
  }
}

void RegularGrid::sendParticles(std::vector<ParticleType> &particles, int &receiver) {
  std::vector<char> buffer;

  for (auto &particle : particles) {
    ParticleSerializationTools::serializeParticle(particle, buffer);
  }

  size_t sizeOfParticleAttributes = sizeof(ParticleAttributes);

  sendDataToNeighbour(buffer, receiver);
}

void RegularGrid::receiveParticles(std::vector<ParticleType> &receivedParticles, int &source) {
  std::vector<char> receiveBuffer;

  receiveDataFromNeighbour(source, receiveBuffer);

  if (!receiveBuffer.empty()) {
    ParticleSerializationTools::deserializeParticleData(receiveBuffer, receivedParticles);
  }
}

void RegularGrid::sendDataToNeighbour(std::vector<char> sendBuffer, const int &neighbour) {
  _sendBuffers.push_back(sendBuffer);

  MPI_Request sendRequest;
  _sendRequests.push_back(sendRequest);

  MPI_Isend(_sendBuffers.back().data(), _sendBuffers.back().size(), MPI_CHAR, neighbour, 0, _communicator,
            &_sendRequests.back());
}

void RegularGrid::receiveDataFromNeighbour(const int &neighbour, std::vector<char> &receiveBuffer) {
  MPI_Status status;
  MPI_Probe(neighbour, 0, _communicator, &status);

  int receiveBufferSize;
  MPI_Get_count(&status, MPI_CHAR, &receiveBufferSize);
  receiveBuffer.resize(receiveBufferSize);

  MPI_Recv(receiveBuffer.data(), receiveBufferSize, MPI_CHAR, neighbour, 0, _communicator, MPI_STATUS_IGNORE);
}

void RegularGrid::waitForSendRequests() {
  std::vector<MPI_Status> sendStates;
  sendStates.resize(_sendRequests.size());
  MPI_Waitall(_sendRequests.size(), _sendRequests.data(), sendStates.data());
  _sendRequests.clear();
  _sendBuffers.clear();
}

bool RegularGrid::isInsideLocalDomain(const std::vector<double> &coordinates) {
  return DomainTools::isInsideDomain(coordinates, _localBoxMin, _localBoxMax);
}
