/**
 * @file RegularGrid.cpp
 * @author J. Körner
 * @date 19.04.2021
 */
#include "RegularGrid.h"

#include <math.h>
#include <algorithm>
#include <functional>
#include <numeric>

#include "DomainTools.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayUtils.h"
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
}  // namespace

RegularGrid::RegularGrid(const int &dimensionCount, const std::vector<double> &globalBoxMin,
              const std::vector<double> &globalBoxMax, const double &cutoffWidth, const double &skinWidth)
              : _dimensionCount(dimensionCount), _cutoffWidth(cutoffWidth), _skinWidth(skinWidth) {

  autopas::AutoPas_MPI_Comm_size(autopas::AUTOPAS_MPI_COMM_WORLD, &_subdomainCount);

  DomainTools::generateDecomposition(_subdomainCount, _dimensionCount, _decomposition);

  initializeMPICommunicator();

  initializeLocalDomain();

  initializeGlobalBox(globalBoxMin, globalBoxMax);

  initializeLocalBox();

  initializeNeighbourIds();

	std::cout << "Neighbours: " << autopas::utils::ArrayUtils::to_string(_neighbourDomainIndices) << std::endl;
	std::cout << "LocalBox: " << autopas::utils::ArrayUtils::to_string(_localBoxMin) << ", " << autopas::utils::ArrayUtils::to_string(_localBoxMax) << std::endl;
}

RegularGrid::~RegularGrid() {}

void RegularGrid::update() { updateLocalBox(); }

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
  autopas::AutoPas_MPI_Cart_create(autopas::AUTOPAS_MPI_COMM_WORLD, _dimensionCount, _decomposition.data(),
                                   periods.data(), true, &_communicator);
  autopas::AutoPas_MPI_Comm_rank(_communicator, &_domainIndex);
}

void RegularGrid::initializeLocalDomain() {
  _domainId.resize(_dimensionCount);
  autopas::AutoPas_MPI_Comm_rank(_communicator, &_domainIndex);

  std::vector<int> periods(_dimensionCount, 1);
  autopas::AutoPas_MPI_Cart_get(_communicator, _dimensionCount, _decomposition.data(), periods.data(),
                                _domainId.data());

  if (_domainId.empty()) {
    _domainId.resize(_dimensionCount, 0);
  }
}

void RegularGrid::initializeLocalBox() {
  _localBoxMin.resize(_dimensionCount);
  _localBoxMax.resize(_dimensionCount);
  updateLocalBox();
  _haloBoxes.resize(4 * _dimensionCount);
  updateHaloBoxes();
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

bool RegularGrid::isInsideLocalDomain(const std::vector<double> &coordinates) {
  return DomainTools::isInsideDomain(coordinates, _localBoxMin, _localBoxMax);
}

bool RegularGrid::isInsideLocalDomain(const std::array<double, 3> &coordinates) {
  std::vector<double> coordinatesVector;
  coordinatesVector.insert(coordinatesVector.begin(), coordinates.begin(), coordinates.end());
  return isInsideLocalDomain(coordinatesVector);
}

void RegularGrid::exchangeHaloParticles(SharedAutoPasContainer &autoPasContainer) {
  int dimensionCount = _localBoxMin.size();
  int neighbourCount = dimensionCount * 2;

  double particlePosition;

  std::array<double, 3> localBoxMin, localBoxMax, leftHaloBoxMin, leftHaloBoxMax, rightHaloBoxMin, rightHaloBoxMax;

  localBoxMin = {_localBoxMin[0], _localBoxMin[1], _localBoxMin[2]};
  localBoxMax = {_localBoxMax[0], _localBoxMax[1], _localBoxMax[2]};

  int leftNeighbour, rightNeighbour;

  int nextDimensionIndex;

  int haloWidth = _cutoffWidth + _skinWidth;

  for (int i = 0; i < dimensionCount && _decomposition[i] > 1; ++i) {
    leftNeighbour = _neighbourDomainIndices[(i * 2) % neighbourCount];
    rightNeighbour = _neighbourDomainIndices[(i * 2 + 1) % neighbourCount];

    std::vector<ParticleType> haloParticles;

    if (leftNeighbour != _domainIndex && rightNeighbour != _domainIndex) {

      std::vector<ParticleType> particlesForLeftNeighbour;
      std::vector<ParticleType> particlesForRightNeighbour;

      std::array<double, 3> periodicShift;
      for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
        double particlePosition = particle->getR()[i];
        if (particlePosition >= _haloBoxes[4 * i] && particlePosition < _haloBoxes[4 * i + 1] ) {
          particlesForLeftNeighbour.push_back(*particle);

          if (_localBoxMin[i] == _globalBoxMin[i]) {
            periodicShift = {0.0, 0.0, 0.0};
            periodicShift[i] = _globalBoxMax[i] - _globalBoxMin[i];
            particlesForLeftNeighbour.back().addR(periodicShift);
          }
        }
      }

      for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
        double particlePosition = particle->getR()[i];
        if (particlePosition > _haloBoxes[4 * i + 2] && particlePosition < _haloBoxes[4 * i + 3]) {
          particlesForRightNeighbour.push_back(*particle);

          if (_localBoxMax[i] == _globalBoxMax[i]) {
            periodicShift = {0.0, 0.0, 0.0};
            periodicShift[i] = - (_globalBoxMax[i] - _globalBoxMin[i]);
            particlesForRightNeighbour.back().addR(periodicShift);
          }
        }
      }

      sendParticles(particlesForLeftNeighbour, leftNeighbour);
      sendParticles(particlesForRightNeighbour, rightNeighbour);

      receiveParticles(haloParticles, leftNeighbour);
      receiveParticles(haloParticles, rightNeighbour);

      waitForSendRequests();
    }

    leftNeighbour = _neighbourDomainIndices[(leftNeighbour + 2) % neighbourCount];
    rightNeighbour = _neighbourDomainIndices[(rightNeighbour + 2) % neighbourCount];

    if (leftNeighbour != _domainIndex && rightNeighbour != _domainIndex) {
      std::vector<ParticleType> particlesForLeftNeighbour;
      std::vector<ParticleType> particlesForRightNeighbour;

      nextDimensionIndex = (i + 1) % dimensionCount;

      std::array<double, 3> periodicShift;
      for (const auto &particle : haloParticles) {
        double particlePosition = particle.getR()[nextDimensionIndex];
        if (particlePosition >= _haloBoxes[4 * nextDimensionIndex] && particlePosition < _haloBoxes[4 * i + 1]) {
          particlesForLeftNeighbour.push_back(particle);

          if (_localBoxMin[nextDimensionIndex] == _globalBoxMin[nextDimensionIndex]) {
            periodicShift = {0.0, 0.0, 0.0};
            periodicShift[nextDimensionIndex] = _globalBoxMax[nextDimensionIndex] - _globalBoxMin[nextDimensionIndex];
            particlesForLeftNeighbour.back().addR(periodicShift);
          }
        }

        if (particlePosition > _haloBoxes[4 * nextDimensionIndex + 2] && particlePosition < _haloBoxes[4 * nextDimensionIndex + 3]) {
          particlesForRightNeighbour.push_back(particle);

          if (_localBoxMax[nextDimensionIndex] == _globalBoxMax[nextDimensionIndex]) {
            periodicShift = {0.0, 0.0, 0.0};
            periodicShift[nextDimensionIndex] = - (_globalBoxMax[nextDimensionIndex] - _globalBoxMin[nextDimensionIndex]);
            particlesForRightNeighbour.back().addR(periodicShift);
          }
        }
      }

      sendParticles(particlesForLeftNeighbour, leftNeighbour);
      sendParticles(particlesForRightNeighbour, rightNeighbour);

      receiveParticles(haloParticles, leftNeighbour);
      receiveParticles(haloParticles, rightNeighbour);

      waitForSendRequests();
    }

    for (auto &particle : haloParticles) {
      autoPasContainer->addOrUpdateHaloParticle(particle);
    }
  }
}

void RegularGrid::exchangeMigratingParticles(SharedAutoPasContainer &autoPasContainer) {
  auto [emigrants, updated] = autoPasContainer->updateContainer();

  int dimensionCount = _localBoxMin.size();
  int neighbourCount = _neighbourDomainIndices.size();

  for (int i = 0; i < dimensionCount; ++i) {
    std::vector<ParticleType> immigrants, migrants;

    int leftNeighbour = _neighbourDomainIndices[(i * 2) % neighbourCount];
    int rightNeighbour = _neighbourDomainIndices[(i * 2 + 1) % neighbourCount];

		std::vector<ParticleType> particlesForLeftNeighbour;
		std::vector<ParticleType> particlesForRightNeighbour;

		std::array<double, 3> periodicShift;

		for (const auto &particle : emigrants) {
			double particlePosition = particle.getR()[i];
			if (particlePosition < _localBoxMin[i]) {
				particlesForLeftNeighbour.push_back(particle);

				if (_localBoxMin[i] == _globalBoxMin[i]) {
					periodicShift = {0.0, 0.0, 0.0};
					periodicShift[i] = _globalBoxMax[i] - _globalBoxMin[i];
					particlesForLeftNeighbour.back().addR(periodicShift);
				}
			}
			else if (particlePosition > _localBoxMax[i]) {
				particlesForRightNeighbour.push_back(particle);

				if (_localBoxMax[i] == _globalBoxMax[i]) {
					periodicShift = {0.0, 0.0, 0.0};
					periodicShift[i] = - (_globalBoxMax[i] - _globalBoxMin[i]);
					particlesForRightNeighbour.back().addR(periodicShift);
				}
			}
		}

		sendParticles(emigrants, leftNeighbour);
		sendParticles(emigrants, rightNeighbour);

		receiveParticles(immigrants, leftNeighbour);
		receiveParticles(immigrants, rightNeighbour);

		for (auto &particle : immigrants) {
			if(isInsideLocalDomain(particle.getR())) {
				if (particle.getID() == 4484) {
					std::cout << "IsInsideLocalDomain" << std::endl;
				}
				autoPasContainer->addParticle(particle);
			}
			else {
				migrants.push_back(particle);
			}
		}

		immigrants.clear();

		waitForSendRequests();
		
    leftNeighbour = _neighbourDomainIndices[(leftNeighbour + 2) % neighbourCount];
    rightNeighbour = _neighbourDomainIndices[(rightNeighbour + 2) % neighbourCount];

    particlesForLeftNeighbour.clear();
    particlesForRightNeighbour.clear();

		for (const auto &particle : emigrants) {
			double particlePosition = particle.getR()[i];
			if (particlePosition < _localBoxMin[i]) {
				particlesForLeftNeighbour.push_back(particle);

				if (_localBoxMin[i] == _globalBoxMin[i]) {
					periodicShift = {0.0, 0.0, 0.0};
					periodicShift[i] = _globalBoxMax[i] - _globalBoxMin[i];
					particlesForLeftNeighbour.back().addR(periodicShift);
				}
			}
			else if (particlePosition > _localBoxMax[i]) {
				particlesForRightNeighbour.push_back(particle);

				if (_localBoxMax[i] == _globalBoxMax[i]) {
					periodicShift = {0.0, 0.0, 0.0};
					periodicShift[i] = - (_globalBoxMax[i] - _globalBoxMin[i]);
					particlesForRightNeighbour.back().addR(periodicShift);
				}
			}
		}

    int nextDimensionIndex = (i + 1) % dimensionCount;

    sendParticles(migrants, leftNeighbour);
    sendParticles(migrants, rightNeighbour);

    receiveParticles(immigrants, leftNeighbour);
    receiveParticles(immigrants, rightNeighbour);

    for (auto &particle : immigrants) {
      if(isInsideLocalDomain(particle.getR())) {
        autoPasContainer->addParticle(particle);
      }
    }

    waitForSendRequests();
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

  autopas::AutoPas_MPI_Request sendRequest;
  _sendRequests.push_back(sendRequest);

  autopas::AutoPas_MPI_Isend(_sendBuffers.back().data(), _sendBuffers.back().size(), autopas::AUTOPAS_MPI_CHAR,
                             neighbour, 0, _communicator, &_sendRequests.back());
}

void RegularGrid::receiveDataFromNeighbour(const int &neighbour, std::vector<char> &receiveBuffer) {
  autopas::AutoPas_MPI_Status status;
  autopas::AutoPas_MPI_Probe(neighbour, 0, _communicator, &status);

  int receiveBufferSize;
  autopas::AutoPas_MPI_Get_count(&status, autopas::AUTOPAS_MPI_CHAR, &receiveBufferSize);
  receiveBuffer.resize(receiveBufferSize);

  autopas::AutoPas_MPI_Recv(receiveBuffer.data(), receiveBufferSize, autopas::AUTOPAS_MPI_CHAR, neighbour, 0,
                            _communicator, AUTOPAS_MPI_STATUS_IGNORE);
}

void RegularGrid::waitForSendRequests() {
  std::vector<autopas::AutoPas_MPI_Status> sendStates;
  sendStates.resize(_sendRequests.size());
  autopas::AutoPas_MPI_Waitall(_sendRequests.size(), _sendRequests.data(), sendStates.data());
  _sendRequests.clear();
  _sendBuffers.clear();
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

void RegularGrid::updateHaloBoxes() {
  for (int dimensionIndex = 0; dimensionIndex < _dimensionCount; ++dimensionIndex) {
    int haloBoxesStartIndex = 4 * dimensionIndex;

    // Halo box values for left neighbour
    _haloBoxes[haloBoxesStartIndex] = _localBoxMin[dimensionIndex] - _skinWidth; 
    _haloBoxes[haloBoxesStartIndex + 1] = _localBoxMin[dimensionIndex] + _cutoffWidth + _skinWidth; 

    // Halo box values for left neighbour
    _haloBoxes[haloBoxesStartIndex + 2] = _localBoxMax[dimensionIndex] - _cutoffWidth - _skinWidth; 
    _haloBoxes[haloBoxesStartIndex + 3] = _localBoxMax[dimensionIndex] + _skinWidth; 
  }
}
