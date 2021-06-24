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

RegularGridDecomposition::RegularGridDecomposition(const int &dimensionCount, const std::vector<double> &globalBoxMin,
                                                   const std::vector<double> &globalBoxMax, const double &cutoffWidth,
                                                   const double &skinWidth)
    : _dimensionCount(dimensionCount), _cutoffWidth(cutoffWidth), _skinWidth(skinWidth) {


#if defined(AUTOPAS_INCLUDE_MPI)
  _mpiIsEnabled = true;
#else
  _mpiIsEnabled = false;
#endif

  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &_subdomainCount);

  if (_subdomainCount == 1){
    _mpiIsEnabled = false;
  }

  if (_mpiIsEnabled) {
    std::cout << "MPI will be used." << std::endl;
  }
  else {
    std::cout << "MPI will not be used." << std::endl;
  }

  DomainTools::generateDecomposition(_subdomainCount, _dimensionCount, _decomposition);

  initializeMPICommunicator();

  initializeLocalDomain();

  initializeGlobalBox(globalBoxMin, globalBoxMax);

  initializeLocalBox();

  initializeNeighbourIds();

  std::cout
    << "DomainIndex: " << _domainIndex << ", "
    << "DomainId: " << autopas::utils::ArrayUtils::to_string(_domainId) << ", "
    << "Neighbours: " << autopas::utils::ArrayUtils::to_string(_neighbourDomainIndices) << ", "
    << "LocalBox: " << autopas::utils::ArrayUtils::to_string(_localBoxMin) << ", "
    << autopas::utils::ArrayUtils::to_string(_localBoxMax) << std::endl;
}

RegularGridDecomposition::~RegularGridDecomposition() {}

void RegularGridDecomposition::update() { updateLocalBox(); }

void RegularGridDecomposition::initializeDecomposition() {
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

void RegularGridDecomposition::initializeMPICommunicator() {
  std::vector<int> periods(_dimensionCount, 1);
  autopas::AutoPas_MPI_Cart_create(AUTOPAS_MPI_COMM_WORLD, _dimensionCount, _decomposition.data(), periods.data(), true,
                                   &_communicator);
  autopas::AutoPas_MPI_Comm_rank(_communicator, &_domainIndex);
}

void RegularGridDecomposition::initializeLocalDomain() {
  _domainId.resize(_dimensionCount);
  autopas::AutoPas_MPI_Comm_rank(_communicator, &_domainIndex);

  std::vector<int> periods(_dimensionCount, 1);
  autopas::AutoPas_MPI_Cart_get(_communicator, _dimensionCount, _decomposition.data(), periods.data(),
                                _domainId.data());

  if (_domainId.empty()) {
    _domainId.resize(_dimensionCount, 0);
  }
}

void RegularGridDecomposition::initializeLocalBox() {
  _localBoxMin.resize(_dimensionCount);
  _localBoxMax.resize(_dimensionCount);
  updateLocalBox();
  _haloBoxes.resize(4 * _dimensionCount);
  updateHaloBoxes();
}

void RegularGridDecomposition::initializeNeighbourIds() {
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

void RegularGridDecomposition::updateLocalBox() {
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

void RegularGridDecomposition::initializeGlobalBox(const std::vector<double> &globalBoxMin,
                                                   const std::vector<double> &globalBoxMax) {
  _globalBoxMin.resize(_dimensionCount);
  _globalBoxMax.resize(_dimensionCount);
  for (int i = 0; i < _dimensionCount; ++i) {
    _globalBoxMin[i] = globalBoxMin[i];
    _globalBoxMax[i] = globalBoxMax[i];
  }
}

bool RegularGridDecomposition::isInsideLocalDomain(const std::vector<double> &coordinates) {
  return DomainTools::isInsideDomain(coordinates, _localBoxMin, _localBoxMax);
}

bool RegularGridDecomposition::isInsideLocalDomain(const std::array<double, 3> &coordinates) {
  std::vector<double> coordinatesVector;
  coordinatesVector.insert(coordinatesVector.begin(), coordinates.begin(), coordinates.end());
  return isInsideLocalDomain(coordinatesVector);
}

void RegularGridDecomposition::exchangeHaloParticles(SharedAutoPasContainer &autoPasContainer) {
  int dimensionCount = _localBoxMin.size();
  int neighbourCount = dimensionCount * 2;

  const std::array<double, 3> globalBoxMin = {_globalBoxMin[0], _globalBoxMin[1], _globalBoxMin[2]};
  const std::array<double, 3> globalBoxMax = {_globalBoxMax[0], _globalBoxMax[1], _globalBoxMax[2]};
  const std::array<double, 3> globalBoxLengths = autopas::utils::ArrayMath::sub(globalBoxMax, globalBoxMin);
  
  for (int i = 0; i < dimensionCount; ++i) {
    int leftNeighbour = _neighbourDomainIndices[(i * 2) % neighbourCount];
    int rightNeighbour = _neighbourDomainIndices[(i * 2 + 1) % neighbourCount];

    std::vector<ParticleType> haloParticles, particlesForLeftNeighbour, particlesForRightNeighbour;

    for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
      std::array<double, 3> particlePosition = particle->getR();
      if (particlePosition[i] >= _haloBoxes[4 * i] && particlePosition[i] < _haloBoxes[4 * i + 1]) {
        particlesForLeftNeighbour.push_back(*particle);
        if (_localBoxMin[i] == _globalBoxMin[i]) {
          particlePosition[i] = particlePosition[i] + globalBoxLengths[i];
          particlesForLeftNeighbour.back().setR(particlePosition);
        }
      }
      else if (particlePosition[i] > _haloBoxes[4 * i + 2] && particlePosition[i] < _haloBoxes[4 * i + 3]) {
        particlesForRightNeighbour.push_back(*particle);
        if (_localBoxMax[i] == _globalBoxMax[i]) {
          particlePosition[i] = particlePosition[i] - globalBoxLengths[i];
          particlesForRightNeighbour.back().setR(particlePosition);
        }
      }
    }

    sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbour, particlesForRightNeighbour, leftNeighbour, rightNeighbour, haloParticles);
    
    for (auto &particle : haloParticles) {
      autoPasContainer->addOrUpdateHaloParticle(particle);
    }

    particlesForLeftNeighbour.clear();
    particlesForRightNeighbour.clear();

    // index of next dimension
    int j = (i + 1) % dimensionCount;

    leftNeighbour = _neighbourDomainIndices[(j * 2) % neighbourCount];
    rightNeighbour = _neighbourDomainIndices[(j * 2 + 1) % neighbourCount];
    
    for (auto &particle : haloParticles) {
      std::array<double, 3> particlePosition = particle.getR();
      if (particlePosition[j] >= _haloBoxes[4 * j] && particlePosition[j] < _haloBoxes[4 * j + 1]) {
        particlesForLeftNeighbour.push_back(particle);
        if (_localBoxMin[j] == _globalBoxMin[j]) {
          particlePosition[j] = particlePosition[j] + globalBoxLengths[j];
          particlesForLeftNeighbour.back().setR(particlePosition);
        }
      }
      else if (particlePosition[j] > _haloBoxes[4 * j + 2] && particlePosition[j] < _haloBoxes[4 * j + 3]) {
        particlesForRightNeighbour.push_back(particle);
        if (_localBoxMax[j] == _globalBoxMax[j]) {
          particlePosition[j] = particlePosition[j] - globalBoxLengths[j];
          particlesForRightNeighbour.back().setR(particlePosition);
        }
      }
    }

    haloParticles.clear();

    sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbour, particlesForRightNeighbour, leftNeighbour, rightNeighbour, haloParticles);

    for (auto &particle : haloParticles) {
      autoPasContainer->addOrUpdateHaloParticle(particle);
    }
  }
}

void RegularGridDecomposition::exchangeMigratingParticles(SharedAutoPasContainer &autoPasContainer) {
  auto [emigrants, updated] = autoPasContainer->updateContainer(false);

  if (updated) {
    const std::array<double, 3> globalBoxMin = {_globalBoxMin[0], _globalBoxMin[1], _globalBoxMin[2]};
    const std::array<double, 3> globalBoxMax = {_globalBoxMax[0], _globalBoxMax[1], _globalBoxMax[2]};
    const std::array<double, 3> globalBoxLengths = autopas::utils::ArrayMath::sub(globalBoxMax, globalBoxMin);

    int neighbourCount = _neighbourDomainIndices.size();

    for (int i = 0; i < _dimensionCount; ++i) {
      std::vector<ParticleType> immigrants, migrants, remainingEmigrants;

      std::vector<ParticleType> particlesForLeftNeighbour;
      std::vector<ParticleType> particlesForRightNeighbour;

      int leftNeighbour = _neighbourDomainIndices[(i * 2) % neighbourCount];
      int rightNeighbour = _neighbourDomainIndices[(i * 2 + 1) % neighbourCount];

      for (const auto &particle : emigrants) {
        std::array<double, 3> position = particle.getR();
        if (position[i] < _localBoxMin[i]){
          particlesForLeftNeighbour.push_back(particle);
          if (_localBoxMin[i] == _globalBoxMin[i]){
            position[i] = std::min(std::nextafter(_globalBoxMax[i], _globalBoxMin[i]), position[i] + globalBoxLengths[i]);
            particlesForLeftNeighbour.back().setR(position);
          }
        }
        else if (position[i] >= _localBoxMax[i]){
          particlesForRightNeighbour.push_back(particle);
          if (_localBoxMax[i] == _globalBoxMax[i]){
            position[i] = std::max(_globalBoxMin[i], position[i] - globalBoxLengths[i]);
            particlesForRightNeighbour.back().setR(position);
          }
        }
        else {
          remainingEmigrants.push_back(particle);
        }
      }
      emigrants = remainingEmigrants;

      sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbour, particlesForRightNeighbour, leftNeighbour, rightNeighbour, immigrants);

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

      leftNeighbour = _neighbourDomainIndices[(j * 2) % neighbourCount];
      rightNeighbour = _neighbourDomainIndices[(j * 2 + 1) % neighbourCount];

      for (const auto &particle : migrants) {
        std::array<double, 3> position = particle.getR();
        if (position[j] < _localBoxMin[j]){
          particlesForLeftNeighbour.push_back(particle);
          if (_localBoxMin[j] == _globalBoxMin[j]){
            position[j] = std::min(std::nextafter(_globalBoxMax[j], _globalBoxMin[j]), position[j] + globalBoxLengths[j]);
            particlesForLeftNeighbour.back().setR(position);
          }
        }
        else if (position[j] >= _localBoxMax[j]){
          particlesForRightNeighbour.push_back(particle);
          if (_localBoxMax[j] == _globalBoxMax[j]){
            position[j] = std::max(_globalBoxMin[j], position[j] - globalBoxLengths[j]);
            particlesForRightNeighbour.back().setR(position);
          }
        }
      }

      sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbour, particlesForRightNeighbour, leftNeighbour, rightNeighbour, immigrants);

      for (auto &particle : immigrants) {
        autoPasContainer->addParticle(particle);
      }

      waitForSendRequests();
    }
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
    ParticleSerializationTools::deserializeParticleData(receiveBuffer, receivedParticles);
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

void RegularGridDecomposition::sendAndReceiveParticlesLeftAndRight(const std::vector<ParticleType>& particlesToLeft, const std::vector<ParticleType>& particlesToRight, const int& leftNeighbour, const int& rightNeighbour, std::vector<ParticleType>& receivedParticles){
  if (_mpiIsEnabled && leftNeighbour != _domainIndex) {
    sendParticles(particlesToLeft, leftNeighbour);
    sendParticles(particlesToRight, rightNeighbour);

    receiveParticles(receivedParticles, leftNeighbour);
    receiveParticles(receivedParticles, rightNeighbour);

    waitForSendRequests();
  }
  else{
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

int RegularGridDecomposition::convertIdToIndex(const std::vector<int> &domainId) {
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

void RegularGridDecomposition::updateHaloBoxes() {
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
