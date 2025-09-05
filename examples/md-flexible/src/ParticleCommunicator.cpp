/**
 * @file ParticleCommunicator.h
 * @author J. Körner
 * @date 28.07.2021
 */
#include "ParticleCommunicator.h"

#include <vector>

#include "ParticleSerializationTools.h"

ParticleCommunicator::ParticleCommunicator(const autopas::AutoPas_MPI_Comm &communicator)
    : _communicator(communicator) {}

void ParticleCommunicator::sendParticles(const std::vector<ParticleType> &particles, const int &receiver) {
  std::vector<char> buffer;

  for (const auto &particle : particles) {
    ParticleSerializationTools::serializeParticle(particle, buffer);
  }

  sendDataToNeighbor(buffer, receiver);
}

void ParticleCommunicator::receiveParticles(std::vector<ParticleType> &receivedParticles, const int &source) {
  std::vector<char> receiveBuffer;

  receiveDataFromNeighbor(source, receiveBuffer);

  if (!receiveBuffer.empty()) {
    ParticleSerializationTools::deserializeParticles(receiveBuffer, receivedParticles);
  }
}

void ParticleCommunicator::sendAndReceiveParticlesLeftAndRight(const std::vector<ParticleType> &particlesToLeft,
                                                               const std::vector<ParticleType> &particlesToRight,
                                                               std::vector<ParticleType> &receivedParticles,
                                                               int leftNeighbor, int rightNeighbor, int dimension) {
  // Create unique tags for each direction: dimension*2 + {0=left, 1=right}
  const int leftTag = dimension * 2;
  const int rightTag = dimension * 2 + 1;

  // Pre-reserve buffers to avoid reallocations during serialization
  const size_t particleSize =
#if MD_FLEXIBLE_MODE == MULTISITE
      200;  // AttributesSize for multisite
#else
      120;  // AttributesSize for single site
#endif

  // Reuse buffers to avoid repeated allocations
  const int leftCount = static_cast<int>(particlesToLeft.size());
  const int rightCount = static_cast<int>(particlesToRight.size());

  // Prepare reusable send buffers
  _reusableLeftSendBuffer.clear();
  _reusableLeftSendBuffer.reserve(sizeof(int) + leftCount * particleSize);
  _reusableLeftSendBuffer.resize(sizeof(int));
  std::memcpy(_reusableLeftSendBuffer.data(), &leftCount, sizeof(int));
  for (const auto &particle : particlesToLeft) {
    ParticleSerializationTools::serializeParticle(particle, _reusableLeftSendBuffer);
  }

  _reusableRightSendBuffer.clear();
  _reusableRightSendBuffer.reserve(sizeof(int) + rightCount * particleSize);
  _reusableRightSendBuffer.resize(sizeof(int));
  std::memcpy(_reusableRightSendBuffer.data(), &rightCount, sizeof(int));
  for (const auto &particle : particlesToRight) {
    ParticleSerializationTools::serializeParticle(particle, _reusableRightSendBuffer);
  }

  // Prepare reusable receive buffers
  const size_t maxParticles = 10000;  // reasonable upper bound
  const size_t maxBufferSize = sizeof(int) + maxParticles * particleSize;
  if (_reusableLeftRecvBuffer.size() < maxBufferSize) {
    _reusableLeftRecvBuffer.resize(maxBufferSize);
  }
  if (_reusableRightRecvBuffer.size() < maxBufferSize) {
    _reusableRightRecvBuffer.resize(maxBufferSize);
  }

  // Use non-blocking sends and receives with unique tags - eliminates MPI_Probe!
  autopas::AutoPas_MPI_Request sendRequests[2];
  autopas::AutoPas_MPI_Request recvRequests[2];

  // Post non-blocking receives first (with unique tags)
  autopas::AutoPas_MPI_Irecv(_reusableLeftRecvBuffer.data(), static_cast<int>(_reusableLeftRecvBuffer.size()),
                             AUTOPAS_MPI_CHAR, leftNeighbor, rightTag, _communicator, &recvRequests[0]);
  autopas::AutoPas_MPI_Irecv(_reusableRightRecvBuffer.data(), static_cast<int>(_reusableRightRecvBuffer.size()),
                             AUTOPAS_MPI_CHAR, rightNeighbor, leftTag, _communicator, &recvRequests[1]);

  // Post non-blocking sends
  autopas::AutoPas_MPI_Isend(_reusableLeftSendBuffer.data(), static_cast<int>(_reusableLeftSendBuffer.size()),
                             AUTOPAS_MPI_CHAR, leftNeighbor, leftTag, _communicator, &sendRequests[0]);
  autopas::AutoPas_MPI_Isend(_reusableRightSendBuffer.data(), static_cast<int>(_reusableRightSendBuffer.size()),
                             AUTOPAS_MPI_CHAR, rightNeighbor, rightTag, _communicator, &sendRequests[1]);

  // Wait for all communications to complete
  autopas::AutoPas_MPI_Waitall(2, sendRequests, AUTOPAS_MPI_STATUS_IGNORE);
  autopas::AutoPas_MPI_Waitall(2, recvRequests, AUTOPAS_MPI_STATUS_IGNORE);

  // Deserialize received particles from left neighbor
  int leftRecvCount;
  std::memcpy(&leftRecvCount, _reusableLeftRecvBuffer.data(), sizeof(int));
  if (leftRecvCount > 0) {
    std::vector<char> leftPayload(_reusableLeftRecvBuffer.begin() + sizeof(int),
                                  _reusableLeftRecvBuffer.begin() + sizeof(int) + leftRecvCount * particleSize);
    ParticleSerializationTools::deserializeParticles(leftPayload, receivedParticles);
  }

  // Deserialize received particles from right neighbor
  int rightRecvCount;
  std::memcpy(&rightRecvCount, _reusableRightRecvBuffer.data(), sizeof(int));
  if (rightRecvCount > 0) {
    std::vector<char> rightPayload(_reusableRightRecvBuffer.begin() + sizeof(int),
                                   _reusableRightRecvBuffer.begin() + sizeof(int) + rightRecvCount * particleSize);
    ParticleSerializationTools::deserializeParticles(rightPayload, receivedParticles);
  }
}

void ParticleCommunicator::waitForSendRequests() {
  std::vector<autopas::AutoPas_MPI_Status> sendStates;
  sendStates.resize(_sendRequests.size());
  autopas::AutoPas_MPI_Waitall(static_cast<int>(_sendRequests.size()), _sendRequests.data(), sendStates.data());
  _sendRequests.clear();
  _sendBuffers.clear();
}

void ParticleCommunicator::sendDataToNeighbor(const std::vector<char> &sendBuffer, const int &neighbour) {
  _sendBuffers.push_back(sendBuffer);

  autopas::AutoPas_MPI_Request sendRequest{};
  _sendRequests.push_back(sendRequest);

  autopas::AutoPas_MPI_Isend(_sendBuffers.back().data(), _sendBuffers.back().size(), AUTOPAS_MPI_CHAR, neighbour, 0,
                             _communicator, &_sendRequests.back());
}

void ParticleCommunicator::receiveDataFromNeighbor(const int &neighbour, std::vector<char> &receiveBuffer) {
  autopas::AutoPas_MPI_Status status;
  autopas::AutoPas_MPI_Probe(neighbour, 0, _communicator, &status);

  int receiveBufferSize = 0;
  autopas::AutoPas_MPI_Get_count(&status, AUTOPAS_MPI_CHAR, &receiveBufferSize);
  receiveBuffer.resize(receiveBufferSize);

  autopas::AutoPas_MPI_Recv(receiveBuffer.data(), receiveBufferSize, AUTOPAS_MPI_CHAR, neighbour, 0, _communicator,
                            AUTOPAS_MPI_STATUS_IGNORE);
}
