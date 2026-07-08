/**
 * @file ParticleCommunicator.h
 * @author J. Körner
 * @date 28.07.2021
 */
#include "ParticleCommunicator.h"

#include <vector>

#include "ParticleSerializationTools.h"

ParticleCommunicator::ParticleCommunicator(const autopas::AutoPas_MPI_Comm &communicator) : _MPIComm(communicator) {}

void ParticleCommunicator::sendParticles(const std::vector<ParticleType> &particles, const int &receiver,
                                         const std::optional<Direction> direction) {
  // Use direction to select a reusable buffer. If no direction, use the "left" buffer arbitrarily.
  const bool useLeftBuffer = direction.has_value() ? (direction.value() == left) : true;
  auto &buffer = useLeftBuffer ? _reusableLeftSendBuffer : _reusableRightSendBuffer;

  // Reserve extra space in reusable buffer if needed.
  buffer.clear();
  buffer.reserve(particles.size() * ParticleSerializationTools::AttributesSize);

  for (const auto &particle : particles) {
    ParticleSerializationTools::serializeParticle(particle, buffer);
  }

  sendDataToNeighbor(buffer, receiver);
}

void ParticleCommunicator::receiveParticles(std::vector<ParticleType> &receivedParticles, const int &source,
                                            const std::optional<Direction> direction) {
  // Use direction to select a reusable buffer. If no direction, use the "left" buffer arbitrarily.
  const bool useLeftBuffer = direction.has_value() ? (direction.value() == left) : true;
  auto &buffer = useLeftBuffer ? _reusableLeftRecvBuffer : _reusableRightRecvBuffer;

  receiveDataFromNeighbor(source, buffer);

  if (not buffer.empty()) {
    ParticleSerializationTools::deserializeParticles(buffer, receivedParticles);
  }
}

void ParticleCommunicator::waitForSendRequests() {
  std::vector<autopas::AutoPas_MPI_Status> sendStates;
  sendStates.resize(_sendRequests.size());
  autopas::AutoPas_MPI_Waitall(static_cast<int>(_sendRequests.size()), _sendRequests.data(), sendStates.data());
  _sendRequests.clear();
  _sendBuffers.clear();
}

void ParticleCommunicator::sendDataToNeighbor(const std::vector<char> &sendBuffer, const int &neighbor) {
  autopas::AutoPas_MPI_Request sendRequest{};
  _sendRequests.push_back(sendRequest);

  autopas::AutoPas_MPI_Isend(_sendBuffers.back().data(), _sendBuffers.back().size(), AUTOPAS_MPI_CHAR, neighbor, 0,
                             _MPIComm, &_sendRequests.back());
}

void ParticleCommunicator::receiveDataFromNeighbor(const int &neighbour, std::vector<char> &receiveBuffer) const {
  autopas::AutoPas_MPI_Status status;
  autopas::AutoPas_MPI_Probe(neighbour, 0, _MPIComm, &status);

  int receiveBufferSize = 0;
  autopas::AutoPas_MPI_Get_count(&status, AUTOPAS_MPI_CHAR, &receiveBufferSize);
  receiveBuffer.resize(receiveBufferSize);

  autopas::AutoPas_MPI_Recv(receiveBuffer.data(), receiveBufferSize, AUTOPAS_MPI_CHAR, neighbour, 0, _MPIComm,
                            AUTOPAS_MPI_STATUS_IGNORE);
}
