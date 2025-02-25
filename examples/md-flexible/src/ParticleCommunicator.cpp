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

  size_t beforeSize = receivedParticles.size();
  int rank;
  autopas::AutoPas_MPI_Comm_rank(_communicator, &rank);

  if (!receiveBuffer.empty()) {
    ParticleSerializationTools::deserializeParticles(receiveBuffer, receivedParticles);
  }
}

void ParticleCommunicator::sendParticlePositions(const std::vector<ParticleType> &particles, const int &receiver) {
#if ZONAL_RECOMMUNICATION_OPT == 0
  sendParticles(particles, receiver);
#else
  std::vector<char> buffer;

  ParticleSerializationTools::serializeParticlePositions(particles, buffer);

  sendDataToNeighbor(buffer, receiver);
#endif
}

void ParticleCommunicator::receiveParticlePositions(std::vector<ParticleType> &receivedParticles, const int &source) {
#if ZONAL_RECOMMUNICATION_OPT == 0
  receiveParticles(receivedParticles, source);
#else
  std::vector<char> receiveBuffer;

  receiveDataFromNeighbor(source, receiveBuffer);

  if (!receiveBuffer.empty()) {
    ParticleSerializationTools::deserializeParticlePositions(receiveBuffer, receivedParticles);
  }
#endif
}

void ParticleCommunicator::sendParticleForces(const std::vector<ParticleType> &particles, const int &receiver) {
#if ZONAL_RECOMMUNICATION_OPT == 0
  sendParticles(particles, receiver);
#else
  std::vector<char> buffer;

  ParticleSerializationTools::serializeParticleForces(particles, buffer);

  sendDataToNeighbor(buffer, receiver);
#endif
}

void ParticleCommunicator::receiveParticleForces(std::vector<ParticleType> &receivedParticles, const int &source) {
#if ZONAL_RECOMMUNICATION_OPT == 0
  receiveParticles(receivedParticles, source);
#else
  std::vector<char> receiveBuffer;

  receiveDataFromNeighbor(source, receiveBuffer);

  if (!receiveBuffer.empty()) {
    ParticleSerializationTools::deserializeParticleForces(receiveBuffer, receivedParticles);
  }
#endif
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
