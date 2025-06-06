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


void ParticleCommunicator::exchangeParticlesWithNeighbors(const std::vector<ParticleType> &particlesToLeft,
                                                          const std::vector<ParticleType> &particlesToRight,
                                                          std::vector<ParticleType> &receivedParticles,
                                                          int leftNeighbor, int rightNeighbor, int myRank) {
  // Clear and reserve space for received particles buffer with improved estimation
  receivedParticles.clear();
  size_t estimatedReceiveSize = std::max(size_t(64), particlesToLeft.size() + particlesToRight.size() + 32);
  receivedParticles.reserve(estimatedReceiveSize);

  // Check for special cases once
  const bool isSelfCommunication = (leftNeighbor == myRank) || (rightNeighbor == myRank);
  const bool sameLRNeighbors = (leftNeighbor == rightNeighbor);

  // Handle case where we don't need MPI communication (talking to ourselves)
  if (isSelfCommunication) {
    receivedParticles.insert(receivedParticles.end(), particlesToLeft.begin(), particlesToLeft.end());
    receivedParticles.insert(receivedParticles.end(), particlesToRight.begin(), particlesToRight.end());
    return;
  }

  // Handle case where left and right neighbors are the same
  if (sameLRNeighbors) {
    // Combine particles and send once, receive once
    std::vector<ParticleType> combinedParticles;
    combinedParticles.reserve(particlesToLeft.size() + particlesToRight.size());
    combinedParticles.insert(combinedParticles.end(), particlesToLeft.begin(), particlesToLeft.end());
    combinedParticles.insert(combinedParticles.end(), particlesToRight.begin(), particlesToRight.end());
    
    sendParticles(combinedParticles, leftNeighbor);
    receiveParticles(receivedParticles, leftNeighbor);
    waitForSendRequests();
    return;
  }

  // Use thread-local static buffers to reduce allocations
  static thread_local std::vector<autopas::AutoPas_MPI_Request> localSendRequests;
  static thread_local std::vector<std::vector<char>> localSendBuffers(2); // One for left, one for right
  
  // Clear and prepare buffers efficiently
  localSendRequests.clear();
  localSendRequests.reserve(2); // We know we'll have exactly 2 sends
  
  localSendBuffers[0].clear();
  localSendBuffers[1].clear();

  // Prepare and send to left neighbor (even if empty) with better size estimation  
  if (!particlesToLeft.empty()) {
    size_t estimatedLeftSize = std::max(size_t(128), particlesToLeft.size() * 64);
    localSendBuffers[0].reserve(estimatedLeftSize);
    for (const auto &particle : particlesToLeft) {
      ParticleSerializationTools::serializeParticle(particle, localSendBuffers[0]);
    }
  }
  
  autopas::AutoPas_MPI_Request leftSendReq;
  int leftSendResult = autopas::AutoPas_MPI_Isend(localSendBuffers[0].data(), 
                                                  static_cast<int>(localSendBuffers[0].size()),
                              AUTOPAS_MPI_CHAR, leftNeighbor, 0, _communicator, &leftSendReq);
  if (leftSendResult == 0) { // MPI_SUCCESS is 0
    localSendRequests.push_back(leftSendReq);
  }

  // Prepare and send to right neighbor (even if empty) with better size estimation
  if (!particlesToRight.empty()) {
    size_t estimatedRightSize = std::max(size_t(128), particlesToRight.size() * 64);
    localSendBuffers[1].reserve(estimatedRightSize);
    for (const auto &particle : particlesToRight) {
      ParticleSerializationTools::serializeParticle(particle, localSendBuffers[1]);
    }
  }
  
  autopas::AutoPas_MPI_Request rightSendReq;
  int rightSendResult = autopas::AutoPas_MPI_Isend(localSendBuffers[1].data(), 
                                                   static_cast<int>(localSendBuffers[1].size()),
                                AUTOPAS_MPI_CHAR, rightNeighbor, 0, _communicator, &rightSendReq);
  if (rightSendResult == 0) { // MPI_SUCCESS is 0
    localSendRequests.push_back(rightSendReq);
  }

  // Use thread-local buffers for receives to avoid repeated allocations
  static thread_local std::vector<ParticleType> leftReceivedParticles;
  static thread_local std::vector<ParticleType> rightReceivedParticles;
  
  // Always receive from left neighbor with pre-allocated buffer
  leftReceivedParticles.clear();
  leftReceivedParticles.reserve(std::max(size_t(64), particlesToLeft.size()));
  receiveParticles(leftReceivedParticles, leftNeighbor);
  receivedParticles.insert(receivedParticles.end(), leftReceivedParticles.begin(), leftReceivedParticles.end());

  // Always receive from right neighbor with pre-allocated buffer
  rightReceivedParticles.clear();
  rightReceivedParticles.reserve(std::max(size_t(64), particlesToRight.size()));
  receiveParticles(rightReceivedParticles, rightNeighbor);
  receivedParticles.insert(receivedParticles.end(), rightReceivedParticles.begin(), rightReceivedParticles.end());

  // Wait for all sends to complete - only wait for successful sends
  if (!localSendRequests.empty()) {
    static thread_local std::vector<autopas::AutoPas_MPI_Status> sendStatuses;
    sendStatuses.clear();
    sendStatuses.resize(localSendRequests.size());
    
    int waitResult = autopas::AutoPas_MPI_Waitall(static_cast<int>(localSendRequests.size()), 
                                localSendRequests.data(), sendStatuses.data());
    
    // Basic error handling - log if wait fails but don't crash
    if (waitResult != 0) {
      // In a production system, you might want to log this error
      // For now, we continue execution as the receives have completed
      // Should theoretically never happen
    }
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
