/**
 * @file ParticleCommunicator.h
 * @author J. Körner
 * @date 28.07.2021
 */
#pragma once

#include <vector>

#include "autopas/utils/WrapMPI.h"
#include "src/TypeDefinitions.h"

/**
 * Provides tools to communicate particles between MPI ranks.
 */
class ParticleCommunicator {
 public:
  /**
   * Constructor
   * @param communicator: The MPI communicator used for communication.
   */
  explicit ParticleCommunicator(const autopas::AutoPas_MPI_Comm &communicator);

  /**
   * Destructor.
   */
  ~ParticleCommunicator() = default;

  /**
   * Sends particles of type ParticleType to a receiver.
   * @param particles The particles to be sent to the receiver.
   * @param receiver The recipient of the particels.
   */
  void sendParticles(const std::vector<ParticleType> &particles, const int &receiver);

  /**
   * Received particles sent by a sender.
   * @param receivedParticles The container where the received particles will be stored.
   * @param source The sender id/rank.
   */
  void receiveParticles(std::vector<ParticleType> &receivedParticles, const int &source);

  /**
   * Waits for all send requests to be finished.
   */
  void waitForSendRequests();

 private:
  /**
   * The MPI communicator used for the communications.
   */
  autopas::AutoPas_MPI_Comm _communicator;

  /**
   * A temporary buffer used for MPI send requests.
   */
  std::vector<autopas::AutoPas_MPI_Request> _sendRequests;

  /**
   * A temporary buffer for data which is sent by MPI_Send.
   */
  std::vector<std::vector<char>> _sendBuffers;

  /**
   * Sends data to a specific neighbour of this domain.
   * @param sendBuffer The buffer which will be sent to the neighbour.
   * @param neighbour The neighbour to which the data will be sent.
   */
  void sendDataToNeighbour(const std::vector<char> &sendBuffer, const int &neighbour);

  /**
   * Received data which has been sent by a specifig neighbour of this domain.
   * @param neighbour The neighbour where the data originates from.
   * @param dataBuffer The buffer where the received data will be stored.
   */
  void receiveDataFromNeighbour(const int &neighbour, std::vector<char> &dataBuffer);
};
