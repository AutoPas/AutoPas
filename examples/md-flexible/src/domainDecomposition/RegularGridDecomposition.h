/**
 * @file RegularGridDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include <list>
#include <memory>

#include "DomainDecomposition.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/WrapMPI.h"
#include "src/TypeDefinitions.h"

/**
 * This class can be used as a domain decomposition which divides the domain in equal sized rectangular subdomains.
 * The number of subdomains is equal to the number of MPI processes available.
 */
class RegularGridDecomposition final : public DomainDecomposition {
 public:
  /**
   * Constructor.
   * @param globalBoxMin: The minimum coordinates of the global domain.
   * @param globalBoxMax: The maximum coordinates of the global domain.
   * @param cutoffWidth: The cutoff width for halo particles.
   * @param skinWidth: The skin width of an autopas container domain.
   */
  RegularGridDecomposition(const std::array<double, 3> &globalBoxMin, const std::array<double, 3> &globalBoxMax,
                           const double &cutoffWidth, const double &skinWidth);

  /**
   * Destructor.
   */
  virtual ~RegularGridDecomposition();

  /**
   * Type for the AutoPas container
   */
  using SharedAutoPasContainer = std::shared_ptr<autopas::AutoPas<ParticleType>>;

  /**
   * Used to update the domain to the current topology.
   * Currently does nothing
   */
  void update() override;

  /**
   * Returns the index of the local domain in the global domain context.
   * @return domain index.
   */
  const int getDomainIndex() override { return _domainIndex; }

  /**
   * Returns the minimum coordinates of global domain.
   * @return bottom left front corner of the global domain.
   */
  const std::array<double, 3> getGlobalBoxMin() override { return _globalBoxMin; }

  /**
   * Returns the maximum coordinates of global domain.
   * @return top right back corner of the global domain.
   */
  const std::array<double, 3> getGlobalBoxMax() override { return _globalBoxMax; }

  /**
   * Returns the minimum coordinates of local domain.
   * @return bottom left front corner of the local domain.
   */
  const std::array<double, 3> getLocalBoxMin() override { return _localBoxMin; }

  /**
   * Returns the maximum coordinates of local domain.
   * @return top right back corner of the local domain.
   */
  const std::array<double, 3> getLocalBoxMax() override { return _localBoxMax; }

  /**
   * Returns the number of domains in each dimension
   * @return vector containing the number of subdomains along each dimension
   */
  const std::array<int, 3> getDecomposition() { return _decomposition; }

  /**
   * Checks if the provided coordinates are located in the local domain.
   * @param coordinates: The coordinates in question.
   * @return true if the coordinates lie inside the local domain, false otherwise.
   */
  bool isInsideLocalDomain(const std::array<double, 3> &coordinates) override;

  /**
   * Exchanges halo particles with all neighbours of the provided AutoPasContainer.
   * @param autoPasContainer The container, where the halo particles originate from.
   */
  void exchangeHaloParticles(SharedAutoPasContainer &autoPasContainer);

  /**
   * Exchanges migrating particles with all neighbours of the provided AutoPasContainer.
   * @param autoPasContainer The container, where the migrating particles originate from.
   */
  void exchangeMigratingParticles(SharedAutoPasContainer &autoPasContainer);

 private:
  /**
   * The number of neighbours of a rectangular domain.
   * This number does not include diagonal neighbours.
   */
  static constexpr int _neighbourCount = 6;

  /**
   * The number of dimensions in the simulation domain.
   */
  static constexpr int _dimensionCount = 3;

  /**
   * Indicates if MPI is enabled and if it will be used.
   * In the case that MPI is enabled, but only one process is being used, this variable will be false.
   */
  bool _mpiIsEnabled;

  /**
   * The number of subdomains in this decomposition.
   */
  int _subdomainCount;

  /**
   * The minimum coordinates of the global domain.
   */
  std::array<double, 3> _globalBoxMin;

  /**
   * The maximum coordinates of the global domain.
   */
  std::array<double, 3> _globalBoxMax;

  /**
   * The decomposition computed depending on the number of subdomains.
   */
  std::array<int, 3> _decomposition;

  /**
   * The MPI communicator containing all processes which own a subdomain in this decomposition.
   */
  autopas::AutoPas_MPI_Comm _communicator;

  /**
   * Stores the domain cutoff width.
   */
  double _cutoffWidth;

  /**
   * Stores the domain skin width.
   */
  double _skinWidth;

  /**
   * The index of the current processor's domain.
   * This also is the rank of the current processor.
   */
  int _domainIndex;

  /**
   * The ID of the current processor's domain.
   */
  std::array<int, 3> _domainId;

  /**
   * The indices of the local domain's neighbours.
   * These correspond to the ranks of the processors which own the neigbour domain.  */
  std::array<int, 6> _neighbourDomainIndices;

  /**
   * The minimum cooridnates of the local domain.
   */
  std::array<double, 3> _localBoxMin;

  /**
   * The maximum cooridnates of the local domain.
   */
  std::array<double, 3> _localBoxMax;

  /**
   * A temporary buffer used for MPI send requests.
   */
  std::vector<autopas::AutoPas_MPI_Request> _sendRequests;

  /**
   * A temporary buffer for data which is sent by MPI_Send.
   */
  std::vector<std::vector<char>> _sendBuffers;

  /**
   * Initializes the decomposition of the domain.
   * This needs to be called before initializeMPICommunicator.
   */
  void initializeDecomposition();

  /**
   * Initializes the MPI communicator.
   * This needs to be called before initializeLocalDomain.
   */
  void initializeMPICommunicator();

  /**
   * Initializes the local domain.
   * This needs to be called before initializeLocalBox.
   */
  void initializeLocalDomain();

  /**
   * Initializes the global domain coordinates.
   */
  void initializeGlobalBox(const std::array<double, 3> &globalBoxMin, const std::array<double, 3> &globalBoxMax);

  /**
   * Initializes the local domain coordinates.
   * This needs to be called after initializeLocalDomain and initialzieGlobalDomain.
   */
  void initializeLocalBox();

  /**
   * Initializes the neighbour ids.
   * This needs to be called after initializeLocalDomain.
   */
  void initializeNeighbourIds();

  /**
   * Updates the local box.
   */
  void updateLocalBox();

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
   * Received data which has been sent by a specifig neighbour of this domain.
   * @param neighbour The neighbour where the data originates from.
   * @param dataBuffer The buffer where the received data will be stored.
   */
  void receiveDataFromNeighbour(const int &neighbour, std::vector<char> &dataBuffer);

  /**
   * Sends data to a specific neighbour of this domain.
   * @param sendBuffer The buffer which will be sent to the neighbour.
   * @param neighbour The neighbour to which the data will be sent.
   */
  void sendDataToNeighbour(std::vector<char> sendBuffer, const int &neighbour);

  /**
   * Sends and also receives particles to and from the left and right neighbours.
   * @param particlesToLeft: Particles which get send to the left neighbour.
   * @param particlesToRight: Particles which get send to the right neighbor.
   * @param leftNeighbour: The left neighbour's index / rank.
   * @param rightNeighbour: The right neighbour's index / rank.
   * @param receivedParticles: Container for the particles received from either neighbour.
   */
  void sendAndReceiveParticlesLeftAndRight(const std::vector<ParticleType> &particlesToLeft,
                                           const std::vector<ParticleType> &particlesToRight, const int &leftNeighbour,
                                           const int &rightNeighbour, std::vector<ParticleType> &receivedParticles);

  /**
   * Waits for all send requests to be finished.
   */
  void waitForSendRequests();

  /**
   * Converts a domain id to the domain index, i.e. rank of the local processor.
   */
  int convertIdToIndex(const std::array<int, 3> &domainIndex);
};
