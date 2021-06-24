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
   * @param argc The argument count passed to the main function.
   * @param argv The argument vector passed to the main function.
   * @param dimensionCount The number of dimensions for this domain decomposition.
   * @param globalBoxMin The minimum coordinates of the global domain.
   * @param globalBoxMax The maximum coordinates of the global domain.
   */
  RegularGridDecomposition(const int &dimensionCount, const std::vector<double> &globalBoxMin,
                           const std::vector<double> &globalBoxMax, const double &cutoffWidth, const double &skinWidth);

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
   */
  const int getDomainIndex() { return _domainIndex; }

  /**
   * Returns the number of dimesnions in the domain decomposition.
   */
  const int getDimensionCount() override { return _dimensionCount; }

  /**
   * Returns the minimum coordinates of global domain.
   */
  const std::vector<double> getGlobalBoxMin() override { return _globalBoxMin; }

  /**
   * Returns the maximum coordinates of global domain.
   */
  const std::vector<double> getGlobalBoxMax() override { return _globalBoxMax; }

  /**
   * Returns the minimum coordinates of local domain.
   */
  const std::vector<double> getLocalBoxMin() override { return _localBoxMin; }

  /**
   * Returns the maximum coordinates of local domain.
   */
  const std::vector<double> getLocalBoxMax() override { return _localBoxMax; }

  /**
   * Returns the number of domains in each dimension
   */
  const std::vector<int> getDecomposition() { return _decomposition; }

  /**
   * Checks if the provided coordinates are located in the local domain.
   */
  bool isInsideLocalDomain(const std::vector<double> &coordinates) override;

  /**
   * Checks if the provided coordinates are located in the local domain.
   * Instead of a vector, the coordinates are of type std::array<double, 3> to be compatible with AutoPas.
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
   * Indicates if MPI is enabled and if it will be used.
   * In the case that MPI is enabled, but only one process is being used, this variable will be false.
   */
  bool _mpiIsEnabled;

  /**
   * The number of dimensions in this decomposition.
   */
  int _dimensionCount;

  /**
   * The number of subdomains in this decomposition.
   */
  int _subdomainCount;

  /**
   * The minimum coordinates of the global domain.
   */
  std::vector<double> _globalBoxMin;

  /**
   * The maximum coordinates of the global domain.
   */
  std::vector<double> _globalBoxMax;

  /**
   * The decomposition computed depending on the number of subdomains.
   */
  std::vector<int> _decomposition;

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
  std::vector<int> _domainId;

  /**
   * The indices of the local domain's neighbours.
   * These correspond to the ranks of the processors which own the neigbour domain.
   */
  std::vector<int> _neighbourDomainIndices;

  /**
   * The minimum cooridnates of the local domain.
   */
  std::vector<double> _localBoxMin;

  /**
   * The maximum cooridnates of the local domain.
   */
  std::vector<double> _localBoxMax;

  /**
   * Stores the maximum and minimum coordinates of all halo boxes.
   * The halo box coordinates differ only at a single index from the local box coordinates. Therfore it is
   * enough to store 2 values for each neighbour.
   * The values for the left neighbour in the nth dimension are starting at index 4 * n.
   * The values for the right neighbour in the nth dimension are starting at index 4 * n + 2.
   */
  std::vector<double> _haloBoxes;

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
  void initializeGlobalBox(const std::vector<double> &globalBoxMin, const std::vector<double> &globalBoxMax);

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
   * Updates the boxes used to identify halo particles.
   */
  void updateHaloBoxes();

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
  void sendAndReceiveParticlesLeftAndRight(const std::vector<ParticleType>& particlesToLeft, const std::vector<ParticleType>& particlesToRight, const int& leftNeighbour, const int& rightNeighbour, std::vector<ParticleType>& receivedParticles);

  /**
   * Waits for all send requests to be finished.
   */
  void waitForSendRequests();

  /**
   * Converts a domain id to the domain index, i.e. rank of the local processor.
   */
  int convertIdToIndex(const std::vector<int> &domainIndex);
};
