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
#include "src/options/BoundaryTypeOption.h"

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
   * @param subdivideDimension: Decides if a dimension will be subdivided.
   * @param cutoffWidth: The cutoff width for halo particles.
   * @param skinWidth: The skin width of an autopas container domain.
   * @param reflWidth: The width of the reflective 'skin' in front of the boundary.
   * @param boundaryConditions: An array of boundary conditions in the x, y, and z directions.
   */
  RegularGridDecomposition(const std::array<double, 3> &globalBoxMin, const std::array<double, 3> &globalBoxMax,
                           const std::array<bool, 3> &subdivideDimension, const double &cutoffWidth,
                           const double &skinWidth, const double &reflWidth, const std::array<autopas::BoundaryTypeOption,3> &boundaryConditions);

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
  int getDomainIndex() const override { return _domainIndex; }

  /**
   * Returns the minimum coordinates of global domain.
   * @return bottom left front corner of the global domain.
   */
  std::array<double, 3> getGlobalBoxMin() const override { return _globalBoxMin; }

  /**
   * Returns the maximum coordinates of global domain.
   * @return top right back corner of the global domain.
   */
  std::array<double, 3> getGlobalBoxMax() const override { return _globalBoxMax; }

  /**
   * Returns the minimum coordinates of local domain.
   * @return bottom left front corner of the local domain.
   */
  std::array<double, 3> getLocalBoxMin() const override { return _localBoxMin; }

  /**
   * Returns the maximum coordinates of local domain.
   * @return top right back corner of the local domain.
   */
  std::array<double, 3> getLocalBoxMax() const override { return _localBoxMax; }

  /**
   * Returns the number of domains in each dimension
   * @return vector containing the number of subdomains along each dimension
   */
  std::array<int, 3> getDecomposition() const { return _decomposition; }

  /**
   * Returns the numnber of subdomains in the decomposition.
   * @return numner of subdomains in the decomposition.
   */
  int getSubdomainCount() const { return _subdomainCount; }

  /**
   * Returns the current processes domain id.
   * @return domain id of the current processor
   */
  const std::array<int, 3> getDomainId() const { return _domainId; }

  /**
   * Returns the number of subdomains in the simulation.
   * @return number of subdomains
   */
  int getNumberOfSubdomains() const;

  /**
   * Checks if the provided coordinates are located in the local domain.
   * @param coordinates: The coordinates in question.
   * @return true if the coordinates lie inside the local domain, false otherwise.
   */
  bool isInsideLocalDomain(const std::array<double, 3> &coordinates) const override;

  /**
   * Calculates and returns the extent of the subdomain with inde subdomainIndex.
   * @param subdomainIndex: The index of the subdomain for which to calculate the extent.
   * @return extent of the subdomain with index subdomainIndex.
   */
  std::array<int, 6> getExtentOfSubdomain(const int subdomainIndex) const;

  /**
   * Exchanges halo particles with all neighbours of the provided AutoPasContainer.
   * @param autoPasContainer: The container, where the halo particles originate from.
   */
  void exchangeHaloParticles(SharedAutoPasContainer &autoPasContainer);

  /**
   * Exchanges migrating particles with all neighbours of the provided AutoPasContainer.
   * @param autoPasContainer: The container, where the migrating particles originate from.
   */
  void exchangeMigratingParticles(SharedAutoPasContainer &autoPasContainer);

  /**
   * Reflects particles within a reflective skin along a boundary.
   * @param autoPasContainer: The container, where the migrating particles originate from.
   */
  void reflectParticlesAtBoundaries(SharedAutoPasContainer &autoPasContainer);

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
  bool _mpiCommunicationNeeded;

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
   * Width of the reflective skin in front of a reflective boundary.
   */
   double _reflWidth;

   /**
    * Boundary condition types.
    */
   std::array<autopas::BoundaryTypeOption, 3> _boundaryType;

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
   * Initializes the reflective skin width.
   */
   void initializeReflWidth(const double &reflWidth);

   /**
    * Initializes the Boundary Condition types.
    */
    void initializeBoundaryConditions(const std::array<autopas::BoundaryTypeOption,3> &boundaryConditions);

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
  void sendAndReceiveParticlesLeftAndRight(std::vector<ParticleType> &particlesToLeft,
                                           std::vector<ParticleType> &particlesToRight, const int &leftNeighbour,
                                           const int &rightNeighbour, std::vector<ParticleType> &receivedParticles);

  /**
   * Waits for all send requests to be finished.
   */
  void waitForSendRequests();

  /**
   * Collects the halo particles for the left neighbour.
   * Halo particle positions will be wrapped around the global domain boundary if necessary.
   * @param autoPasContainer: The autopas container which owns the potential halo particles.
   * @param direction: The direction along which the neighbour is located.
   * @param haloParticles: The container the identified halo particles are gathered in to.
   */
  void collectHaloParticlesForLeftNeighbour(SharedAutoPasContainer &autoPasContainer, const size_t &direction,
                                            std::vector<ParticleType> &haloParticles);

  /**
   * Collects the halo particles for the right neighbour.
   * Halo particle positions will be wrapped around the global domain boundary if necessary.
   * @param autoPasContainer: The autopas container which owns the potential halo particles.
   * @param direction: The direction along which the neighbour is located.
   * @param haloParticles: The container the identified halo particles are gathered in to.
   */
  void collectHaloParticlesForRightNeighbour(SharedAutoPasContainer &autoPasContainer, const size_t &direction,
                                             std::vector<ParticleType> &haloParticles);

  /**
   * Categorizes the provided particles as particles for the left or the right neighbour and adds them to the respective
   * output vector. Particle positions will be wrapped around the global domain boundary if necessary.
   * @param particles: The particles which need to be categorized.
   * @param direction: The index of the dimension along which the left and right neighbour lie.
   * @param leftNeighbourParticles: Contains the particles for the left neighbour after function execution.
   * @param rightNeighbourParticles: Contains the particles for the right neighbour after function execution.
   * @param uncategorizedParticles: Contains particles which could neither be assigned to the left nor the right
   * neighbour.
   */
  void categorizeParticlesIntoLeftAndRightNeighbour(const std::vector<ParticleType> &particles, const size_t &direction,
                                                    std::vector<ParticleType> &leftNeighbourParticles,
                                                    std::vector<ParticleType> &rightNeighbourParticles,
                                                    std::vector<ParticleType> &uncategorizedParticles);
};
