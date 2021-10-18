/**
 * @file RegularGridDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#if defined(AUTOPAS_ENABLE_ALLLBL)
#include <ALL.hpp>
#endif
#include <memory>

#include "DomainDecomposition.h"
#include "autopas/utils/WrapMPI.h"
#include "src/TypeDefinitions.h"
#include "src/configuration/MDFlexConfig.h"

/**
 * This class can be used as a domain decomposition which divides the domain in equal sized rectangular subdomains.
 * The number of subdomains is equal to the number of MPI processes available.
 */
class RegularGridDecomposition final : public DomainDecomposition {
 public:
  /**
   * Constructor.
   * @param configuration: The configuration for definig the decomposition properties
   */
  RegularGridDecomposition(const MDFlexConfig &configuration);

  /**
   * Destructor.
   */
  ~RegularGridDecomposition() override;

  /**
   * Used to update the domain to the current topology.
   * Handles the diffuse load balancing by resizing the domains according to their work done.
   * Currently the metric used as work is the timins created by the simulation timer. Although they are of type 'long'
   * the update function use type 'double' because the work will be implicitly converted to 'double' during load
   * balancing anyway.
   * @param work: The work performed in the AutoPas container.
   */
  void update(const double &work) override;

  /**
   * Returns the index of the local domain in the global domain context.
   * @return domain index.
   */
  [[nodiscard]] int getDomainIndex() const override { return _domainIndex; }

  /**
   * Returns the minimum coordinates of global domain.
   * @return bottom left front corner of the global domain.
   */
  [[nodiscard]] std::array<double, 3> getGlobalBoxMin() const override { return _globalBoxMin; }

  /**
   * Returns the maximum coordinates of global domain.
   * @return top right back corner of the global domain.
   */
  [[nodiscard]] std::array<double, 3> getGlobalBoxMax() const override { return _globalBoxMax; }

  /**
   * Returns the minimum coordinates of local domain.
   * @return bottom left front corner of the local domain.
   */
  [[nodiscard]] std::array<double, 3> getLocalBoxMin() const override { return _localBoxMin; }

  /**
   * Returns the maximum coordinates of local domain.
   * @return top right back corner of the local domain.
   */
  [[nodiscard]] std::array<double, 3> getLocalBoxMax() const override { return _localBoxMax; }

  /**
   * Returns the number of domains in each dimension
   * @return vector containing the number of subdomains along each dimension
   */
  [[nodiscard]] std::array<int, 3> getDecomposition() const { return _decomposition; }

  /**
   * Returns the numnber of subdomains in the decomposition.
   * @return numner of subdomains in the decomposition.
   */
  [[nodiscard]] int getSubdomainCount() const { return _subdomainCount; }

  /**
   * Returns the current processes domain id.
   * @return domain id of the current processor
   */
  [[nodiscard]] const std::array<int, 3> getDomainId() const { return _domainId; }

  /**
   * Returns the number of subdomains in the simulation.
   * @return number of subdomains
   */
  [[nodiscard]] int getNumberOfSubdomains() const;

  /**
   * Checks if the provided coordinates are located in the local domain.
   * @param coordinates: The coordinates in question.
   * @return true if the coordinates lie inside the local domain, false otherwise.
   */
  [[nodiscard]] bool isInsideLocalDomain(const std::array<double, 3> &coordinates) const override;

  /**
   * Calculates and returns the extent of the subdomain with inde subdomainIndex.
   * @param subdomainIndex: The index of the subdomain for which to calculate the extent.
   * @return extent of the subdomain with index subdomainIndex.
   */
  [[nodiscard]] std::array<int, 6> getExtentOfSubdomain(const int subdomainIndex) const;

  /**
   * Exchanges halo particles with all neighbors of the provided AutoPasContainer.
   * @param autoPasContainer: The container, where the halo particles originate from.
   */
  void exchangeHaloParticles(SharedAutoPasContainer &autoPasContainer);

  /**
   * Exchanges migrating particles with all neighbors of the provided AutoPasContainer.
   * @param autoPasContainer: The container, where the migrating particles originate from.
   * @param emigrants: The emigrating particles send to neighbors.
   */
  void exchangeMigratingParticles(SharedAutoPasContainer &autoPasContainer, std::vector<ParticleType> &emigrants);

 private:
  /**
   * The number of neighbors of a rectangular domain.
   * This number does not include diagonal neighbors.
   */
  static constexpr int _neighborCount = 6;

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
   * Contains the planar communicators along each dimension where the current process is a part of.
   * A planar communicator contains all processes with the same coordinate in a single dimension.
   */
  std::array<autopas::AutoPas_MPI_Comm, 3> _planarCommunicators;

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
   * The indices of the local domain's neighbors.
   * These correspond to the ranks of the processors which own the neigbour domain.  */
  std::array<int, 6> _neighborDomainIndices;

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
   * Defines which load balancer will be used.
   */
  LoadBalancerOption _loadBalancerOption;

#if defined(AUTOPAS_ENABLE_ALLLBL)
  /**
   * The ALL load balancer used for diffuse load balancing
   * We cannot use a shared pointer here, because whenn the load balancer is deleted, it calls MPI_Comm_free after
   * we call MPI_Finalize().
   */
  std::unique_ptr<ALL::ALL<double, double>> _allLoadBalancer;
#endif

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
   * Initializes the neighbor ids.
   * This needs to be called after initializeLocalDomain.
   */
  void initializeNeighborIds();

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
   * Received data which has been sent by a specifig neighbor of this domain.
   * @param neighbor The neighbor where the data originates from.
   * @param dataBuffer The buffer where the received data will be stored.
   */
  void receiveDataFromNeighbor(const int &neighbor, std::vector<char> &dataBuffer);

  /**
   * Sends data to a specific neighbor of this domain.
   * @param sendBuffer The buffer which will be sent to the neighbor.
   * @param neighbor The neighbor to which the data will be sent.
   */
  void sendDataToNeighbor(std::vector<char> sendBuffer, const int &neighbor);

  /**
   * Sends and also receives particles to and from the left and right neighbors.
   * @param particlesToLeft: Particles which get send to the left neighbor.
   * @param particlesToRight: Particles which get send to the right neighbor.
   * @param leftNeighbor: The left neighbor's index / rank.
   * @param rightNeighbor: The right neighbor's index / rank.
   * @param receivedParticles: Container for the particles received from either neighbor.
   */
  void sendAndReceiveParticlesLeftAndRight(std::vector<ParticleType> &particlesToLeft,
                                           std::vector<ParticleType> &particlesToRight, const int &leftNeighbor,
                                           const int &rightNeighbor, std::vector<ParticleType> &receivedParticles);

  /**
   * Waits for all send requests to be finished.
   */
  void waitForSendRequests();

  /**
   * Collects the halo particles for the left neighbor.
   * Halo particle positions will be wrapped around the global domain boundary if necessary.
   * @param autoPasContainer: The autopas container which owns the potential halo particles.
   * @param direction: The direction along which the neighbor is located.
   * @param haloParticles: The container the identified halo particles are gathered in to.
   */
  void collectHaloParticlesForLeftNeighbor(SharedAutoPasContainer &autoPasContainer, const size_t &direction,
                                           std::vector<ParticleType> &haloParticles);

  /**
   * Collects the halo particles for the right neighbor.
   * Halo particle positions will be wrapped around the global domain boundary if necessary.
   * @param autoPasContainer: The autopas container which owns the potential halo particles.
   * @param direction: The direction along which the neighbor is located.
   * @param haloParticles: The container the identified halo particles are gathered in to.
   */
  void collectHaloParticlesForRightNeighbor(SharedAutoPasContainer &autoPasContainer, const size_t &direction,
                                            std::vector<ParticleType> &haloParticles);

  /**
   * Categorizes the provided particles as particles for the left or the right neighbor and adds them to the respective
   * output vector. Particle positions will be wrapped around the global domain boundary if necessary.
   * @param particles: The particles which need to be categorized.
   * @param direction: The index of the dimension along which the left and right neighbor lie.
   * @param leftNeighborParticles: Contains the particles for the left neighbor after function execution.
   * @param rightNeighborParticles: Contains the particles for the right neighbor after function execution.
   * @param uncategorizedParticles: Contains particles which could neither be assigned to the left nor the right
   * neighbor.
   */
  void categorizeParticlesIntoLeftAndRightNeighbor(const std::vector<ParticleType> &particles, const size_t &direction,
                                                   std::vector<ParticleType> &leftNeighborParticles,
                                                   std::vector<ParticleType> &rightNeighborParticles,
                                                   std::vector<ParticleType> &uncategorizedParticles);

  /**
   * Balances the subdomains of the grid decomposition using the inverted pressure balancing algorithm.
   * @param work: The work performed by the process owning this sudomain.
   */
  void balanceWithInvertedPressureLoadBalancer(const double &work);

#if defined(AUTOPAS_ENABLE_ALLLBL)
  /**
   * Balances the subdomains of the grid decomposition using the ALL load balancer.
   * @param work: The work performed by the process owning this sudomain.
   */
  void balanceWithAllLoadBalancer(const double &work);
#endif
};
