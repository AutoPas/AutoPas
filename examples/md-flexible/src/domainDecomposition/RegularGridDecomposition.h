/**
 * @file RegularGridDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#if defined(MD_FLEXIBLE_ENABLE_ALLLBL)
#include <ALL.hpp>
#endif
#include <memory>

#include "DomainDecomposition.h"
#include "autopas/utils/WrapMPI.h"
#include "src/TypeDefinitions.h"
#include "src/configuration/MDFlexConfig.h"
#include "src/options/BoundaryTypeOption.h"

namespace {
/**
 * Sixth Root of Two precomputed here, as it is used a lot in reflecting boundaries.
 */
const double sixthRootOfTwo = std::pow(2., 1. / 6.);
}  // namespace

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
  explicit RegularGridDecomposition(const MDFlexConfig &configuration);

  /**
   * Destructor.
   */
  ~RegularGridDecomposition() override;

  /**
   * Used to update the domain to the current topology.
   * Handles the diffuse load balancing by resizing the domains according to their work done.
   * Currently the metric used as work is the timings created by the simulation timer. Although they are of type 'long'
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
  [[nodiscard]] const std::array<double, 3> &getGlobalBoxMin() const override { return _globalBoxMin; }

  /**
   * Returns the maximum coordinates of global domain.
   * @return top right back corner of the global domain.
   */
  [[nodiscard]] const std::array<double, 3> &getGlobalBoxMax() const override { return _globalBoxMax; }

  /**
   * Returns the minimum coordinates of local domain.
   * @return bottom left front corner of the local domain.
   */
  [[nodiscard]] const std::array<double, 3> &getLocalBoxMin() const override { return _localBoxMin; }

  /**
   * Returns the maximum coordinates of local domain.
   * @return top right back corner of the local domain.
   */
  [[nodiscard]] const std::array<double, 3> &getLocalBoxMax() const override { return _localBoxMax; }

  /**
   * Returns the number of domains in each dimension
   * @return vector containing the number of subdomains along each dimension
   */
  [[nodiscard]] const std::array<int, 3> &getDecomposition() const { return _decomposition; }

  /**
   * Returns the numnber of subdomains (=ranks) in the decomposition.
   * @return numner of subdomains in the decomposition.
   */
  [[nodiscard]] int getNumberOfSubdomains() const { return _subdomainCount; }

  /**
   * Returns the current processes domain id.
   * @return domain id of the current processor
   */
  [[nodiscard]] const std::array<int, 3> &getDomainId() const { return _domainId; }

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
  [[nodiscard]] std::array<int, 6> getExtentOfSubdomain(int subdomainIndex) const;

  /**
   * Exchanges halo particles with all neighbors of the provided AutoPasContainer.
   * @param autoPasContainer: The container, where the halo particles originate from.
   */
  void exchangeHaloParticles(AutoPasType &autoPasContainer);

  /**
   * Exchanges migrating particles with all neighbors of the provided AutoPasContainer.
   * @param autoPasContainer: The container, where the migrating particles originate from.
   * @param emigrants: The emigrating particles to send to neighbors.
   */
  void exchangeMigratingParticles(AutoPasType &autoPasContainer, std::vector<ParticleType> &emigrants);

  /**
   * Reflects particles within a reflective skin along the inside of a boundary.
   *
   * Particle reflection occurs by interacting the particle with particle mirrored onto the other side of the boundary.
   * Iteraction occurs using the AoS variant of the chosen functor. Particle reflection only occurs if the particle
   * would experience a repulsive effect (i.e. is within the 6th root of sigma from the border).
   *
   * @param autoPasContainer: The container, where the migrating particles originate from.
   * @param PPL: Particle Properties Library (needed to get particle's sigma)
   */
  void reflectParticlesAtBoundaries(AutoPasType &autoPasContainer, ParticlePropertiesLibraryType &PPL);

  /**
   * Getter for the communicator that encompasses the whole decomposition.
   * @return
   */
  autopas::AutoPas_MPI_Comm getCommunicator() const;

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
   * Defines which load balancer will be used.
   */
  LoadBalancerOption _loadBalancerOption;

  /**
   * Stores the domain cutoff width.
   */
  double _cutoffWidth;

  /**
   * Stores the domain skin width.
   */
  double _skinWidth;
  /**
   * Stores the rebuild Frequency.
   */
  unsigned int _rebuildFrequency;
  /**
   * The greatest distance from a reflective boundary at which a particle might experience reflection.
   */
  double _maxReflectiveSkin;
  /**
   * The minimum coordinates of the global domain.
   */
  std::array<double, _dimensionCount> _globalBoxMin;

  /**
   * The maximum coordinates of the global domain.
   */
  std::array<double, _dimensionCount> _globalBoxMax;

  /**
   * Boundary condition types of all dimensions.
   */
  std::array<options::BoundaryTypeOption, _dimensionCount> _boundaryType;

  /**
   * Indicates if MPI is enabled and if it will be used.
   * In the case that MPI is enabled, but only one process is being used, this variable will be false.
   */
  bool _mpiCommunicationNeeded;

  /**
   * The number of subdomains (=ranks) in this decomposition.
   */
  int _subdomainCount{};

  /**
   * The decomposition computed depending on the number of subdomains.
   * Number of ranks per dimension.
   */
  std::array<int, 3> _decomposition{};

  /**
   * Indicator to MPI to view all communication dimensions as periodic.
   * @note For usage in MPI functions, the const needs to be casted away.
   */
  std::vector<int> _periods{_dimensionCount, 1};

  /**
   * The MPI communicator containing all processes which own a subdomain in this decomposition.
   */
  autopas::AutoPas_MPI_Comm _communicator{};

  /**
   * The index of the current processor's domain in _communicator.
   * This also is the rank of the current processor.
   */
  int _domainIndex{};

  /**
   * The 3D ID of the current processor's domain.
   */
  std::array<int, 3> _domainId{};

  /**
   * Contains the planar communicators along each dimension where the current process is a part of.
   * A planar communicator contains all processes with the same coordinate in a single dimension.
   */
  std::array<autopas::AutoPas_MPI_Comm, _dimensionCount> _planarCommunicators{};

  /**
   * The indices of the local domain's neighbors.
   * These correspond to the ranks of the processors which own the neighbor domain.
   */
  std::array<int, _neighborCount> _neighborDomainIndices{};

  /**
   * The minimum coordinates of the local domain.
   */
  std::array<double, _dimensionCount> _localBoxMin{};

  /**
   * The maximum coordinates of the local domain.
   */
  std::array<double, _dimensionCount> _localBoxMax{};

  /**
   * Buffer for particles to be sent to the left neighbor.
   * @note This buffer exists solely to avoid reallocations.
   */
  std::vector<ParticleType> _particlesForLeftNeighbor{};
  /**
   * Buffer for particles to be sent to the right neighbor.
   * @note This buffer exists solely to avoid reallocations.
   */
  std::vector<ParticleType> _particlesForRightNeighbor{};

  /**
   * Buffer for communicating particles.
   * @note This buffer exists solely to avoid reallocations.
   */
  std::vector<ParticleType> _receivedParticlesBuffer{};

  /**
   * Buffer for halo particles that have to be inserted locally.
   */
  std::vector<ParticleType> _haloParticles{};

#if defined(MD_FLEXIBLE_ENABLE_ALLLBL)
  /**
   * The ALL load balancer used for diffuse load balancing
   * We cannot use a shared pointer here, because whenn the load balancer is deleted, it calls MPI_Comm_free after
   * we call MPI_Finalize().
   */
  std::unique_ptr<ALL::ALL<double, double>> _allLoadBalancer;
#endif

  /**
   * Initializes the MPI communicator.
   * This needs to be called before initializeLocalDomain.
   */
  void initializeMPICommunicator();

  /**
   * Initialize _planarCommunicators and domainID.
   * This needs to be called before initializeLocalBox.
   */
  void initializeLocalDomain();

  /**
   * Initializes the local domain coordinates, aka _localBoxMin and _localBoxMax.
   * This needs to be called after initializeLocalDomain.
   */
  void initializeLocalBox();

  /**
   * Initializes _neighborDomainIndices.
   * This needs to be called after initializeLocalDomain.
   */
  void initializeNeighborIndices();

  /**
   * Sends and also receives particles to and from the left and right neighbours.
   * @param particlesToLeft: Particles which get send to the left neighbour.
   * @param particlesToRight: Particles which get send to the right neighbor.
   * @param receivedParticlesBuffer In-out parameter where the received particles will be placed in.
   * @param leftNeighbor: The left neighbor's index / rank.
   * @param rightNeighbor: The right neighbor's index / rank.
   */
  void sendAndReceiveParticlesLeftAndRight(const std::vector<ParticleType> &particlesToLeft,
                                           const std::vector<ParticleType> &particlesToRight,
                                           std::vector<ParticleType> &receivedParticlesBuffer, int leftNeighbor,
                                           int rightNeighbor);

  /**
   * Helper function to reduce code duplication between collectHaloParticlesForLeftNeighbor() and
   * collectHaloParticlesForRightNeighbor().
   * @param autoPasContainer
   * @param direction Along which dimension should we collect.
   * @param boxMin lower left front corner of the region that is collected.
   * @param boxMax upper right rear  corner of the region that is collected.
   * @param atGlobalBoundary Indicator if we are at a global boundary and periodics are necessary.
   * @param wrapAroundDistance Distance (incl sign for direction) that the particle needs to be relocated.
   * @param haloParticlesBuffer in-out buffer for haloParticles#
   */
  void collectHaloParticlesAux(AutoPasType &autoPasContainer, size_t direction,
                               const std::array<double, _dimensionCount> &boxMin,
                               const std::array<double, _dimensionCount> &boxMax, bool atGlobalBoundary,
                               double wrapAroundDistance, std::vector<ParticleType> &haloParticlesBuffer);
  /**
   * Collects the halo particles for the left neighbour.
   * Halo particle positions will be wrapped around the global domain boundary if necessary.
   * @param autoPasContainer: The autopas container which owns the potential halo particles.
   * @param direction: The direction along which the neighbor is located.
   * @param particlesForLeftNeighborBuffer in-out buffer for haloParticles.
   */
  void collectHaloParticlesForLeftNeighbor(AutoPasType &autoPasContainer, size_t direction,
                                           std::vector<ParticleType> &particlesForLeftNeighborBuffer);

  /**
   * Collects the halo particles for the right neighbor.
   * Halo particle positions will be wrapped around the global domain boundary if necessary.
   * @param autoPasContainer: The autopas container which owns the potential halo particles.
   * @param direction: The direction along which the neighbor is located.
   * @param particlesForRightNeighborBuffer in-out buffer for haloParticles.
   */
  void collectHaloParticlesForRightNeighbor(AutoPasType &autoPasContainer, size_t direction,
                                            std::vector<ParticleType> &particlesForRightNeighborBuffer);

  /**
   * Categorizes the provided particles as particles for the left or the right neighbor and adds them to the respective
   * output vector. Particle positions will be wrapped around the global domain boundary if necessary.
   * @param particles: The particles which need to be categorized.
   * @param direction: The index of the dimension along which the left and right neighbor lie.
   * @return a tuple consisting of:
   *    leftNeighborParticles: Contains the particles for the left neighbor after function execution.
   *    rightNeighborParticles: Contains the particles for the right neighbor after function execution.
   *    uncategorizedParticles: Contains particles which could neither be assigned to the left nor the right neighbor.
   */
  std::tuple<std::vector<ParticleType>, std::vector<ParticleType>, std::vector<ParticleType>>
  categorizeParticlesIntoLeftAndRightNeighbor(const std::vector<ParticleType> &particles, size_t direction);

  /**
   * Balances the subdomains of the grid decomposition using the inverted pressure balancing algorithm.
   * @param work: The work performed by the process owning this sudomain.
   */
  void balanceWithInvertedPressureLoadBalancer(double work);

  /**
   * Balances the subdomains of the grid decomposition using the ALL load balancer.
   * @param work: The work performed by the process owning this sudomain.
   *
   * @note If md-flexible is compiled without ALL this function throws an exception.
   */
  void balanceWithAllLoadBalancer(const double &work);
};
