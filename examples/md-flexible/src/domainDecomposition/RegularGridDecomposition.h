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
#include "DomainTools.h"
#include "src/ParticleSerializationTools.h"

/**
 * This class can be used as a domain decomposition which divides the domain in equal sized rectangular subdomains.
 * The number of subdomains is equal to the number of MPI processes available.
 */
template <class ParticleClass>
class RegularGridDecomposition final : public DomainDecomposition {
 public:
  /**
   * Constructor.
   * @param globalBoxMin: The minimum coordinates of the global domain.
   * @param globalBoxMax: The maximum coordinates of the global domain.
   * @param subdivideDimension: Decides if a dimension will be subdivided.
   * @param cutoffWidth: The cutoff width for halo particles.
   * @param skinWidth: The skin width of an autopas container domain.
   * @param boundaryConditions: An array of boundary conditions in the x, y, and z directions.
   */
  RegularGridDecomposition(const std::array<double, 3> &globalBoxMin, const std::array<double, 3> &globalBoxMax,
                           const std::array<bool, 3> &subdivideDimension, double cutoffWidth, double skinWidth,
                           const std::array<options::BoundaryTypeOption, 3> &boundaryConditions)
    : _cutoffWidth(cutoffWidth),
    _skinWidth(skinWidth),
    _globalBoxMin(globalBoxMin),
    _globalBoxMax(globalBoxMax),
    _boundaryType(boundaryConditions) {
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &_subdomainCount);

  int rank;
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);

#if defined(AUTOPAS_INCLUDE_MPI)
  _mpiCommunicationNeeded = true;
#else
  _mpiCommunicationNeeded = false;
#endif

  if (_subdomainCount == 1) {
    _mpiCommunicationNeeded = false;
  }

  DomainTools::generateDecomposition(_subdomainCount, subdivideDimension, _decomposition);

  initializeMPICommunicator();

  initializeLocalDomain();

  initializeLocalBox();

  initializeNeighbourIds();
}

  /**
   * Destructor.
   */
  virtual ~RegularGridDecomposition() {}

  /**
   * Type for the AutoPas container
   */
  using SharedAutoPasContainer = std::shared_ptr<autopas::AutoPas<ParticleClass>>;

  /**
   * Used to update the domain to the current topology.
   * Currently does nothing
   */
  void update() override  { updateLocalBox(); }

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
  int getNumberOfSubdomains() const {
    return std::accumulate(_decomposition.begin(), _decomposition.end(), 1, std::multiplies<int>());
  }

  /**
   * Checks if the provided coordinates are located in the local domain.
   * @param coordinates: The coordinates in question.
   * @return true if the coordinates lie inside the local domain, false otherwise.
   */
  bool isInsideLocalDomain(const std::array<double, 3> &coordinates) const override {
    return DomainTools::isInsideDomain(coordinates, _localBoxMin, _localBoxMax);
  }

  /**
   * Calculates and returns the extent of the subdomain with inde subdomainIndex.
   * @param subdomainIndex: The index of the subdomain for which to calculate the extent.
   * @return extent of the subdomain with index subdomainIndex.
   */
  std::array<int, 6> getExtentOfSubdomain(const int subdomainIndex) const {
    return DomainTools::getExtentOfSubdomain(subdomainIndex, _decomposition);
  }

  /**
   * Exchanges halo particles with all neighbours of the provided AutoPasContainer.
   * @param autoPasContainer: The container, where the halo particles originate from.
   */
  void exchangeHaloParticles(SharedAutoPasContainer &autoPasContainer) {
    for (int i = 0; i < _dimensionCount; ++i) {
      std::vector<ParticleClass> haloParticles{};
      for (int j = i; j < _dimensionCount; ++j) {
        const size_t dimensionIndex = j % _dimensionCount;

        // completely bypass Halo particle exchange in this dimension if boundaries in this direction are not periodic
        // *and* if both local boundaries are the global boundaries in this dimension
        if (_boundaryType[dimensionIndex] != options::BoundaryTypeOption::periodic and
            _localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex] and
            _localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex])
          continue;

        std::vector<ParticleClass> particlesForLeftNeighbour{};
        std::vector<ParticleClass> particlesForRightNeighbour{};

        if (j == i) {
          collectHaloParticlesForLeftNeighbour(autoPasContainer, dimensionIndex, particlesForLeftNeighbour);
          collectHaloParticlesForRightNeighbour(autoPasContainer, dimensionIndex, particlesForRightNeighbour);
        }

        double leftHaloMin = _localBoxMin[dimensionIndex] - _skinWidth;
        double leftHaloMax = _localBoxMin[dimensionIndex] + _cutoffWidth + _skinWidth;
        double rightHaloMin = _localBoxMax[dimensionIndex] - _cutoffWidth - _skinWidth;
        double rightHaloMax = _localBoxMax[dimensionIndex] + _skinWidth;

        for (const auto &particle : haloParticles) {
          auto position = particle.getR();

          // check left boundary is a periodic (global) boundary
          if (_boundaryType[dimensionIndex] == options::BoundaryTypeOption::periodic and
              _localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex]) {
            if (position[dimensionIndex] >= leftHaloMin and position[dimensionIndex] < leftHaloMax) {
              particlesForLeftNeighbour.push_back(particle);

              // Apply boundary condition
              if (_localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex]) {
                position[dimensionIndex] =
                    position[dimensionIndex] + (_globalBoxMax[dimensionIndex] - _globalBoxMin[dimensionIndex]);
                particlesForLeftNeighbour.back().setR(position);
              }
            }
          }

          // check right boundary is a periodic (global) boundary
          if (_boundaryType[dimensionIndex] == options::BoundaryTypeOption::periodic and
              _localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
            if (position[dimensionIndex] >= rightHaloMin and position[dimensionIndex] < rightHaloMax) {
              particlesForRightNeighbour.push_back(particle);

              // Apply boundary condition
              if (_localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
                position[dimensionIndex] =
                    position[dimensionIndex] - (_globalBoxMax[dimensionIndex] - _globalBoxMin[dimensionIndex]);
                particlesForRightNeighbour.back().setR(position);
              }
            }
          }
        }

        // See documentation for _neighbourDomainIndices to explain the indexing
        int leftNeighbour = _neighbourDomainIndices[(dimensionIndex * 2) % _neighbourCount];
        int rightNeighbour = _neighbourDomainIndices[(dimensionIndex * 2 + 1) % _neighbourCount];
        sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbour, particlesForRightNeighbour, leftNeighbour,
                                            rightNeighbour, haloParticles);
      }

      for (const auto &particle : haloParticles) {
        autoPasContainer->addHaloParticle(particle);
      }
    }
  }

  /**
   * Exchanges migrating particles with all neighbours of the provided AutoPasContainer.
   * @param autoPasContainer: The container, where the migrating particles originate from.
   */
  void exchangeMigratingParticles(SharedAutoPasContainer &autoPasContainer) {
    auto emigrants = autoPasContainer->updateContainer();

    const std::array<double, _dimensionCount> globalBoxLength =
        autopas::utils::ArrayMath::sub(_globalBoxMax, _globalBoxMin);

    for (int i = 0; i < _dimensionCount; ++i) {
      for (int j = i; j < _dimensionCount; ++j) {
        const size_t dimensionIndex = j % _dimensionCount;
        // completely bypass particle exchange in this dimension if boundaries in this direction are not periodic *and*
        // if both local boundaries are the global boundaries in this dimension
        if (_boundaryType[dimensionIndex] != options::BoundaryTypeOption::periodic and
            _localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex] and
            _localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex])
          continue;

        std::vector<ParticleClass> immigrants, remainingEmigrants;
        std::vector<ParticleClass> particlesForLeftNeighbour;
        std::vector<ParticleClass> particlesForRightNeighbour;

        // See documentation for _neighbourDomainIndices to explain the indexing
        int leftNeighbour = _neighbourDomainIndices[(dimensionIndex * 2) % _neighbourCount];
        int rightNeighbour = _neighbourDomainIndices[(dimensionIndex * 2 + 1) % _neighbourCount];

        for (const auto &particle : emigrants) {
          auto position = particle.getR();
          // check left boundary is a periodic (global) boundary
          if (_boundaryType[dimensionIndex] == options::BoundaryTypeOption::periodic and
              _localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex]) {
            if (position[dimensionIndex] < _localBoxMin[dimensionIndex]) {
              particlesForLeftNeighbour.push_back(particle);

              // Apply boundary condition
              if (_localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex]) {
                position[dimensionIndex] =
                    std::min(std::nextafter(_globalBoxMax[dimensionIndex], _globalBoxMin[dimensionIndex]),
                             position[dimensionIndex] + globalBoxLength[dimensionIndex]);
                particlesForLeftNeighbour.back().setR(position);
              }
            }
          }
          // check right boundary is a periodic (global) boundary
          if (_boundaryType[dimensionIndex] == options::BoundaryTypeOption::periodic and
              _localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
            if (position[dimensionIndex] >= _localBoxMax[dimensionIndex]) {
              particlesForRightNeighbour.push_back(particle);

              // Apply boundary condition
              if (_localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
                position[dimensionIndex] =
                    std::max(_globalBoxMin[dimensionIndex], position[dimensionIndex] - globalBoxLength[dimensionIndex]);
                particlesForRightNeighbour.back().setR(position);
              }
            } else {
              remainingEmigrants.push_back(particle);
            }
          }
        }
        emigrants = remainingEmigrants;

        sendAndReceiveParticlesLeftAndRight(particlesForLeftNeighbour, particlesForRightNeighbour, leftNeighbour,
                                            rightNeighbour, immigrants);

        for (const auto &particle : immigrants) {
          if (isInsideLocalDomain(particle.getR())) {
            autoPasContainer->addParticle(particle);
          } else {
            emigrants.push_back(particle);
          }
        }

        immigrants.clear();
      }
    }
  }

  /**
   * Reflects particles within a reflective skin along the inside of a boundary.
   * @param autoPasContainer: The container, where the migrating particles originate from.
   */
  void reflectParticlesAtBoundaries(SharedAutoPasContainer &autoPasContainer) {
    std::array<double, _dimensionCount> reflSkinMin{}, reflSkinMax{};

    for (int dimensionIndex = 0; dimensionIndex < _dimensionCount; ++dimensionIndex) {
      // skip if boundary is not reflective
      if (_boundaryType[dimensionIndex] != options::BoundaryTypeOption::reflective) continue;

      auto reflect = [&](bool isUpper) {
        for (auto p = autoPasContainer->getRegionIterator(reflSkinMin, reflSkinMax, autopas::IteratorBehavior::owned);
             p.isValid(); ++p) {
          auto vel = p->getV();
          // reverse velocity in dimension if towards boundary
          if ((vel[dimensionIndex] < 0) xor isUpper) {
            vel[dimensionIndex] *= -1;
          }
          p->setV(vel);
        }
      };

      // apply if we are at a global boundary on lower end of the dimension
      if (_localBoxMin[dimensionIndex] == _globalBoxMin[dimensionIndex]) {
        reflSkinMin = _globalBoxMin;
        reflSkinMax = _globalBoxMax;
        reflSkinMax[dimensionIndex] = _globalBoxMin[dimensionIndex] + autoPasContainer->getVerletSkin() / 2;

        reflect(false);
      }
      // apply if we are at a global boundary on upper end of the dimension
      if (_localBoxMax[dimensionIndex] == _globalBoxMax[dimensionIndex]) {
        reflSkinMin = _globalBoxMin;
        reflSkinMax = _globalBoxMax;
        reflSkinMin[dimensionIndex] = _globalBoxMax[dimensionIndex] - autoPasContainer->getVerletSkin() / 2;

        reflect(true);
      }
    }
  }

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
   * Boundary condition types.
   */
  std::array<options::BoundaryTypeOption, 3> _boundaryType;

  /**
   * A temporary buffer used for MPI send requests.
   */
  std::vector<autopas::AutoPas_MPI_Request> _sendRequests;

  /**
   * A temporary buffer for data which is sent by MPI_Send.
   */
  std::vector<std::vector<char>> _sendBuffers;

  /**
   * Initializes the MPI communicator.
   * This needs to be called before initializeLocalDomain.
   */
  void initializeMPICommunicator() {
    std::vector<int> periods(3, 1);
    autopas::AutoPas_MPI_Cart_create(AUTOPAS_MPI_COMM_WORLD, 3, _decomposition.data(), periods.data(), true,
                                     &_communicator);
    autopas::AutoPas_MPI_Comm_rank(_communicator, &_domainIndex);
  }

  /**
   * Initializes the local domain.
   * This needs to be called before initializeLocalBox.
   */
  void initializeLocalDomain() {
    _domainId = {0, 0, 0};
    autopas::AutoPas_MPI_Comm_rank(_communicator, &_domainIndex);

    std::vector<int> periods(3, 1);
    autopas::AutoPas_MPI_Cart_get(_communicator, 3, _decomposition.data(), periods.data(), _domainId.data());
  }

  /**
   * Initializes the local domain coordinates.
   * This needs to be called after initializeLocalDomain and initialzieGlobalDomain.
   */
  void initializeLocalBox() { updateLocalBox(); }

  /**
   * Initializes the neighbour ids.
   * This needs to be called after initializeLocalDomain.
   */
  void initializeNeighbourIds() {
    for (int i = 0; i < 3; ++i) {
      auto neighbourIndex = i * 2;
      auto preceedingNeighbourId = _domainId;
      preceedingNeighbourId[i] = (--preceedingNeighbourId[i] + _decomposition[i]) % _decomposition[i];
      _neighbourDomainIndices[neighbourIndex] = DomainTools::convertIdToIndex(preceedingNeighbourId, _decomposition);

      ++neighbourIndex;
      auto succeedingNeighbourId = _domainId;
      succeedingNeighbourId[i] = (++succeedingNeighbourId[i] + _decomposition[i]) % _decomposition[i];
      _neighbourDomainIndices[neighbourIndex] = DomainTools::convertIdToIndex(succeedingNeighbourId, _decomposition);
    }
  }

  /**
   * Updates the local box.
   */
  void updateLocalBox() {
    for (int i = 0; i < 3; ++i) {
      double localBoxWidth = (_globalBoxMax[i] - _globalBoxMin[i]) / static_cast<double>(_decomposition[i]);

      _localBoxMin[i] = _domainId[i] * localBoxWidth + _globalBoxMin[i];
      _localBoxMax[i] = (_domainId[i] + 1) * localBoxWidth + _globalBoxMin[i];

      if (_domainId[i] == 0) {
        _localBoxMin[i] = _globalBoxMin[i];
      } else if (_domainId[i] == _decomposition[i] - 1) {
        _localBoxMax[i] = _globalBoxMax[i];
      }
    }
  }

  /**
   * Sends particles of type ParticleClass to a receiver.
   * @param particles The particles to be sent to the receiver.
   * @param receiver The recipient of the particels.
   */
  void sendParticles(const std::vector<ParticleClass> &particles, const int &receiver) {
    std::vector<char> buffer;

    for (const auto &particle : particles) {
      ParticleSerializationTools::serializeParticle(particle, buffer);
    }

    sendDataToNeighbour(buffer, receiver);
  }

  /**
   * Received particles sent by a sender.
   * @param receivedParticles The container where the received particles will be stored.
   * @param source The sender id/rank.
   */
  void receiveParticles(std::vector<ParticleClass> &receivedParticles, const int &source) {
    std::vector<char> receiveBuffer;

    receiveDataFromNeighbour(source, receiveBuffer);

    if (!receiveBuffer.empty()) {
      ParticleSerializationTools::deserializeParticles(receiveBuffer, receivedParticles);
    }
  }

  /**
   * Received data which has been sent by a specifig neighbour of this domain.
   * @param neighbour The neighbour where the data originates from.
   * @param dataBuffer The buffer where the received data will be stored.
   */
  void receiveDataFromNeighbour(const int &neighbour, std::vector<char> &receiveBuffer) {
    autopas::AutoPas_MPI_Status status;
    autopas::AutoPas_MPI_Probe(neighbour, 0, _communicator, &status);

    int receiveBufferSize;
    autopas::AutoPas_MPI_Get_count(&status, AUTOPAS_MPI_CHAR, &receiveBufferSize);
    receiveBuffer.resize(receiveBufferSize);

    autopas::AutoPas_MPI_Recv(receiveBuffer.data(), receiveBufferSize, AUTOPAS_MPI_CHAR, neighbour, 0, _communicator,
                              AUTOPAS_MPI_STATUS_IGNORE);
  }

  /**
   * Sends data to a specific neighbour of this domain.
   * @param sendBuffer The buffer which will be sent to the neighbour.
   * @param neighbour The neighbour to which the data will be sent.
   */
  void sendDataToNeighbour(std::vector<char> sendBuffer, const int &neighbour) {
    _sendBuffers.push_back(sendBuffer);

    autopas::AutoPas_MPI_Request sendRequest;
    _sendRequests.push_back(sendRequest);

    autopas::AutoPas_MPI_Isend(_sendBuffers.back().data(), _sendBuffers.back().size(), AUTOPAS_MPI_CHAR, neighbour, 0,
                               _communicator, &_sendRequests.back());
  }

  /**
   * Sends and also receives particles to and from the left and right neighbours.
   * @param particlesToLeft: Particles which get send to the left neighbour.
   * @param particlesToRight: Particles which get send to the right neighbor.
   * @param leftNeighbour: The left neighbour's index / rank.
   * @param rightNeighbour: The right neighbour's index / rank.
   * @param receivedParticles: Container for the particles received from either neighbour.
   */
  void sendAndReceiveParticlesLeftAndRight(std::vector<ParticleClass> &particlesToLeft,
                                           std::vector<ParticleClass> &particlesToRight, const int &leftNeighbour,
                                           const int &rightNeighbour, std::vector<ParticleClass> &receivedParticles) {
    if (_mpiCommunicationNeeded and leftNeighbour != _domainIndex) {
      sendParticles(particlesToLeft, leftNeighbour);
      sendParticles(particlesToRight, rightNeighbour);

      receiveParticles(receivedParticles, leftNeighbour);
      receiveParticles(receivedParticles, rightNeighbour);

      waitForSendRequests();
    } else {
      receivedParticles.insert(receivedParticles.end(), particlesToLeft.begin(), particlesToLeft.end());
      receivedParticles.insert(receivedParticles.end(), particlesToRight.begin(), particlesToRight.end());
    }
  }

  /**
   * Waits for all send requests to be finished.
   */
  void waitForSendRequests() {
    std::vector<autopas::AutoPas_MPI_Status> sendStates;
    sendStates.resize(_sendRequests.size());
    autopas::AutoPas_MPI_Waitall(_sendRequests.size(), _sendRequests.data(), sendStates.data());
    _sendRequests.clear();
    _sendBuffers.clear();
  }

  /**
   * Collects the halo particles for the left neighbour.
   * Halo particle positions will be wrapped around the global domain boundary if necessary.
   * @param autoPasContainer: The autopas container which owns the potential halo particles.
   * @param direction: The direction along which the neighbour is located.
   * @param haloParticles: The container the identified halo particles are gathered in to.
   */
  void collectHaloParticlesForLeftNeighbour(SharedAutoPasContainer &autoPasContainer, const size_t &direction,
                                            std::vector<ParticleClass> &haloParticles) {
    std::array<double, _dimensionCount> boxMin, boxMax;

    // Calculate halo box for left neighbour
    for (int i = 0; i < _dimensionCount; ++i) {
      boxMin[i] = _localBoxMin[i] - _skinWidth;
      boxMax[i] = _skinWidth + _localBoxMax[i];
    }
    boxMax[direction] = _localBoxMin[direction] + _cutoffWidth + _skinWidth;

    // Collect the halo particles for the left neighbour
    for (auto particle = autoPasContainer->getRegionIterator(boxMin, boxMax, autopas::IteratorBehavior::owned);
         particle.isValid(); ++particle) {
      std::array<double, _dimensionCount> position = particle->getR();
      haloParticles.push_back(*particle);

      // Apply boundary condition
      if (_localBoxMin[direction] == _globalBoxMin[direction]) {
        position[direction] = position[direction] + (_globalBoxMax[direction] - _globalBoxMin[direction]);
        haloParticles.back().setR(position);
      }
    }
  }

  /**
   * Collects the halo particles for the right neighbour.
   * Halo particle positions will be wrapped around the global domain boundary if necessary.
   * @param autoPasContainer: The autopas container which owns the potential halo particles.
   * @param direction: The direction along which the neighbour is located.
   * @param haloParticles: The container the identified halo particles are gathered in to.
   */
  void collectHaloParticlesForRightNeighbour(SharedAutoPasContainer &autoPasContainer, const size_t &direction,
                                             std::vector<ParticleClass> &haloParticles) {
    std::array<double, _dimensionCount> boxMin, boxMax;

    // Calculate left halo box of right neighbour
    for (int i = 0; i < _dimensionCount; ++i) {
      boxMin[i] = _localBoxMin[i] - _skinWidth;
      boxMax[i] = _localBoxMax[i] + _skinWidth;
    }
    boxMin[direction] = _localBoxMax[direction] - _cutoffWidth - _skinWidth;

    // Collect the halo particles for the right neighbour
    for (auto particle = autoPasContainer->getRegionIterator(boxMin, boxMax, autopas::IteratorBehavior::owned);
         particle.isValid(); ++particle) {
      std::array<double, _dimensionCount> position = particle->getR();
      haloParticles.push_back(*particle);

      // Apply boundary condition
      if (_localBoxMax[direction] == _globalBoxMax[direction]) {
        position[direction] = position[direction] - (_globalBoxMax[direction] - _globalBoxMin[direction]);
        haloParticles.back().setR(position);
      }
    }
  }

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
  void categorizeParticlesIntoLeftAndRightNeighbour(const std::vector<ParticleClass> &particles, const size_t &direction,
                                                    std::vector<ParticleClass> &leftNeighbourParticles,
                                                    std::vector<ParticleClass> &rightNeighbourParticles,
                                                    std::vector<ParticleClass> &uncategorizedParticles) {
    for (const auto &particle : particles) {
      std::array<double, _dimensionCount> position = particle.getR();
      if (position[direction] < _localBoxMin[direction]) {
        leftNeighbourParticles.push_back(particle);

        // Apply boundary condition
        if (_localBoxMin[direction] == _globalBoxMin[direction]) {
          position[direction] = std::min(std::nextafter(_globalBoxMax[direction], _globalBoxMin[direction]),
                                         position[direction] + _globalBoxMax[direction] - _globalBoxMin[direction]);
          leftNeighbourParticles.back().setR(position);
        }
      } else if (position[direction] >= _localBoxMax[direction]) {
        rightNeighbourParticles.push_back(particle);

        // Apply boundary condition
        if (_localBoxMax[direction] == _globalBoxMax[direction]) {
          position[direction] = std::max(_globalBoxMin[direction],
                                         position[direction] - (_globalBoxMax[direction] - _globalBoxMin[direction]));
          rightNeighbourParticles.back().setR(position);
        }
      } else {
        uncategorizedParticles.push_back(particle);
      }
    }
  }
};
