/**
 * @file ClusterTower.h
 * @author humig
 * @date 27.07.19
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/ParticleDeletedObserver.h"
#include "autopas/containers/verletClusterLists/Cluster.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/markParticleAsDeleted.h"

namespace autopas::internal {

/**
 * This class represents one tower for clusters in the VerletClusterLists container.
 *
 * A ClusterTower contains multiple clusters that are stacked on top (z-direction) of each other. It saves all particles
 * in a FullParticleCell, provides methods to generate and work on the clusters contained, and handles the generation of
 * dummy particles to make sure that each cluster is full.
 *
 * Only the last cluster of the ClusterTower is filled up with dummy particles, since all others are guaranteed to
 * already be full.
 *
 * To use this container:
 *  1. Use addParticle() to add all particles you want.
 *  2. Then call generateClusters(). This copies the last particle as often as it is necessary to fill up the last
 * cluster. (maximum clusterSize-1 times).
 *  3. Generate your neighbor lists somehow.
 *  4. Call setDummyValues() to replace the copies of the last particle made in generateClusters() with
 * dummies.
 *
 * If you want to add more particles after calling generateClusters(), definitely make sure to call clear() before
 * calling addParticle() again, since doing otherwise will mess up the dummy particles and actual particles will likely
 * get lost.
 *
 * @tparam Particle_T
 * @tparam clusterSize
 */
template <class Particle_T>
class ClusterTower : public FullParticleCell<Particle_T> {
 public:
  /**
   * Type that holds or refers to the actual particles.
   */
  using StorageType = typename FullParticleCell<Particle_T>::StorageType;

  /**
   * Dummy constructor.
   */
  ClusterTower() : _clusterSize(0) {}

  /**
   * Constructor
   * @param clusterSize of all clusters in this tower.
   */
  explicit ClusterTower(size_t clusterSize) : _clusterSize(clusterSize) {}

  CellType getParticleCellTypeAsEnum() override { return CellType::ClusterTower; }

  /**
   * Clears all particles from the tower and resets it to be ready for new particles.
   */
  void clear() override {
    _clusters.clear();
    FullParticleCell<Particle_T>::clear();
    _firstOwnedCluster = _clusters.end();
    _firstTailHaloCluster = _clusters.end();
    _numDummyParticles = 0;
  }

  /**
   * Generates the clusters for the particles in this cluster tower.
   *
   * Copies the last particle as often as necessary to fill up the last cluster. This makes sure that iteration over
   * clusters already works after this, while the bounding box of the last cluster is also not messed up by dummy
   * particles. This is necessary for rebuilding the neighbor lists.
   *
   * @return Returns the number of clusters in the tower.
   */
  size_t generateClusters() {
    if (getNumActualParticles() == 0) {
      _clusters.resize(0, Cluster<Particle_T>(nullptr, _clusterSize));
      _firstOwnedCluster = _clusters.end();
      _firstTailHaloCluster = _clusters.end();
      return 0;
    }
    this->sortByDim(2);

    // if the number of particles is divisible by the cluster size this is 0
    const auto sizeLastCluster = this->_particles.size() % _clusterSize;
    _numDummyParticles = sizeLastCluster == 0 ? 0 : _clusterSize - sizeLastCluster;

    // fill the last cluster with dummy copies of the last particle
    auto lastParticle = this->_particles[this->_particles.size() - 1];
    markParticleAsDeleted(lastParticle);
    for (size_t i = 0; i < _numDummyParticles; i++) {
      this->addParticle(lastParticle);
    }

    // Mark start of the different clusters by adding pointers to _particles
    const size_t numClusters = this->_particles.size() / _clusterSize;
    _clusters.resize(numClusters, Cluster<Particle_T>(nullptr, _clusterSize));
    bool foundFirstOwnedCluster = false;
    bool foundFirstTailHaloCluster = false;

    const auto isAnyOfClusterOwned = [&](const auto &cluster) {
      for (size_t i = 0; i < _clusterSize; ++i) {
        if (cluster[i].isOwned()) {
          return true;
        }
      }
      return false;
    };

    for (size_t index = 0; index < numClusters; index++) {
      _clusters[index].reset(&(this->_particles[_clusterSize * index]));

      const bool clusterContainsOwnedParticles =
          not foundFirstTailHaloCluster and isAnyOfClusterOwned(_clusters[index]);

      // identify the range of owned clusters
      if (not foundFirstOwnedCluster and clusterContainsOwnedParticles) {
        _firstOwnedCluster = _clusters.begin() + index;
        foundFirstOwnedCluster = true;
      }
      if (not foundFirstTailHaloCluster and foundFirstOwnedCluster and not clusterContainsOwnedParticles) {
        _firstTailHaloCluster = _clusters.begin() + index;
        foundFirstTailHaloCluster = true;
      }
    }
    // for safety, make sure that worst case (e.g. no particles) the iterators point to reasonable positions
    if (not foundFirstOwnedCluster) {
      _firstOwnedCluster = _clusters.end();
    }
    if (not foundFirstTailHaloCluster) {
      _firstTailHaloCluster = _clusters.end();
    }

    return getNumClusters();
  }

  /**
   * Replaces the copies of the last particle made in generateClusters() with dummies.
   * Dummy particles have ID std::numeric_limits<size_t>::max().
   *
   * @param dummyStartX The x-coordinate for all dummies.
   * @param dummyDistZ The distance in z-direction that all generated dummies will have from each other.
   */
  void setDummyValues(double dummyStartX, double dummyDistZ) {
    auto &lastCluster = getCluster(getNumClusters() - 1);
    for (size_t index = 1; index <= _numDummyParticles; index++) {
      auto &dummy = lastCluster[_clusterSize - index];
      dummy = lastCluster[0];  // use first Particle in last cluster as dummy particle!
      dummy.setOwnershipState(OwnershipState::dummy);
      dummy.setR({
        static_cast<Particle_T::ParticleSoAFloatPrecision>(dummyStartX),
        0,
        static_cast<Particle_T::ParticleSoAFloatPrecision>(dummyDistZ) * static_cast<Particle_T::ParticleSoAFloatPrecision>(index)
      });
      dummy.setID(std::numeric_limits<size_t>::max());
    }
  }

  /**
   * More or less inverse operation of setDummyValues().
   * It sets the positions of the dummy particles to the position of the last actual particle in the tower.
   */
  void setDummyParticlesToLastActualParticle() {
    if (_numDummyParticles > 0) {
      auto &lastCluster = getCluster(getNumClusters() - 1);
      auto lastActualParticle = lastCluster[_clusterSize - _numDummyParticles - 1];
      for (size_t index = 1; index <= _numDummyParticles; index++) {
        lastCluster[_clusterSize - index] = lastActualParticle;
      }
    }
  }

  /**
   * Loads the particles into the SoA stored in this tower and generates the SoAView for each cluster.
   * @tparam Functor The type of the functor to use.
   * @param functor The functor to use for loading the particles into the SoA.
   */
  template <class Functor>
  void loadSoA(Functor *functor) {
    functor->SoALoader(*this, this->_particleSoABuffer, 0, /*skipSoAResize*/ false);
    for (size_t index = 0; index < getNumClusters(); index++) {
      auto &cluster = getCluster(index);
      cluster.setSoAView({&(this->_particleSoABuffer), index * _clusterSize, (index + 1) * _clusterSize});
    }
  }

  /**
   * Extracts the SoA into the particles/clusters.
   * @tparam Functor The type of the functor to use.
   * @param functor The functor to use for extracting the SoA into the particles/clusters.
   */
  template <class Functor>
  void extractSoA(Functor *functor) {
    functor->SoAExtractor(*this, this->_particleSoABuffer, 0);
  }

  /**
   * Returns a rvalue reference to a std::vector containing all particles of this tower that are not dummies.
   *
   * clear() has to called afterwards!
   * @return
   */
  std::vector<Particle_T> &&collectAllActualParticles() {
    if (not this->_particles.empty()) {
      // Workaround to remove requirement of default constructible particles.
      // This function will always only shrink the array, particles are not actually inserted.
      this->_particles.resize(getNumActualParticles(), this->_particles[0]);
    }
    return std::move(this->_particles);
  }

  /**
   * Collect all particles that should not be in this tower based on their current position.
   * The particles are deleted from the tower.
   *
   * @note also deletes dummies.
   *
   * @param boxMin
   * @param boxMax
   * @return Vector of particles that should be stored somewhere else.
   */
  std::vector<Particle_T> collectOutOfBoundsParticles(const std::array<double, 3> &boxMin,
                                                      const std::array<double, 3> &boxMax) {
    // make sure to get rid of all dummies
    deleteDummyParticles();
    // move all particles that are not in the tower to the back of the storage
    const auto firstOutOfBoundsParticleIter =
        std::partition(this->_particles.begin(), this->_particles.end(),
                       [&](const auto &p) { return utils::inBox(p.getR(), boxMin, boxMax); });

    // copy out of bounds particles
    std::vector<Particle_T> outOfBoundsParticles(firstOutOfBoundsParticleIter, this->_particles.end());
    // shrink the particle storage so all out of bounds particles are cut away
    this->_particles.resize(std::distance(this->_particles.begin(), firstOutOfBoundsParticleIter),
                            *this->_particles.begin());
    return outOfBoundsParticles;
  }

  /**
   * Returns the number of dummy particles in the tower that were inserted in the tower to fill the last cluster.
   * @return.
   */
  [[nodiscard]] size_t getNumTailDummyParticles() const { return _numDummyParticles; }

  /**
   * Get the number of all particles stored in this tower (owned, halo and dummy).
   * @return number of particles stored in this tower (owned, halo and dummy).
   */
  [[nodiscard]] size_t size() const override { return this->_particles.size(); }

  /**
   * Get the number of all particles saved in the tower without tailing dummies that are used to fill up clusters
   * (owned + halo + dummies excluding tailing dummies).
   * @return
   */
  [[nodiscard]] unsigned long getNumActualParticles() const {
    return this->_particles.size() - getNumTailDummyParticles();
  }

  /**
   * Returns the number of clusters in the tower.
   * @return the number of clusters in the tower.
   */
  [[nodiscard]] size_t getNumClusters() const { return _clusters.size(); }

  /**
   * Returns a reference to the std::vector holding the clusters of this container.
   * @return a reference to the std::vector holding the clusters of this container.
   */
  [[nodiscard]] auto &getClusters() { return _clusters; }

  /**
   * Returns an iterator to the first cluster that contains at least one owned particle.
   * @return
   */
  [[nodiscard]] const typename std::vector<autopas::internal::Cluster<Particle_T>>::iterator &getFirstOwnedCluster()
      const {
    return _firstOwnedCluster;
  }

  /**
   * Returns the index of the first cluster that contains at least one owned particle.
   * @return
   */
  [[nodiscard]] size_t getFirstOwnedClusterIndex() const {
    return std::distance(_clusters.cbegin(), decltype(_clusters.cbegin()){_firstOwnedCluster});
  }

  /**
   * Returns an iterator to the first particle after the owned clusters, that contains no owned particles anymore.
   * @return
   */
  [[nodiscard]] const typename std::vector<autopas::internal::Cluster<Particle_T>>::iterator &getFirstTailHaloCluster()
      const {
    return _firstTailHaloCluster;
  }

  /**
   * Returns the index of the first particle after the owned clusters, that contains no owned particles anymore.
   * @return
   */
  [[nodiscard]] size_t getFirstTailHaloClusterIndex() const {
    return std::distance(_clusters.cbegin(), decltype(_clusters.cbegin()){_firstTailHaloCluster});
  }

  /**
   * Returns the cluster at position index.
   * @param index The index of the cluster to return.
   * @return the cluster at position index.
   */
  [[nodiscard]] auto &getCluster(size_t index) { return _clusters[index]; }

  /**
   * @copydoc getCluster(size_t)
   */
  [[nodiscard]] auto &getCluster(size_t index) const { return _clusters[index]; }

  [[nodiscard]] bool isEmpty() const override { return getNumActualParticles() == 0; }

  void deleteDummyParticles() override {
    // call super function to do the actual delete
    FullParticleCell<Particle_T>::deleteDummyParticles();
    _numDummyParticles = 0;
  }

  void deleteByIndex(size_t index) override {
    /// @note The implementation of this function prevents a regionIterator to make sorted assumptions of particles
    /// inside a cell! supporting this would mean that the deleted particle should be swapped to the end of the valid
    /// particles. See also https://github.com/AutoPas/AutoPas/issues/435

    // swap particle that should be deleted to end of actual particles.
    std::swap(this->_particles[index], this->_particles[getNumActualParticles() - 1]);
    if (getNumTailDummyParticles() != 0) {
      // swap particle that should be deleted (now at end of actual particles) with last dummy particle.
      std::swap(this->_particles[getNumActualParticles() - 1], this->_particles[this->_particles.size() - 1]);
    }

    this->_particles.pop_back();

    if (_particleDeletionObserver) {
      _particleDeletionObserver->notifyParticleDeleted();
    }
  }

  void setCellLength(std::array<double, 3> &) override {
    autopas::utils::ExceptionHandler::exception("ClusterTower::setCellLength(): Not supported!");
  }

  [[nodiscard]] std::array<double, 3> getCellLength() const override {
    autopas::utils::ExceptionHandler::exception("ClusterTower::getCellLength(): Not supported!");
    return {0, 0, 0};
  }

  /**
   * Set cluster size.
   * @param clusterSize
   */
  void setClusterSize(size_t clusterSize) { _clusterSize = clusterSize; }

  /**
   * Get cluster size
   * @return
   */
  [[nodiscard]] size_t getClusterSize() const { return _clusterSize; }

  /**
   * Set the ParticleDeletionObserver, which is called, when a particle is deleted.
   * @param observer
   */
  void setParticleDeletionObserver(internal::ParticleDeletedObserver *observer) {
    _particleDeletionObserver = observer;
  };

  /**
   * Get reference to internal particle vector.
   * @return
   */
  StorageType &particleVector() { return this->_particles; }

 private:
  /**
   * The number of particles in a full cluster.
   */
  size_t _clusterSize;

  /**
   * The clusters that are contained in this tower.
   */
  std::vector<Cluster<Particle_T>> _clusters;

  /**
   * Iterator pointing to the first cluster that contains at least one owned particle.
   * This cluster might consist of owned, halo, and, if the number of particles in the tower is smaller
   * than the cluster size, dummy particles.
   */
  typename decltype(_clusters)::iterator _firstOwnedCluster{_clusters.end()};
  /**
   * Iterator pointing to the first cluster after the owned clusters that contains no owned particles anymore.
   * This cluster might consist of halo and dummy particles.
   */
  typename decltype(_clusters)::iterator _firstTailHaloCluster{_clusters.end()};

  /**
   * The number of dummy particles in this tower.
   */
  size_t _numDummyParticles{};

  internal::ParticleDeletedObserver *_particleDeletionObserver{nullptr};
};
}  // namespace autopas::internal
