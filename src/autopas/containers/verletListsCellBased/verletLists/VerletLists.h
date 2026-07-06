/**
 * @file VerletLists.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include "VerletListHelpers.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLListIterationTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"

namespace autopas {

/**
 * Verlet Lists container.
 * The VerletLists class uses neighborhood lists to calculate pairwise
 * interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * Cells are created using a cell size of at least cutoff + skin radius.
 * @tparam Particle_T
 */
template <class Particle_T>
class VerletLists : public VerletListsLinkedBase<Particle_T> {
  /**
   * Type of the Particle.
   */
  using ParticleType = Particle_T;
  /**
   * Type of the ParticleCell used by the underlying linked cells.
   */
  using ParticleCellType = FullParticleCell<Particle_T>;

 public:
  /**
   * Enum that specifies how the verlet lists should be build
   */
  enum BuildVerletListType {
    /**
     * Build it using AoS
     */
    VerletAoS,
    /**
     * Build it using AoS
     */
    VerletSoA,
  };

  /**
   * Constructor of the VerletLists class.
   * The neighbor lists are build using a search radius of cutoff + skin.
   * @param boxMin The lower corner of the domain.
   * @param boxMax The upper corner of the domain.
   * @param cutoff The cutoff radius of the interaction.
   * @param skin The skin radius per timestep.
   * @param buildVerletListType Specifies how the verlet list should be build, see BuildVerletListType
   * @param cellSizeFactor cell size factor ralative to cutoff
   */
  VerletLists(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double cutoff,
              const double skin, const BuildVerletListType buildVerletListType = BuildVerletListType::VerletSoA,
              const double cellSizeFactor = 1.0)
      : VerletListsLinkedBase<Particle_T>(boxMin, boxMax, cutoff, skin, cellSizeFactor),
        _buildVerletListType(buildVerletListType) {}

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::verletLists; }

  void computeInteractions(TraversalInterface *traversal) override {
    // Check if traversal is allowed for this container and give it the data it needs.
    auto *verletTraversalInterface = dynamic_cast<VLTraversalInterface<ParticleCellType> *>(traversal);
    if (verletTraversalInterface) {
      verletTraversalInterface->setCellsAndNeighborLists(this->_linkedCells.getCells(), _neighborList,
                                                         _particleToIndex);
    } else {
      utils::ExceptionHandler::exception("trying to use a traversal of wrong type in VerletLists::computeInteractions");
    }

    traversal->initTraversal();
    traversal->traverseParticles();
    traversal->endTraversal();
  }

  /**
   * Returns the flat CRS neighbor list.
   * Offsets and indices are valid after the most recent rebuildNeighborLists() call.
   * @return the CRS neighbor list
   */
  const typename VerletListHelpers<Particle_T>::NeighborListCRS &getNeighborList() const { return _neighborList; }

  /**
   * Returns the particle-pointer-to-SoA-index map.
   * Built once per rebuild alongside the CRS list.
   * @return the particle index map
   */
  const std::unordered_map<const Particle_T *, size_t> &getParticleIndex() const { return _particleToIndex; }

  /**
   * Rebuilds the neighbor lists, marks them valid and resets the internal counter.
   * Builds the CRS directly — no separate AoS→SoA conversion pass needed.
   * @note This function will be called in computeInteractions()!
   * @param traversal
   */
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    this->_verletBuiltNewton3 = traversal->getUseNewton3();
    updateNeighborLists(traversal->getUseNewton3());
    // the neighbor list is now valid
    this->_neighborListIsValid.store(true, std::memory_order_relaxed);
  }

 private:
  /**
   * Builds the particle-pointer -> SoA-index map and returns the particle count N.
   * Iterates over owned + halo + dummy particles in the same order as the SoA loader, so index i here matches row i in
   * the SoA buffer.
   * @return Total number of particles (owned + halo + dummy).
   */
  size_t buildParticleIndex() {
    _particleToIndex.clear();
    // Reserve generously to avoid rehashing (particles + halo estimate).
    _particleToIndex.reserve(this->_linkedCells.size() * 2);
    size_t idx = 0;
    for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter, ++idx) {
      _particleToIndex[&(*iter)] = idx;
    }
    return idx;
  }

  /**
   * Converts the temporary per-particle ragged lists produced by the generator
   * functor into the flat CRS representation stored in _neighborList.
   *
   * Two passes:
   *  1. Prefix-sum over list sizes -> fills _neighborList.offsets.
   *  2. Copy each inner vector into the matching slice of _neighborList.indices.
   *
   * @param tempLists  Ragged list: tempLists[i] holds SoA indices of particle i's neighbors.
   */
  void buildCRS(const std::vector<std::vector<size_t>> &tempLists) {
    const size_t numParticles = tempLists.size();
    // Step 1: compute offsets via prefix sum
    _neighborList.offsets.resize(numParticles + 1);
    _neighborList.offsets[0] = 0;
    for (size_t i = 0; i < numParticles; ++i) {
      _neighborList.offsets[i + 1] = _neighborList.offsets[i] + tempLists[i].size();
    }

    // Step 2: copy neighbor indices into the flat array
    const size_t totalNeighbors = _neighborList.offsets[numParticles];
    _neighborList.indices.resize(totalNeighbors);
    for (size_t i = 0; i < numParticles; ++i) {
      std::copy(tempLists[i].begin(), tempLists[i].end(),
                _neighborList.indices.begin() + static_cast<std::ptrdiff_t>(_neighborList.offsets[i]));
    }

    AutoPasLog(DEBUG, "VerletLists::buildCRS: {} particles, {} total neighbors, avg list size {:.2f}", numParticles,
               totalNeighbors,
               numParticles > 0 ? static_cast<double>(totalNeighbors) / static_cast<double>(numParticles) : 0.0);
  }

  /**
   * Rebuilds _particleToIndex and _neighborList from scratch.
   *
   * Sequence:
   *  1. buildParticleIndex()         — assign a stable index to every particle.
   *  2. Allocate temporary ragged lists (one per particle).
   *  3. Run LCC08Traversal with VerletListGeneratorFunctor — populates tempLists.
   *  4. buildCRS()                   — prefix-sum + copy into flat CRS.
   *
   * @param useNewton3  Whether the force traversal will use Newton's third law.
   */
  virtual void updateNeighborLists(bool useNewton3) {
    const size_t N = buildParticleIndex();

    // Temporary ragged lists written by the generator functor.
    std::vector<std::vector<size_t>> tempLists(N);

    typename VerletListHelpers<Particle_T>::VerletListGeneratorFunctor f(tempLists, _particleToIndex,
                                                                         this->getCutoff() + this->getVerletSkin());

    /// @todo autotune build traversal
    DataLayoutOption dataLayout;
    if (_buildVerletListType == BuildVerletListType::VerletAoS) {
      dataLayout = DataLayoutOption::aos;
    } else if (_buildVerletListType == BuildVerletListType::VerletSoA) {
      dataLayout = DataLayoutOption::soa;
    } else {
      utils::ExceptionHandler::exception("VerletLists::updateNeighborLists(): unsupported BuildVerletListType: {}",
                                         static_cast<int>(_buildVerletListType));
    }
    auto traversal =
        LCC08Traversal<ParticleCellType, typename VerletListHelpers<Particle_T>::VerletListGeneratorFunctor>(
            this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), f, this->getInteractionLength(),
            this->_linkedCells.getCellBlock().getCellLength(), dataLayout, useNewton3);
    this->_linkedCells.computeInteractions(&traversal);

    buildCRS(tempLists);
  }

  /**
   * Mapping of every particle pointer to its dense SoA index.
   * Built once per rebuild in buildParticleIndex() before the traversal.
   * Shared with VerletListGeneratorFunctor during list construction and with the traversal for the AoS force
   * computations.
   */
  std::unordered_map<const Particle_T *, size_t> _particleToIndex;

  /**
   * Flat CRS neighbor list.
   * Both AoS and SoA traversal paths read from this neighbor structure.
   */
  typename VerletListHelpers<Particle_T>::NeighborListCRS _neighborList;

  /**
   * Specifies which data layout is used when building the neighbor lists.
   */
  BuildVerletListType _buildVerletListType;
};

}  // namespace autopas
