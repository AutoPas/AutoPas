/**
 * @file VerletLists.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include "VerletListHelpers.h"
#include "autopas/baseFunctors/InteractionListGeneratorFunctor.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLListIterationTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * Verlet Lists container.
 * The VerletLists class uses neighborhood lists to calculate pairwise or triwise interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of the interaction.
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
   * Enum that specifies how the verlet lists should be built
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
   * The neighbor lists are built using a search radius of cutoff + skin.
   * @param boxMin The lower corner of the domain.
   * @param boxMax The upper corner of the domain.
   * @param cutoff The cutoff radius of the interaction.
   * @param skin The skin radius per timestep.
   * @param buildVerletListType Specifies how the verlet list should be built, see BuildVerletListType
   * @param cellSizeFactor cell size factor relative to cutoff
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
      verletTraversalInterface->setCellsAndNeighborLists(this->_linkedCells.getCells(), _neighborList, _particleToIndex,
                                                         _indexToParticle, _neighborPairsList);
    } else {
      utils::ExceptionHandler::exception(
          "VerletLists::computeInteractions(): Trying to use a traversal of wrong type.");
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
  const VerletListHelpers<Particle_T>::NeighborListCRS &getNeighborList() const { return _neighborList; }

  /**
   * Returns the flat CRS neighbor pairs list.
   * Offsets and pairs are valid after the most recent updateNeighborPairsList() call.
   * @return the CRS neighbor pairs list
   */
  const VerletListHelpers<Particle_T>::NeighborPairsListCRS &getNeighborPairsList() const { return _neighborPairsList; }

  /**
   * Build the pair neighbor list if necessary without fully rebuilding the other neighbor lists.
   * @param traversal
   */
  void prepareForTraversal(TraversalInterface *traversal) override {
    if (traversal->getTraversalType() == TraversalOption::vl_pair_list_iteration) {
      if (not _pairListIsValid) {
        this->updateNeighborPairsList(this->_verletBuiltNewton3);
      }
    }
  }

  /**
   * Returns the particle-pointer-to-SoA-index map.
   * Built once per rebuild alongside the CRS list.
   * @return the particle index map
   */
  const std::unordered_map<const Particle_T *, size_t> &getParticleIndex() const { return _particleToIndex; }

  /**
   * Rebuilds the neighbor lists, marks them valid and resets the internal counter.
   * Builds the CRS directly — no separate AoS→SoA conversion pass needed.
   * @note This function will be called in computeInteractions()
   * @param traversal
   */
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    this->_verletBuiltNewton3 = traversal->getUseNewton3();
    const auto buildWithN3 = traversal->getUseNewton3();

    // Check for triwise traversals
    switch (traversal->getTraversalType()) {
      // Standard pairwise traversal
      case TraversalOption::vl_list_iteration: {
        this->updateNeighborLists(buildWithN3);
        if (buildWithN3) {
          this->modifyNeighborListsForTriwiseTraversal(TraversalOption::vl_list_iteration);
        }
        break;
      }
      case TraversalOption::vl_list_intersection: {
        this->updateNeighborLists(buildWithN3);
        this->modifyNeighborListsForTriwiseTraversal(TraversalOption::vl_list_intersection);
        break;
      }
      case TraversalOption::vl_pair_list_iteration: {
        // build 3Body Verlet lists through VLIteration traversal
        this->updateNeighborPairsList(buildWithN3);
        break;
      }
      // Default builds normal neighbor lists.
      default: {
        this->updateNeighborLists(buildWithN3);
      }
    }

    // the neighbor list is now valid
    this->_neighborListIsValid.store(true, std::memory_order_relaxed);
    _pairListIsValid = false;
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
    _particleToIndex.reserve(this->_linkedCells.size());
    _indexToParticle.clear();
    _indexToParticle.reserve(this->_linkedCells.size());
    size_t idx = 0;
    for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter, ++idx) {
      _particleToIndex[&(*iter)] = idx;
      _indexToParticle.push_back(&*iter);
    }
    return idx;
  }

  /**
   * Rebuilds _particleToIndex and _neighborList from scratch.
   *
   * Dispatches to the single-pass path when only one thread is available (lower overhead, no atomics) and to the
   * two-pass lock-free path when multiple threads are active (eliminates false sharing and malloc contention).
   *
   * @param useNewton3  Whether the force traversal will use Newton's third law.
   */
  virtual void updateNeighborLists(bool useNewton3) {
    const size_t N = buildParticleIndex();
    const double interactionLength = this->getInteractionLength();

    DataLayoutOption dataLayout;
    if (_buildVerletListType == BuildVerletListType::VerletAoS) {
      dataLayout = DataLayoutOption::aos;
    } else if (_buildVerletListType == BuildVerletListType::VerletSoA) {
      dataLayout = DataLayoutOption::soa;
    } else {
      utils::ExceptionHandler::exception("VerletLists::updateNeighborLists(): unsupported BuildVerletListType: {}",
                                         static_cast<int>(_buildVerletListType));
    }

    if (autopas_get_max_threads() == 1) {
      updateNeighborListsSingleThread(N, interactionLength, dataLayout, useNewton3);
    } else {
      updateNeighborListsMultiThread(N, interactionLength, dataLayout, useNewton3);
    }
  }

  /**
   * Single-threaded rebuild: One traversal with VerletListGeneratorFunctor writing into per-particle
   * std::vector<size_t>, followed by a serial prefix-sum + copy into the flat CRS.
   */
  void updateNeighborListsSingleThread(size_t N, double interactionLength, DataLayoutOption dataLayout,
                                       bool useNewton3) {
    std::vector<std::vector<size_t>> tempLists(N);

    typename VerletListHelpers<Particle_T>::CRSNeighborListPolicy policy(tempLists, _particleToIndex);

    InteractionListGeneratorFunctor<Particle_T, typename VerletListHelpers<Particle_T>::CRSNeighborListPolicy> f(
        policy, interactionLength, useNewton3);

    auto traversal = LCC08Traversal<
        ParticleCellType,
        InteractionListGeneratorFunctor<Particle_T, typename VerletListHelpers<Particle_T>::CRSNeighborListPolicy>>(
        this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), f, this->getInteractionLength(),
        this->_linkedCells.getCellBlock().getCellLength(), dataLayout, useNewton3);
    this->_linkedCells.computeInteractions(&traversal);

    // Prefix-sum + copy:
    _neighborList.offsets.resize(N + 1);
    _neighborList.offsets[0] = 0;
    for (size_t i = 0; i < N; ++i) {
      _neighborList.offsets[i + 1] = _neighborList.offsets[i] + tempLists[i].size();
    }
    const size_t totalNeighbors = _neighborList.offsets[N];
    _neighborList.indices.resize(totalNeighbors);
    for (size_t i = 0; i < N; ++i) {
      std::copy(tempLists[i].begin(), tempLists[i].end(),
                _neighborList.indices.begin() + static_cast<std::ptrdiff_t>(_neighborList.offsets[i]));
    }

    AutoPasLog(DEBUG, "VerletLists::updateNeighborLists (1T): {} particles, {} neighbors, avg {:.2f}", N,
               totalNeighbors, N > 0 ? static_cast<double>(totalNeighbors) / static_cast<double>(N) : 0.0);
  }

  /**
   * Multithreaded rebuild:
   *
   * Pass 1 (parallel, VerletListCounterFunctor):
   *   Count neighbors per particle into cache-line-padded atomics.
   *   No heap allocation, no false sharing between adjacent particles.
   *
   * Prefix sum (serial, O(N)):
   *   Compute CRS offsets; allocate the flat indices array.
   *
   * Pass 2 (parallel, VerletListFillerFunctor):
   *   Write neighbor indices directly into the pre-allocated CRS slice
   *   via per-particle atomic fetch-add fill cursors (also padded).
   *
   * Optimal when autopas_get_max_threads() > 1.
   */
  void updateNeighborListsMultiThread(size_t N, double interactionLength, DataLayoutOption dataLayout,
                                      bool useNewton3) {
    using PaddedAtomic = VerletListHelpers<Particle_T>::VerletListCounterFunctor::PaddedAtomic;

    // Pass 1: Count neighbors per particle
    std::vector<PaddedAtomic> counts(N);
    {
      typename VerletListHelpers<Particle_T>::VerletListCounterFunctor counter(counts, _particleToIndex,
                                                                               interactionLength);
      auto traversal =
          LCC08Traversal<ParticleCellType, typename VerletListHelpers<Particle_T>::VerletListCounterFunctor>(
              this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), counter, this->getInteractionLength(),
              this->_linkedCells.getCellBlock().getCellLength(), dataLayout, useNewton3);
      this->_linkedCells.computeInteractions(&traversal);
    }

    // Prefix sum: counts -> CRS offsets, allocate flat indices
    _neighborList.offsets.resize(N + 1);
    _neighborList.offsets[0] = 0;
    for (size_t i = 0; i < N; ++i) {
      _neighborList.offsets[i + 1] = _neighborList.offsets[i] + counts[i].value.load(std::memory_order_relaxed);
    }
    const size_t totalNeighbors = _neighborList.offsets[N];
    _neighborList.indices.resize(totalNeighbors);

    AutoPasLog(DEBUG, "VerletLists::updateNeighborLists (MT): {} particles, {} neighbors, avg {:.2f}", N,
               totalNeighbors, N > 0 ? static_cast<double>(totalNeighbors) / static_cast<double>(N) : 0.0);

    // Pass 2: Fill CRS indices directly
    // Reuse the PaddedAtomic array as fill cursors, seeding each with offsets[i].
    for (size_t i = 0; i < N; ++i) {
      counts[i].value.store(_neighborList.offsets[i], std::memory_order_relaxed);
    }
    {
      typename VerletListHelpers<Particle_T>::VerletListFillerFunctor filler(_neighborList, counts, _particleToIndex,
                                                                             interactionLength);
      auto traversal =
          LCC08Traversal<ParticleCellType, typename VerletListHelpers<Particle_T>::VerletListFillerFunctor>(
              this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), filler, this->getInteractionLength(),
              this->_linkedCells.getCellBlock().getCellLength(), dataLayout, useNewton3);
      this->_linkedCells.computeInteractions(&traversal);
    }
  }

  /**
   * Modifies neighbor lists for specialized triwise traversals:
   * - vl_list_iteration (Newton3 on): Re-orients (halo -> owned) edges to (owned -> halo) edges so owned
   *   particles own all pairwise connections to their halo neighbors, and clears halo particle neighbor lists.
   * - vl_list_intersection: Generates halo-halo edges between halo particles that share an owned neighbor.
   *
   * @param traversalOption The triwise traversal type.
   */
  void modifyNeighborListsForTriwiseTraversal(TraversalOption traversalOption) {
    using namespace utils::ArrayMath::literals;

    const size_t N = _neighborList.size();
    std::vector<std::vector<size_t>> tempLists(N);
    for (size_t i = 0; i < N; ++i) {
      tempLists[i].assign(_neighborList.begin(i), _neighborList.end(i));
    }

    if (traversalOption == TraversalOption::vl_list_iteration) {
      // Re-orient halo -> owned edges to owned -> halo
      for (size_t i = 0; i < N; ++i) {
        if (not _indexToParticle[i]->isHalo()) {
          continue;
        }

        for (const auto &neighbor : _neighborList.getNeighbors(i)) {
          if (not _indexToParticle[neighbor]->isOwned()) {
            // skip halo neighbors, we only want to find owned neighbors
            continue;
          }
          tempLists[neighbor].push_back(i);
        }
        tempLists[i].clear();  // clear halo particle's list, we don't need it anymore
      }
    }

    if (traversalOption == TraversalOption::vl_list_intersection) {
      // Connect halo-halo pairs sharing an owned particle
      const double interactionLength = this->getInteractionLength();
      const double interactionLengthSquared = interactionLength * interactionLength;
      for (size_t i = 0; i < N; ++i) {
        if (not _indexToParticle[i]->isOwned()) {
          continue;
        }
        std::vector<size_t> haloNeighbors;
        for (const size_t neighborIdx : _neighborList.getNeighbors(i)) {
          if (not _indexToParticle[neighborIdx]->isOwned()) {
            haloNeighbors.push_back(neighborIdx);
          }
        }
        const size_t numHalos = haloNeighbors.size();
        for (size_t j = 0; j < numHalos; ++j) {
          const size_t h1 = haloNeighbors[j];
          const auto &pos1 = _indexToParticle[h1]->getR();
          for (size_t k = j + 1; k < numHalos; ++k) {
            const size_t h2 = haloNeighbors[k];
            const auto &pos2 = _indexToParticle[h2]->getR();
            const auto dist = pos1 - pos2;
            if (utils::ArrayMath::dot(dist, dist) < interactionLengthSquared) {
              tempLists[h1].push_back(h2);
              tempLists[h2].push_back(h1);
            }
          }
        }
      }
    }

    // Rebuild the CRS neighbor list
    AUTOPAS_OPENMP(parallel for schedule(dynamic))
    for (size_t i = 0; i < N; ++i) {
      std::ranges::sort(tempLists[i]);
      auto [first, last] = std::ranges::unique(tempLists[i]);
      tempLists[i].erase(first, last);
    }

    _neighborList.offsets.resize(N + 1);
    _neighborList.offsets[0] = 0;
    for (size_t i = 0; i < N; ++i) {
      _neighborList.offsets[i + 1] = _neighborList.offsets[i] + tempLists[i].size();
    }
    _neighborList.indices.resize(_neighborList.offsets[N]);
    for (size_t i = 0; i < N; ++i) {
      std::copy(tempLists[i].begin(), tempLists[i].end(),
                _neighborList.indices.begin() + static_cast<std::ptrdiff_t>(_neighborList.offsets[i]));
    }
  }

  /**
   * Rebuilds _neighborPairsList from scratch.
   *
   * @param useNewton3 Whether the traversal will use Newton's third law.
   */
  void updateNeighborPairsList(bool useNewton3) {
    updateNeighborLists(useNewton3);
    // this->addSharedHaloNeighbors(useNewton3);
    const size_t N = _neighborList.size();
    const double interactionLength = this->getInteractionLength();

    DataLayoutOption dataLayout = DataLayoutOption::aos;
    if (_buildVerletListType == BuildVerletListType::VerletSoA) {
      // there are no SoA 3-body traversals, so we print out a warning and use AoS Layout instead
      AutoPasLog(WARN, "Pair Verlet Lists can currently only be built with AoS DataLayout, using that instead!");
    }

    if (autopas_get_max_threads() == 1) {
      updateNeighborPairsListSingleThread(N, interactionLength, dataLayout, useNewton3);
    } else {
      updateNeighborPairsListMultiThread(N, interactionLength, dataLayout, useNewton3);
    }
    _pairListIsValid = true;
  }

  /**
   * Single-threaded rebuild of neighbor pairs list.
   */
  void updateNeighborPairsListSingleThread(size_t N, double interactionLength, DataLayoutOption dataLayout,
                                           bool useNewton3) {
    std::vector<std::vector<std::pair<size_t, size_t>>> tempLists(N);
    typename VerletListHelpers<Particle_T>::PairVerletListGeneratorFunctor f(tempLists, _particleToIndex,
                                                                             interactionLength);
    auto traversal = VLListIterationTraversal<ParticleCellType,
                                              typename VerletListHelpers<Particle_T>::PairVerletListGeneratorFunctor>(
        f, dataLayout, useNewton3);
    this->computeInteractions(&traversal);

    _neighborPairsList.offsets.resize(N + 1);
    _neighborPairsList.offsets[0] = 0;
    for (size_t i = 0; i < N; ++i) {
      _neighborPairsList.offsets[i + 1] = _neighborPairsList.offsets[i] + tempLists[i].size();
    }
    const size_t totalPairs = _neighborPairsList.offsets[N];
    _neighborPairsList.pairs.resize(totalPairs);
    for (size_t i = 0; i < N; ++i) {
      std::copy(tempLists[i].begin(), tempLists[i].end(),
                _neighborPairsList.pairs.begin() + static_cast<std::ptrdiff_t>(_neighborPairsList.offsets[i]));
    }

    AutoPasLog(DEBUG, "VerletLists::updateNeighborPairsList (1T): {} particles, {} pairs, avg {:.2f}", N, totalPairs,
               N > 0 ? static_cast<double>(totalPairs) / static_cast<double>(N) : 0.0);
  }

  /**
   * Multi-threaded rebuild of neighbor pairs list.
   */
  void updateNeighborPairsListMultiThread(size_t N, double interactionLength, DataLayoutOption dataLayout,
                                          bool useNewton3) {
    using PaddedAtomic = typename VerletListHelpers<Particle_T>::VerletListCounterFunctor::PaddedAtomic;

    // Pass 1: Count neighbor pairs per particle
    std::vector<PaddedAtomic> counts(N);
    {
      typename VerletListHelpers<Particle_T>::PairVerletListCounterFunctor counter(counts, _particleToIndex,
                                                                                   interactionLength);
      auto traversal = VLListIterationTraversal<ParticleCellType,
                                                typename VerletListHelpers<Particle_T>::PairVerletListCounterFunctor>(
          counter, dataLayout, useNewton3);
      this->computeInteractions(&traversal);
    }

    // Prefix sum: counts -> CRS offsets, allocate flat pairs
    _neighborPairsList.offsets.resize(N + 1);
    _neighborPairsList.offsets[0] = 0;
    for (size_t i = 0; i < N; ++i) {
      _neighborPairsList.offsets[i + 1] =
          _neighborPairsList.offsets[i] + counts[i].value.load(std::memory_order_relaxed);
    }
    const size_t totalPairs = _neighborPairsList.offsets[N];
    _neighborPairsList.pairs.resize(totalPairs);

    AutoPasLog(DEBUG, "VerletLists::updateNeighborPairsList (MT): {} particles, {} pairs, avg {:.2f}", N, totalPairs,
               N > 0 ? static_cast<double>(totalPairs) / static_cast<double>(N) : 0.0);

    // Pass 2: Fill neighbor pairs directly
    for (size_t i = 0; i < N; ++i) {
      counts[i].value.store(_neighborPairsList.offsets[i], std::memory_order_relaxed);
    }
    {
      typename VerletListHelpers<Particle_T>::PairVerletListFillerFunctor filler(_neighborPairsList, counts,
                                                                                 _particleToIndex, interactionLength);
      auto traversal = VLListIterationTraversal<ParticleCellType,
                                                typename VerletListHelpers<Particle_T>::PairVerletListFillerFunctor>(
          filler, dataLayout, useNewton3);
      this->computeInteractions(&traversal);
    }
  }

  /**
   * Mapping of every particle pointer to its dense SoA index.
   * Built once per rebuild in buildParticleIndex() before the traversal.
   * Shared with VerletListGeneratorFunctor during list construction and with the traversal for the AoS force
   * computations.
   */
  std::unordered_map<const Particle_T *, size_t> _particleToIndex;

  /**
   * Flat array mapping SoA index i -> pointer to particle i.
   * Built in initTraversal() in the same cell/particle iteration order as VerletLists::buildParticleIndex().
   * Used during AoS force computations to resolve CRS neighbor indices to particles.
   */
  std::vector<ParticleType *> _indexToParticle;

  /**
   * Flat CRS neighbor list.
   * Both AoS and SoA traversal paths read from this neighbor structure.
   */
  VerletListHelpers<Particle_T>::NeighborListCRS _neighborList;

  /**
   * Flat CRS Neighbor Pairs List. (To find triplets.)
   */
  typename VerletListHelpers<Particle_T>::NeighborPairsListCRS _neighborPairsList;

  /**
   * Shows if the SoA neighbor list is currently valid.
   */
  bool _soaListIsValid{false};

  /**
   * Shows if the pair list for triwise interactions is currently valid.
   */
  bool _pairListIsValid{false};

  /**
   * Specifies which data layout is used when building the neighbor lists.
   */
  BuildVerletListType _buildVerletListType;
};

}  // namespace autopas
