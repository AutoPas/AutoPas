/**
 * @file VerletLists.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include "VerletListHelpers.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticSelectorMacros.h"

namespace autopas {

/**
 * Verlet Lists container.
 * The VerletLists class uses neighborhood lists to calculate pairwise
 * interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * Cells are created using a cell size of at least cutoff + skin radius.
 * @note This container supports the [blackbox mode](@ref md_BlackBoxMode).
 * @note This class does NOT work with RMM cells and is not intended to!
 * @tparam Particle
 * @todo deleting particles should also invalidate the verlet lists - should be
 * implemented somehow
 */
template <class Particle>
class VerletLists
    : public VerletListsLinkedBase<Particle, typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                   typename VerletListHelpers<Particle>::SoAArraysType> {
  typedef VerletListHelpers<Particle> verlet_internal;
  typedef FullParticleCell<Particle> ParticleCell;
  typedef typename VerletListHelpers<Particle>::SoAArraysType SoAArraysType;
  typedef typename VerletListHelpers<Particle>::VerletListParticleCellType LinkedParticleCell;

 private:
  static const std::vector<TraversalOptions>& VLApplicableTraversals() {
    /// @todo: implement some traversals for this
    static const std::vector<TraversalOptions> v{};
    return v;
  }

 public:
  /**
   * Enum that specifies how the verlet lists should be build
   */
  enum BuildVerletListType {
    VerletAoS,  /// build it using AoS
    VerletSoA   /// build it using SoA
  };

  /**
   * Constructor of the VerletLists class.
   * The neighbor lists are build using a search radius of cutoff + skin.
   * The rebuildFrequency should be chosen, s.t. the particles do not move more
   * than a distance of skin/2 between two rebuilds of the lists.
   * @param boxMin The lower corner of the domain
   * @param boxMax The upper corner of the domain
   * @param cutoff The cutoff radius of the interaction
   * @param skin The skin radius
   * @param rebuildFrequency Specifies after how many pair-wise traversals the
   * neighbor lists are to be rebuild. A frequency of 1 means that they are
   * always rebuild, 10 means they are rebuild after 10 traversals.
   * @param buildVerletListType Specifies how the verlet list should be build, see BuildVerletListType
   * @param blackBoxMode Indicates whether the [blackbox mode](@ref md_BlackBoxMode) shall be used.
   */
  VerletLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
              const double skin, const unsigned int rebuildFrequency = 1,
              const BuildVerletListType buildVerletListType = BuildVerletListType::VerletSoA,
              const bool blackBoxMode = false)
      : VerletListsLinkedBase<Particle, LinkedParticleCell, SoAArraysType>(
            boxMin, boxMax, cutoff, skin, rebuildFrequency, allVLApplicableTraversals(), blackBoxMode),
        _soaListIsValid(false),
        _soa(),
        _buildVerletListType(buildVerletListType) {}

  /**
   * Lists all traversal options applicable for the Verlet Lists container.
   * @return Vector of all applicable traversal options.
   */
  static const std::vector<TraversalOptions>& allVLApplicableTraversals() {
    // @FIXME This is a workaround because this container does not yet use traversals like it should
    static const std::vector<TraversalOptions> v{TraversalOptions::dummyTraversal};
    return v;
  }

  ContainerOptions getContainerType() override { return ContainerOptions::verletLists; }

  /**
   * @copydoc LinkedCells::iteratePairwiseAoS
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwiseAoS(ParticleFunctor* f, Traversal* traversal, bool useNewton3 = true) {
    if (this->needsRebuild(useNewton3)) {  // if we need to rebuild the list, we should rebuild it!
      rebuildVerletLists(useNewton3);
    }
    this->iterateVerletListsAoS(f, useNewton3);
    // we iterated, so increase traversal counter
    this->_traversalsSinceLastRebuild++;

    if (this->_blackBoxMode) {
      iterateBoundaryPartsAoS(f, useNewton3);
    }
  }

  /**
   * @copydoc LinkedCells::iteratePairwiseSoA
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwiseSoA(ParticleFunctor* f, Traversal* traversal, bool useNewton3 = true) {
    if (this->needsRebuild(useNewton3)) {
      rebuildVerletLists(useNewton3);
      generateSoAListFromAoSVerletLists();
    } else if (not _soaListIsValid) {
      generateSoAListFromAoSVerletLists();
    }
    iterateVerletListsSoA(f, useNewton3);
    this->_traversalsSinceLastRebuild++;

    if (this->_blackBoxMode) {
      iterateBoundaryPartsSoA(f, useNewton3);
    }
  }

  /**
   * Get the actual neighbour list.
   * @return the neighbour list
   */
  typename verlet_internal::AoS_verletlist_storage_type& getVerletListsAoS() { return _aosNeighborLists; }

  /**
   * Checks whether the neighbor lists are valid.
   * A neighbor list is valid if all pairs of particles whose interaction should
   * be calculated are represented in the neighbor lists.
   * @param useNewton3 specified whether newton 3 should be used
   * @return whether the list is valid
   * @note This check involves pair-wise interaction checks and is thus
   * relatively costly.
   */
  bool checkNeighborListsAreValid(bool useNewton3 = true) {
    // if a particle was added or deleted, ... the list is definitely invalid
    if (not this->_neighborListIsValid) {
      return false;
    }
    // if a particle moved more than skin/2 outside of its cell the list is
    // invalid
    if (this->isContainerUpdateNeeded()) {
      return false;
    }

    // particles can also simply be very close already:
    typename verlet_internal::template VerletListValidityCheckerFunctor<LinkedParticleCell> validityCheckerFunctor(
        _aosNeighborLists, ((this->getCutoff() - this->_skin) * (this->getCutoff() - this->_skin)));

    auto traversal =
        C08Traversal<LinkedParticleCell,
                     typename verlet_internal::template VerletListValidityCheckerFunctor<LinkedParticleCell>, false,
                     true>(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &validityCheckerFunctor);
    this->_linkedCells.iteratePairwiseAoS(&validityCheckerFunctor, &traversal, useNewton3);

    return validityCheckerFunctor.neighborlistsAreValid();
  }

  TraversalSelector<ParticleCell> generateTraversalSelector(std::vector<TraversalOptions> traversalOptions) override {
    //    std::vector<TraversalOptions> allowedAndApplicable;
    //
    //    std::sort(traversalOptions.begin(), traversalOptions.end());
    //    std::set_intersection(this->_applicableTraversals.begin(), this->_applicableTraversals.end(),
    //    traversalOptions.begin(),
    //                          traversalOptions.end(), std::back_inserter(allowedAndApplicable));
    /// @todo dummyTraversal is a workaround because this container does not yet use traversals like it should
    return TraversalSelector<ParticleCell>({0, 0, 0}, {dummyTraversal});
  }

 protected:
  /**
   * Rebuilds the verlet lists, marks them valid and resets the internal counter.
   * @note This function will be called in iteratePairwiseAoS() and iteratePairwiseSoA() appropriately!
   * @param useNewton3
   */
  void rebuildVerletLists(bool useNewton3 = true) {
    this->_verletBuiltNewton3 = useNewton3;
    this->updateVerletListsAoS(useNewton3);
    // the neighbor list is now valid
    this->_neighborListIsValid = true;
    this->_traversalsSinceLastRebuild = 0;
  }

  /**
   * update the verlet lists for AoS usage
   * @param useNewton3 Specifies whether newton 3 should be used or not.
   */
  virtual void updateVerletListsAoS(bool useNewton3) {
    updateIdMapAoS();
    typename verlet_internal::VerletListGeneratorFunctor f(_aosNeighborLists, this->getCutoff());

    /// @todo autotune traversal
    switch (_buildVerletListType) {
      case BuildVerletListType::VerletAoS: {
        if (this->_blackBoxMode) {
          AUTOPAS_WITH_STATIC_BOOL(useNewton3, {
            auto traversal =
                C08Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor, false,
                             c_useNewton3, inner>(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
            this->_linkedCells.iteratePairwiseAoS(&f, &traversal);
          })
        } else {
          AUTOPAS_WITH_STATIC_BOOL(useNewton3, {
            auto traversal =
                C08Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor, false,
                             c_useNewton3>(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
            this->_linkedCells.iteratePairwiseAoS(&f, &traversal);
          })
        }
        break;
      }
      case BuildVerletListType::VerletSoA: {
        if (this->_blackBoxMode) {
          AUTOPAS_WITH_STATIC_BOOL(useNewton3, {
            auto traversal =
                C08Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor, true,
                             c_useNewton3, inner>(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
            this->_linkedCells.iteratePairwiseSoA(&f, &traversal);
          })
        } else {
          AUTOPAS_WITH_STATIC_BOOL(useNewton3, {
            auto traversal =
                C08Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor, true,
                             c_useNewton3>(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
            this->_linkedCells.iteratePairwiseSoA(&f, &traversal);
          })
        }
        break;
      }
      default:
        utils::ExceptionHandler::exception("VerletLists::updateVerletListsAoS(): unsupported BuildVerletListType: {}",
                                           _buildVerletListType);
        break;
    }
    _soaListIsValid = false;
  }

  /**
   * iterate over the verlet lists using the AoS traversal
   * @tparam ParticleFunctor
   * @param f
   * @param useNewton3
   */
  template <class ParticleFunctor>
  void iterateVerletListsAoS(ParticleFunctor* f, const bool useNewton3) {
    /// @todo optimize iterateVerletListsAoS, e.g. by using openmp-capable traversals

#if defined(AUTOPAS_OPENMP)
    if (not useNewton3) {
      size_t buckets = _aosNeighborLists.bucket_count();
      /// @todo find a sensible chunk size
#pragma omp parallel for schedule(dynamic)
      for (size_t b = 0; b < buckets; b++) {
        auto endIter = _aosNeighborLists.end(b);
        for (auto it = _aosNeighborLists.begin(b); it != endIter; ++it) {
          Particle& i = *(it->first);
          for (auto j_ptr : it->second) {
            Particle& j = *j_ptr;
            f->AoSFunctor(i, j, false);
          }
        }
      }
    } else
#endif
    {
      for (auto& list : _aosNeighborLists) {
        Particle& i = *list.first;
        for (auto j_ptr : list.second) {
          Particle& j = *j_ptr;
          f->AoSFunctor(i, j, useNewton3);
        }
      }
    }
  }

  /**
   * Iterate over the boundary dependent parts of the domain using the SoA traversal.
   * @tparam ParticleFunctor
   * @param f
   * @param useNewton3
   */
  template <class ParticleFunctor>
  void iterateBoundaryPartsSoA(ParticleFunctor* f, const bool useNewton3) {
    // load SoA's
    loadBoundarySoA(f);

    if (useNewton3) {
      C08Traversal<ParticleCell, ParticleFunctor, /*useSoA*/ true, /*useNewton3*/ true, BlackBoxTraversalOption::outer>
          traversal(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), f);
      traversal.traverseCellPairs(this->_linkedCells.getCells());
    } else {
      C08Traversal<ParticleCell, ParticleFunctor, /*useSoA*/ true, /*useNewton3*/ false, BlackBoxTraversalOption::outer>
          traversal(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), f);
      traversal.traverseCellPairs(this->_linkedCells.getCells());
    }

    // extract SoA's
    extractBoundarySoA(f);
  }

  /**
   * Iterate over the boundary dependent parts of the domain using the SoA traversal.
   * @tparam ParticleFunctor
   * @param f
   * @param useNewton3
   */
  template <class ParticleFunctor>
  void iterateBoundaryPartsAoS(ParticleFunctor* f, const bool useNewton3) {
    if (useNewton3) {
      C08Traversal<ParticleCell, ParticleFunctor, /*useSoA*/ false, /*useNewton3*/ true, BlackBoxTraversalOption::outer>
          traversal(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), f);
      traversal.traverseCellPairs(this->_linkedCells.getCells());
    } else {
      C08Traversal<ParticleCell, ParticleFunctor, /*useSoA*/ false, /*useNewton3*/ false,
                   BlackBoxTraversalOption::outer>
          traversal(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), f);
      traversal.traverseCellPairs(this->_linkedCells.getCells());
    }
  }

  /**
   * Iterate over the verlet lists using the SoA traversal
   * @tparam ParticleFunctor
   * @param f
   * @param useNewton3
   */
  template <class ParticleFunctor>
  void iterateVerletListsSoA(ParticleFunctor* f, const bool useNewton3) {
    /// @todo optimize iterateVerletListsSoA, e.g. by using traversals with openmp possibilities

    // load data from cells into soa
    loadVerletSoA(f);

    /// @todo here you can (sort of) use traversals, by modifying iFrom and iTo.
    const size_t iFrom = 0;
    const size_t iTo = _soaNeighborLists.size();

#if defined(AUTOPAS_OPENMP)
    if (not useNewton3) {
      /// @todo find a sensible chunk size
      const size_t chunkSize = std::max((iTo - iFrom) / (omp_get_max_threads() * 10), 1ul);
#pragma omp parallel for schedule(dynamic, chunkSize)
      for (size_t i = iFrom; i < iTo; i++) {
        f->SoAFunctor(_soa, _soaNeighborLists, i, i + 1, useNewton3);
      }
    } else
#endif
    {
      // iterate over SoA
      f->SoAFunctor(_soa, _soaNeighborLists, iFrom, iTo, useNewton3);
    }

    // extract SoA
    extractVerletSoA(f);
  }

  /**
   * Update the AoS id maps.
   * The Id Map is used to map the id of a particle to the actual particle
   */
  void updateIdMapAoS() {
    /// @todo: potentially adapt to _blackBox -- only consider inner parts -- not necessary, but might be useful
    _aosNeighborLists.clear();
    // DON'T simply parallelize this loop!!! this needs modifications if you
    // want to parallelize it!
    for (auto iter = this->begin(); iter.isValid(); ++iter) {
      // create the verlet list entries for all particles
      _aosNeighborLists[&(*iter)];
    }
  }

  /**
   * Loops over relevant parts of the boundary.
   * @tparam LoopBody
   * @param dims
   * @param loopBody
   */
  template <class LoopBody>
  void boundaryRelevantTraversal(const std::array<size_t, 3>& dims, LoopBody loopBody) {
#if defined(AUTOPAS_OPENMP)
#pragma omp parallel
#endif
    {
      // lower z
#if defined(AUTOPAS_OPENMP)
#pragma omp for collapse(3) nowait
#endif
      for (size_t z = 0; z < 3; ++z) {
        for (size_t y = 0; y < dims[1]; ++y) {
          for (size_t x = 0; x < dims[0]; ++x) {
            loopBody(x, y, z);
          }
        }
      }
      // upper z
#if defined(AUTOPAS_OPENMP)
#pragma omp for collapse(3) nowait
#endif
      for (size_t z = dims[2] - 3; z < dims[2]; ++z) {
        for (size_t y = 0; y < dims[1]; ++y) {
          for (size_t x = 0; x < dims[0]; ++x) {
            loopBody(x, y, z);
          }
        }
      }

      // lower y
#if defined(AUTOPAS_OPENMP)
#pragma omp for collapse(3) nowait
#endif
      for (size_t z = 3; z < dims[2] - 3; ++z) {
        for (size_t y = 0; y < 3; ++y) {
          for (size_t x = 0; x < dims[0]; ++x) {
            loopBody(x, y, z);
          }
        }
      }
      // upper y
#if defined(AUTOPAS_OPENMP)
#pragma omp for collapse(3) nowait
#endif
      for (size_t z = 3; z < dims[2] - 3; ++z) {
        for (size_t y = dims[1] - 3; y < dims[1]; ++y) {
          for (size_t x = 0; x < dims[0]; ++x) {
            loopBody(x, y, z);
          }
        }
      }

      // lower x
#if defined(AUTOPAS_OPENMP)
#pragma omp for collapse(3) nowait
#endif
      for (size_t z = 3; z < dims[2] - 3; ++z) {
        for (size_t y = 3; y < dims[1] - 3; ++y) {
          for (size_t x = 0; x < 3; ++x) {
            loopBody(x, y, z);
          }
        }
      }
      // upper x
#if defined(AUTOPAS_OPENMP)
#pragma omp for collapse(3) nowait
#endif
      for (size_t z = 3; z < dims[2] - 3; ++z) {
        for (size_t y = 3; y < dims[1] - 3; ++y) {
          for (size_t x = dims[0] - 3; x < dims[0]; ++x) {
            loopBody(x, y, z);
          }
        }
      }
    }
  }

  /**
   * Load the SoA's for all relevant cells of the boundary traversal using
   * functor.SoALoader(...)
   * @tparam ParticleFunctor The type of the functor.
   * @param functor The SoALoader method of this functor is used.
   */
  template <class ParticleFunctor>
  void loadBoundarySoA(ParticleFunctor* functor) {
    auto& cells = this->_linkedCells.getCells();
    auto& dims = this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo();
    boundaryRelevantTraversal(dims, [&](size_t x, size_t y, size_t z) {
      auto& cell = cells[utils::ThreeDimensionalMapping::threeToOneD(x, y, z, dims)];
      functor->SoALoader(cell, cell._particleSoABuffer, 0);
    });
  }

  /**
   * Extracts the particle information from the cell SoA's using
   * functor.SoAExtractor(...)
   * @tparam ParticleFunctor The type of the functor.
   * @param functor The SoAExtractor method of this functor is used.
   */
  template <class ParticleFunctor>
  void extractBoundarySoA(ParticleFunctor* functor) {
    auto& cells = this->_linkedCells.getCells();
    auto& dims = this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo();
    boundaryRelevantTraversal(dims, [&](size_t x, size_t y, size_t z) {
      auto& cell = cells[utils::ThreeDimensionalMapping::threeToOneD(x, y, z, dims)];
      functor->SoAExtractor(cell, cell._particleSoABuffer, 0);
    });
  }

  /**
   * Load the particle information from the cell and store it in the global SoA
   * using functor.SoALoader(...)
   * @tparam ParticleFunctor The type of the functor.
   * @param functor The SoALoader method of this functor is used.
   */
  template <class ParticleFunctor>
  void loadVerletSoA(ParticleFunctor* functor) {
    size_t offset = 0;
    /// @todo adapt to _blackBoxMode
    for (auto& cell : this->_linkedCells.getCells()) {
      functor->SoALoader(cell, _soa, offset);
      offset += cell.numParticles();
    }
  }

  /**
   * Extracts the particle information from the global SoA using
   * functor.SoAExtractor(...)
   * @tparam ParticleFunctor The type of the functor.
   * @param functor The SoAExtractor method of this functor is used.
   */
  template <class ParticleFunctor>
  void extractVerletSoA(ParticleFunctor* functor) {
    size_t offset = 0;
    /// @todo adapt to _blackBoxMode
    for (auto& cell : this->_linkedCells.getCells()) {
      functor->SoAExtractor(cell, _soa, offset);
      offset += cell.numParticles();
    }
  }

  /**
   * Converts the verlet list stored for AoS usage into one for SoA usage
   */
  void generateSoAListFromAoSVerletLists() {
    /// @todo adapt to _blackBoxMode
    /// @todo openmp parallelization?

    // resize the list to the size of the aos neighborlist
    _soaNeighborLists.resize(_aosNeighborLists.size());
    // clear the aos 2 soa map
    _aos2soaMap.clear();

    _aos2soaMap.reserve(_aosNeighborLists.size());
    size_t i = 0;
    for (auto iter = this->begin(); iter.isValid(); ++iter, ++i) {
      // set the map
      _aos2soaMap[&(*iter)] = i;
    }
    size_t accumulatedListSize = 0;
    for (auto& aosList : _aosNeighborLists) {
      accumulatedListSize += aosList.second.size();
      size_t i_id = _aos2soaMap[aosList.first];
      // each soa neighbor list should be of the same size as for aos
      _soaNeighborLists[i_id].resize(aosList.second.size());
      size_t j = 0;
      for (auto neighbor : aosList.second) {
        _soaNeighborLists[i_id][j] = _aos2soaMap.at(neighbor);
        j++;
      }
    }
    AutoPasLog(debug,
               "VerletLists::generateSoAListFromAoSVerletLists: average verlet list "
               "size is {}",
               static_cast<double>(accumulatedListSize) / _aosNeighborLists.size());
    _soaListIsValid = true;
  }

 private:
  /// verlet lists.
  typename verlet_internal::AoS_verletlist_storage_type _aosNeighborLists;

  /// map converting from the aos type index (Particle *) to the soa type index
  /// (continuous, size_t)
  std::unordered_map<Particle*, size_t> _aos2soaMap;

  /// verlet list for SoA:
  std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> _soaNeighborLists;

  // specifies if the SoA neighbor list is currently valid
  bool _soaListIsValid;

  /// global SoA of verlet lists
  SoA<typename Particle::SoAArraysType> _soa;

  /// specifies how the verlet lists are build
  BuildVerletListType _buildVerletListType;
};

}  // namespace autopas
