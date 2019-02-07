/**
 * @file VerletListsCells.h
 * @author nguyen
 * @date 30.08.18
 */

#pragma once

#include "VerletListsCellsHelpers.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/traversals/C01Traversal.h"
#include "autopas/containers/linkedCells/traversals/C08Traversal.h"
#include "autopas/containers/linkedCells/traversals/C18Traversal.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

/**
 * Linked Cells with Verlet Lists container.
 * The VerletListsCells class uses neighborhood lists for each cell
 * to calculate pairwise interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * Cells are created using a cell size of at least cutoff + skin radius.
 * @tparam Particle
 */
template <class Particle>
class VerletListsCells
    : public VerletListsLinkedBase<Particle, typename VerletListsCellsHelpers<Particle>::VerletListParticleCellType> {
  typedef VerletListsCellsHelpers<Particle> verlet_internal;
  typedef FullParticleCell<Particle> ParticleCell;
  typedef typename VerletListsCellsHelpers<Particle>::VerletListParticleCellType LinkedParticleCell;

 private:
  const std::vector<TraversalOptions>& VLCApplicableTraversals() {
    static const std::vector<TraversalOptions> v{TraversalOptions::slicedVerlet, TraversalOptions::c18Verlet,
                                                 TraversalOptions::c01Verlet};
    return v;
  }

 public:
  /**
   * Constructor of the VerletListsCells class.
   * The neighbor lists are build using a search radius of cutoff + skin.
   * The rebuildFrequency should be chosen, s.t. the particles do not move more
   * than a distance of skin/2 between two rebuilds of the lists.
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param cutoff the cutoff radius of the interaction
   * @param skin the skin radius
   * @param rebuildFrequency specifies after how many pair-wise traversals the
   * neighbor lists are to be rebuild. A frequency of 1 means that they are
   * always rebuild, 10 means they are rebuild after 10 traversals
   * @param buildTraversal the traversal used to build the verletlists
   */
  VerletListsCells(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                   const TraversalOptions buildTraversal, const double skin = 0,
                   const unsigned int rebuildFrequency = 1)
      : VerletListsLinkedBase<Particle, LinkedParticleCell>(boxMin, boxMax, cutoff, skin, rebuildFrequency,
                                                            VLCApplicableTraversals(),
                                                            /*blackbox not yet supported*/ false),
        _buildTraversal(buildTraversal) {}

  ContainerOptions getContainerType() override { return ContainerOptions::verletListsCells; }

  /**
   * Function to iterate over all pairs of particles.
   * This function only handles short-range interactions.
   * @tparam the type of ParticleFunctor
   * @tparam Traversal
   * @param f functor that describes the pair-potential
   * @param traversal the traversal that will be used
   * @param useNewton3 whether newton 3 optimization should be used
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwiseAoS(ParticleFunctor* f, Traversal* traversal, bool useNewton3 = true) {
    if (this->needsRebuild(useNewton3)) {
      updateVerletLists(useNewton3);
    }

    if (useNewton3) {
      if (auto vTraversal =
              dynamic_cast<autopas::VerletListsCellsTraversal<Particle, ParticleFunctor, true>*>(traversal))
        vTraversal->traverseCellVerlet(_neighborLists);
      else
        autopas::utils::ExceptionHandler::exception("wrong type of traversal in VerletListCells.h. TraversalID: {}",
                                                    traversal->getTraversalType());
    } else {
      if (auto vTraversal =
              dynamic_cast<autopas::VerletListsCellsTraversal<Particle, ParticleFunctor, false>*>(traversal))
        vTraversal->traverseCellVerlet(_neighborLists);
      else
        autopas::utils::ExceptionHandler::exception("wrong type of traversal in VerletListCells.h. TraversalID: {}",
                                                    traversal->getTraversalType());
    }

    // we iterated, so increase traversal counter
    this->_traversalsSinceLastRebuild++;
  }

  /**
   * Dummy function. (Uses AoS instead)
   * @tparam the type of ParticleFunctor
   * @tparam Traversal
   * @param f functor that describes the pair-potential
   * @param traversal the traversal that will be used
   * @param useNewton3 whether newton 3 optimization should be used
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwiseSoA(ParticleFunctor* f, Traversal* traversal, bool useNewton3 = true) {
    iteratePairwiseAoS(f, traversal, useNewton3);
  }

  TraversalSelector<ParticleCell> generateTraversalSelector(std::vector<TraversalOptions> traversalOptions) override {
    std::vector<TraversalOptions> allowedAndApplicable;

    std::sort(traversalOptions.begin(), traversalOptions.end());
    std::set_intersection(this->_applicableTraversals.begin(), this->_applicableTraversals.end(),
                          traversalOptions.begin(), traversalOptions.end(), std::back_inserter(allowedAndApplicable));

    return TraversalSelector<ParticleCell>(this->getCellsPerDimension(), allowedAndApplicable);
  }

  /**
   * Get the neighbors list of a particle.
   * @param particle
   * @param useNewton3
   * @return the neighbor list of the particle
   */
  std::vector<Particle*>& getVerletList(Particle* particle, bool useNewton3 = true) {
    if (this->needsRebuild(useNewton3)) {
      updateVerletLists(useNewton3);
    }
    auto indices = _cellMap[particle];
    return _neighborLists[indices.first][indices.second].second;
  }

 protected:
  /**
   * Update the verlet lists.
   * @param useNewton3 use newton3?
   */
  void updateVerletLists(bool useNewton3) {
    this->_verletBuiltNewton3 = useNewton3;

    // create a Verlet Lists for each cell
    _neighborLists.clear();
    auto& cells = this->_linkedCells.getCells();
    size_t cellsSize = cells.size();
    _neighborLists.resize(cellsSize);
    for (size_t cellIndex = 0; cellIndex < cellsSize; ++cellIndex) {
      size_t i = 0;
      for (auto iter = cells[cellIndex].begin(); iter.isValid(); ++iter, ++i) {
        Particle* particle = &*iter;
        _neighborLists[cellIndex].push_back(
            std::pair<Particle*, std::vector<Particle*>>(particle, std::vector<Particle*>()));
        _cellMap[particle] = std::pair<size_t, size_t>(cellIndex, i);
      }
    }

    typename verlet_internal::VerletListGeneratorFunctor f(_neighborLists, _cellMap, this->getCutoff());

    switch (_buildTraversal) {
      case c08: {
        if (useNewton3) {
          auto traversal =
              C08Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor, false, true>(
                  this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
          this->_linkedCells.iteratePairwiseAoS(&f, &traversal);
        } else {
          auto traversal =
              C08Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor, false, false>(
                  this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
          this->_linkedCells.iteratePairwiseAoS(&f, &traversal);
        }
        break;
      }
      case c18: {
        if (useNewton3) {
          auto traversal =
              C18Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor, false, true>(
                  this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
          this->_linkedCells.iteratePairwiseAoS(&f, &traversal);
        } else {
          auto traversal =
              C18Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor, false, false>(
                  this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
          this->_linkedCells.iteratePairwiseAoS(&f, &traversal);
        }
        break;
      }
      case c01: {
        if (not useNewton3) {
          auto traversal =
              C01Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor, false, false>(
                  this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
          this->_linkedCells.iteratePairwiseAoS(&f, &traversal);
        } else {
          utils::ExceptionHandler::exception("VerletListsCells::updateVerletLists(): c01 does not support newton3");
        }
        break;
      }
      default:
        utils::ExceptionHandler::exception("VerletListsCells::updateVerletLists(): unsupported Traversal: {}",
                                           _buildTraversal);
        break;
    }

    // the neighbor list is now valid
    this->_neighborListIsValid = true;
    this->_traversalsSinceLastRebuild = 0;
  }

 private:
  /// verlet lists for each particle for each cell
  typename verlet_internal::VerletList_storage_type _neighborLists;

  /// mapping each particle to its corresponding cell and position in this cell
  std::unordered_map<Particle*, std::pair<size_t, size_t>> _cellMap;

  // the traversal used to build the verletlists
  TraversalOptions _buildTraversal;
};

}  // namespace autopas
