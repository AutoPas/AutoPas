/**
 * @file VerletListsCells.h
 * @author nguyen
 * @date 30.08.18
 */

#pragma once

#include "VerletListsCellsHelpers.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/traversals/C01Traversal.h"
#include "autopas/containers/linkedCells/traversals/C08Traversal.h"
#include "autopas/containers/linkedCells/traversals/C18Traversal.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VerletListsCellsTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/TraversalOption.h"
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

 public:
  /**
   * Constructor of the VerletListsCells class.
   * The neighbor lists are build using a search radius of cutoff + skin.
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param cutoff the cutoff radius of the interaction
   * @param skin the skin radius
   * @param buildTraversal the traversal used to build the verletlists
   * @param cellSizeFactor cell size factor relative to cutoff
   */
  VerletListsCells(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                   const TraversalOption buildTraversal, const double skin = 0, const double cellSizeFactor = 1.0)
      : VerletListsLinkedBase<Particle, LinkedParticleCell>(
            boxMin, boxMax, cutoff, skin, compatibleTraversals::allVLCCompatibleTraversals(), cellSizeFactor),
        _buildTraversal(buildTraversal) {}

  ContainerOption getContainerType() const override { return ContainerOption::verletListsCells; }

  void iteratePairwise(TraversalInterface *traversal) override {
    AutoPasLog(debug, "Using traversal {}.", utils::StringUtils::to_string(traversal->getTraversalType()));

    // Check if traversal is allowed for this container and give it the data it needs.
    auto vTraversal = dynamic_cast<autopas::VerletListsCellsTraversal<Particle> *>(traversal);
    if (vTraversal) {
      vTraversal->setVerletList(_neighborLists);
    } else {
      autopas::utils::ExceptionHandler::exception("wrong type of traversal in VerletListCells.h. TraversalID: {}",
                                                  traversal->getTraversalType());
    }

    traversal->initTraversal();
    traversal->traverseParticlePairs();
    traversal->endTraversal();
  }

  /**
   * Get the neighbors list of a particle.
   * @param particle
   * @return the neighbor list of the particle
   */
  const std::vector<Particle *> &getVerletList(const Particle *particle) const {
    const auto indices = _cellMap.at(const_cast<Particle *>(particle));
    return _neighborLists.at(indices.first).at(indices.second).second;
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    bool useNewton3 = traversal->getUseNewton3();
    this->_verletBuiltNewton3 = useNewton3;

    // create a Verlet Lists for each cell
    _neighborLists.clear();
    auto &cells = this->_linkedCells.getCells();
    size_t cellsSize = cells.size();
    _neighborLists.resize(cellsSize);
    for (size_t cellIndex = 0; cellIndex < cellsSize; ++cellIndex) {
      size_t i = 0;
      for (auto iter = cells[cellIndex].begin(); iter.isValid(); ++iter, ++i) {
        Particle *particle = &*iter;
        _neighborLists[cellIndex].push_back(
            std::pair<Particle *, std::vector<Particle *>>(particle, std::vector<Particle *>()));
        _cellMap[particle] = std::pair<size_t, size_t>(cellIndex, i);
      }
    }

    typename verlet_internal::VerletListGeneratorFunctor f(_neighborLists, _cellMap,
                                                           this->getCutoff() + this->getSkin());

    // FIXME: Clang compiler bug makes this necessary
    switch (static_cast<TraversalOption>(_buildTraversal)) {
//    switch (_buildTraversal) {
      case TraversalOption::c08: {
        autopas::utils::withStaticBool(useNewton3, [&](auto n3) {
          auto buildTraversal = C08Traversal<LinkedParticleCell, decltype(f), DataLayoutOption::aos, n3>(
              this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getInteractionLength(),
              this->_linkedCells.getCellBlock().getCellLength());
          this->_linkedCells.iteratePairwise(&buildTraversal);
        });
        break;
      }
      case TraversalOption::c18: {
        autopas::utils::withStaticBool(useNewton3, [&](auto n3) {
          auto buildTraversal = C18Traversal<LinkedParticleCell, decltype(f), DataLayoutOption::aos, n3>(
              this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getInteractionLength(),
              this->_linkedCells.getCellBlock().getCellLength());
          this->_linkedCells.iteratePairwise(&buildTraversal);
        });
        break;
      }
      case TraversalOption::c01: {
        if (useNewton3) {
          utils::ExceptionHandler::exception("VerletListsCells::updateVerletLists(): c01 does not support newton3");
        } else {
          auto buildTraversal = C01Traversal<LinkedParticleCell, decltype(f), DataLayoutOption::aos, false>(
              this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getInteractionLength(),
              this->_linkedCells.getCellBlock().getCellLength());
          this->_linkedCells.iteratePairwise(&buildTraversal);
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
  }

  /**
   * Return the cell length of the underlying linked cells structure, normally needed only for unit tests.
   * @return
   */
  const std::array<double, 3> &getCellLength() { return this->_linkedCells.getCellBlock().getCellLength(); }

 private:
  /// verlet lists for each particle for each cell
  typename verlet_internal::VerletList_storage_type _neighborLists;

  /// mapping each particle to its corresponding cell and position in this cell
  std::unordered_map<Particle *, std::pair<size_t, size_t>> _cellMap;

  // the traversal used to build the verletlists
  TraversalOption _buildTraversal;
};

}  // namespace autopas
