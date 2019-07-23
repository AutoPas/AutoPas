/**
 * @file VerletClusterCells.h
 * @author jspahl
 * @date 25.3.19
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CellBorderAndFlagManager.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/cellPairTraversals/VerletClusterTraversalInterface.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/CudaDeviceVector.h"

namespace autopas {

/**
 * Particles are divided into clusters.
 * The VerletClusterCells class uses neighborhood lists for each cluster pair
 * to calculate pairwise interactions.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * @tparam Particle
 */
template <class Particle>
class VerletClusterCells : public ParticleContainer<Particle, FullParticleCell<Particle>> {
 public:
  /**
   * Constructor of the VerletClusterLists class.
   * The neighbor lists are build using a estimated density.
   * The box is divided into cuboids with roughly the
   * same side length. The rebuildFrequency should be chosen, s.t. the particles do
   * not move more than a distance of skin/2 between two rebuilds of the lists.
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param cutoff the cutoff radius of the interaction
   * @param skin the skin radius
   * @param clusterSize size of clusters
   */
  VerletClusterCells(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff,
                     double skin = 0, int clusterSize = 32)
      : ParticleContainer<Particle, FullParticleCell<Particle>>(boxMin, boxMax, cutoff, skin),
        _boxMinWithHalo(ArrayMath::subScalar(boxMin, cutoff + skin)),
        _boxMaxWithHalo(ArrayMath::addScalar(boxMax, cutoff + skin)),
        _clusterSize(clusterSize),
        _isValid(false) {
    this->_cells.resize(1);
    _dummyStarts = {0};
  }

  /**
   * Lists all traversal options applicable for the Verlet Lists container.
   * @return Vector of all applicable traversal options.
   */
  static const std::vector<TraversalOption> &allVCLApplicableTraversals() {
    // traversal not used but prevents usage of newton3
    static const std::vector<TraversalOption> v{TraversalOption::verletClusterCellsTraversal};
    return v;
  }

  ContainerOption getContainerType() override { return ContainerOption::verletClusterCells; }

  /**
   * Function to iterate over all pairs of particles.
   * This function only handles short-range interactions.
   * @tparam the type of ParticleFunctor
   * @tparam Traversal
   * @param f functor that describes the pair-potential
   * @param traversal not used
   * @param useNewton3 whether newton 3 optimization should be used
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwise(ParticleFunctor *f, Traversal *traversal) {
    auto *traversalInterface = dynamic_cast<VerletClusterTraversalInterface<FullParticleCell<Particle>> *>(traversal);
    if (!traversalInterface) {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in VerletClusterCells::iteratePairwise");
    }

    traversalInterface->setVerletListPointer(&_clusterSize, &_neighborCellIds, &_neighborMatrixDim, &_neighborMatrix);

    if (traversalInterface->getSignature() != _lastTraversalSig or (not _isValid)) {
      if (!_isValid) {
        rebuild();
      }
      traversalInterface->rebuildVerlet(_cellsPerDim, this->_cells, _boundingBoxes,
                                        std::ceil(this->getInteractionLength() * _gridSideLengthReciprocal),
                                        this->getInteractionLength());
      _lastTraversalSig = traversalInterface->getSignature();
    }

    traversal->initTraversal(this->_cells);
    traversalInterface->traverseCellPairs(this->_cells);
    traversal->endTraversal(this->_cells);
  }

  /**
   * @copydoc VerletLists::addParticle()
   */
  void addParticle(Particle &p) override {
    if (autopas::utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax())) {
      _isValid = false;
      this->_cells[0].resize(_dummyStarts[0]);
      p.setOwned(true);
      // add particle somewhere, because lists will be rebuild anyways
      this->_cells[0].addParticle(p);
      ++_dummyStarts[0];
    } else {
      utils::ExceptionHandler::exception(
          "VerletCluster: trying to add particle that is not inside the bounding box.\n" + p.toString());
    }
  }

  /**
   * @copydoc VerletLists::addHaloParticle()
   */
  void addHaloParticle(Particle &haloParticle) override {
    _isValid = false;
    if (autopas::utils::notInBox(haloParticle.getR(), this->getBoxMin(), this->getBoxMax())) {
      if (autopas::utils::inBox(haloParticle.getR(), _boxMinWithHalo, _boxMaxWithHalo)) {
        this->_cells[0].resize(_dummyStarts[0]);
        haloParticle.setOwned(false);
        // add particle somewhere, because lists will be rebuild anyways
        this->_cells[0].addParticle(haloParticle);
        ++_dummyStarts[0];
      } else {
        utils::ExceptionHandler::exception(
            "VerletCluster: trying to add halo particle that is not near the bounding box.\n" +
            haloParticle.toString());
      }
    } else {
      utils::ExceptionHandler::exception(
          "VerletCluster: trying to add halo particle that is inside the bounding box.\n" + haloParticle.toString());
    }
  }

  /**
   * Update a halo particle of the container with the given haloParticle.
   * @param haloParticle Particle to be updated.
   * @return Returns true if the particle was updated, false if no particle could be found.
   */
  bool updateHaloParticle(Particle &haloParticle) override {
    int xt = (int)((haloParticle.getR()[0] - _boxMinWithHalo[0]) * _gridSideLengthReciprocal);
    int yt = (int)((haloParticle.getR()[1] - _boxMinWithHalo[1]) * _gridSideLengthReciprocal);
    for (int x = -1; x <= 1; ++x) {
      for (int y = -1; y <= 1; ++y) {
        size_t index = xt + yt * _cellsPerDim[0];
        if (0 <= index && index < this->_cells.size()) {
          Particle range = haloParticle;
          auto pos = haloParticle.getR();
          pos[2] -= this->getSkin();
          range.setR(pos);
          auto lower =
              std::lower_bound(this->_cells[index]._particles.begin(), this->_cells[index]._particles.end(), range,
                               [](const Particle &a, const Particle &b) { return a.getR()[2] < b.getR()[2]; });

          pos = haloParticle.getR();
          pos[2] += this->getSkin();
          range.setR(pos);
          auto upper =
              std::upper_bound(this->_cells[index]._particles.begin(), this->_cells[index]._particles.end(), range,
                               [](const Particle &a, const Particle &b) { return a.getR()[2] < b.getR()[2]; });

          for (; lower != upper; ++lower) {
            if (haloParticle.getID() == lower->getID()) {
              *lower = haloParticle;
              return true;
            }
          }
        }
      }
    }
    return false;
  }

  /**
   * Rebuilds the neighbor lists.
   * @param traversal The used traversal.
   */
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    auto *traversalInterface = dynamic_cast<VerletClusterTraversalInterface<FullParticleCell<Particle>> *>(traversal);
    if (!traversalInterface) {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in VerletClusterCells::iteratePairwise");
    }
    if (!_isValid) {
      rebuild();
    }
    traversalInterface->setVerletListPointer(&_clusterSize, &_neighborCellIds, &_neighborMatrixDim, &_neighborMatrix);

    traversalInterface->rebuildVerlet(_cellsPerDim, this->_cells, _boundingBoxes,
                                      std::ceil(this->getInteractionLength() * _gridSideLengthReciprocal),
                                      this->getInteractionLength());
    _lastTraversalSig = traversalInterface->getSignature();
  }

  /**
   * @copydoc VerletLists::deleteHaloParticles
   */
  void deleteHaloParticles() override {
    _isValid = false;

    for (size_t i = 0; i < this->_cells.size(); ++i) {
    	std::remove_if(this->_cells[i]._particles.begin(), this->_cells[i]._particles.end(), [](Particle& p){return not p.isOwned();});
    }
  }

  /**
   * @copydoc VerletLists::updateContainer()
   */
  std::vector<Particle> updateContainer() override {
    AutoPasLog(debug, "updating container");
    _isValid = false;

    return rebuild();
  }

  bool isContainerUpdateNeeded() override {
    if (not _isValid) {
      return true;
    }
    for (size_t i = 0; i < this->_cells.size(); ++i) {
      size_t pid = 0;
      const size_t end = (_boundingBoxes[i].size() > 0) ? _boundingBoxes[i].size() - 1 : 0;

      for (size_t cid = 0; cid < end; ++cid)
        for (unsigned int pic = 0; pic < _clusterSize; ++pic) {
          if (not particleInSkinOfBox(_boundingBoxes[i][cid], this->_cells[i][pid])) {
            return true;
          }
          ++pid;
        }
      for (unsigned int pic = 0; pic < _clusterSize && pid < _dummyStarts[i]; ++pic) {
        if (not particleInSkinOfBox(_boundingBoxes[i][_boundingBoxes[i].size() - 1], this->_cells[i][pid])) {
          return true;
        }
        ++pid;
      }
    }

    return false;
  }

  TraversalSelectorInfo getTraversalSelectorInfo() override {
    return TraversalSelectorInfo(_cellsPerDim, this->getInteractionLength());
  }

  ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(new internal::ParticleIterator<Particle, FullParticleCell<Particle>>(
        &this->_cells, 0, &_cellBorderFlagManager, behavior));
  }

  ParticleIteratorWrapper<Particle> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                      const std::array<double, 3> &higherCorner,
                                                      IteratorBehavior behavior = IteratorBehavior::haloAndOwned,
                                                      bool incSearchRegion = false) override {
    std::vector<size_t> cellsOfInterest(this->_cells.size());
    std::iota(cellsOfInterest.begin(), cellsOfInterest.end(), 0);

    return ParticleIteratorWrapper<Particle>(new internal::RegionParticleIterator<Particle, FullParticleCell<Particle>>(
        &this->_cells, lowerCorner, higherCorner, cellsOfInterest, &_cellBorderFlagManager, behavior));
  }

  /**
   * Deletes all Dummy Particles in the container
   */
  void deleteDummyParticles() override {
    for (size_t i = 0; i < this->_cells.size(); ++i) {
      this->_cells[i].resize(_dummyStarts[i]);
    }
    _isValid = false;
  }

 protected:
  /**
   * Recalculate grids and clusters,
   * build verlet lists and pad clusters.
   */
  std::vector<Particle> rebuild() {
    deleteDummyParticles();
    _boundingBoxes.clear();
    // get the dimensions and volumes of the box
    std::array<double, 3> boxSize{};
    double volume = 1.0;

    for (int d = 0; d < 3; ++d) {
      boxSize[d] = _boxMaxWithHalo[d] - _boxMinWithHalo[d];
      volume *= boxSize[d];
    }

    // get all particles and clear clusters
    std::vector<Particle> invalidParticles;
    std::vector<Particle> outsideParticles;

    for (size_t i = 0; i < this->_cells.size(); ++i) {
      for (auto &p : this->_cells[i]._particles) {
        if (utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax())) {
          invalidParticles.push_back(p);
        } else {
          if (p.isOwned()) {
            outsideParticles.push_back(p);
          } else {
            invalidParticles.push_back(p);
          }
        }
      }
      this->_cells[i].clear();
    }
    // estimate particle density
    double density = (std::max(1.0, (double)invalidParticles.size())) / volume;

    // guess optimal grid side length
    _gridSideLength = std::cbrt(((double)_clusterSize) / density);
    _gridSideLengthReciprocal = 1 / _gridSideLength;

    // get cells per dimension
    size_t sizeGrid = 1;
    for (int d = 0; d < 2; d++) {
      _cellsPerDim[d] = static_cast<size_t>(std::ceil(boxSize[d] * _gridSideLengthReciprocal));
      sizeGrid *= _cellsPerDim[d];
    }
    _cellsPerDim[2] = static_cast<size_t>(1);

    // resize to number of grids
    if (this->_cells.size() < sizeGrid) {
      this->_cells.resize(sizeGrid);
    }
    _dummyStarts.clear();
    _dummyStarts.resize(sizeGrid);
    _boundingBoxes.resize(sizeGrid);

    // put particles into grid cells
    for (size_t i = 0; i < invalidParticles.size(); ++i) {
      int index =
          (int)((invalidParticles[i].getR()[0] - _boxMinWithHalo[0]) * _gridSideLengthReciprocal) +
          (int)((invalidParticles[i].getR()[1] - _boxMinWithHalo[1]) * _gridSideLengthReciprocal) * _cellsPerDim[0];
      this->_cells[index].addParticle(invalidParticles[i]);
    }

    // sort by last dimension
#if defined(AUTOPAS_OPENMP)
#pragma omp parallel for schedule(guided)
#endif
    for (size_t i = 0; i < sizeGrid; ++i) {
      this->_cells[i].sortByDim(2);
      const auto numParticles = this->_cells[i].numParticles();
      _dummyStarts[i] = numParticles;
      unsigned int numDummys = _clusterSize;
      if (numParticles > 0) {
        numDummys -= (numParticles % _clusterSize);
      }

      for (unsigned int j = 0; j < numDummys; ++j) {
        Particle dummyParticle = Particle();
        dummyParticle.setR({_boxMaxWithHalo[0] + 8 * this->getInteractionLength() + static_cast<double>(i),
                            _boxMaxWithHalo[1] + 8 * this->getInteractionLength() + static_cast<double>(j),
                            _boxMaxWithHalo[2] + 8 * this->getInteractionLength()});
        dummyParticle.setID(ULONG_MAX);
        this->_cells[i].addParticle(dummyParticle);
      }
    }

#if defined(AUTOPAS_OPENMP)
#pragma omp parallel for schedule(guided)
#endif
    for (size_t i = 0; i < sizeGrid; ++i) {
      const size_t nClusters = this->_cells[i].numParticles() / _clusterSize;

      _boundingBoxes[i].resize(nClusters, {_boxMaxWithHalo[0], _boxMaxWithHalo[1], _boxMaxWithHalo[2],
                                           _boxMinWithHalo[0], _boxMinWithHalo[1], _boxMinWithHalo[2]});

      for (size_t cid = 0; cid < nClusters; ++cid)
        for (size_t pid = cid * _clusterSize; pid < _dummyStarts[i]; ++pid) {
          expandBoundingBox(_boundingBoxes[i][cid], this->_cells[i][pid]);
        }
    }

    _isValid = true;
    return outsideParticles;
  }

 private:
  /**
   * Expands a bounding Box such the Particle is in it
   * @param box
   * @param p
   */
  void expandBoundingBox(std::array<double, 6> &box, Particle &p) {
    for (int i = 0; i < 3; ++i) {
      box[i] = std::min(box[i], p.getR()[i]);
      box[3 + i] = std::max(box[3 + i], p.getR()[i]);
    }
  }

  /**
   * Checks if particle is within skin of bounding box
   * @param box
   * @param p
   *    */
  bool particleInSkinOfBox(std::array<double, 6> &box, Particle &p) {
    for (int i = 0; i < 3; ++i) {
      if (box[0 + i] - this->getSkin() > p.getR()[i] || box[3 + i] + this->getSkin() < p.getR()[i]) return false;
    }
    return true;
  }
  std::array<double, 3> _boxMinWithHalo, _boxMaxWithHalo;

  /// indices where dummy particles in the cells start
  std::vector<size_t> _dummyStarts;

  // number of particles in a cluster
  unsigned int _clusterSize;

  // id of neighbor clusters of a clusters
  std::vector<std::vector<std::vector<std::pair<size_t, size_t>>>> _neighborCellIds;

  size_t _neighborMatrixDim;
  utils::CudaDeviceVector<unsigned int> _neighborMatrix;

  // bounding boxes of all clusters (xmin,ymin,zmin,xmax,ymax,zmax)
  std::vector<std::vector<std::array<double, 6>>> _boundingBoxes;

  // side length of xy-grid and reciprocal
  double _gridSideLength;
  double _gridSideLengthReciprocal;

  // dimensions of grid
  std::array<size_t, 3> _cellsPerDim;

  // specifies if the neighbor list is currently valid
  bool _isValid;

  /// Signature of the last Traversal to trigger rebuild when a new one is used
  std::tuple<TraversalOption, DataLayoutOption, bool> _lastTraversalSig;

  class VerletClusterCellsCellBorderAndFlagManager : public internal::CellBorderAndFlagManager {
   public:
    bool cellCanContainHaloParticles(size_t index1d) const override { return true; }

    bool cellCanContainOwnedParticles(size_t index1d) const override { return true; }

  } _cellBorderFlagManager;
};

namespace internal {
/**
 * ParticleIterator class to access particles inside a VerletClusterCells container efficently.
 * The particles can be accessed using "iterator->" or "*iterator". The next
 * particle using the ++operator, e.g. "++iterator".
 * @tparam Particle Type of the particle that is accessed.
 */
template <class Particle>
class VerletClusterCellsParticleIterator : public ParticleIteratorInterfaceImpl<Particle> {
};
} // namespace internal
}  // namespace autopas
