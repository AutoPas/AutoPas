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
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/cellPairTraversals/VerletClusterTraversalInterface.h"
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
  using ParticleFloatType = typename Particle::ParticleFloatingPointType;

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
   * @param rebuildFrequency specifies after how many pair-wise traversals the
   * neighbor lists are to be rebuild. A frequency of 1 means that they are
   * always rebuild, 10 means they are rebuild after 10 traversals.
   * @param clusterSize size of clusters
   */
  VerletClusterCells(const std::array<ParticleFloatType, 3> boxMin, const std::array<ParticleFloatType, 3> boxMax,
                     ParticleFloatType cutoff, ParticleFloatType skin = 0, unsigned int rebuildFrequency = 1,
                     int clusterSize = 32)
      : ParticleContainer<Particle, FullParticleCell<Particle>>(boxMin, boxMax, cutoff + skin,
                                                                allVCLApplicableTraversals()),
        _firstHaloClusterId(1),
        _clusterSize(clusterSize),
        _traversalsSinceLastRebuild(UINT_MAX),
        _rebuildFrequency(rebuildFrequency),
        _isValid(false),
        _cellBorderFlagManager(&_firstHaloClusterId) {
    this->_cells.resize(1);
    _dummyStarts[0] = 0;
  }

  /**
   * Lists all traversal options applicable for the Verlet Lists container.
   * @return Vector of all applicable traversal options.
   */
  static const std::vector<TraversalOption>& allVCLApplicableTraversals() {
    // traversal not used but prevents usage of newton3
    static const std::vector<TraversalOption> v{TraversalOption::verletClusterCellsTraversal};
    return v;
  }

  std::vector<TraversalOption> getAllTraversals() override { return allVCLApplicableTraversals(); }

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
  void iteratePairwise(ParticleFunctor* f, Traversal* traversal, bool useNewton3 = false) {
    auto* traversalInterface = dynamic_cast<VerletClusterTraversalInterface<FullParticleCell<Particle>>*>(traversal);
    if (!traversalInterface) {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in VerletClusterCells::iteratePairwise");
    }

    traversalInterface->setVerletListPointer(&_clusterSize, &_neighborCellIds, &_neighborMatrixDim, &_neighborMatrix);
    if (needsRebuild()) {
      this->rebuild();
      traversalInterface->rebuildVerlet(_cellsPerDim, this->_cells, _boundingBoxes, this->getCutoff());
      _lastTraversalSig = traversalInterface->getSignature();
    }
    if (traversalInterface->getSignature() != _lastTraversalSig) {
      traversalInterface->rebuildVerlet(_cellsPerDim, this->_cells, _boundingBoxes, this->getCutoff());
      _lastTraversalSig = traversalInterface->getSignature();
    }

    traversal->initTraversal(this->_cells);
    traversalInterface->traverseCellPairs(this->_cells);
    traversal->endTraversal(this->_cells);

    // we iterated, so increase traversal counter
    _traversalsSinceLastRebuild++;
  }

  /**
   * @copydoc VerletLists::addParticle()
   */
  void addParticle(Particle& p) override {
    if (autopas::utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax())) {
      _isValid = false;
      this->_cells[0].resize(_dummyStarts[0]);
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
  void addHaloParticle(Particle& haloParticle) override {
    _isValid = false;
    if (autopas::utils::notInBox(haloParticle.getR(), this->getBoxMin(), this->getBoxMax())) {
      _haloInsertQueue.push_back(haloParticle);
    } else {
      utils::ExceptionHandler::exception(
          "VerletCluster: trying to add halo particle that is inside the bounding box.\n" + haloParticle.toString());
    }
  }

  /**
   * @copydoc VerletLists::deleteHaloParticles
   */
  void deleteHaloParticles() override { this->_cells.resize(_firstHaloClusterId); }

  /**
   * @copydoc VerletLists::updateContainer()
   */
  void updateContainer() override {
    AutoPasLog(debug, "updating container");
    _isValid = false;
    rebuild();
  }

  bool isContainerUpdateNeeded() override {
    if (not _isValid) return true;

    size_t i = 0;
    for (i = 0; i < _firstHaloClusterId; ++i) {
      auto dit = _dummyStarts.find(i);
      unsigned int end = _clusterSize;
      if (dit != _dummyStarts.end()) {
        end = dit->second;
      }
      for (unsigned int j = 0; j < end; ++j) {
        if (autopas::utils::notInBox(this->_cells[i]._particles[j].getR(), this->getBoxMin(), this->getBoxMax()))
          return true;
      }
    }
    for (; i < this->_cells.size(); ++i) {
      for (auto& p : this->_cells[i]._particles) {
        if (autopas::utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax())) return true;
      }
    }
    return false;
  }

  TraversalSelector<FullParticleCell<Particle>> generateTraversalSelector() override {
    return TraversalSelector<FullParticleCell<Particle>>(_cellsPerDim);
  }

  /**
   * Specifies whether the neighbor lists need to be rebuild.
   * @return true if the neighbor lists need to be rebuild, false otherwise
   */
  bool needsRebuild() {
    AutoPasLog(debug, "VerletLists: neighborlist is valid: {}", _isValid);
    // if the neighbor list is NOT valid or we have not rebuild for _rebuildFrequency steps
    return (not _isValid) or (_traversalsSinceLastRebuild >= _rebuildFrequency);
  }

  ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(new internal::ParticleIterator<Particle, FullParticleCell<Particle>>(
        &this->_cells, 0, &_cellBorderFlagManager, behavior));
  }

  ParticleIteratorWrapper<Particle> getRegionIterator(std::array<ParticleFloatType, 3> lowerCorner,
                                                      std::array<ParticleFloatType, 3> higherCorner,
                                                      IteratorBehavior behavior = IteratorBehavior::haloAndOwned,
                                                      bool incSearchRegion = false) override {
    std::vector<size_t> cellsOfInterest;
    for (size_t i = 0; i < _boundingBoxes.size(); ++i) {
      if (boxesOverlap(lowerCorner, higherCorner, _boundingBoxes[i])) cellsOfInterest.push_back(i);
    }
    return ParticleIteratorWrapper<Particle>(new internal::RegionParticleIterator<Particle, FullParticleCell<Particle>>(
        &this->_cells, lowerCorner, higherCorner, cellsOfInterest, &_cellBorderFlagManager, behavior));
  }

  /**
   * Deletes all Dummy Particles in the container
   */
  void deleteDummyParticles() override {
    for (auto& pair : _dummyStarts) {
      this->_cells[pair.first].resize(pair.second);
    }
    _isValid = false;
  }

 protected:
  /**
   * Recalculate grids and clusters,
   * build verlet lists and pad clusters.
   */
  void rebuild() {
    deleteDummyParticles();
    _boundingBoxes.clear();
    // get the dimensions and volumes of the box
    std::array<ParticleFloatType, 3> boxSize{};
    ParticleFloatType volume = 1.0;

    for (int d = 0; d < 3; ++d) {
      boxSize[d] = this->getBoxMax()[d] - this->getBoxMin()[d];
      volume *= boxSize[d];
    }

    // get all particles and clear clusters
    std::vector<Particle> invalidParticles;

    for (size_t i = 0; i < this->_cells.size(); ++i) {
      for (auto& p : this->_cells[i]._particles) {
        if (utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax())) {
          invalidParticles.push_back(p);

        } else {
          _haloInsertQueue.push_back(p);
        }
      }
      this->_cells[i].clear();
    }
    // estimate particle density
    ParticleFloatType density = (ParticleFloatType)invalidParticles.size() / volume;

    // guess optimal grid side length
    _gridSideLength = std::cbrt((ParticleFloatType)_clusterSize / density);

    _gridSideLengthReciprocal = 1 / _gridSideLength;

    // get cells per dimension
    size_t sizeGrid = 1;
    for (int d = 0; d < 2; d++) {
      _cellsPerDim[d] = static_cast<size_t>(std::ceil(boxSize[d] * _gridSideLengthReciprocal));

      sizeGrid *= _cellsPerDim[d];
    }
    _cellsPerDim[2] = static_cast<size_t>(1);

    // resize to number of grids
    this->_cells.resize(sizeGrid);
    _dummyStarts.clear();

    // put particles into grid cells
    for (size_t i = 0; i < invalidParticles.size(); ++i) {
      int index =
          (int)((invalidParticles[i].getR()[0] - this->getBoxMin()[0]) * _gridSideLengthReciprocal) +
          (int)((invalidParticles[i].getR()[1] - this->getBoxMin()[1]) * _gridSideLengthReciprocal) * _cellsPerDim[0];
      this->_cells[index].addParticle(invalidParticles[i]);
    }

    size_t numOwnClusters = 0;
    size_t numHaloClusters = 0;

    // sort by last dimension

    for (size_t i = 0; i < this->_cells.size(); ++i) {
      this->_cells[i].sortByDim(2);
      size_t numParticles = this->_cells[i].numParticles();
      if (numParticles)
        numOwnClusters += ((numParticles + _clusterSize - 1) / _clusterSize);
      else
        ++numOwnClusters;
    }

    _firstHaloClusterId = numOwnClusters;

    size_t sizeHaloGrid = 6;  // one halo cell per side of the domain to further cut later
    this->_cells.resize(numOwnClusters + sizeHaloGrid);

    // put halo particles into grid cells
    for (auto& haloParticle : _haloInsertQueue) {
      int index = _firstHaloClusterId;
      for (int i = 0; i < 3; ++i) {
        if (haloParticle.getR()[i] < this->getBoxMin()[i]) {
          index += 2 * i;
          break;
        } else if (haloParticle.getR()[i] > this->getBoxMax()[i]) {
          index += 2 * i + 1;
          break;
        }
      }
      this->_cells[index].addParticle(haloParticle);
    }

    // sort halo cells

    for (int i = 0; i < 2; ++i) {
      this->_cells[_firstHaloClusterId + i].sortByDim(i + 1);
      this->_cells[_firstHaloClusterId + i + 1].sortByDim(i + 1);
    }
    this->_cells[_firstHaloClusterId + 4].sortByDim(0);
    this->_cells[_firstHaloClusterId + 5].sortByDim(0);

    // calculate number of required halo cells
    for (size_t i = 0; i < sizeHaloGrid; ++i) {
      size_t numParticles = this->_cells[_firstHaloClusterId + i].numParticles();
      if (numParticles)
        numHaloClusters += ((numParticles + _clusterSize - 1) / _clusterSize);
      else
        ++numHaloClusters;
    }

    this->_cells.resize(numOwnClusters + numHaloClusters);
    _boundingBoxes.resize(this->_cells.size(),
                          {this->getBoxMax()[0] + 8 * this->getCutoff(), this->getBoxMax()[1] + 8 * this->getCutoff(),
                           this->getBoxMax()[2] + 8 * this->getCutoff(), this->getBoxMin()[0] - 8 * this->getCutoff(),
                           this->getBoxMin()[1] - 8 * this->getCutoff(), this->getBoxMin()[2] - 8 * this->getCutoff()});
    splitZ(0, sizeGrid);
    splitZ(_firstHaloClusterId, _firstHaloClusterId + sizeHaloGrid);

    _haloInsertQueue.clear();
    _traversalsSinceLastRebuild = 0;
    _isValid = true;
  }

 private:
  void splitZ(size_t start, size_t end) {
    size_t cid = end;
    for (size_t i = start; i < end; ++i) {
      size_t dummyStart = this->_cells[i].numParticles() % _clusterSize;
      // no cell with only dummy particles except when grid already empty
      if (dummyStart == 0 && !(this->_cells[i].numParticles() == 0)) {
        dummyStart = _clusterSize;
      }
      for (size_t clusterStart = dummyStart; clusterStart < this->_cells[i].numParticles();
           clusterStart += _clusterSize) {
        for (size_t pid = 0; pid < _clusterSize; ++pid) {
          Particle& p = this->_cells[i]._particles[clusterStart + pid];
          expandBoundingBox(_boundingBoxes[cid], p);
          this->_cells[cid].addParticle(p);
        }
        ++cid;
      }
      this->_cells[i].resize(dummyStart);
      for (auto& p : this->_cells[i]._particles) {
        expandBoundingBox(_boundingBoxes[i], p);
      }
      _dummyStarts[i] = dummyStart;
      for (size_t pid = dummyStart; pid < _clusterSize; ++pid) {
        // add dummy Particles with ID ULONG_MAX
        Particle dummyParticle = Particle({this->getBoxMax()[0] + 8 * this->getCutoff() +
                                               static_cast<typename Particle::ParticleFloatingPointType>(i),
                                           this->getBoxMax()[1] + 8 * this->getCutoff() +
                                               static_cast<typename Particle::ParticleFloatingPointType>(pid),
                                           this->getBoxMax()[2] + 8 * this->getCutoff()},
                                          {0., 0., 0.}, ULONG_MAX);
        this->_cells[i].addParticle(dummyParticle);
      }
    }
  }
  /**
   * Expands a bounding Box such the Particle is in it
   * @param box
   * @param p
   */
  void expandBoundingBox(std::array<typename Particle::ParticleFloatingPointType, 6>& box, Particle& p) {
    for (int i = 0; i < 3; ++i) {
      box[i] = std::min(box[i], p.getR()[i]);
      box[3 + i] = std::max(box[3 + i], p.getR()[i]);
    }
  }

  /**
   * Returns true if the two boxes are within distance
   * @param box1
   * @param box2
   * @param distance betwwen the boxes to return true
   * @return true if the boxes are overlapping
   */
  inline bool boxesOverlap(std::array<ParticleFloatType, 3>& box1lowerCorner,
                           std::array<ParticleFloatType, 3>& box1higherCorner,
                           std::array<typename Particle::ParticleFloatingPointType, 6>& box2) {
    for (int i = 0; i < 3; ++i) {
      if (box1lowerCorner[0 + i] > box2[3 + i] || box1higherCorner[i] < box2[0 + i]) return false;
    }
    return true;
  }

  /// first Cluster ID to be a Halo Cluster
  size_t _firstHaloClusterId;

  /// indices where dummy particles in the cells start
  std::unordered_map<size_t, size_t> _dummyStarts;

  /// Stores Halo Particles to be inserted
  std::vector<Particle> _haloInsertQueue;

  // number of particles in a cluster
  unsigned int _clusterSize;

  // id of neighbor clusters of a clusters
  std::vector<std::vector<size_t>> _neighborCellIds;

  size_t _neighborMatrixDim;
  utils::CudaDeviceVector<unsigned int> _neighborMatrix;

  // bounding boxes of all clusters (xmin,ymin,zmin,xmax,ymax,zmax) including skin and cutoff
  std::vector<std::array<typename Particle::ParticleFloatingPointType, 6>> _boundingBoxes;

  // side length of xy-grid and reciprocal
  ParticleFloatType _gridSideLength;
  ParticleFloatType _gridSideLengthReciprocal;

  // dimensions of grid
  std::array<size_t, 3> _cellsPerDim;

  /// how many pairwise traversals have been done since the last traversal
  unsigned int _traversalsSinceLastRebuild;

  /// specifies after how many pairwise traversals the neighbor list is to be
  /// rebuild
  unsigned int _rebuildFrequency;

  // specifies if the neighbor list is currently valid
  bool _isValid;

  /// Signature of the last Traversal to trigger rebuild when a new one is used
  std::tuple<TraversalOption, DataLayoutOption, bool> _lastTraversalSig;

  class VerletClusterCellsCellBorderAndFlagManager : public CellBorderAndFlagManager {
   public:
    VerletClusterCellsCellBorderAndFlagManager(size_t* firstHaloClusterId) { _firstHaloClusterId = firstHaloClusterId; }
    bool isHaloCell(size_t index1d) const override { return *_firstHaloClusterId <= index1d; }

    bool isOwningCell(size_t index1d) const override { return index1d < *_firstHaloClusterId; }

   private:
    size_t* _firstHaloClusterId;
  } _cellBorderFlagManager;
};

}  // namespace autopas
