/**
 * @file SlicedBlkBasedTraversal.h
 * @date 09 Apr 2020
 * @author henkel
 */

#pragma once

#include <algorithm>
#include <queue>
#include <vector>

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/utils/DataLayoutConverter.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the 3 dimensional slice (cuboid) traversal.
 *
 * The traversal cuts the simulation domain in one cuboid per thread. To achieve this
 * each dimension of the 3D-simulation domain is cut by the cubic root of threads.
 * The cuboids are assigned from the longest dimension towards the shortest.
 * The cuboids at the borders of the simulation domain are adapted to include possible
 * cut off due the flooring of the cubic root. Each cuboid (cellblock) is seperated
 * into 27 smaller cuboids (subcube), their size depending on the overlap length.
 * Each subcube forms a overlap-region with the corresponding subcubes from other cuboids.
 * Each overlap-region for the corresponding subcubes is locked by a thread accessing it,
 * leaving the largest subcube in the middle of the domain, with no necessary locking
 * for the corresponding thread. The lock is lifted as soon as the thread is finished
 * processing their subcube of the overlap-region.
 *
 * Improvements could include:
 * The queuing of locked overlap-regions and continuing to the next. The calculation
 * of threads necessary and unused. Allow optimization to make cuboids larger to shrink
 * the size of the outer subcubes and increase the size of the core subcube to optimal
 * relations.
 *
 * @tparam ParticleCell The type of cells.
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class SlicedBlkBasedTraversal : public CellPairTraversal<ParticleCell> {
 public:
  /**
   * Constructor of the sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  explicit SlicedBlkBasedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                   const double interactionLength, const std::array<double, 3> &cellLength)
      : CellPairTraversal<ParticleCell>(dims),
        _overlap{},
        _dimsPerLength{},
        _interactionLength(interactionLength),
        _cellLength(cellLength),
        _overlapAxis(),
        _cellBlockDimensions{},
        locks(),
        _dataLayoutConverter(pairwiseFunctor) {
    init(dims);
  }

  /**
   * Checks if the traversal is applicable to the current state of the domain.
   * @return true if the traversal can be applied.
   */
  bool isApplicable() const override {
    const bool atLeast27Cells = _cellBlockDimensions[0].front() > 2 and _cellBlockDimensions[1].front() > 2 and
                                _cellBlockDimensions[2].front() > 2;
    return not(dataLayout == DataLayoutOption::cuda) and atLeast27Cells;
  }

  /**
   * Load Data Layouts required for this Traversal if cells have been set through setCellsToTraverse().
   */
  void initTraversal() override {
    if (this->_cells) {
      auto &cells = *(this->_cells);
#ifdef AUTOPAS_OPENMP
      /// @todo find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
      for (size_t i = 0; i < cells.size(); ++i) {
        _dataLayoutConverter.loadDataLayout(cells[i]);
      }
    }
  }

  /**
   * Write Data to AoS if cells have been set through setCellsToTraverse().
   */
  void endTraversal() override {
    if (this->_cells) {
      auto &cells = *(this->_cells);
#ifdef AUTOPAS_OPENMP
      /// @todo find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
      for (size_t i = 0; i < cells.size(); ++i) {
        _dataLayoutConverter.storeDataLayout(cells[i]);
      }
    }
  }

 protected:
  /**
   * Resets the cell structure of the traversal.
   * @param dims
   */
  void init(const std::array<unsigned long, 3> &dims);

  /**
   * The main traversal of the sliced traversal.
   *
   * @copydetails C01BasedTraversal::c01Traversal()
   *
   * @tparam allCells Defines whether or not to iterate over all cells with the loop body given as argument. By default
   * (allCells=false) it will not iterate over all cells and instead skip the last few cells, because they will be
   * covered by the base step. If you plan to use the default base step of the traversal on this function, use
   * allCells=false, if you plan to just iterate over all cells, e.g., to iterate over verlet lists saved within the
   * cells, use allCells=true. For the sliced step if allCells is false, iteration will not occur over the last layer of
   * cells (for _overlap=1) (in x, y and z direction).
   */
  template <bool allCells = false, typename LoopBody>
  inline void slicedBlkTraversal(LoopBody &&loopBody);

  /**
   * Overlap of interacting cells. Array allows asymmetric cell sizes.
   */
  std::array<unsigned long, 3> _overlap;

 private:
  /**
   * Store ids of dimensions ordered by number of cells per dimensions.
   */
  std::array<int, 3> _dimsPerLength;

  /**
   * Interaction length (cutoff + skin).
   */
  double _interactionLength;

  /**
   * Cell length in CellBlock3D.
   */
  std::array<double, 3> _cellLength;

  /**
   * Overlap of interacting cells along the longest axis.
   */
  std::array<unsigned long, 3> _overlapAxis;

  /**
   * The number of cells per sliced dimension
   */
  std::array<std::vector<unsigned long>, 3> _cellBlockDimensions;
  std::vector<AutoPasLock> locks;

  /**
   * Data Layout Converter to be used with this traversal.
   */
  utils::DataLayoutConverter<PairwiseFunctor, dataLayout> _dataLayoutConverter;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void SlicedBlkBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::init(
    const std::array<unsigned long, 3> &dims) {
  using array3D = std::array<unsigned long, 3>;

  for (unsigned int d = 0; d < 3; d++) {
    _overlap[d] = std::ceil(_interactionLength / _cellLength[d]);
  }

  // order dimensions by descending length
  auto minMaxElem = std::minmax_element(this->_cellsPerDimension.begin(), this->_cellsPerDimension.end());
  _dimsPerLength[0] = (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.second);
  _dimsPerLength[2] = (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.first);
  _dimsPerLength[1] = 3 - (_dimsPerLength[0] + _dimsPerLength[2]);

  // order overlapAxis by descending length of dimensions
  _overlapAxis[0] = _overlap[_dimsPerLength[0]];
  _overlapAxis[1] = _overlap[_dimsPerLength[1]];
  _overlapAxis[2] = _overlap[_dimsPerLength[2]];

  // we need to calculate the size of the standard cube only once by one axis
  // if we want different shapes of boxes or rectangles, we need to adapt the 2D or 3D vectors
  auto numSlices = (size_t)autopas_get_max_threads();
  array3D numCellBlocks3D = {0, 0, 0};

  // calculate cubic root of the number of slices to find the amount each dimension needs to be split into
  auto numSlicesCqrt = std::cbrt(numSlices);

  numSlicesCqrt = floor(numSlicesCqrt);
  numCellBlocks3D[0] = numSlicesCqrt;
  numCellBlocks3D[1] = numSlicesCqrt;
  numCellBlocks3D[2] = numSlicesCqrt;

  // clear _cellBlockDimensions for each dimension
  _cellBlockDimensions[0].clear();
  _cellBlockDimensions[1].clear();
  _cellBlockDimensions[2].clear();

  // Insert the number of cells per dimensional slice by cellblock.
  // E.g. each cellblock is 3 cells long (in the longest dimension)
  //      then _cellBlockDimensions[0] would be [3,3,3,...] afterwards.
  for (int j = 0; j < 3; ++j) {
    _cellBlockDimensions[j].resize(numCellBlocks3D[j]);
    std::fill(_cellBlockDimensions[j].begin(), _cellBlockDimensions[j].end(),
              static_cast<unsigned long>(this->_cellsPerDimension[_dimsPerLength[j]] / numCellBlocks3D[j]));
  }

  // Calculate the rest of the cells per dimension, which where cutoff by possible floor of division of dimension.
  // Accounts for possible floor of the numSlicesCqrt
  array3D rest;
  for (int k = 0; k < 3; ++k) {
    rest[k] = this->_cellsPerDimension[_dimsPerLength[k]] - _cellBlockDimensions[k].size() * numSlicesCqrt;
    // Expand the cellBlocks bordering a dimensional end with rest, by the rest.
    if (rest[k] != 0) {
      _cellBlockDimensions[k].back() + rest[k];
    }
  }

  // Each cellblock should be at least 3x3x3 cells big,
  //   checking the first and last element of each dimension is enough. SECOND CHECK
  for (auto &cellBlockDimension : _cellBlockDimensions) {
    if (cellBlockDimension.front() < (unsigned long)3) {
      return;
    }
    if (cellBlockDimension.back() < (unsigned long)3) {
      return;
    }
  }

  // Decreasing last _cellBlockDimensions by _overlapAxis to account for the way we handle base cells?
  // No, do not do something like _sliceThickness.back() -= _overlapLongestAxis;
  // This is handled per subBlock.

  // We actually need only at maximum of 8 locks for only a few cellblocks at worst case at the same time.
  // But each CellBlock can have 27 (26) different lock locations.
  locks.resize(27 * (_cellBlockDimensions[0].size() + _cellBlockDimensions[1].size() + _cellBlockDimensions[2].size()));
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
template <bool allCells, typename LoopBody>
void SlicedBlkBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::slicedBlkTraversal(
    LoopBody &&loopBody) {
  using array3D = std::array<unsigned long, 3>;
  using subBlock = std::array<std::array<unsigned long, 2>, 3>;
  using subBlockArray = std::array<subBlock, 27>;

  // fill a vector with the starting coordinates lower left front corner {0,0,0} for each different cellblock
  //  and the upper right back corner {2,2,2}. So we can build a spanning vector for the cuboid.
  // We know each of those starting coordinates form a cubic mesh. -> Simple iteration is enough
  std::vector<std::array<array3D, 2>> _cellBlocks;  //+spanning vector of the cube as 2nd array

  _cellBlocks.reserve(_cellBlockDimensions[0].size() * _cellBlockDimensions[1].size() * _cellBlockDimensions[2].size());
  unsigned long accumulator_firstDim = 0;
  for (auto &it : _cellBlockDimensions[0]) {
    unsigned long accumulator_secondDim = 0;
    for (auto &jt : _cellBlockDimensions[1]) {
      unsigned long accumulator_thirdDim = 0;
      for (auto &kt : _cellBlockDimensions[2]) {
        _cellBlocks.emplace_back(std::array<array3D, 2>{
            array3D  // The 000 Edge of the Cube
            {accumulator_firstDim, accumulator_secondDim, accumulator_thirdDim},
            array3D{// The 222 Edge of the Cube
                    accumulator_firstDim + it, accumulator_secondDim + jt, accumulator_thirdDim + kt}});
        accumulator_thirdDim += kt;
      }

      accumulator_secondDim += jt;
    }

    accumulator_firstDim += it;
  }
  auto _threads = _cellBlocks.size();

#ifdef AUTOPAS_OPENMP
  // although every thread gets exactly one iteration (=cellblock) this is faster than a normal parallel region
#pragma omp parallel for schedule(static, 1) num_threads(_threads)
#endif

  for (int m = 0; m < _threads; ++m) {
    auto cellblock = _cellBlocks[m];

    /**
     * Splitting each cellblock = slice into its subdomains = subBlocks
     *            220      		     221      		    222
     *           / |                / |                / |
     *       120   |            121   |            122   |
     *      / |   210          / |   211          / |   212
     *  020   |  / |       021   |  / |       022   |  / |
     *   |   110   |        |   111   |        |   112   |
     *   |  / |   200       |  / |   201       |  / |   202
     *  010   |  /         011   |  /         012   |  /
     *   |   100            |   101            |   102
     *   |  /               |  /               |  /
     *  000................001................002
     *
     * Above model shows the numbering of each sub_block, be aware that they can be different in cellsize.
     * And the _cellBlocks bordering the end of a dimension, will have different sizes for subBlocks bordering a
     *      dimension end.
     * Only, but very often, if the simulation-dimension will not be perfectly splittable into #thread-amount of
     *      equally sized cubes.
     */
    // advanced subBlocks which carry knowledge about their location in space (i,j,k)
    // subBlocks accessing [x][0-2][0-1]
    subBlockArray subBlocks;

    // loop over the dimensions and build the subBlocks
    // Attention: subBlocks include overlap and are therefore bigger than the cellblock they are part of!
    int subBlockNumber = 0;
    for (unsigned long i = 0; i < 3; ++i) {
      for (unsigned long j = 0; j < 3; ++j) {
        for (unsigned long k = 0; k < 3; ++k) {
          // iterate over k first -> use for longest dimension
          if (k == 0) {
            subBlocks[subBlockNumber][0][0] = cellblock[0][0];
          } else if (k == 1) {
            subBlocks[subBlockNumber][0][0] = cellblock[0][0] + _overlapAxis[0];
          } else if (k == 2) {
            subBlocks[subBlockNumber][0][0] = cellblock[0][0] + cellblock[0][1] - _overlapAxis[0];
          }
          subBlocks[subBlockNumber][0][1] = k;

          // second longest dimension
          if (j == 0) {
            subBlocks[subBlockNumber][1][0] = cellblock[1][0];
          } else if (j == 1) {
            subBlocks[subBlockNumber][1][0] = cellblock[1][0] + _overlapAxis[1];
          } else if (j == 2) {
            subBlocks[subBlockNumber][1][0] = cellblock[1][0] + cellblock[1][1] - _overlapAxis[1];
          }
          subBlocks[subBlockNumber][1][1] = j;

          if (i == 0) {
            subBlocks[subBlockNumber][2][0] = cellblock[2][0];
          } else if (i == 1) {
            subBlocks[subBlockNumber][2][0] = cellblock[2][0] + _overlapAxis[2];
          } else if (i == 2) {
            subBlocks[subBlockNumber][2][0] = cellblock[2][0] + cellblock[2][1] - _overlapAxis[2];
          }
          subBlocks[subBlockNumber][2][1] = i;

          // one block finished building
          subBlockNumber++;
        }
      }
    }

    // std::queue<subBlock> subBlockQueue;      // For using openmp testlock
    // iterate over subBlocks
    for (auto &subBlock : subBlocks) {
      // Calculate the 000 cell of the subBlock to lock --> Will be lock for all the subBlocks
      // This cell is not in the Cellblock unless it is the 111 cell
      array3D subBlockSpanForLock;
      for (int i = 0; i < 3; ++i) {
        subBlockSpanForLock[i] = subBlock[i][0];
        if (subBlock[i][1] == 0) {
          subBlockSpanForLock[i] -= _overlapAxis[i];
        } else if (subBlock[i][1] == 2) {
          subBlockSpanForLock[i] += _overlapAxis[i];
        }
      }

      // Check if sub_block is locked || Lock first cell of the subBlock corresponding for the whole subBlock
      locks[subBlockSpanForLock[0]].lock();
      // if not -> queue.push(subBlock);

      // Calculate the actual subblock without overlap for iteration.
      // Reminder: We know we can iterate through it, because we locked the subBlock with overlap to other cellblocks
      //           and in our own cellblock there is no other thread writing into.
      array3D subblockStartIteration;
      for (int i = 0; i < 3; ++i) {
        subblockStartIteration[i] = subBlock[i][0];
        if (subBlock[i][1] == 0) {
          subblockStartIteration[i] += _overlapAxis[i];
        } else if (subBlock[i][1] == 2) {
          subblockStartIteration[i] -= _overlapAxis[i];
        }
      }
      // >1 Cell is locked in sub_block -> goto next sub_block and put this into the queue NOT POSSIBLE because _locked
      // = private loop only over the part of the subBlock, which is actually part of the cellBlock
      array3D subBlockEndIteration;
      for (int i = 0; i < 3; ++i) {
        subBlockEndIteration[i] = subBlock[i][0];
        if (subBlock[i][1] == 1) {
          subBlockEndIteration[i] -= _overlapAxis[i];
        } else if (subBlock[i][1] == 2) {
          subBlockEndIteration[i] += _overlapAxis[i];
        }
      }
      // loop over the cells in the subBlock
      for (int j = subblockStartIteration[0]; j < subBlockEndIteration[0]; ++j) {
        for (int k = subblockStartIteration[0]; k < subBlockEndIteration[0]; ++k) {
          for (int l = subblockStartIteration[0]; l < subBlockEndIteration[0]; ++l) {
            loopBody(j, k, l);
          }
        }
      }
      // finished? -> unlock sub_block + DO NOT QUEUE
      locks[subBlockSpanForLock[0]].unlock();
    }

    // iterate over queue if queue != empty
    // same as above, but waiting until empty
  }

}  // namespace autopas
}  // namespace autopas
