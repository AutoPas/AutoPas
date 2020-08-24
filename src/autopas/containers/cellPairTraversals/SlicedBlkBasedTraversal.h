/**
 * @file SlicedBlkBasedTraversal.h
 * @date 09 Apr 2020
 * @author henkel
 */

#pragma once

#include <algorithm>
#include <map>
#include <queue>
#include <vector>

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/utils/DataLayoutConverter.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

class array;
namespace autopas {

/**
 * This class provides the 3 dimensional slice (cuboid) traversal.
 *
 * The traversal cuts the simulation domain in one cuboid per thread. To achieve this
 * each dimension of the 3D-simulation domain is seperately cut by the cubic root of threads.
 * The cuboids at the borders of the simulation domain are adapted to include possible
 * cut off due the flooring of the cubic root. Each cuboid (block) is seperated
 * into 27 smaller cuboids (sub-block), their size depending on the overlap length.
 * Each block overlaps with the other blocks through the sub-blocks at their borders.
 * Therefore each sub-block besides those at the start, the end and the actual middle sub-block,
 * are sub-blocks of multiple blocks.
 * Each sub-block has its own unique lock. Each step needs to lock always 8 sub-blocks.
 * To avoid deadlocks, the locks are tested first (whereas the testing is locked by a single
 * master-lock) and if all could be set, the thread proceeds. If not successfull, it aborts
 * the current iteration and tries to iterate through the next sub-block and tries the current one later.
 *
 * The largest sub-block is the middle (111) sub-block in a block and increases the larger the
 * simulation domain becomes and the less threads can be used.
 *
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
   * Constructor of the sliced blk traversal.
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
        _cellBlockDimensions{},
        locks(),
        _overlapAxis({0, 0, 0}),
        _dataLayoutConverter(pairwiseFunctor) {
    init(dims);
  }

  /**
   * Checks if the traversal is applicable to the current state of the domain.
   *
   * @return true if the traversal can be applied.
   */
  bool isApplicable() const override {
    const bool atLeast27CellsForEachCellBlock = _cellBlockDimensions[0].front() >= _overlapAxis[0] * 2 + 1 and
                                                _cellBlockDimensions[1].front() >= _overlapAxis[1] * 2 + 1 and
                                                _cellBlockDimensions[2].front() >= _overlapAxis[2] * 2 + 1 and
                                                _cellBlockDimensions[0].back() >= _overlapAxis[0] * 2 + 1 and
                                                _cellBlockDimensions[1].back() >= _overlapAxis[1] * 2 + 1 and
                                                _cellBlockDimensions[2].back() >= _overlapAxis[2] * 2 + 1;
    const bool atLeast3CellsOnEachAxis = this->_cellsPerDimension[0] >= _overlapAxis[0] * 2 + 1 and
                                         this->_cellsPerDimension[1] >= _overlapAxis[1] * 2 + 1 and
                                         this->_cellsPerDimension[2] >= _overlapAxis[2] * 2 + 1;
    return not(dataLayout == DataLayoutOption::cuda) and atLeast27CellsForEachCellBlock and atLeast3CellsOnEachAxis;
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
   *
   * @tparam allCells Defines whether or not to iterate over all cells with the loop body given as argument. By default
   * (allCells=false) it will not iterate over all cells and instead skip the last few cells, because they will be
   * covered by the base step. If you plan to use the default base step of the traversal on this function, use
   * allCells=false, if you plan to just iterate over all cells, e.g., to iterate over verlet lists saved within the
   * cells, use allCells=true. For the sliced step if allCells is false, the last layer will not be iterated over,
   * but its corresponding sub-blocks will be locked.
   *
   * @param dims The dimensions of the simulation domain in cells.
   */
  void init(const std::array<unsigned long, 3> &dims);

  /**
   * allCells Defines whether or not to iterate over all cells with the loop body given as argument. By default
   * (allCells=false) it will not iterate over all cells and instead skip the last few cells, because they will be
   * covered by the base step. If you plan to use the default base step of the traversal on this function, use
   * allCells=false, if you plan to just iterate over all cells, e.g., to iterate over verlet lists saved within the
   * cells, use allCells=true. For the sliced step if allCells is false, the last layer will not be iterated over,
   * but its corresponding sub-blocks will be locked.
   */
  bool allCells = false;

  /**
   * The main traversal of the sliced blk traversal.
   *
   * @copydetails C01BasedTraversal::c01Traversal()
   */
  template <typename LoopBody>
  inline void slicedBlkTraversal(LoopBody and loopBody);

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
   * The number of cells per dimension.
   */
  std::array<unsigned long, 3> _dims;

  /**
   * The number of cells per sliced dimension.
   */
  std::array<std::vector<unsigned long>, 3> _cellBlockDimensions;
  std::vector<AutoPasLock> locks;

  /**
   * The Masterlock to avoid deadlocks while threads are locking sub-blocks.
   */
  AutoPasLock _masterlock;

  /**
   * The normal number of cells by dimension for each Block.
   * Does not apply to Blocks at the end of a dimensional border unless dimension sizes are a natural cubic root.
   */
  std::array<unsigned long, 3> _numCellsInBlock;

  /**
   * Data Layout Converter to be used with this traversal.
   */
  utils::DataLayoutConverter<PairwiseFunctor, dataLayout> _dataLayoutConverter;

  /**
   * A vector holding the two opposing corners of each block as array of coordinates in an array.
   */
  std::vector<std::array<std::array<unsigned long, 3>, 2>> blocks;

  /**
   * A vector holding all subblocks.
   * Each vector entry x corresponds to block at blocks[x]. Each entry holds an array of 27 sub blocks.
   * Each sub block has a starting corner and unique order identity (000-222) assigned as follows:
   *            220      	    221    	       222
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
   */
  std::vector<std::array<std::array<std::array<unsigned long, 2>, 3>, 27>> _allSubBlocks;

  /**
   * A map converting consecutive ints to the subBlock Order ints e.g. 000->0, 111-14, 201->11, ...
   */
  std::map<int, int> _subBlockOrderToSubBlockIndex;

  /**
   * A function mapping a vector of the corresponding subBlocks to the corresponding lock ids.
   */
  int lockToSubBlocks(int x, int y, int z) {
    // 0 => my lock and 1 => neighbour lock
    if (x == 0 and  y == 0 and  z == 0) return 0;
    else if (x == 1 and  y == 0 and  z == 0) return 100;
    else if (x == 0 and  y == 1 and  z == 0) return 10;
    else if (x == 0 and  y == 0 and  z == 1) return 1;

    else if (x == 1 and  y == 1 and  z == 0) return 110;
    else if (x == 1 and  y == 0 and  z == 1) return 101;
    else if (x == 0 and  y == 1 and  z == 1) return 11;
    else if (x == 1 and  y == 1 and  z == 1) return 111;
    else {
      AutoPasLog(error, "Requested not a SubBlock for Locking.");
    }
  }

  /**
   * An array holding the locks in order of iteration over the subblocks of a block of cells. Is the same for all
   * blocks/ threads.
   */
  std::array<std::array<unsigned long, 3>, 8> _cellBlockTraverseOrderByLocks;

  /**
   * A map converting the cellblock Coordinates (relative Coordinates iterated by Cellblock) to their index
   */
  std::map<std::array<unsigned long, 3>, unsigned long> _cellBlocksToIndex;
  /**
   * A map converting a cellblock index to their relative cellblock coordinates.
   */
  std::map<unsigned long, std::array<unsigned long, 3>> _indexToCellBlock;

  /**
   * Testing the locking of all neighbour sub-blocks of the given one, including the given one.
   *
   * @param blockNumber The number of the block/ thread which wants the sub-block to lock.
   * @param order The block unique order id of the sub-block. e.g. 000-222 as an array. = {0,0,0}-{2,2,2}
   * @return A boolean corresponding true if the locks could been set, or false otherwise.
   */
  bool locking2x2x2sub_blocks(unsigned long blockNumber, std::array<unsigned long, 3> order) {
    // IF LOCK POSSIBLE, LOCK & SET MASTER LOCK & return true
    // IF LOCK IMPOSSIBLE UNLOCK locked Locks & SET MASTER LOCK & return false
    bool LocksThatCanBeSet[8] = {false, false, false, false, false, false, false, false};
    int no_locksSet = 0;
    if(_masterlock.testlock()) {
      for (unsigned long x = 0; x < 2; ++x) {
        for (unsigned long y = 0; y < 2; ++y) {
          for (unsigned long z = 0; z < 2; ++z) {
            LocksThatCanBeSet[no_locksSet] =
                locks[uniqueLockCalculation(blockNumber, {x + order[0], y + order[1], z + order[2]})].testlock();
            no_locksSet++;
          }
        }
      }

      if (LocksThatCanBeSet[0] and  LocksThatCanBeSet[1] and  LocksThatCanBeSet[2] and  LocksThatCanBeSet[3] and
          LocksThatCanBeSet[4] and  LocksThatCanBeSet[5] and  LocksThatCanBeSet[6] and  LocksThatCanBeSet[7]) {
        _masterlock.unlock();
        return true;
      } else {
        no_locksSet = 0;
        for (unsigned long x = 0; x < 2; ++x) {
          for (unsigned long y = 0; y < 2; ++y) {
            for (unsigned long z = 0; z < 2; ++z) {
              if (LocksThatCanBeSet[no_locksSet]) {
                locks[uniqueLockCalculation(blockNumber, {x + order[0], y + order[1], z + order[2]})].unlock();
              }
              no_locksSet++;
            }
          }
        }
        _masterlock.unlock();
        return false;
      }
    }
  }

  /**
   * Unlocks all neigbor sub-blocks and the current sub-block.
   * No need to test these, as we always want to unlock after we finished a sub-block.
   * @param blockNumber The number of the block/ thread which wants the sub-block to lock.
   * @param order The block unique order id of the sub-block. e.g. 000-222 as an array. = {0,0,0}-{2,2,2}
   */
  void unlocking2x2x2sub_blocks(unsigned long blockNumber, std::array<unsigned long, 3> order) {
    for (unsigned long x = 0; x < 2; ++x) {
      for (unsigned long y = 0; y < 2; ++y) {
        for (unsigned long z = 0; z < 2; ++z) {
          locks[uniqueLockCalculation(blockNumber, {x + order[0], y + order[1], z + order[2]})].unlock();
        }
      }
    }
  }

  /**
   * Calculates a unique id of the sub-block for locking
   *
   * @param blockNumber The number of the block, which the sub-block is part of
   * @param order The order of the sub-block, which to lock.
   * @return The number for locking.
   */
  unsigned long uniqueLockCalculation(unsigned long blockNumber, std::array<unsigned long, 3> order) {
    std::array<std::array<unsigned long, 3>, 2> currentBlock = blocks[blockNumber];
    auto startingEdge = currentBlock[0];
    auto endingEdge = currentBlock[1];
    unsigned long coordinates[3] = {0, 0, 0};
    coordinates[2] = blockNumber % _cellBlockDimensions[2].size();
    coordinates[1] = ((blockNumber - coordinates[2]) / _cellBlockDimensions[2].size()) % _cellBlockDimensions[1].size();
    coordinates[0] = ((blockNumber - coordinates[2] - coordinates[1] * _cellBlockDimensions[2].size()) /
                      (_cellBlockDimensions[2].size() * _cellBlockDimensions[1].size()));

    return ((order[2] + coordinates[0] * 2) * (_cellBlockDimensions[2].size() * 2 + 1) +
            (order[1] + coordinates[1] * 2)) *
               (_cellBlockDimensions[1].size() * 2 + 1) +
           order[0] + coordinates[2] * 2;
  }
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void SlicedBlkBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::init(
    const std::array<unsigned long, 3> &dims) {
  using array3D = std::array<unsigned long, 3>;
  using subBlock = std::array<std::array<unsigned long, 2>, 3>;
  using subBlocksSingleCellblock = std::array<subBlock, 27>;
  using subBlocksAllCellblocks = std::vector<subBlocksSingleCellblock>;

  for (unsigned int d = 0; d < 3; d++) {
    _overlap[d] = std::ceil(_interactionLength / _cellLength[d]);
  }
  _overlapAxis = _overlap;

  // init dimensions by descending length of dimensions
  for (int l = 0; l < 3; ++l) {
    _dims[l] = dims[l];
  }

  auto max_threads = (size_t)autopas_get_max_threads();
  _numCellsInBlock = {0ul, 0ul, 0ul};

  // calculate cubic root of the number of slices to find the amount each dimension needs to be split into
  auto numSlicesCqrt = std::cbrt(max_threads);
  numSlicesCqrt = std::floor(numSlicesCqrt);
  for (int n = 0; n < 3; ++n) {
    _numCellsInBlock[n] = numSlicesCqrt;
  }

  // clear _cellBlockDimensions for each dimension
  _cellBlockDimensions[0].clear();
  _cellBlockDimensions[1].clear();
  _cellBlockDimensions[2].clear();

  unsigned long min_slice_length[3] = {0, 0, 0};
  unsigned long slice_length[3] = {0, 0, 0};
  for (int axis = 0; axis < 3; ++axis) {
    min_slice_length[axis] = _overlapAxis[axis] * 2 + 1;
  }

  array3D no_possible_slices = _numCellsInBlock;
  array3D no_slices = _numCellsInBlock;
  for (int axis = 0; axis < 3; ++axis) {
    no_possible_slices[axis] = std::floor(this->_cellsPerDimension[axis] / min_slice_length[axis]);
  }

  for (int axis = 0; axis < 3; ++axis) {
    if (_numCellsInBlock[axis] > no_possible_slices[axis]) {
      no_slices[axis] = no_possible_slices[axis];
    }
    slice_length[axis] = std::floor(this->_cellsPerDimension[axis] / no_slices[axis]);
  }

  unsigned long current_slice_length[3] = {0, 0, 0};

  for (int axis = 0; axis < 3; ++axis) {
    _cellBlockDimensions[axis].resize(no_slices[axis]);
    for (int i = 0; i < no_slices[axis]; ++i) {
      _cellBlockDimensions[axis][i] = slice_length[axis];
      current_slice_length[axis] += slice_length[axis];
    }
    // rest calculation
    if (current_slice_length[axis] < this->_cellsPerDimension[axis]) {
      _cellBlockDimensions[axis].back() += (this->_cellsPerDimension[axis] - current_slice_length[axis]);
    } else if (current_slice_length[axis] > this->_cellsPerDimension[axis]) {
      AutoPasLog(error, "Negative overlap or other errors might have occured.");
    }
    if (allCells) {
      // This is necessary because of the overlap handling.
      _cellBlockDimensions[axis].back() += _overlapAxis[axis];
    }
  }

  // Each block should be at least overlap[axis]*2+1 in each dimension long.
  //   Checking the first and last element of each dimension is enough. (Might be the same. That is ok.)
  for (int axis = 0; axis < 3; ++axis) {
    if (_cellBlockDimensions[axis].front() < (_overlapAxis[axis] * 2 + 1)) {
      return;
    }
    if (_cellBlockDimensions[axis].back() < (_overlapAxis[axis] * 2 + 1)) {
      return;
    }
  }

  // Mapping SubBlock Id to SubBlock Index
  _subBlockOrderToSubBlockIndex = {{0, 0},    {100, 1},  {200, 2},  {10, 3},   {110, 4},  {210, 5},  {20, 6},
                                   {120, 7},  {220, 8},  {1, 9},    {101, 10}, {201, 11}, {11, 12},  {111, 13},
                                   {211, 14}, {21, 15},  {121, 16}, {221, 17}, {2, 18},   {102, 19}, {202, 20},
                                   {12, 21},  {112, 22}, {212, 23}, {22, 24},  {122, 25}, {222, 26}};

  // SubBlock Traverse Order
  _cellBlockTraverseOrderByLocks[0] = {0, 0, 0};
  _cellBlockTraverseOrderByLocks[1] = {0, 0, 1};
  _cellBlockTraverseOrderByLocks[2] = {0, 1, 1};
  _cellBlockTraverseOrderByLocks[3] = {0, 1, 0};

  _cellBlockTraverseOrderByLocks[4] = {1, 1, 0};
  _cellBlockTraverseOrderByLocks[5] = {1, 0, 0};
  _cellBlockTraverseOrderByLocks[6] = {1, 0, 1};
  _cellBlockTraverseOrderByLocks[7] = {1, 1, 1};

  // We need 8 = locks per block + 1 for all sub-blocks bordering the end of the simulation domain
  locks.resize((_cellBlockDimensions[0].size() * 2 + 1) * (_cellBlockDimensions[1].size() * 2 + 1) *
               (_cellBlockDimensions[2].size() * 2 + 1));

  /**
   * CELLBLOCK GENERATION:
   *
   *  We fill a vector with the starting coordinates lower left front corner {0,0,0} for each different cellblock
   *  and the upper right back corner {2,2,2}. So we can build a spanning vector for the block.
   *  We know each of those starting coordinates form a cubic mesh. -> Simple iteration is enough
   *  Please consider that the Spanning Vector MUST have the same direction as the c08 cell handling direction!!!
   */

  AutoPasLog(debug, "_cellBlockDimensions Size: " + std::to_string(_cellBlockDimensions[0].size()) + " " +
                        std::to_string(_cellBlockDimensions[1].size()) + " " +
                        std::to_string(_cellBlockDimensions[2].size()) + " ");

  blocks.resize(_cellBlockDimensions[0].size() * _cellBlockDimensions[1].size() * _cellBlockDimensions[2].size());

  AutoPasLog(debug, "Amount of Blocks: " + std::to_string(blocks.size()));

  unsigned long acc_x = 0ul;
  unsigned long acc_y = 0ul;
  unsigned long acc_z = 0ul;
  unsigned long cellBlockIterator = 0ul;
  std::array<unsigned long, 3> cellBlockOrder = {0, 0, 0};
  unsigned long x = 0ul;
  unsigned long y = 0ul;
  unsigned long z = 0ul;
  unsigned long max_cellblocks =
      _cellBlockDimensions[0].size() * _cellBlockDimensions[1].size() * _cellBlockDimensions[2].size();

  // each step builds the previous cellblock spanning vector + the next cellblock starting point
  for (unsigned long xit = 0; xit < _cellBlockDimensions[0].size(); ++xit) {
    acc_y = 0;
    for (unsigned long yit = 0; yit < _cellBlockDimensions[1].size(); ++yit) {
      acc_z = 0;
      for (unsigned long zit = 0; zit < _cellBlockDimensions[2].size(); ++zit) {
        // for parallelization we calculate the cellBlockIterator from the block coordinates
        // cellBlockIterator = utils::ThreeDimensionalMapping::threeToOneD(xit, yit, zit, max_cellblocks);
        cellBlockIterator = (zit * _cellBlockDimensions[1].size() + yit) * _cellBlockDimensions[0].size() + xit;

        x = acc_x + xit;
        y = acc_y + yit;
        z = acc_z + zit;

        blocks[cellBlockIterator][0][0] = x;
        blocks[cellBlockIterator][0][1] = y;
        blocks[cellBlockIterator][0][2] = z;

        if (x + _cellBlockDimensions[0][xit] < _dims[0] - _overlapAxis[0]) {
          blocks[cellBlockIterator][1][0] = x + _cellBlockDimensions[0][xit];
        } else {
          blocks[cellBlockIterator][1][0] = x + _cellBlockDimensions[0][xit] - _overlapAxis[0];
        }

        if (y + _cellBlockDimensions[1][yit] < _dims[1] - _overlapAxis[1]) {
          blocks[cellBlockIterator][1][1] = y + _cellBlockDimensions[1][yit];
        } else {
          blocks[cellBlockIterator][1][1] = y + _cellBlockDimensions[1][yit] - _overlapAxis[1];
        }

        if (z + _cellBlockDimensions[2][zit] < _dims[2] - _overlapAxis[2]) {
          blocks[cellBlockIterator][1][2] = z + _cellBlockDimensions[2][zit];
        } else {
          blocks[cellBlockIterator][1][2] = z + _cellBlockDimensions[2][zit] - _overlapAxis[2];
        }

        cellBlockOrder = {xit, yit, zit};
        _cellBlocksToIndex[cellBlockOrder] = cellBlockIterator;
        _indexToCellBlock[cellBlockIterator] = cellBlockOrder;

        // cellBlockIterator++;
        acc_z += _cellBlockDimensions[2][zit];
      }
      acc_y += _cellBlockDimensions[1][yit];
    }
    acc_x += _cellBlockDimensions[0][xit];
  }

  /**
   * Splitting each cellblock = slice into its subdomains = subBlocks
   *            220      	   221    	      222
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
   * Above model shows the numbering of each sub_block, be aware that they can be different in cellsizes.
   * And the blocks bordering the end of a dimension, will have different sizes for subBlocks bordering a
   * dimension end. Very often, if the simulation-dimension will not be perfectly splittable into the
   * thread-amount of equally sized cubes.
   *
   * Each block including a 0 or a 2 is also a sub-block of another block. Unless the block borders the domain.
   * Each sub-block saves its order and its starting coordinates. The end coordinates can be calculated by order
   * and overlap.
   *
   * Generation of Subblocks always is in the Order: 000 - 100 - 200 - 010 - 110 - 210 - 020 - 120 - 220 - 001 - ...
   */

  _allSubBlocks.resize(blocks.size());
  auto _threads = blocks.size();
  if (_threads != max_threads) {
    AutoPasLog(debug, "Sliced Bulk traversal only using {} threads because the number of cells is too small.",
               _threads);
  }

  for (unsigned long n = 0; n < _threads; ++n) {
    auto &block = blocks[n];

    // subBlocks accessing: [0-27][0-2][0-1] -->
    //  [subBlockNumber][dimension][0: coordinates of starting cell of subblock, 1: 000 to 222 subblock description]
    subBlocksSingleCellblock _subBlocksSingleCellBlock;

    // loop over the dimensions and build the subBlocks
    int subBlockNumber = 0;
    for (unsigned long xAxis = 0; xAxis < 3; ++xAxis) {
      for (unsigned long yAxis = 0; yAxis < 3; ++yAxis) {
        for (unsigned long zAxis = 0; zAxis < 3; ++zAxis) {
          subBlock currentSubBlock;
          if (zAxis == 0) {
            currentSubBlock[0][0] = block[0][0];
          } else if (zAxis == 1) {
            currentSubBlock[0][0] = block[0][0] + this->_overlapAxis[0];
          } else if (zAxis == 2) {
            currentSubBlock[0][0] = block[1][0] + 1;  // is in the next block
          }
          currentSubBlock[0][1] = zAxis;

          if (yAxis == 0) {
            currentSubBlock[1][0] = block[0][1];
          } else if (yAxis == 1) {
            currentSubBlock[1][0] = block[0][1] + this->_overlapAxis[1];
          } else if (yAxis == 2) {
            currentSubBlock[1][0] = block[1][1] + 1;  // is in the next block
          }
          currentSubBlock[1][1] = yAxis;

          if (xAxis == 0) {
            currentSubBlock[2][0] = block[0][2];
          } else if (xAxis == 1) {
            currentSubBlock[2][0] = block[0][2] + this->_overlapAxis[2];
          } else if (xAxis == 2) {
            currentSubBlock[2][0] = block[1][2] + 1;  // is in the next block
          }
          currentSubBlock[2][1] = xAxis;
          _subBlocksSingleCellBlock[subBlockNumber] = currentSubBlock;

          // one block finished building
          subBlockNumber++;
        }
      }
    }
    // Make this by reference, if the compiler does not so automatically
    _allSubBlocks[n] = _subBlocksSingleCellBlock;
    /** SubBlock Order in _allSubBlocks:
     * 000 -> 100 -> 200 ->
     * 010 -> 110 -> 210 ->
     * 020 -> 120 -> 220 ->
     * 001 -> 101 -> 201 ->
     * 011 -> 111 -> 211 ->
     * 021 -> 121 -> 221 ->
     * 002 -> 102 -> 202 ->
     * 012 -> 112 -> 212 ->
     * 022 -> 122 -> 222
     */
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
template <typename LoopBody>
void SlicedBlkBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::slicedBlkTraversal(
    LoopBody and loopBody) {
  using array3D = std::array<unsigned long, 3>;
  using subBlock = std::array<std::array<unsigned long, 2>, 3>;
  using subBlocksOfSingleBlock = std::array<subBlock, 27>;
  using subBlocksAllCellblocks = std::vector<subBlocksOfSingleBlock>;
  auto _threads = blocks.size();

#ifdef AUTOPAS_OPENMP
// although every thread gets exactly one iteration (=cellblock) this is faster than a normal parallel region
#pragma omp parallel for schedule(static, 1) num_threads(_threads)
#endif

  for (unsigned long m = 0; m < _threads; ++m) {
    auto &block = blocks[m];

    subBlocksOfSingleBlock &currentSubBlocksInThisCellBlock = _allSubBlocks[m];

    std::queue<int> traversalOrderQueue[8];
    for (int n = 0; n < 8; ++n) {
      traversalOrderQueue->push(n);
    }

    // Iterate over own sub-blocks (excluding all with 2)
    while (!traversalOrderQueue->empty()) {
      // todo: maybe make an indicator if a broken traversalOrder was used e.g. >1000 Iterations per thread
      //  or tell scheduler to give thread less priority after x iterations of this
      int traversalOrder = traversalOrderQueue->front();
      traversalOrderQueue->pop();
      auto &order = _cellBlockTraverseOrderByLocks[traversalOrder];
      if (!locking2x2x2sub_blocks(m, order)) {
        traversalOrderQueue->push(traversalOrder);
      } else {
        int subBlockId = lockToSubBlocks(order[0], order[1], order[2]);
        subBlock currentSubBlock = currentSubBlocksInThisCellBlock[_subBlockOrderToSubBlockIndex[subBlockId]];

        // Create the last point of sub-block for spanning vector
        array3D _subBlockEndIteration;
        for (unsigned long i = 0; i < 3; ++i) {
          if (currentSubBlock[i][1] == 0) {
            _subBlockEndIteration[i] = currentSubBlock[i][0] + _overlapAxis[i] - 1;
          } else if (currentSubBlock[i][1] == 1) {
            _subBlockEndIteration[i] = block[1][i];
          }
        }

        // loop over the cells in the subBlock
        for (unsigned long j = currentSubBlock[0][0]; j <= _subBlockEndIteration[0] and j < _dims[0] - _overlapAxis[0];
             ++j) {
          for (unsigned long k = currentSubBlock[1][0];
               k <= _subBlockEndIteration[1] and k < _dims[1] - _overlapAxis[1]; ++k) {
            for (unsigned long l = currentSubBlock[2][0];
                 l <= _subBlockEndIteration[2] and l < _dims[2] - _overlapAxis[2]; ++l) {
              loopBody(j, k, l);
            }
          }
        }

        unlocking2x2x2sub_blocks(m, order);
      }
    }
  }
}

}  // namespace autopas
