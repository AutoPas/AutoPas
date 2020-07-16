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
 * each dimension of the 3D-simulation domain is cut by the cubic root of threads.
 * The cuboids are assigned from the longest dimension towards the shortest.
 * The cuboids at the borders of the simulation domain are adapted to include possible
 * cut off due the flooring of the cubic root. Each cuboid (cellblock) is seperated
 * into 27 smaller cuboids (subcube), their size depending on the overlap length.
 * Each subcube forms an overlap-region with the corresponding subcubes from other cuboids.
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
        _overlapAxis({0, 0, 0}),
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
    const bool atLeast27CellsForEachCellBlock = _cellBlockDimensions[0].front() > _overlapAxis[0] * 2 and
                                                _cellBlockDimensions[1].front() > _overlapAxis[1] * 2 and
                                                _cellBlockDimensions[2].front() > _overlapAxis[2] * 2 and
                                                _cellBlockDimensions[0].back() >= _overlapAxis[0] * 3 and
                                                _cellBlockDimensions[1].back() >= _overlapAxis[1] * 3 and
                                                _cellBlockDimensions[2].back() >= _overlapAxis[2] * 3;
    // @todo: if allcells == false -> _cellBlockDimensions[x].back >= _overlapAxis[x] * 4
    return not(dataLayout == DataLayoutOption::cuda) and atLeast27CellsForEachCellBlock;
  }

  /**
   * allCells Defines whether or not to iterate over all cells with the loop body given as argument. By default
   * (allCells=false) it will not iterate over all cells and instead skip the last few cells, because they will be
   * covered by the base step. If you plan to use the default base step of the traversal on this function, use
   * allCells=false, if you plan to just iterate over all cells, e.g., to iterate over verlet lists saved within the
   * cells, use allCells=true. For the sliced step if allCells is false, iteration will not occur over the last layer of
   * cells (for _overlap=1) (in x, y and z direction).
   */
  bool allCells = false;

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
   * The main traversal of the sliced blk traversal.
   *
   * @copydetails C01BasedTraversal::c01Traversal()
   */
  template <typename LoopBody>
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
   * The number of cells per dimension
   */
  std::array<unsigned long, 3> _dims;

  /**
   * The number of cells per sliced dimension
   */
  std::array<std::vector<unsigned long>, 3> _cellBlockDimensions;
  std::vector<AutoPasLock> locks;

  /**
   * The normal number of cells by dimension for each Block.
   * Does not apply to Blocks at the end of a dimensional border unless dimension sizes are a natural cubic root.
   */
  std::array<unsigned long, 3> _numCellsInCellBlock;

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
   * Each sub block has a starting corner and unique order identity (000-222).
   */
  std::vector<std::array<std::array<std::array<unsigned long, 2>, 3>, 27>> _allSubBlocks;

  /**
   * A map converting consecutive ints to the subBlock Order ints e.g. 000, 111, 201, ...
   */
  std::map<int, int> _subBlockBlockCoordinatesToSubBlockIndex;

  /**
   * A function mapping a vector of the corresponding subBlocks to the corresponding lock ids
   */
  std::vector<int> locksToSubBlocks(int x, int y, int z) {
    // 0 => my lock and 1 => neighbour lock
    // x,y,z == 0 does not include 111 as this does not need a lock!
    if (x == 0 && y == 0 && z == 0) return {000, 001, 100, 101, 010, 011, 110};
    if (x == 1 && y == 0 && z == 0) return {002, 102, 012, 112};
    if (x == 0 && y == 1 && z == 0) return {200, 201, 210, 211};
    if (x == 0 && y == 0 && z == 1) return {020, 021, 120, 121};

    if (x == 1 && y == 1 && z == 0) return {202, 212};
    if (x == 1 && y == 0 && z == 1) return {022, 122};
    if (x == 0 && y == 1 && z == 1) return {220, 221};
    if (x == 1 && y == 1 && z == 1) return {222};
  }

  /**
   * An array holding the locks in order of iteration over the subblocks of a block of cells. Is the same for all blocks.
   */
  std::array<std::array<int, 3>, 8> _cellBlockTraverseOrderByLocks;

  /**
   * A map converting the cellblock Coordinates (relative Coordinates iterated by Cellblock) to their index
   */
  std::map<std::array<unsigned long, 3>, unsigned long> _cellBlocksToIndex;
  /**
   * A map converting a cellblock index to their relative cellblock coordinates.
   */
  std::map<unsigned long, std::array<unsigned long, 3>> _indexToCellBlock;


  /**
   * Function returning the closest cellblock startingCell for a given cell
   * @tparam cell An array consisting of the x, y and z coordinates of the given cell.
   */
  std::array<unsigned long, 3> cellToCellblockStartingCell(std::array<unsigned long, 3> cell) {
    // calculating div for each dimension of the cell and the container
    std::array<ldiv_t, 3> d{};
    for (int j = 0; j < 3; ++j) {
      d[j] = (std::ldiv(cell[j], _numCellsInCellBlock[j]));
    }

    std::array<unsigned long, 3> returnvalue = {cell[0], cell[1], cell[2]};
    // check if last cellblock
    for (int i = 0; i < 3; ++i) {
      if (_cellBlockDimensions[i].size() == d[i].quot) {
        // calculate last cellblock starting point
        returnvalue[i] = _cellBlockDimensions[i].size() * _numCellsInCellBlock[i];
      }
    }
    // return std::array holding the StartingCell Coordinates
    return {returnvalue[0] - d[0].rem, returnvalue[1] - d[1].rem, returnvalue[2] - d[2].rem};
  }

  std::string debugHelperFunctionIntArrayToString(std::array<std::vector<unsigned long>, 3> _cellBlockDimensions) {
    std::string str;
    for (int j = 0; j < 3; ++j) {
      str += std::to_string(j) + ": ";
      for (auto i : _cellBlockDimensions[j]) {
        str += std::to_string(i) + " ";
      }
    }
    return str;
  }
  std::string debugHelperFunctionIntArrayToString(std::array<unsigned long, 3> intArray) {
    std::string str;
    for (auto i : intArray) {
      str += std::to_string(i) + " ";
    }
    return str;
  }
  std::string debugHelperFunctionIntArrayToString(std::array<int, 3> intArray) {
    std::string str;
    for (auto i : intArray) {
      str += std::to_string(i) + " ";
    }
    return str;
  }
  std::string debugHelperFunctionIntArrayToString(std::array<double, 3> intArray) {
    std::string str;
    for (auto i : intArray) {
      str += std::to_string(i) + " ";
    }
    return str;
  }
  std::string debugHelperFunctionSubBlockArrayToString(std::array<std::array<unsigned long, 2>, 3> subBlock) {
    std::string str;
    str += " d1: " + std::to_string(subBlock[0][0]);
    str += " d2: " + std::to_string(subBlock[1][0]);
    str += " d3: " + std::to_string(subBlock[2][0]);
    return str;
  }
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void SlicedBlkBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::init(
    const std::array<unsigned long, 3> &dims) {

  AutoPasLog(debug, "Init");
  using array3D = std::array<unsigned long, 3>;
  using subBlock = std::array<std::array<unsigned long, 2>, 3>;
  using subBlocksSingleCellblock = std::array<subBlock, 27>;
  using subBlocksAllCellblocks = std::vector<subBlocksSingleCellblock>;

  for (unsigned int d = 0; d < 3; d++) {
    _overlap[d] = std::ceil(_interactionLength / _cellLength[d]);
  }
  AutoPasLog(debug, "_cellLength: " + debugHelperFunctionIntArrayToString(_cellLength));

  // order overlapAxis by descending length of dimensions
  for (int i = 0; i < 3; ++i) {
    _overlapAxis[i] = _overlap[i];
  }
  AutoPasLog(debug, "_overlapAxis: " + debugHelperFunctionIntArrayToString(_overlapAxis));

  // init dimensions by descending length of dimensions
  for (int l = 0; l < 3; ++l) {
    _dims[l] = dims[l];
  }
  AutoPasLog(debug, "_dims: " + debugHelperFunctionIntArrayToString(_dims));

  // we need to calculate the size of the standard cube only once by one axis
  // if we want different shapes of boxes or rectangles, we need to adapt the 2D or 3D vectors
  auto numSlices = (size_t)autopas_get_max_threads();
  _numCellsInCellBlock = {0, 0, 0};

  // calculate cubic root of the number of slices to find the amount each dimension needs to be split into
  auto numSlicesCqrt = std::cbrt(numSlices);
  numSlicesCqrt = floor(numSlicesCqrt);
  for (int n = 0; n < 3; ++n) {
    _numCellsInCellBlock[n] = numSlicesCqrt;
  }
  AutoPasLog(debug, "_numCellsInCellBlock: " + debugHelperFunctionIntArrayToString(_numCellsInCellBlock));

  // clear _cellBlockDimensions for each dimension
  _cellBlockDimensions[0].clear();
  _cellBlockDimensions[1].clear();
  _cellBlockDimensions[2].clear();

  // Insert the number of cells per dimensional slice by cellblock.
  // E.g. each cellblock has length 3 in x axis direction
  //      then _cellBlockDimensions[0] = [3,3,3,...] afterwards.
  for (int j = 0; j < 3; ++j) {
    _cellBlockDimensions[j].resize(_numCellsInCellBlock[j]);
    std::fill(_cellBlockDimensions[j].begin(), _cellBlockDimensions[j].end(),
              static_cast<unsigned long>(this->_cellsPerDimension[j] / _numCellsInCellBlock[j]));
  }

  // Calculate the rest of the cells per dimension, which where cutoff by possible floor of division of dimension.
  // Accounts for possible floor of the numSlicesCqrt
  array3D rest;
  for (int k = 0; k < 3; ++k) {
    rest[k] = this->_cellsPerDimension[k] - (_cellBlockDimensions[k].size() * numSlicesCqrt);
    // Expand the cellBlocks bordering a dimensional end with rest, by the rest.
    if (rest[k] != 0) {
      _cellBlockDimensions[k].back() + rest[k];
    }
    // TODO: If rest[k] + _cellBlockDimensions[k].back - overlapAxis[k] < 3
    // then add rest[k] + cBD[k].back to _cellBlockDimensions[k].secondToLast and delete back element
    // right now it would not be applicable (the comparison towards the .back of the vector in applicable function)
  }
  AutoPasLog(debug, "_cellBlockDimensions: " + debugHelperFunctionIntArrayToString(_cellBlockDimensions));

  // Each cellblock should be at least 3x3x3 cells big,
  //   checking the first and last element of each dimension is enough.
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

  // We need only a maximum of 3 locks per cellblocks at worst case at the same time.
  locks.resize(_cellBlockDimensions[0].size() * _cellBlockDimensions[1].size() * _cellBlockDimensions[2].size() * 3);

  _subBlockBlockCoordinatesToSubBlockIndex = {
      {000, 0},  {100, 1},  {200, 2},  {010, 3},  {110, 4},  {210, 5},  {020, 6},  {120, 7},  {220, 8},
      {001, 9},  {101, 10}, {201, 11}, {011, 12}, {111, 13}, {211, 14}, {021, 15}, {121, 16}, {221, 17},
      {002, 18}, {102, 19}, {202, 20}, {012, 21}, {112, 22}, {212, 23}, {022, 24}, {122, 25}, {222, 26}};

  // INIT locking steps, minimizing lock changes per step to max. 1 change
  _cellBlockTraverseOrderByLocks[0] = {0, 0, 0};
  _cellBlockTraverseOrderByLocks[1] = {1, 0, 0};
  _cellBlockTraverseOrderByLocks[2] = {1, 1, 0};
  _cellBlockTraverseOrderByLocks[3] = {0, 1, 0};

  _cellBlockTraverseOrderByLocks[4] = {0, 1, 1};
  _cellBlockTraverseOrderByLocks[5] = {0, 0, 1};
  _cellBlockTraverseOrderByLocks[6] = {1, 0, 1};
  _cellBlockTraverseOrderByLocks[7] = {1, 1, 1};

//  for (int m = 0; m < 2; ++m) {
//    for (int i = 0; i < 2; ++i) {
//      for (int j = 0; j < 2; ++j) {
//        _cellBlockTraverseOrderByLocks[(m-1+i*2+j*4)+1] = {m,i,j};
//      }
//    }
//  }

  AutoPasLog(debug, "_cellBlockDimensions Size: " + std::to_string(_cellBlockDimensions[0].size()) + " " +
                    std::to_string(_cellBlockDimensions[1].size()) + " " +
                    std::to_string(_cellBlockDimensions[2].size()) + " ");
  blocks.resize(_cellBlockDimensions[0].size() * _cellBlockDimensions[1].size() * _cellBlockDimensions[2].size());
  AutoPasLog(debug, "CellBlocks Array Size reserved: " +
                    std::to_string(_cellBlockDimensions[0].size() * _cellBlockDimensions[1].size() *
                                   _cellBlockDimensions[2].size()));



  // fill a vector with the starting coordinates lower left front corner {0,0,0} for each different cellblock
  // and the upper right back corner {2,2,2}. So we can build a spanning vector for the block.
  // We know each of those starting coordinates form a cubic mesh. -> Simple iteration is enough
  // Please consider that the Spanning Vector MUST have the same direction as the c08 cell handling direction
  // TODO: Take ^^ into account and remove the second condition in the for loops around loopbody

  unsigned long acc_x = 0ul;
  unsigned long acc_y = 0ul;
  unsigned long acc_z = 0ul;
  unsigned long cellBlockIterator = 0ul;
  std::array<unsigned long, 3> cellBlockOrder = {0,0,0};
  unsigned long x = 0ul;
  unsigned long y = 0ul;
  unsigned long z = 0ul;
  unsigned long max_cellblocks = _cellBlockDimensions[0].size() * _cellBlockDimensions[1].size() * _cellBlockDimensions[2].size();

  // each step builds the previous cellblock spanning vector + the next cellblock starting point
  for (unsigned long xit = 0; xit < _cellBlockDimensions[0].size(); ++xit) {

    acc_y = 0;
    for (unsigned long yit = 0; yit < _cellBlockDimensions[1].size(); ++yit) {

      acc_z = 0;
      for (unsigned long zit = 0; zit < _cellBlockDimensions[2].size(); ++zit) {
        // for parallelization we calculate the cellBlockIterator from the block coordinates
        cellBlockIterator = zit + yit * _cellBlockDimensions[2].size() + xit * _cellBlockDimensions[1].size();

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

        cellBlockOrder = {xit,yit,zit};
        _cellBlocksToIndex[cellBlockOrder] = cellBlockIterator;
        _indexToCellBlock[cellBlockIterator] = cellBlockOrder;

        // Debug Output
        AutoPasLog(debug, "BLOCK " + std::to_string(cellBlockIterator) + " START: " +
                          std::to_string(blocks[cellBlockIterator][0][0]) + " " +
                          std::to_string(blocks[cellBlockIterator][0][1]) + " " +
                          std::to_string(blocks[cellBlockIterator][0][2]) +
                          "| END: " + std::to_string(blocks[cellBlockIterator][1][0]) + " " +
                          std::to_string(blocks[cellBlockIterator][1][1]) + " " +
                          std::to_string(blocks[cellBlockIterator][1][2]) + " ");

        // cellBlockIterator++;
        acc_z += _cellBlockDimensions[2][zit];
      }
      acc_y += _cellBlockDimensions[1][yit];
    }
    acc_x += _cellBlockDimensions[0][xit];
  }
  AutoPasLog(debug, "Finished Accumulators: " + std::to_string(acc_x) + " " +
                    std::to_string(acc_y) + " " + std::to_string(acc_z) + " ");

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
   * dimension end. Only, but very often, if the simulation-dimension will not be perfectly splittable into the
   * thread-amount of equally sized cubes.
   *
   * Generation of Subblocks always is in the Order: 000 - 100 - 200 - 010 - 110 - 210 - 020 - 120 - 220 - 001 - ...
   */

  _allSubBlocks.resize(blocks.size());
  auto _threads = blocks.size();

#ifdef AUTOPAS_OPENMP
  // although every thread gets exactly one iteration (=cellblock) this is faster than a normal parallel region
#pragma omp parallel for schedule(static, 1) num_threads(_threads)
#endif

  for (unsigned long n = 0; n < _threads; ++n) {
    auto &block = blocks[n];

    // subBlocks accessing: [0-27][0-2][0-1] -->
    //  [subBlockNumber][dimension][0: coordinates of starting cell of subblock, 1: 000 to 222 subblock description]
    subBlocksSingleCellblock _subBlocksSingleCellBlock;
    AutoPasLog(debug, "SubBlock creation for cellBlock " + std::to_string(n) + " started");

    // loop over the dimensions and build the subBlocks

    int subBlockNumber = 0;
    for (unsigned long i = 0; i < 3; ++i) {
      for (unsigned long j = 0; j < 3; ++j) {
        for (unsigned long k = 0; k < 3; ++k) {
          subBlock currentSubBlock;
          if (k == 0) {
            currentSubBlock[0][0] = block[0][0];
          } else if (k == 1) {
            currentSubBlock[0][0] = block[0][0] + this->_overlapAxis[0];
          } else if (k == 2) {
            currentSubBlock[0][0] = block[1][0] - this->_overlapAxis[0] + 1;
          }
          currentSubBlock[0][1] = k;

          // second longest dimension
          if (j == 0) {
            currentSubBlock[1][0] = block[0][1];
          } else if (j == 1) {
            currentSubBlock[1][0] = block[0][1] + this->_overlapAxis[1];
          } else if (j == 2) {
            currentSubBlock[1][0] = block[1][1] - this->_overlapAxis[1] + 1;
          }
          currentSubBlock[1][1] = j;

          if (i == 0) {
            currentSubBlock[2][0] = block[0][2];
          } else if (i == 1) {
            currentSubBlock[2][0] = block[0][2] + this->_overlapAxis[2];
          } else if (i == 2) {
            currentSubBlock[2][0] = block[1][2] - this->_overlapAxis[2] + 1;
          }
          currentSubBlock[2][1] = i;
          _subBlocksSingleCellBlock[subBlockNumber] = currentSubBlock;

            AutoPasLog(debug, "SB:" + std::to_string(subBlockNumber) +
                              " | Dim0: " + std::to_string(_subBlocksSingleCellBlock[subBlockNumber][0][0]) +
                              " Dim1: " + std::to_string(_subBlocksSingleCellBlock[subBlockNumber][1][0]) +
                              " Dim2: " + std::to_string(_subBlocksSingleCellBlock[subBlockNumber][2][0]) +
                              " | Order: " + std::to_string(_subBlocksSingleCellBlock[subBlockNumber][0][1]) +
                              std::to_string(_subBlocksSingleCellBlock[subBlockNumber][1][1]) +
                              std::to_string(_subBlocksSingleCellBlock[subBlockNumber][2][1]));
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
    LoopBody &&loopBody) {
  using array3D = std::array<unsigned long, 3>;
  using subBlock = std::array<std::array<unsigned long, 2>, 3>;
  using subBlocksSingleCellblock = std::array<subBlock, 27>;
  using subBlocksAllCellblocks = std::vector<subBlocksSingleCellblock>;

  auto _threads = blocks.size();

  for (unsigned long m = 0; m < _threads; ++m) {
    auto &block = blocks[m];

    AutoPasLog(debug, "Iteration for BLOCK " + std::to_string(m) + " started");

    // Calculate my_locks
    unsigned long my_locks[3] = {};
    unsigned long nb_locks[3] = {};

    for (unsigned long n = 0ul; n < 3ul; ++n) {                    // if the block exists, their own locks do too
      my_locks[n] = m + n;
    }

    array3D cellBlockRelativeCoordinates = _indexToCellBlock[m];
    for (unsigned long axis = 0ul; axis < 3ul; ++axis) {
      cellBlockRelativeCoordinates[axis]++;
      AutoPasLog(debug, "_cellBlocksToIndex.count(cellBlockRelativeCoordinates): " + std::to_string(_cellBlocksToIndex.count(cellBlockRelativeCoordinates)));
      if (_cellBlocksToIndex[cellBlockRelativeCoordinates]) {             // check if the key exists
        nb_locks[axis] = _cellBlocksToIndex[cellBlockRelativeCoordinates] + axis;
      } else {
        nb_locks[axis] = 0ul - 1ul;
      }

      cellBlockRelativeCoordinates[axis]--;
    }

    // Iterate over 111 --> no locks needed
    subBlock middleSubBlock = _allSubBlocks[m][13];
    array3D _middlesubBlockEnd;
    for (unsigned long i = 0; i < 3; ++i) {
      _middlesubBlockEnd[i] = middleSubBlock[i][0] + block[1][i] - block[0][i] - 2 * _overlapAxis[i];
    }

    for (unsigned long j = middleSubBlock[0][0]; j <= _middlesubBlockEnd[0]; ++j) {
      for (unsigned long k = middleSubBlock[1][0]; k <= _middlesubBlockEnd[1]; ++k) {
        for (unsigned long l = middleSubBlock[2][0]; l <= _middlesubBlockEnd[2]; ++l) {
          loopBody(j,k,l);
        }
      }
    }

    // lock corresponding locks
    for (int traversalOrder = 0; traversalOrder < 8; ++traversalOrder) {  // Iterate over traverse array

      for (int i = 0; i < 3; ++i) {                                       // Iterate over locks
        if (_cellBlockTraverseOrderByLocks[traversalOrder][i] == 0) {     // my lock or neigbour lock?
          if (my_locks[i] != 0 - 1) {                                     // no need to lock at borders
            locks[my_locks[i]].lock();                                    // lock necessary my locks
          }
        } else {                        // == _cellBlockTraverseOrderByLocks[traversalOrder][i] == 1
          if (nb_locks[i] != 0 - 1) {   // no need to lock at borders
            locks[nb_locks[i]].lock();  // lock necessary nb locks
          }
        }
      }
      subBlocksSingleCellblock &currentSubBlocksInThisCellBlock = _allSubBlocks[m];

      // Iterate over corrsponding subblocks
      // call subBlockLoop for subBlock
      subBlock currentSubBlock;
      auto &order = _cellBlockTraverseOrderByLocks[traversalOrder];
//        AutoPasLog(debug, "Order: " + std::to_string(order[0]) + std::to_string(order[1]) + std::to_string(order[2]));
        for (auto &lockedSubBlocks : locksToSubBlocks(order[0], order[1], order[2])) {
          currentSubBlock = currentSubBlocksInThisCellBlock[_subBlockBlockCoordinatesToSubBlockIndex[lockedSubBlocks]];
          // build the opposite corner from the subBlock to allow iteration through the cells
          // for 0 add overlap, for 1 add cellblocklength - overlap, for 2 add overlap

          AutoPasLog(debug, "LockedSubBlocks: " + std::to_string(lockedSubBlocks) + " | " + "currentSubBlock: "
                                + std::to_string(currentSubBlock[0][1]) + " " + std::to_string(currentSubBlock[1][1]) + " " + std::to_string(currentSubBlock[2][1]));

          array3D _subBlockEndIteration;
          for (unsigned long i = 0; i < 3; ++i) {
            if (currentSubBlock[i][1] == 0 || currentSubBlock[i][1] == 2) {
              _subBlockEndIteration[i] = currentSubBlock[i][0] + _overlapAxis[i] - 1;
            } else {
              // e.g. subBlockEnd 001 = subBlockStart from 001 + cellBlockLength = _cellBlockDimensions[i][m] =
              // cellBlockEnd - cellBlockStart
              //      - 000 Length - 002 Length = 2*overLap[i] -1 to stay in the subBlock and compensate for overlap > 1
              _subBlockEndIteration[i] =
                  currentSubBlock[i][0] + block[1][i] - block[0][i] - 2 * _overlapAxis[i];
            }
          }

          // loop over the cells in the subBlock
          for (unsigned long j = currentSubBlock[0][0];
               j <= _subBlockEndIteration[0] and j < _dims[0] - _overlapAxis[0]; ++j) {
            for (unsigned long k = currentSubBlock[1][0];
                 k <= _subBlockEndIteration[1] and k < _dims[1] - _overlapAxis[1]; ++k) {
              for (unsigned long l = currentSubBlock[2][0];
                   l <= _subBlockEndIteration[2] and l < _dims[2] - _overlapAxis[2]; ++l) {
                loopBody(j,k,l);
              }
            }
          }
        }


      // Unlock corresponding locks
      for (int i = 0; i < 3; ++i) {                                     // Iterate over locks
        if (_cellBlockTraverseOrderByLocks[traversalOrder][i] == 0) {   // my lock or neigbour lock?
          if (my_locks[i] != 0 - 1) {                                   // no need to unlock at borders
            locks[my_locks[i]].unlock();                                              // lock necessary my locks
          }
        } else {                          // == _cellBlockTraverseOrderByLocks[traversalOrder][i] == 1
          if (nb_locks[i] != 0 - 1) {     // no need to unlock at borders
            locks[nb_locks[i]].unlock();  // lock necessary nb locks
          }
        }
      }
    }

    // unlock corresponding SINGLE lock (not unlocking all to keep subblock execution order compared to other threads)

    AutoPasLog(debug, "CELLBLOCK: " + std::to_string(m) + " finished!");
  }
}

}  // namespace autopas
