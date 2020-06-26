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
        _overlapAxis({0,0,0}),
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
   * The normal number of cells by dimension for each CellBlock.
   * Does not apply to CellBlocks at the end of a dimensional border unless dimension sizes are a natural cubic root.
   */
  std::array<unsigned long, 3> _numCellsInCellBlock;

  /**
   * Data Layout Converter to be used with this traversal.
   */
  utils::DataLayoutConverter<PairwiseFunctor, dataLayout> _dataLayoutConverter;

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
      for (auto i: _cellBlockDimensions[j]) {
        str += std::to_string(i) + " ";
      }
    }
    return str;
  }
  std::string debugHelperFunctionIntArrayToString(std::array<unsigned long, 3> intArray) {
    std::string str;
    for (auto i: intArray) {
      str += std::to_string(i) + " ";
    }
    return str;
  }
  std::string debugHelperFunctionIntArrayToString(std::array<int, 3> intArray) {
    std::string str;
    for (auto i: intArray) {
      str += std::to_string(i) + " ";
    }
    return str;
  }
  std::string debugHelperFunctionIntArrayToString(std::array<double, 3> intArray) {
    std::string str;
    for (auto i: intArray) {
      str += std::to_string(i) + " ";
    }
    return str;
  }

  std::string debugHelperFunctionSubBlockArrayToString(std::array<std::array<unsigned long, 2>, 3> subBlock) {
    std::string str;
    str += " d1: " +std::to_string(subBlock[0][0]);
    str += " d2: " +std::to_string(subBlock[1][0]);
    str += " d3: " +std::to_string(subBlock[2][0]);
    return str;
  }
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void SlicedBlkBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::init(
    const std::array<unsigned long, 3> &dims) {
  using array3D = std::array<unsigned long, 3>;
  AutoPasLog(debug, "Init?");

  for (unsigned int d = 0; d < 3; d++) {
    _overlap[d] = std::ceil(_interactionLength / _cellLength[d]);
  }
  AutoPasLog(debug, "_cellLength: " + debugHelperFunctionIntArrayToString(_cellLength));

  // order dimensions by descending length
  auto minMaxElem = std::minmax_element(this->_cellsPerDimension.begin(), this->_cellsPerDimension.end());
  _dimsPerLength[0] = (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.second);
  _dimsPerLength[2] = (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.first);
  _dimsPerLength[1] = 3 - (_dimsPerLength[0] + _dimsPerLength[2]);

  // order overlapAxis by descending length of dimensions
  _overlapAxis[0] = _overlap[_dimsPerLength[0]];
  _overlapAxis[1] = _overlap[_dimsPerLength[1]];
  _overlapAxis[2] = _overlap[_dimsPerLength[2]];

  AutoPasLog(debug, "_overlapAxis: " + debugHelperFunctionIntArrayToString(_overlapAxis));

  // init dimensions by descending length of dimensions
  _dims[0] = dims[_dimsPerLength[0]];
  _dims[1] = dims[_dimsPerLength[1]];
  _dims[2] = dims[_dimsPerLength[2]];

  AutoPasLog(debug, "_dimsPerLength: " + debugHelperFunctionIntArrayToString(_dimsPerLength));
  AutoPasLog(debug, "_dims: " + debugHelperFunctionIntArrayToString(_dims));

  // we need to calculate the size of the standard cube only once by one axis
  // if we want different shapes of boxes or rectangles, we need to adapt the 2D or 3D vectors
  auto numSlices = (size_t)autopas_get_max_threads();
  _numCellsInCellBlock = {0, 0, 0};

  // calculate cubic root of the number of slices to find the amount each dimension needs to be split into
  auto numSlicesCqrt = std::cbrt(numSlices);

  numSlicesCqrt = floor(numSlicesCqrt);
  _numCellsInCellBlock[0] = numSlicesCqrt;
  _numCellsInCellBlock[1] = numSlicesCqrt;
  _numCellsInCellBlock[2] = numSlicesCqrt;
  AutoPasLog(debug, "_numCellsInCellBlock: " + debugHelperFunctionIntArrayToString(_numCellsInCellBlock));

  // clear _cellBlockDimensions for each dimension
  _cellBlockDimensions[0].clear();
  _cellBlockDimensions[1].clear();
  _cellBlockDimensions[2].clear();

  // Insert the number of cells per dimensional slice by cellblock.
  // E.g. each cellblock is 3 cells long (in the longest dimension)
  //      then _cellBlockDimensions[0] would be [3,3,3,...] afterwards.
  for (int j = 0; j < 3; ++j) {
    _cellBlockDimensions[j].resize(_numCellsInCellBlock[j]);
    std::fill(_cellBlockDimensions[j].begin(), _cellBlockDimensions[j].end(),
              static_cast<unsigned long>(this->_cellsPerDimension[_dimsPerLength[j]] / _numCellsInCellBlock[j]));
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

  // We actually need only at maximum of 8 locks for only a few cellblocks at worst case at the same time.
  // But each CellBlock can have 27 (26) different lock locations.
  locks.resize(27 * (_cellBlockDimensions[0].size() + _cellBlockDimensions[1].size() + _cellBlockDimensions[2].size()));
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
template <typename LoopBody>
void SlicedBlkBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::slicedBlkTraversal(
    LoopBody &&loopBody) {
  using array3D = std::array<unsigned long, 3>;
  using subBlock = std::array<std::array<unsigned long, 2>, 3>;
  using subBlocksSingleCellblock = std::array<subBlock, 27>;
  using subBlocksAllCellblocks = std::vector<subBlocksSingleCellblock>;
  // fill a vector with the starting coordinates lower left front corner {0,0,0} for each different cellblock
  // and the upper right back corner {2,2,2}. So we can build a spanning vector for the cuboid.
  // We know each of those starting coordinates form a cubic mesh. -> Simple iteration is enough
  std::vector<std::array<array3D, 2>> _cellBlocks;  //+spanning vector of the cube as 2nd array

  AutoPasLog(debug, "_cellBlockDimensions: " + std::to_string(_cellBlockDimensions[0].size()) + " " +
                                             std::to_string(_cellBlockDimensions[1].size()) + " " +
                                             std::to_string(_cellBlockDimensions[2].size()) + " " );
  _cellBlocks.resize(_cellBlockDimensions[0].size() * _cellBlockDimensions[1].size() * _cellBlockDimensions[2].size());
  AutoPasLog(debug, "CellBlocks Array Size reserved: " +
                        std::to_string(_cellBlockDimensions[0].size() * _cellBlockDimensions[1].size() * _cellBlockDimensions[2].size()));

  unsigned long accumulator_firstDim = 0;
  unsigned long accumulator_secondDim;
  unsigned long accumulator_thirdDim;
  unsigned long cellBlockIterator = 0;

  AutoPasLog(debug, "_cellBlockDimensions at[0]: " +
                    std::to_string(_cellBlockDimensions[0][0]) + " " +
                    std::to_string(_cellBlockDimensions[1][0]) + " " +
                    std::to_string(_cellBlockDimensions[2][0]) + " " );

  AutoPasLog(debug, "CellBlockCreation Started:");
  for (auto &it : _cellBlockDimensions[0]) {
    accumulator_secondDim = 0;
    for (auto &jt : _cellBlockDimensions[1]) {
      accumulator_thirdDim = 0;
      for (auto &kt : _cellBlockDimensions[2]) {
        // TODO: FOR ALLCELLS we should subtract from _cellBlockDimensions[i].end the overlap[i]


//        _cellBlocks.push_back({
//            {accumulator_firstDim, accumulator_secondDim, accumulator_thirdDim}, // The 000 Edge of the Cube
//            {accumulator_firstDim + it, accumulator_secondDim + jt, accumulator_thirdDim + kt} // The 222 Edge of the Cube
//        });

        _cellBlocks[cellBlockIterator][0][0] = (unsigned long) accumulator_firstDim;
        _cellBlocks[cellBlockIterator][0][1] = (unsigned long) accumulator_secondDim;
        _cellBlocks[cellBlockIterator][0][2] = (unsigned long) accumulator_thirdDim;
        _cellBlocks[cellBlockIterator][1][0] = (unsigned long) accumulator_firstDim + it;
        _cellBlocks[cellBlockIterator][1][1] = (unsigned long) accumulator_secondDim + jt;
        _cellBlocks[cellBlockIterator][1][2] = (unsigned long) accumulator_thirdDim + kt;

        AutoPasLog(debug, "Cellblock START: " +
                          std::to_string(_cellBlocks[cellBlockIterator][0][0]) + " " +
                          std::to_string(_cellBlocks[cellBlockIterator][0][1]) + " " +
                          std::to_string(_cellBlocks[cellBlockIterator][0][2]) + " " );
        AutoPasLog(debug, "Cellblock END: " +
                          std::to_string(_cellBlocks[cellBlockIterator][1][0]) + " " +
                          std::to_string(_cellBlocks[cellBlockIterator][1][1]) + " " +
                          std::to_string(_cellBlocks[cellBlockIterator][1][2]) + " " );


        accumulator_thirdDim += kt;
        cellBlockIterator++;
      }

      accumulator_secondDim += jt;
    }

    accumulator_firstDim += it;
  }
  AutoPasLog(debug, "Finished Accumulators: " +
                    std::to_string(accumulator_firstDim) + " " +
                    std::to_string(accumulator_secondDim) + " " +
                    std::to_string(accumulator_thirdDim) + " " );

  auto _threads = _cellBlocks.size();

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
   * Above model shows the numbering of each sub_block, be aware that they can be different in cellsize.
   * And the _cellBlocks bordering the end of a dimension, will have different sizes for subBlocks bordering a
   * dimension end. Only, but very often, if the simulation-dimension will not be perfectly splittable into
   * thread-amount of equally sized cubes.
   *
   * Iteration of Subblocks always is in the Order: 000 - 100 - 200 - 010 - 110 - 210 - 020 - 120 - 220 - 001 - ...
   */

  subBlocksAllCellblocks _allSubBlocks;
  _allSubBlocks.resize(_cellBlocks.size());

#ifdef AUTOPAS_OPENMP
  // although every thread gets exactly one iteration (=cellblock) this is faster than a normal parallel region
#pragma omp parallel for schedule(static, 1) num_threads(_threads)
#endif

  for (unsigned long n = 0; n < _threads; ++n) {
    auto cellblock = _cellBlocks[n];

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
          // iterate over k first -> use for longest dimension
          if (k == 0) {
            currentSubBlock[0][0] = cellblock[0][0];
          } else if (k == 1) {
            currentSubBlock[0][0] = cellblock[0][0] + this->_overlapAxis[0];
          } else if (k == 2) {
            currentSubBlock[0][0] = cellblock[1][0] - this->_overlapAxis[0] + 1;
          }
          currentSubBlock[0][1] = k;

          // second longest dimension
          if (j == 0) {
            currentSubBlock[1][0] = cellblock[0][1];
          } else if (j == 1) {
            currentSubBlock[1][0] = cellblock[0][1] + this->_overlapAxis[1];
          } else if (j == 2) {
            currentSubBlock[1][0] = cellblock[1][1] - this->_overlapAxis[1] + 1;
          }
          currentSubBlock[1][1] = j;

          if (i == 0) {
            currentSubBlock[2][0] = cellblock[0][2];
          } else if (i == 1) {
            currentSubBlock[2][0] = cellblock[0][2] + this->_overlapAxis[2];
          } else if (i == 2) {
            currentSubBlock[2][0] = cellblock[1][2] - this->_overlapAxis[2] + 1;
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
  }

  for (unsigned long m = 0; m < _threads; ++m) {
    auto cellblock = _cellBlocks[m];

    AutoPasLog(debug, "SubBlock iteration for cellblock " + std::to_string(m) + " started");
    // Iterate over the subblocks in following order:
    /**
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


    for (auto &subBlock : _allSubBlocks[m]) {

      // Idea: Calculate the cell of the subBlock to lock -->
      //      Will be the lock for all corresponding subBlocks in other cellBlocks
      // Detail: All Subblocks at a 0 edge will subtract the overlap to reach a subBlock with 2
      //      Making all corresponding subBlocks in neighbouring cellBlocks lock the same subBlock starting cell.
      array3D _cellToLockForSubBlock;
      for (unsigned long i = 0; i < 3; ++i) {
        if (subBlock[i][1] == 0) {
          _cellToLockForSubBlock[i] = subBlock[i][0] - _overlapAxis[i];
        } else {
          _cellToLockForSubBlock[i] = subBlock[i][0];
        }
      }
      AutoPasLog(debug, "_cellToLockForSubBlock: " + debugHelperFunctionIntArrayToString(_cellToLockForSubBlock));

      // The cellblocks on the starting borders of the container would become negative,
      // but because they are unsigned they match to rather large indexes. This is not a problem, it is still unique.
      // Unless long is too short of an index.
      const unsigned long index =
          (_cellToLockForSubBlock[2] * _dims[1] + _cellToLockForSubBlock[1]) * _dims[0] + _cellToLockForSubBlock[0];
      // Check if sub_block is locked || Lock first cell of the subBlock neighbours
      // locks[index].lock();

      // build the opposite corner from the subBlock to allow iteration through the cells
        // for 0 add overlap, for 1 add cellblocklength - overlap, for 2 add overlap
      array3D _subBlockEndIteration;
      for (unsigned long i = 0; i < 3; ++i) {
        if (subBlock[i][1] == 0 || subBlock[i][1] == 2) {
          _subBlockEndIteration[i] = subBlock[i][0] + _overlapAxis[i] - 1;
        } else {
          // e.g. subBlockEnd 001 = subBlockStart from 001 + cellBlockLength = _cellBlockDimensions[i][m] = cellBlockEnd - cellBlockStart
          //      - 000 Length - 002 Length = 2*overLap[i] -1 to stay in the subBlock and compensate for overlap > 1
          _subBlockEndIteration[i] = subBlock[i][0] + cellblock[1][i] - cellblock[0][i] - 2*_overlapAxis[i];
        }
      }
      AutoPasLog(debug, "SubBlockLock " + std::to_string(index) +
                        " Dim0: " + std::to_string(subBlock[0][0]) +
                        " Dim1: " + std::to_string(subBlock[1][0]) +
                        " Dim2: " + std::to_string(subBlock[2][0]) +
                        " | Order: " + std::to_string(subBlock[0][1]) +
                        std::to_string(subBlock[1][1]) +
                        std::to_string(subBlock[2][1]));
       AutoPasLog(debug, "loop from: " + std::to_string(subBlock[0][0]) + ":" + std::to_string(_subBlockEndIteration[0]) +
                            " " + std::to_string(subBlock[1][0]) + ":" + std::to_string(_subBlockEndIteration[1]) +
                            " " + std::to_string(subBlock[2][0]) + ":" + std::to_string(_subBlockEndIteration[2]));

       // loop over the cells in the subBlock
      for (unsigned long j = subBlock[0][0]; j <= _subBlockEndIteration[0] and j < _dims[0] - _overlapAxis[0]; ++j) {
        for (unsigned long k = subBlock[1][0]; k <= _subBlockEndIteration[1] and k < _dims[1] - _overlapAxis[1]; ++k) {
          for (unsigned long l = subBlock[2][0]; l <= _subBlockEndIteration[2] and l < _dims[2] - _overlapAxis[2]; ++l ) {
            // AutoPasLog(debug, "loopBody: " + std::to_string(j) + " " + std::to_string(k) + " "+ std::to_string(l) + " ");
            array3D idArray = {};
            idArray[_dimsPerLength[0]] = j;
            idArray[_dimsPerLength[1]] = k;
            idArray[_dimsPerLength[2]] = l;
            // AutoPasLog(debug, "idArray: " + std::to_string(idArray[0]) + " " + std::to_string(idArray[1]) + " "+ std::to_string(idArray[2]) + " ");
            loopBody(idArray[0], idArray[1], idArray[2]);
          }
        }
      }
      AutoPasLog(debug, "SubBlock Id: " + std::to_string(subBlock[0][1]) +
                        std::to_string(subBlock[1][1]) +
                        std::to_string(subBlock[2][1]));
      // finished? -> unlock sub_block + DO NOT QUEUE
      // locks[index].unlock();
    }
    AutoPasLog(debug, "CELLBLOCK: " + std::to_string(m) + " finished!");
  }
}

}  // namespace autopas
