/**
 * @file C04hcp.h
 * @author sabrinakrallmann
 * @date 30.03.2020
 */

#pragma once

using namespace std;

namespace autopas{



/**
 * This class provides the c04 hcp traversal.
 *
 * The traversal uses the c04 base step performed on every single cell. Since
 * these steps overlap a domain coloring with four colors is applied.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>

class C04HCP : public C08BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
                     public LinkedCellTraversalInterface<ParticleCell> {

  public:

  /**
   * Constructor of c04hcp
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length.
   * @param cellLength cell length.
   */
  C04HCP(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
         const double interactionLength, const std::array<double, 3> &cellLength)
          : C08BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                     interactionLength, cellLength),
            _cellHandler(pairwiseFunctor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap),
            _end(utils::ArrayMath::subScalar(utils::ArrayUtils::static_cast_array<long>(this->_cellsPerDimension), 1l)){
  }

  void traverseParticlePairs() override;

  TraversalOption getTraversalType() const override { return TraversalOption::c04HCP; }

  DataLayoutOption getDataLayout() const override { return dataLayout; }

  bool getUseNewton3() const override { return useNewton3; }

    bool isApplicable() const override { //TODO: Individual Restrictions
      if (dataLayout == DataLayoutOption::cuda) {
        return false;
      }
      const double minLength = *std::min_element(this->_cellLength.cbegin(), this->_cellLength.cend());
      const unsigned long minDim = *std::min_element(this->_cellsPerDimension.cbegin(), this->_cellsPerDimension.cend());

      return minLength >= this->_interactionLength and minDim > 3;
    }

  private:

  void traverseSingleColor(std::vector<ParticleCell> &cells, int color);

  void processBasePack6(std::vector<ParticleCell> &cells, const std::array<long, 3> &base3DIndex);

  C08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> _cellHandler;

  const std::array<long, 3> _end;

};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void C04HCP<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::processBasePack6(
  std::vector<ParticleCell> &cells, const std::array<long, 3> &base3DIndex) {

  using utils::ThreeDimensionalMapping::threeToOneD;
  std::array<long, 3> index;
  const std::array<long, 3> signedDims = utils::ArrayUtils::static_cast_array<long>(this->_cellsPerDimension);

  for(long z = 0; z < 3; ++z){ //go through the six cells
    for(long x = 0; x < 2; ++x){
      index[0] = base3DIndex[0] + x;
      index[1] = base3DIndex[1];
      index[2] = base3DIndex[2] + z;

      bool isIn = true;
      for (int d = 0; d < 3; ++d){
        isIn &= (index[d] >= 0l) and (index[d] <= (_end[d] - this->_overlap[d])); //prevent using overlapping cells
      }

      if(isIn){ //skip cells outside radius
        const unsigned long ulIndex = threeToOneD(index, signedDims);
        _cellHandler.processBaseCell(cells, ulIndex);
      }

    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void C04HCP<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseParticlePairs() {
  auto &cells = *(this->_cells);
#if defined(AUTOPAS_OPENMP)
#pragma omp parallel
#endif
  {
    for (int color = 0; color < 4; ++color) {
      traverseSingleColor(cells, color);

#if defined(AUTOPAS_OPENMP)
      if (color < 3) {
#pragma omp barrier
      }
#endif
    }
  }  // close parallel region
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void C04HCP<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseSingleColor(
        std::vector<ParticleCell> &cells, int color) {

  // determine a starting point of one of the grids
  std::array<long, 3> startOfThisColor{};

  // only need to shift z-dimension, because in each z dimension the colors behave like the 0-color in the z = 0 dimension
  switch (color) {
    case 0:
      startOfThisColor = {0l, 0l, 0l};
      break;
    case 1:
      startOfThisColor = {-4l, 0l, 1l};
      break;
    case 2:
      startOfThisColor = {-4l, 0l, -2l};
      break;
    case 3:
      startOfThisColor = {-2l, 0l, -1l};
      break;
  }


  // to fix compiler complaints about perfectly nested loop.
  long endX = _end[0];
  long startY = startOfThisColor[1], endY = _end[1];
  long startZ = startOfThisColor[2], endZ = _end[2];

  std::array<long, 3> startsOfX{};
  startsOfX[0] = startOfThisColor[0];
  startsOfX[1] = startOfThisColor[0] -4;
  startsOfX[2] = startOfThisColor[0] -2;

  //iterate over cartesian grid first time
#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(dynamic, 1) collapse(2) nowait
#endif
  for (long z = startZ; z < endZ; z += 4) {
    for (long y = startY; y < endY; y+=2) {
      for (long x = startsOfX[(z-startZ) % 12 / 4]; x < endX; x += 6) { // color starts every 6th column again
        const std::array<long, 3> base3DIndex = {x, y, z};
        processBasePack6(cells, base3DIndex);
      }
    }
  }

  startY++;
  startsOfX[0] +=3;
  startsOfX[1] +=3;
  startsOfX[2] +=3;

  //iterate over cartesian grid after shift
#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(dynamic, 1) collapse(2) nowait
#endif
  for (long z = startZ; z < endZ; z += 4) {
    for (long y = startY; y < endY; y+=2) {
      for (long x = startsOfX[(z-startZ) % 12 / 4]; x < endX; x += 6) { // color starts every 6th column again
        const std::array<long, 3> base3DIndex = {x, y, z};
        processBasePack6(cells, base3DIndex);
      }
    }
  }

}


}
