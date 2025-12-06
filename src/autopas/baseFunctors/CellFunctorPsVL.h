/**
 * @file CellFunctorPsVL.h
 *
 * @date 06.12.2025
 * @author Lars Doll
 */

#pragma once

#include "autopas/cells/SortedCellView.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas::internal {
/**
 * A cell functor. This functor is built from the normal Functor of the template
 * type ParticleFunctor. It is an internal object to handle interactions between
 * two cells of particles.
 * @tparam ParticleCell
 * @tparam ParticleFunctor the functor which
 * @tparam bidirectional if no newton3 is used processCellPair(cell1, cell2) should also handle processCellPair(cell2,
 * cell1)
 */
template <class ParticleCell, class ParticleFunctor, bool bidirectional = true>
class CellFunctorPsVL {
 public:
  /**
   * The constructor of CellFunctor.
   * @param f The ParticleFunctor which should be used for the interaction.
   * @param sortingCutoff This parameter indicates the maximal distance the sorted particles are to interact. This
   * parameter is only relevant for optimization (sorting). This parameter normally should be the cutoff, for building
   * verlet lists, this should be cutoff+skin.
   * @param dataLayout The data layout to be used.
   * @param useNewton3 Parameter to specify whether newton3 is used or not.
   * @param orientationLists
   */
  explicit CellFunctorPsVL(ParticleFunctor *f, const double sortingCutoff, DataLayoutOption dataLayout, bool useNewton3,
    const std::vector<std::vector<SortedCellView<typename ParticleCell::ParticleType>>> &orientationLists)
      : _functor(f), _sortingCutoff(sortingCutoff), _dataLayout(dataLayout), _useNewton3(useNewton3), _orientationList(orientationLists) {}

  /**
   * Process the interactions inside one cell.
   * @param cell All pairwise interactions of particles inside this cell are calculated.
   */
  void processCell(ParticleCell &cell);

  /**
   * Process the interactions between the particles of cell1 with particles of cell2.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2. If no parameter or {0, 0, 0} is
   * given, sorting will be disabled.
   */
  void processCellPair(ParticleCell &cell1, ParticleCell &cell2,
                       const std::array<double, 3> &sortingDirection = {0., 0., 0.});

  /**
   * Getter
   * @return
   */
  [[nodiscard]] DataLayoutOption getDataLayout() const { return _dataLayout; }

  /**
   * Getter
   * @return
   */
  [[nodiscard]] bool getNewton3() const { return _useNewton3; }

  /**
   * Getter
   * @return
   */
  [[nodiscard]] bool getBidirectional() const { return bidirectional; }

 private:
  /**
   * Applies the functor to all particle pairs exploiting newtons third law of motion.
   * There is only one version of this function as newton3 is always allowed to be applied inside of a cell.
   * The value of _useNewton3 defines whether or whether not to apply the aos version functor in a newton3 fashion or
   * not:
   * - if _useNewton3 is true: the aos functor will be applied once for each pair (only i,j), passing newton3=true.
   * - if _useNewton3 is false: the aos functor will be applied twice for each pair (i,j and j,i), passing
   * newton3=false.
   * @param cell
   */
  void processCellAoS(ParticleCell &cell);

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * exploiting newtons third law of motion.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairAoSN3(ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection);

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * without exploiting newtons third law of motion.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairAoSNoN3(ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection);

  void processCellPairSoAN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellPairSoANoN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellSoAN3(ParticleCell &cell);

  void processCellSoANoN3(ParticleCell &cell);

  ParticleFunctor *_functor;

  const double _sortingCutoff;

  DataLayoutOption _dataLayout;

  bool _useNewton3;

  std::vector<std::vector<SortedCellView<typename ParticleCell::ParticleType>>>& _orientationList;

};

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCell(ParticleCell &cell) {

}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellPair(
    ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection) {

}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellAoS(ParticleCell &cell) {

}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellPairAoSN3(
    ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection) {

}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellPairAoSNoN3(
    ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection) {
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellPairSoAN3(ParticleCell &cell1,
                                                                                     ParticleCell &cell2) {
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellPairSoANoN3(ParticleCell &cell1,
                                                                                       ParticleCell &cell2) {
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellSoAN3(ParticleCell &cell) {

}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellSoANoN3(ParticleCell &cell) {

}
}  // namespace autopas::internal
