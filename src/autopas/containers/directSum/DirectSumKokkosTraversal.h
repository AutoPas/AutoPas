/**
 * @file DirectSumKokkosTraversal.h
 * @author M. Geitner
 * @date 24.06.19
 */

#pragma once

#include <vector>

#include "DirectSumTraversalInterface.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ParticleCell.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/pairwiseFunctors/KokkosCellFunctor.h"
#include "autopas/particles/KokkosParticle.h"
#include "autopas/utils/KokkosDataLayoutConverter.h"
#ifdef AUTOPAS_KOKKOS
#include <Kokkos_Core.hpp>
#include "autopas/pairwiseFunctors/KokkosLJFunctor.h"
#endif

namespace autopas {

/**KOKKOS_FLOAT
 * This sum defines the traversal typically used by the DirectSum container.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class DirectSumKokkosTraversal : public CellPairTraversal<ParticleCell, DataLayout, useNewton3>,
                                 public DirectSumTraversalInterface<ParticleCell> {
  /**
   * Constructor for the DirectSum traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
 public:
  DirectSumKokkosTraversal(PairwiseFunctor *pairwiseFunctor)
      : CellPairTraversal<ParticleCell, DataLayout, useNewton3>({2, 1, 1}),
        _cellFunctor(internal::KokkosCellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, DataLayout, useNewton3>
                (pairwiseFunctor)),
        _dataLayoutConverter(pairwiseFunctor) {}

  TraversalOption getTraversalType() const override { return TraversalOption::kokkosDirectSumTraversal; }

  bool isApplicable() const override {
    return true;//std::is_same<typename ParticleCell::ParticleType, KokkosParticle>::value;
  }
  /**
   * @copydoc LinkedCellTraversalInterface::traverseCellPairs()
   * @note This function expects a vector of exactly two cells. First cell is the main region, second is halo.
   */
  void traverseCellPairs(std::vector<ParticleCell> &cells) override;

  void initTraversal(std::vector<ParticleCell> &cells) override {
#ifdef AUTOPAS_KOKKOS
      //std::cout << "DataLayoutOption:  " << utils::StringUtils::to_string(DataLayout) << "\n";
    // nothing to do here, is aos
      switch (DataLayout) {
          case DataLayoutOption::kokkos: {

          //std::cout << "init started" << "\n";
          //init structure

          for (unsigned int c = 0; c < cells.size(); c++) {
          _dataLayoutConverter.storeDataLayout(cells[c]);
          }

        }

          break;
          default:;
              //std::cout << "not soa" << utils::StringUtils::to_string(DataLayout) << "\n";
      }
      //auto buffer = cells[0]._particleKokkosBuffer;
      //std::cout << buffer.extent(0) << ", " << buffer.extent(1) << ", " << buffer.extent(2) << "\n";
      //buffer = cells[1]._particleKokkosBuffer;
      //std::cout << buffer.extent(0) << ", " << buffer.extent(1) << ", " << buffer.extent(2) << "\n";
      //std::cout << "-------------init end-------------\n";
      /*for (unsigned int c = 0; c < cells.size(); c++) {
          for(unsigned int i = 0; i  < cells[c]._particles.size(); i++){
              std::cout << cells[c]._particles[i].toString() << "\n";
          }
      }
      std::cout << "-------------------end 1 --------------------------------------\n";
       */
#endif
  }
  void endTraversal(std::vector<ParticleCell> &cells) override {
#ifdef AUTOPAS_KOKKOS
    // nothing to do here
    //std::cout << "endTraversal" << "\n";
      switch (DataLayout) {
          case DataLayoutOption::kokkos: {
              //copy data to particle
              for (unsigned int c = 0; c < cells.size(); c++) {
                  _dataLayoutConverter.loadDataLayout(cells[c]);
              }
              break;
          }
          default:;
      }
      /*
      for (unsigned int c = 0; c < cells.size(); c++) {
          for(unsigned int i = 0; i  < cells[c]._particles.size(); i++){
              std::cout << cells[c]._particles[i].toString() << "\n";
          }
      }
       */
#endif
  }


 private:
  /**
   * KokkosCellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::KokkosCellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, DataLayout, useNewton3>
      _cellFunctor;

  utils::KokkosDataLayoutConverter<PairwiseFunctor, DataLayout> _dataLayoutConverter;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
void DirectSumKokkosTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseCellPairs(
    std::vector<ParticleCell> &cells) {
  // Assume cell[0] is the main domain and cell[1] is the halo

#ifdef AUTOPAS_KOKKOS

  _cellFunctor.processCell(cells[0]);
  _cellFunctor.processCellPair(cells[0], cells[1]);
#endif
}


}  // namespace autopas