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
    return std::is_same<typename ParticleCell::ParticleType, KokkosParticle>::value;
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
          /*
              //std::cout << "Num Particles in Cell: " << cells[c].numParticles() << "\n";

              FloatMatrix3Type buffer("_buffer", cells[c].numParticles(), 3, KOKKOS_DIM);

              FloatMatrix3Type::HostMirror h_matrix = Kokkos::create_mirror_view(buffer);
              for (unsigned int x = 0; x < cells[c]._particles.size(); x++) {
                  //position
                  auto arr_r = cells[c]._particles[x].getR();
                  auto arr_f = cells[c]._particles[x].getF();
                  auto arr_v = cells[c]._particles[x].getV();
                  for (unsigned int i = 0; i < KOKKOS_DIM; i++) {
                      h_matrix(x, 0, i) = arr_r[i];//position
                      h_matrix(x, 1, i) = arr_f[i];//force
                      h_matrix(x, 2, i) = arr_v[i];//velocity
                  }
              }
              Kokkos::deep_copy(buffer, h_matrix);
              //std::cout << buffer.extent(0) << ", " << buffer.extent(1) << ", " << buffer.extent(2) << "\n";
              //Test structure
              //Kokkos::parallel_for(cell._particles.size(), KOKKOS_LAMBDA(const int i){
              //cell._particleKokkosBuffer(0, 0, 0) += 1;
              //});
              cells[c]._particleKokkosBuffer =  buffer;*/
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
#endif
  }
  void endTraversal(std::vector<ParticleCell> &cells) override {
#ifdef AUTOPAS_KOKKOS
    // nothing to do here
    //std::cout << "endTraversal" << "\n";
      if(DataLayout == DataLayoutOption::kokkos){
          //copy data to particle
          for (unsigned int c = 0; c < cells.size(); c++) {
              _dataLayoutConverter.loadDataLayout(cells[c]);
          }
      }
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