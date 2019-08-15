/**
*
*@file CBasedKokkosTraversal.h
*@author M. Geitner
*@date 16.07.2019
*
*
*/



#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/KokkosDataLayoutConverter.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using base steps based on cell coloring.
 * This class is designed for using kokkos. if no kokkos is used, please use the base class CBasedTraversal.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 * @tparam collapseDepth Set the depth of loop collapsion for OpenMP. Loop variables from outer to inner loop: z,y,x
 */
    template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3,
            int collapseDepth = 3>
    class CBasedKokkosTraversal : public CellPairTraversal<ParticleCell, dataLayout, useNewton3> {
    protected:
        /**
         * Constructor of the CBasedKokkosTraversal.
         * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
         * y and z direction.
         * @param pairwiseFunctor The functor that defines the interaction of two particles.
         * @param cutoff Cutoff radius.
         * @param cellLength cell length.
         */
        explicit CBasedKokkosTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                 const double cutoff, const std::array<double, 3> &cellLength)
                : CellPairTraversal<ParticleCell, dataLayout, useNewton3>(dims),
                  _cutoff(cutoff),
                  _cellLength(cellLength),
                  _dataLayoutConverter(pairwiseFunctor) {
          for (unsigned int d = 0; d < 3; d++) {
            _overlap[d] = std::ceil(_cutoff / _cellLength[d]);
          }
        }

        /**
         * Destructor of CBasedKokkosTraversal.
         */
        ~CBasedKokkosTraversal() override = default;

    public:
        /**
         * load Data Layouts required for this Traversal.
         * @param cells where the data should be loaded
         */
        void initTraversal(std::vector<ParticleCell> &cells) override {
            if(dataLayout == DataLayoutOption::kokkos){
                //copy data to particle
                /*check particle attributes
                for (unsigned int c = 0; c < cells.size(); c++) {
                    for(Particle p:cells[c]._particles){
                        std::cout << p.toString() << "\n";
                    }
                }
                 */
                //std::cout << "--------------------\n";
                //std::chrono::high_resolution_clock::time_point startCalc, stopCalc;

                //startCalc = std::chrono::high_resolution_clock::now();
                for (unsigned int c = 0; c < cells.size(); c++) {
                    _dataLayoutConverter.storeDataLayout(cells[c]);
                }
                //stopCalc = std::chrono::high_resolution_clock::now();

                //auto durationCalc = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
                //auto durationCalcSec = durationCalc * 1e-6;
                //std::cout << "Time init: " << durationCalc << " \u03bcs (" << durationCalcSec << "s)" << std::endl;
            }
            /*for (unsigned int c = 0; c < cells.size(); c++) {
                for(Particle p:cells[c]._particles){
                    std::cout << p.toString() << "\n";
                }
            }*/
        }

        /**
         * write Data to AoS.
         * @param cells for which the data should be written back
         */
        void endTraversal(std::vector<ParticleCell> &cells) override {
            if(dataLayout == DataLayoutOption::kokkos){
                //copy data to particle
                //std::chrono::high_resolution_clock::time_point startCalc, stopCalc;

                //startCalc = std::chrono::high_resolution_clock::now();
                for (unsigned int c = 0; c < cells.size(); c++) {
                    _dataLayoutConverter.loadDataLayout(cells[c]);
                }
                //stopCalc = std::chrono::high_resolution_clock::now();

                //auto durationCalc = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
                //auto durationCalcSec = durationCalc * 1e-6;
                //std::cout << "Time exit: " << durationCalc << " \u03bcs (" << durationCalcSec << "s)" << std::endl;
            }
            //check particle attributes
            /*for (unsigned int c = 0; c < cells.size(); c++) {
                for(unsigned int i = 0; i  < cells[c]._particles.size(); i++){
                    std::cout << cells[c]._particles[i].toString() << "\n";
                }
            }*/
        }

    protected:
        /**
         * The main traversal of the CTraversal.
         * @tparam LoopBody type of the loop body
         * @param loopBody The body of the loop as a function. Normally a lambda function, that takes as as parameters
         * (x,y,z). If you need additional input from outside, please use captures (by reference).
         * @param end 3D index until interactions are processed (exclusive)
         * @param stride dimension of stride (depends on coloring)
         * @param offset initial offset
         */
        template <typename LoopBody>
        inline void cTraversal(LoopBody &&loopBody, const std::array<unsigned long, 3> &end,
                               const std::array<unsigned long, 3> &stride,
                               const std::array<unsigned long, 3> &offset = {0ul, 0ul, 0ul});

        /**
         * cutoff radius.
         */
        const double _cutoff;

        /**
         * cell length in CellBlock3D.
         */
        const std::array<double, 3> _cellLength;

        /**
         * overlap of interacting cells. Array allows asymmetric cell sizes.
         */
        std::array<unsigned long, 3> _overlap;

    private:
        /**
         * Data Layout Converter to be used with this traversal
         */
        utils::KokkosDataLayoutConverter<PairwiseFunctor, dataLayout> _dataLayoutConverter;
    };

    template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3, int collapseDepth>
    template <typename LoopBody>
    inline void CBasedKokkosTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, collapseDepth>::cTraversal(
            LoopBody &&loopBody, const std::array<unsigned long, 3> &end, const std::array<unsigned long, 3> &stride,
            const std::array<unsigned long, 3> &offset) {

      {
        //std::cout << "startc08BasedKokkos\n";
        const unsigned long numColors = stride[0] * stride[1] * stride[2];
        for (unsigned long col = 0; col < numColors; ++col) {
          std::array<unsigned long, 3> startWithoutOffset(utils::ThreeDimensionalMapping::oneToThreeD(col, stride));
          std::array<unsigned long, 3> start(ArrayMath::add(startWithoutOffset, offset));

          // intel compiler demands following:
          const unsigned long start_x = start[0], start_y = start[1], start_z = start[2];
          const unsigned long end_x = end[0], end_y = end[1], end_z = end[2];
          const unsigned long stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];
#ifdef AUTOPAS_KOKKOS



          int iterationsZ = 1 + (end_z - 1 - start_z)/stride_z;
          if(end_z < start_z || iterationsZ < 0) iterationsZ = 0;

          unsigned int iterationsY = 1 + (end_y - 1 - start_y)/stride_y;
          if(end_y < start_y || iterationsY < 0) iterationsY = 0;

            unsigned int iterationsX = 1 + (end_x - 1 - start_x)/stride_x;
            if(end_x < start_x || iterationsX < 0) iterationsX = 0;


/*
          Kokkos::parallel_for(team_policy( iterationsZ, Kokkos::AUTO), KOKKOS_LAMBDA ( const member_type &teamMember){
              const int i = teamMember.league_rank();
              //std::cout << "League Size: " << teamMember.team_size() <<"\n";
              Kokkos::parallel_for(iterationsY, [&] (const int l) {
                for (unsigned long x = start_x; x < end_x; x += stride_x) {
                  // Don't exchange order of execution (x must be last!), it would break other code
                  loopBody(x, start_y + l * stride_y, start_z + i * stride_z);
                }
              });
          });
*/
/*
            typedef Kokkos::TeamPolicy<> team_policy;
            typedef Kokkos::TeamPolicy<>::member_type member_type;
            Kokkos::parallel_for(team_policy( iterationsZ, Kokkos::AUTO), KOKKOS_LAMBDA ( const member_type &teamMember){
                const int i = teamMember.league_rank();
                //std::cout << "League Size: " << teamMember.team_size() <<"\n";
                Kokkos::parallel_for(Kokkos::TeamVectorRange(teamMember, iterationsY), [&] (const int l) {
                    for (unsigned long x = start_x; x < end_x; x += stride_x) {
                        // Don't exchange order of execution (x must be last!), it would break other code
                        loopBody(x, start_y + l * stride_y, start_z + i * stride_z);
                    }
                });
            });
*/

            //Test version, currently decreases perormance in OpenMP Kokkos, perhaps useful for CUDA?
            //best result

            Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {iterationsZ, iterationsY, iterationsX}),
                    KOKKOS_LAMBDA(int z, int y, int x){
                        //std::cout << x << ", " << y << ", " << z << "\n";
                        // Don't exchange order of execution (x must be last!), it would break other code
                        loopBody(start_x + x * stride_x, start_y + y * stride_y, start_z + z * stride_z);
            });




          //std::cout << start_z <<" | " << end_z << " | " << stride_z << " | " << iterationsZ << "\n";
            /*Kokkos::parallel_for(range_policy(0, iterationsZ), KOKKOS_LAMBDA(const int i){
              for (unsigned long y = start_y; y < end_y; y += stride_y) {
                for (unsigned long x = start_x; x < end_x; x += stride_x) {
                  // Don't exchange order of execution (x must be last!), it would break other code
                  loopBody(x, y, start_z + i * stride_z);
                }
              }
            }
          );*/
          //std::cout << "loop end\n";
#else
          //if kokkos is diabled
          for (unsigned long z = start_z; z < end_z; z += stride_z) {
              for (unsigned long y = start_y; y < end_y; y += stride_y) {
                for (unsigned long x = start_x; x < end_x; x += stride_x) {
                  // Don't exchange order of execution (x must be last!), it would break other code
                  loopBody(x, y, z);
                }
              }
            }
#endif
        }
      }
    }

}  // namespace autopas
