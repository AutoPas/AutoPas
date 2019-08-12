/**
 *
 * @author M. Geitner
 * @file C08KokkosTraversal.h
 * @date 11.07.19
 *
 */


#pragma once
#include "C08KokkosCellHandler.h"
#include "autopas/containers/cellPairTraversals/C08BasedKokkosTraversal.h"
#include "autopas/utils/KokkosDataLayoutConverter.h"
namespace autopas {

/**
 * This class provides the base for traversals using the c08 base step.
 *
 * The traversal is defined in the function c08Traversal and uses 8 colors, such that interactions between the base
 * cell and all adjacent cells with greater ID in each direction are safe, even when using newton3 optimizations.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
    template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
    class C08KokkosTraversal :  public C08BasedKokkosTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>,
                                public LinkedCellTraversalInterface<ParticleCell>{
    public:
        /**
         * Constructor of the c08 traversal.
         * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
         * y and z direction.
         * @param pairwiseFunctor The functor that defines the interaction of two particles.
         * @param cutoff Cutoff radius.
         * @param cellLength cell length.
         */
        C08KokkosTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                   const double cutoff = 1.0, const std::array<double, 3> &cellLength = {1.0, 1.0, 1.0})
                : C08BasedKokkosTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>(dims, pairwiseFunctor, cutoff, cellLength),
                     _cellHandler(pairwiseFunctor, this->_cellsPerDimension, cutoff, cellLength, this->_overlap),
                     _functor(pairwiseFunctor){}

        TraversalOption getTraversalType() const override {
          return TraversalOption::kokkosc08;
        }

        bool isApplicable() const override {
          return std::is_same<typename ParticleCell::ParticleType, KokkosParticle>::value;
        }


        void traverseCellPairs(std::vector<ParticleCell> &cells) override;
    private:


    protected:

    private:
        //pairwise functor
        /**
         * Pairs for processBaseCell().
         */
        std::vector<int> _cellOffsets;

        /***
         * KokkosCellHandler for c08
         */
        C08KokkosCellHandler<ParticleCell, PairwiseFunctor, DataLayout, useNewton3> _cellHandler;

        /**
         * Pairwise Functor to be used
         */
        PairwiseFunctor *_functor;



    };

    template<class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
    void C08KokkosTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::traverseCellPairs(
            std::vector<ParticleCell> &cells) {
      this->c08Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
          unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
          _cellHandler.processBaseCell(cells, baseIndex);
      });
    }


}  // namespace autopas
