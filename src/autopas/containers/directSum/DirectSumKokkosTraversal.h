/**
 * @file DirectSumKokkosTraversal.h
 * @author M. Geitner
 * @date 04/18/2019
 */

#pragma once

/*
#ifndef AUTOPAS_DIRECTSUMKOKKOSTRAVERSAL_H
#define AUTOPAS_DIRECTSUMKOKKOSTRAVERSAL_H

#endif //AUTOPAS_DIRECTSUMKOKKOSTRAVERSAL_H
*/
#include <vector>
#include <autopas/cells/ParticleCell.h>
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "DirectSumTraversal.h"


namespace autopas {

    template<class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
    class DirectSumKokkosTraversal
            : public CellPairTraversal<ParticleCell>, public DirectSumTraversalInterface<ParticleCell> {


    public
        :
        DirectSumKokkosTraversal(PairwiseFunctor *pairwiseFunctor) :
                CellPairTraversal<ParticleCell>({2, 1, 1}),
                _cellFunctor(
                        CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, useSoA, useNewton3>(
                                pairwiseFunctor)) {}


        TraversalOption getTraversalType() override { return TraversalOption::directSumKokkosTraversal; }

        void traverseCellPairs(std::vector<ParticleCell> &cells) override;

    private:
        /**
         * CellFunctor to be used for the traversal defining the interaction between two cells.
         */
        CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, useSoA, useNewton3> _cellFunctor;
    };


    template<class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
    void DirectSumKokkosTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::traverseCellPairs(
            std::vector<ParticleCell> &cells) {
        // Assume cell[0] is the main domain and cell[1] is the halo
        _cellFunctor.processCell(cells[0]);
        _cellFunctor.processCellPair(cells[0], cells[1]);
    }

}