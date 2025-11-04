/**
 * @file KokkosDsNaiveParallelTraversal.h
 * @date 31. October 2025
 * @author Luis Gall
 */

#pragma once

#include "autopas/containers/TraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

template <class Functor>
class KokkosDsNaiveParallelTraversal : public TraversalInterface {

public:

    explicit KokkosDsNaiveParallelTraversal(Functor *functor, DataLayoutOption dataLayout, bool useNewton3)
        : TraversalInterface(dataLayout, useNewton3), _functor{functor} {}

    [[nodiscard]] TraversalOption getTraversalType() const final { return TraversalOption::kokkos_ds_naive_parallel; }

    [[nodiscard]] bool isApplicable() const final { 
        // TODO
        return true;
    }

    void initTraversal() final {
        // TODO
    }

    void traverseParticles() final {
        // TODO
    }

    void endTraversal() final {
        // TODO
    }

private:

    const Functor *_functor;

};

}