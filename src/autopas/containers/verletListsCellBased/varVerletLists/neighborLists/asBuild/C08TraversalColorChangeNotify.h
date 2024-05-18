/**
 * @file C08TraversalColorChangeNotify.h
 * @author humig
 * @date 09.07.19
 */

#pragma once

#include "ColorChangeObserver.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"

namespace autopas {

/**
 * @copydoc LCC08Traversal
 *
 * It furthermore notifies an observer when color changes during the traversal happen.
 */
template <class ParticleCell, class PairwiseFunctor>
class C08TraversalColorChangeNotify : public LCC08Traversal<ParticleCell, PairwiseFunctor> {
 public:
  /**
   * Constructor of the traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param observer The observer to notify when a color change happens during the traversal.
   * @param dataLayout The data layout with which this traversal should be initialised.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  C08TraversalColorChangeNotify(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                double interactionLength, const std::array<double, 3> &cellLength,
                                ColorChangeObserver *observer, DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3),
        LCC08Traversal<ParticleCell, PairwiseFunctor>(dims, pairwiseFunctor, interactionLength, cellLength, dataLayout,
                                                      useNewton3),
        _observer(observer) {}

 protected:
  /**
   * Notifies the observer.
   * @param newColor The new color to notify the observer about.
   */
  void notifyColorChange(unsigned long newColor) override {
    if (_observer) _observer->receiveColorChange(newColor);
  }

 private:
  /**
   * The observer of the traversal.
   */
  ColorChangeObserver *_observer;
};

}  // namespace autopas