/**
 * @file C08TraversalColorChangeNotify.h
 * @author humig
 * @date 09.07.19
 */

#pragma once

#include "../../../../linkedCells/traversals/C08Traversal.h"
#include "ColorChangeObserver.h"

namespace autopas {

/**
 * @copydoc C08Traversal
 *
 * It furthermore notifies an observer when color changes during the traversal happen.
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
class C08TraversalColorChangeNotify : public C08Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> {
 public:
  /**
   * Constructor of the traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param observer The observer to notify when a color change happens during the traversal.
   */
  C08TraversalColorChangeNotify(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                ColorChangeObserver *observer)
      : C08Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor),
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