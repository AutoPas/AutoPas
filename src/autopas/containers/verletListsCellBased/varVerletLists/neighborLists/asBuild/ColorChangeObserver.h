/**
 * @file ColorChangeObserver.h
 * @author humig
 * @date 09.07.19
 */

#pragma once

namespace autopas {

/**
 * Observer interface that specifies handling color changes.
 */
class ColorChangeObserver {
 public:
  /**
   * Gets called when the color changes during the observed traversal.
   * @param newColor The new color that the traversal handles now.
   */
  virtual void receiveColorChange(uint64_t newColor) = 0;
};

}  // namespace autopas