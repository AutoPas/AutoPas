/**
 * @file OctreeInnerNode.h
 *
 * @author Johannes Spies
 * @date 18.05.2021
 */

#pragma once

namespace autopas {

/**
 * Used to index children inside the octree along the axis.
 */
enum Direction {
  NEG_X = (1 << 0),
  POS_X = (1 << 1),
  NEG_Y = (1 << 2),
  POS_Y = (1 << 3),
  NEG_Z = (1 << 4),
  POS_Z = (1 << 5),
};

/**
 * Convert configuration pair (axis,sign) to a Direction.
 *
 * axis=0 means x axis
 * axis=1 means y axis
 * axis=2 means z axis
 *
 * sign=false means negative
 * sign=true means positive
 *
 * For instance, (1,true) should point towards the positive y direction.
 *
 * @param axis The axis index.
 * @param sign The sign of the direction.
 * @return The bit using the Direction enum.
 */
static inline Direction getDirection(unsigned axis, bool sign) {
  if (axis > 2) {
    throw std::runtime_error("[OctreeDirection] Axis must be in range [0,2]");
  }

  // Find the bit that corresponds to the (axis,sign) configuration.
  // The bits are look like the following
  // <- msb                     lsb ->
  // z axis    | y axis    | x axis
  // [pos neg] | [pos neg] | [pos neg]
  unsigned signShift = sign ? 1 : 0;
  unsigned axisShift = 2 * axis;
  return static_cast<Direction>((1u << axisShift) << signShift);
}

/**
 * Convert a direction to a flat index.
 * @param d The direction to convert.
 * @return An integer in the range [0,7].
 */
static inline unsigned directionToIndex(Direction d) {
  // TODO(johannes): Find a less jumpy way...
  switch (d) {
    case NEG_X:
      return 0;
    case POS_X:
      return 1;
    case NEG_Y:
      return 2;
    case POS_Y:
      return 3;
    case NEG_Z:
      return 4;
    case POS_Z:
      return 5;
    default:
      throw std::runtime_error("[OctreeDirection] Direction must be a valid direction from the enum");
  }
}

}  // namespace autopas
