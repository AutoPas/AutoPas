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

enum Face {
  O = 0,  // omega/unknown
  L = 1,
  R = 2,
  D = 3,
  U = 4,
  B = 5,
  F = 6,
};

typedef int unsigned Any;

static constexpr Any buildEdge(Face f1, Face f2) {
  assert(f1 != f2);  // TODO: How can I make this a static assert?
  return (f1 << 3) | f2;
}

static constexpr Any buildVertex(Face f1, Face f2, Face f3) {
  assert((f1 != f2) && (f2 != f3));  // TODO: How can I make this a static assert?
  return (f1 << 6) | (f2 << 3) | f3;
}

enum Edge {
  OO = 0,  // omega/unknown
  LD = buildEdge(L, D),
  LU = buildEdge(L, U),
  LB = buildEdge(L, B),
  LF = buildEdge(L, F),
  RD = buildEdge(R, D),
  RU = buildEdge(R, U),
  RB = buildEdge(R, B),
  RF = buildEdge(R, F),
  DB = buildEdge(D, B),
  DF = buildEdge(D, F),
  UB = buildEdge(U, B),
  UF = buildEdge(U, F),
};

enum Vertex {
  OOO = 0,  // omega/unknown
  LDB = buildVertex(L, D, B),
  LDF = buildVertex(L, D, F),
  LUB = buildVertex(L, U, B),
  LUF = buildVertex(L, U, F),
  RDB = buildVertex(R, D, B),
  RDF = buildVertex(R, D, F),
  RUB = buildVertex(R, U, B),
  RUF = buildVertex(R, U, F),
};

inline Face *getFaces() {
  static Face table[] = {L, R, D, U, B, F};
  return table;
}

inline Edge *getEdges() {
  static Edge table[] = {LD, LU, LB, LF, RD, RU, RB, RF, DB, DF, UB, UF};
  return table;
}

inline Vertex *VERTICES() {
  static Vertex table[] = {LDB, LDF, LUB, LUF, RDB, RDF, RUB, RUF, OOO};
  return table;
}

typedef Vertex Octant;

inline bool ADJ(Any direction, Vertex octant) {
  static std::array<std::array<bool, 8>, 1 << 9> table;

  // TODO: Is this actually initialized static??
  table[L] = {true, true, true, true, false, false, false, false};
  table[R] = {false, false, false, false, true, true, true, true};
  table[D] = {true, true, false, false, true, true, false, false};
  table[U] = {false, false, true, true, false, false, true, true};
  table[B] = {true, false, true, false, true, false, true, false};
  table[F] = {false, true, false, true, false, true, false, true};
  table[LD] = {true, true, false, false, false, false, false, false};
  table[LU] = {false, false, true, true, false, false, false, false};
  table[LB] = {true, false, true, false, false, false, false, false};
  table[LF] = {false, true, false, true, false, false, false, false};
  table[RD] = {false, false, false, false, true, true, false, false};
  table[RU] = {false, false, false, false, false, false, true, true};
  table[RB] = {false, false, false, false, true, false, true, false};
  table[RF] = {false, false, false, false, false, true, false, true};
  table[DB] = {true, false, false, false, true, false, false, false};
  table[DF] = {false, true, false, false, false, true, false, false};
  table[UB] = {false, false, true, false, false, false, true, false};
  table[UF] = {false, false, false, true, false, false, false, true};
  table[LDB] = {true, false, false, false, false, false, false, false};
  table[LDF] = {false, true, false, false, false, false, false, false};
  table[LUB] = {false, false, true, false, false, false, false, false};
  table[LUF] = {false, false, false, true, false, false, false, false};
  table[RDB] = {false, false, false, false, true, false, false, false};
  table[RDF] = {false, false, false, false, false, true, false, false};
  table[RUB] = {false, false, false, false, false, false, true, false};
  table[RUF] = {false, false, false, false, false, false, false, true};

  return table[direction][octant];
}

inline Octant REFLECT(Any direction, Octant octant) {
  static std::array<std::array<Octant, 8>, 1 << 9> table;

  table[L] = {RDB, RDF, RUB, RUF, LDB, LDF, LUB, LUF};
  table[R] = {RDB, RDF, RUB, RUF, LDB, LDF, LUB, LUF};
  table[D] = {LUB, LUF, LDB, LDF, RUB, RUF, RDB, RDF};
  table[U] = {LUB, LUF, LDB, LDF, RUB, RUF, RDB, RDF};
  table[B] = {LDF, LDB, LUF, LUB, RDF, RDB, RUF, RUB};
  table[F] = {LDF, LDB, LUF, LUB, RDF, RDB, RUF, RUB};
  table[LD] = {RUB, RUF, RDB, RDF, LUB, LUF, LDB, LDF};
  table[LU] = {RUB, RUF, RDB, RDF, LUB, LUF, LDB, LDF};
  table[LB] = {RDF, RDB, RUF, RUB, LDF, LDB, LUF, LUB};
  table[LF] = {RDF, RDB, RUF, RUB, LDF, LDB, LUF, LUB};
  table[RD] = {RUB, RUF, RDB, RDF, LUB, LUF, LDB, LDF};
  table[RU] = {RUB, RUF, RDB, RDF, LUB, LUF, LDB, LDF};
  table[RB] = {RDF, RDB, RUF, RUB, LDF, LDB, LUF, LUB};
  table[RF] = {RDF, RDB, RUF, RUB, LDF, LDB, LUF, LUB};
  table[DB] = {LUF, LUB, LDF, LDB, RUF, RUB, RDF, RDB};
  table[DF] = {LUF, LUB, LDF, LDB, RUF, RUB, RDF, RDB};
  table[UB] = {LUF, LUB, LDF, LDB, RUF, RUB, RDF, RDB};
  table[UF] = {LUF, LUB, LDF, LDB, RUF, RUB, RDF, RDB};
  table[LDB] = {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB};
  table[LDF] = {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB};
  table[LUB] = {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB};
  table[LUF] = {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB};
  table[RDB] = {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB};
  table[RDF] = {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB};
  table[RUB] = {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB};
  table[RUF] = {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB};

  return table[direction][octant];
}

inline Face COMMON_FACE(Any direction, Vertex octant) {
  static std::array<std::array<Face, 8>, 1 << 9> table;

  table[LD] = {O, O, L, L, D, D, O, O};
  table[LU] = {L, L, O, O, O, O, U, U};
  table[LB] = {O, L, O, L, B, O, B, O};
  table[LF] = {L, O, L, O, O, F, O, F};
  table[RD] = {D, D, O, O, O, O, R, R};
  table[RU] = {O, O, U, U, R, R, O, O};
  table[RB] = {B, O, B, O, O, R, O, R};
  table[RF] = {O, F, O, F, R, O, R, O};
  table[DB] = {O, D, B, O, O, D, B, O};
  table[DF] = {D, O, O, F, D, O, O, F};
  table[UB] = {B, O, O, U, B, O, O, U};
  table[UF] = {O, F, U, O, O, F, U, O};
  table[LDB] = {O, O, O, L, O, D, B, O};
  table[LDF] = {O, O, L, O, D, O, O, F};
  table[LUB] = {O, L, O, O, B, O, O, U};
  table[LUF] = {L, O, O, O, O, F, U, O};
  table[RDB] = {O, D, B, O, O, O, O, R};
  table[RDF] = {D, O, O, F, O, O, R, O};
  table[RUB] = {B, O, O, U, O, R, O, O};
  table[RUF] = {O, F, U, O, R, O, O, O};

  return table[direction][octant];
}

inline Edge COMMON_EDGE(Any direction, Vertex octant) {
  static std::array<std::array<Edge, 8>, 1 << 9> table;

  table[LDB] = {OO, LD, LB, OO, DB, OO, OO, OO};
  table[LDF] = {LD, OO, OO, LF, OO, DF, OO, OO};
  table[LUB] = {LB, OO, OO, LU, OO, OO, UB, OO};
  table[LUF] = {OO, LF, LU, OO, OO, OO, OO, UF};
  table[RDB] = {DB, OO, OO, OO, OO, RD, RB, OO};
  table[RDF] = {OO, DF, OO, OO, RD, OO, OO, RF};
  table[RUB] = {OO, OO, UB, OO, RB, OO, OO, RU};
  table[RUF] = {OO, OO, OO, UF, OO, RF, RU, OO};

  return table[direction][octant];
}

}  // namespace autopas
