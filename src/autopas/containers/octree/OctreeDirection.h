/**
 * @file OctreeDirection.h
 *
 * @author Johannes Spies
 * @date 18.05.2021
 */

#pragma once

namespace autopas::octree {

/**
 * A datatype that is wide enough to hold faces, edges or vertices.
 */
using Any = int unsigned;

/**
 * This enum can be used to index the faces of a cube including an "invalid" face.
 */
enum Face : Any {
  O = 0,  // omega/unknown
  L = 1,
  R = 2,
  D = 3,
  U = 4,
  B = 5,
  F = 6,
};

/**
 * Create a bitfield for an edge given by the template parameters.
 *
 * @tparam f1
 * @tparam f2
 * @return
 */
template <Face f1, Face f2>
static constexpr Any buildEdge() {
  static_assert(f1 != f2, "Faces must be different");
  return (f1 << 3) | f2;
}

/**
 * Create a bitfield for a vertex given by the template parameters.
 *
 * @tparam f1
 * @tparam f2
 * @tparam f3
 * @return
 */
template <Face f1, Face f2, Face f3>
static constexpr Any buildVertex() {
  static_assert((f1 != f2) and (f2 != f3), "Faces must be different");
  return (f1 << 6) | (f2 << 3) | f3;
}

/**
 * This enum can be used to index all edges of a cube including an "invalid" edge.
 */
enum Edge : Any {
  OO = 0,  // omega/unknown
  LD = buildEdge<L, D>(),
  LU = buildEdge<L, U>(),
  LB = buildEdge<L, B>(),
  LF = buildEdge<L, F>(),
  RD = buildEdge<R, D>(),
  RU = buildEdge<R, U>(),
  RB = buildEdge<R, B>(),
  RF = buildEdge<R, F>(),
  DB = buildEdge<D, B>(),
  DF = buildEdge<D, F>(),
  UB = buildEdge<U, B>(),
  UF = buildEdge<U, F>(),
};

/**
 * This enum can be used to index all vertices of a cube including an "invalid" vertex.
 */
enum Vertex : Any {
  OOO = 0,  // omega/unknown
  LDB = buildVertex<L, D, B>(),
  LDF = buildVertex<L, D, F>(),
  LUB = buildVertex<L, U, B>(),
  LUF = buildVertex<L, U, F>(),
  RDB = buildVertex<R, D, B>(),
  RDF = buildVertex<R, D, F>(),
  RUB = buildVertex<R, U, B>(),
  RUF = buildVertex<R, U, F>(),
};

/**
 * This namespace contains all tables that can be used publicly.
 */
namespace Tables {
/**
 * All available faces for a cube. The "omega" face `O` is excluded.
 */
constexpr static std::array<Face, 6> faces = {L, R, D, U, B, F};

/**
 * All available edges for a cube. The "omega" edge `OO` is excluded.
 */
constexpr static std::array<Edge, 12> edges = {LD, LU, LB, LF, RD, RU, RB, RF, DB, DF, UB, UF};

/**
 * All available vertices for a cube. The "omega" vertex `OOO` is excluded.
 */
constexpr static std::array<Vertex, 8> vertices = {LDB, LDF, LUB, LUF, RDB, RDF, RUB, RUF};
}  // namespace Tables

/**
 * Map an arbitrary vertex to a flat index.
 *
 * @param vertex An element from the Vertex enum
 * @return A flat index in the range of 0 to 7 for any valid vertex. For invalid input -1 is returned.
 */
inline int vertexToIndex(Vertex vertex) {
  // @todo This is very slow and could be sped up.
  //   Mentioned in https://github.com/AutoPas/AutoPas/issues/623
  int result = -1;
  for (int i = 0; i < 8; ++i) {
    if (vertex == Tables::vertices[i]) {
      result = i;
      break;
    }
  }
  return result;
}

/**
 * A vertex is also capable of specifying an arbitrary octant in 3D.
 */
using Octant = Vertex;

/**
 * Check if f is a face.
 *
 * @param f The parameter to check
 * @return true iff f is in the list returned from getFaces()
 */
template <typename T>
inline bool isFace(T f) {
  return std::find(Tables::faces.begin(), Tables::faces.end(), f) != Tables::faces.end();
}

/**
 * Check if e is an edge.
 *
 * @tparam T
 * @param e The parameter to check
 * @return true iff e is in the list returned from getEdges()
 */
template <typename T>
inline bool isEdge(T e) {
  return std::find(Tables::edges.begin(), Tables::edges.end(), e) != Tables::edges.end();
}

/**
 * Check if v is a vertex.
 *
 * @tparam T
 * @param v The parameter to check
 * @return true iff v is in the list returned from VERTICES()
 */
template <typename T>
inline bool isVertex(T v) {
  return std::find(Tables::vertices.begin(), Tables::vertices.end(), v) != Tables::vertices.end();
}

/**
 * This namespace contains the tables that are anonymous within this file. These tables should only be used by their
 * corresponding functions and are therefore private.
 */
namespace {
/**
 * A `constexpr` for creating the table for `ADJ()`.
 *
 * @return A LUT
 */
constexpr std::array<std::array<bool, 8>, 1 << 9> createADJTable() {
  std::array<std::array<bool, 8>, 1 << 9> table{};
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
  return table;
}

/**
 * A `constexpr` for creating the table for `REFLECT()`.
 *
 * @return A LUT
 */
constexpr std::array<std::array<Octant, 8>, 1 << 9> createREFLECTTable() {
  std::array<std::array<Octant, 8>, 1 << 9> table{};
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
  return table;
}

/**
 * A `constexpr` for creating the table for `COMMON_FACE()`.
 *
 * @return A LUT
 */
constexpr std::array<std::array<Face, 8>, 1 << 9> createCOMMONFACETable() {
  std::array<std::array<Face, 8>, 1 << 9> table{};
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
  return table;
}

/**
 * A `constexpr` for creating the table for `COMMON_EDGE()`.
 *
 * @return A LUT
 */
constexpr std::array<std::array<Edge, 8>, 1 << 9> createCOMMONEDGETable() {
  std::array<std::array<Edge, 8>, 1 << 9> table{};
  table[LDB] = {OO, LD, LB, OO, DB, OO, OO, OO};
  table[LDF] = {LD, OO, OO, LF, OO, DF, OO, OO};
  table[LUB] = {LB, OO, OO, LU, OO, OO, UB, OO};
  table[LUF] = {OO, LF, LU, OO, OO, OO, OO, UF};
  table[RDB] = {DB, OO, OO, OO, OO, RD, RB, OO};
  table[RDF] = {OO, DF, OO, OO, RD, OO, OO, RF};
  table[RUB] = {OO, OO, UB, OO, RB, OO, OO, RU};
  table[RUF] = {OO, OO, OO, UF, OO, RF, RU, OO};
  return table;
}

/**
 * A `constexpr` for creating the table for `getOppositeDirection()`.
 *
 * @return A LUT
 */
constexpr std::array<Any, 1 << 9> createOppositeDirectionTable() {
  std::array<Any, 1 << 9> table{};
  table[L] = R;
  table[R] = L;
  table[D] = U;
  table[U] = D;
  table[B] = F;
  table[F] = B;
  table[LD] = RU;
  table[LU] = RD;
  table[LB] = RF;
  table[LF] = RB;
  table[RD] = LU;
  table[RU] = LD;
  table[RB] = LF;
  table[RF] = LB;
  table[DB] = UF;
  table[DF] = UB;
  table[UB] = DF;
  table[UF] = DB;
  table[LDB] = RUF;
  table[LDF] = RUB;
  table[LUB] = RDF;
  table[LUF] = RDB;
  table[RDB] = LUF;
  table[RDF] = LUB;
  table[RUB] = LDF;
  table[RUF] = LDB;
  return table;
}

/**
 * A `constexpr` for creating the table for `getAllowedDirections()`.
 *
 * @return A LUT
 */
constexpr std::array<std::array<Octant, 5>, 1 << 9> createAllowedDirectionsTable() {
  std::array<std::array<Octant, 5>, 1 << 9> table{};
  table[L] = {LDB, LDF, LUB, LUF, OOO};
  table[R] = {RDB, RDF, RUB, RUF, OOO};
  table[D] = {LDB, LDF, RDB, RDF, OOO};
  table[U] = {LUB, LUF, RUB, RUF, OOO};
  table[B] = {LDB, LUB, RDB, RUB, OOO};
  table[F] = {LDF, LUF, RDF, RUF, OOO};
  table[LD] = {LDB, LDF, OOO, OOO, OOO};
  table[LU] = {LUB, LUF, OOO, OOO, OOO};
  table[LB] = {LDB, LUB, OOO, OOO, OOO};
  table[LF] = {LDF, LUF, OOO, OOO, OOO};
  table[RD] = {RDB, RDF, OOO, OOO, OOO};
  table[RU] = {RUB, RUF, OOO, OOO, OOO};
  table[RB] = {RDB, RUB, OOO, OOO, OOO};
  table[RF] = {RDF, RUF, OOO, OOO, OOO};
  table[DB] = {LDB, RDB, OOO, OOO, OOO};
  table[DF] = {LDF, RDF, OOO, OOO, OOO};
  table[UB] = {LUB, RUB, OOO, OOO, OOO};
  table[UF] = {LUF, RUF, OOO, OOO, OOO};
  table[LDB] = {LDB, OOO, OOO, OOO, OOO};
  table[LDF] = {LDF, OOO, OOO, OOO, OOO};
  table[LUB] = {LUB, OOO, OOO, OOO, OOO};
  table[LUF] = {LUF, OOO, OOO, OOO, OOO};
  table[RDB] = {RDB, OOO, OOO, OOO, OOO};
  table[RDF] = {RDF, OOO, OOO, OOO, OOO};
  table[RUB] = {RUB, OOO, OOO, OOO, OOO};
  table[RUF] = {RUF, OOO, OOO, OOO, OOO};
  return table;
}

/**
 * A LUT containing the entries for `ADJ()`.
 */
constexpr std::array<std::array<bool, 8>, 1 << 9> adjTable = createADJTable();

/**
 * A LUT containing the entries for `REFLECT()`.
 */
constexpr std::array<std::array<Octant, 8>, 1 << 9> reflectTable = createREFLECTTable();

/**
 * A LUT containing the entries for `COMMON_FACE()`.
 */
constexpr std::array<std::array<Face, 8>, 1 << 9> commonFaceTable = createCOMMONFACETable();

/**
 * A LUT containing the entries for `COMMON_EDGE()`.
 */
constexpr std::array<std::array<Edge, 8>, 1 << 9> commonEdgeTable = createCOMMONEDGETable();

/**
 * A LUT containing the entries for `getOppositeDirection()`.
 */
constexpr std::array<Any, 1 << 9> oppositeDirectionTable = createOppositeDirectionTable();

/**
 * A LUT containing the entries for `getAllowedDirections()`.
 */
constexpr std::array<std::array<Octant, 5>, 1 << 9> allowedDirectionsTable = createAllowedDirectionsTable();
}  // namespace

/**
 * This function implements a LUT obtained from the Samet paper:
 * "ADJ(I,O) is true if and only if octant O is adjacent to the Ith face, edge, or vertex of O's containing block"
 *
 * @param direction The direction I to search the adjacent neighbor
 * @param octant The octant O from the paper
 * @return true if the neighbor is adjacent, false otherwise
 */
inline bool ADJ(Any direction, Vertex octant) {
  // Check the argument preconditions
  if ((not isFace(direction)) and (not isEdge(direction)) and (not isVertex(direction))) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid direction");
  }

  if (not isVertex(octant)) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid octant");
  }

  int flatOctant = vertexToIndex(octant);
  bool result = adjTable[direction][flatOctant];
  return result;
}

/**
 * This function implements a LUT obtained from the Samet paper:
 * "REFLECT(I,O) yields the SONTYPE value of the block of equal size (not necessarily a brother) that shares the Ith
 * face, edge, or vertex of a block having SONTYPE value O"
 *
 * @param direction The direction I to search the reflected neighbor
 * @param octant The octant O
 * @return The octant resulting from the reflection
 */
inline Octant REFLECT(Any direction, Octant octant) {
  // Check the argument preconditions
  if ((not isFace(direction)) and (not isEdge(direction)) and (not isVertex(direction))) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid direction");
  }

  if (not isVertex(octant)) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid octant");
  }

  int flatOctant = vertexToIndex(octant);
  Octant result = reflectTable[direction][flatOctant];
  return result;
}

/**
 * This function implements a LUT obtained from the Samet paper:
 * "COMMON_FACE(I,O) yields the type on the face (i.e., label), of O's containing block, that is common to octant O and
 * its neighbor in the Ith direction."
 *
 * @param direction The direction I
 * @param octant The octant O
 * @return The face that is shared between the neighbors
 */
inline Face COMMON_FACE(Any direction, Vertex octant) {
  // Check the argument preconditions
  if ((not isEdge(direction)) and (not isVertex(direction))) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid direction");
  }

  if (not isVertex(octant)) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid octant");
  }

  int flatOctant = vertexToIndex(octant);
  Face result = commonFaceTable[direction][flatOctant];
  return result;
}

/**
 * This function implements a LUT obtained from the Samet paper:
 * "COMMON_EDGE(I,O) yields the type of the edge (i.e., label), of O's containing block, that is common to octant O and
 * its neighbor in the Ith direction."
 *
 * @param direction The direction I
 * @param octant The octant O
 * @return The edge that is shared between the neighbors
 */
inline Edge COMMON_EDGE(Any direction, Vertex octant) {
  // Check the argument preconditions
  if (not isVertex(direction)) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid direction");
  }

  if (not isVertex(octant)) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid octant");
  }

  int flatOctant = vertexToIndex(octant);
  Edge result = commonEdgeTable[direction][flatOctant];

  if (not(isEdge(result) or result == OO)) {
    throw std::runtime_error("[OctreeDirection.h] Invalid output");
  }

  return result;
}

/**
 * Convert any direction to a direction that is directly opposing the given direction.
 *
 * @param direction Any direction (Face, Edge or Vertex)
 * @return A direction that is opposing the given direction
 */
inline Any getOppositeDirection(Any direction) {
  using namespace autopas;

  if ((not isFace(direction)) and (not isEdge(direction)) and (not isVertex(direction))) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid direction");
  }

  Any result = oppositeDirectionTable[direction];

  if (not(isFace(result) or isEdge(result) or isVertex(result))) {
    throw std::runtime_error("[OctreeDirection.h] Invalid output");
  }

  return result;
}

/**
 * Get a list of octants that are along the given direction.
 *
 * @param along A direction the returned vertices should be admissible to
 * @return A list of octants that fits the given direction
 */
inline std::vector<Octant> getAllowedDirections(Any along) {
  using namespace autopas;

  if ((not isFace(along)) and (not isEdge(along)) and (not isVertex(along))) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid direction");
  }

  auto resultWithOOO = allowedDirectionsTable[along];
  std::vector<Octant> result;
  // At most, there should be 4 elements in the result vector.
  result.reserve(4);
  for (auto v : resultWithOOO) {
    if (v == OOO) break;
    result.push_back(v);
  }

  // Check post-conditions
  for (auto v : result) {
    if (not isVertex(v)) {
      throw std::runtime_error("[OctreeDirection.h] Result contains an illegal vertex.");
    }
  }

  return result;
}

}  // namespace autopas::octree
