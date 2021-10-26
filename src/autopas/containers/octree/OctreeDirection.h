/**
 * @file OctreeDirection.h
 *
 * @author Johannes Spies
 * @date 18.05.2021
 */

#pragma once

namespace autopas {

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
  static_assert((f1 != f2) && (f2 != f3), "Faces must be different");
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

namespace Faces {
/**
 * All available faces for a cube. The "omega" face `O` is excluded.
 */
constexpr static std::array<Face, 6> table = {L, R, D, U, B, F};
}  // namespace Faces

/**
 * Get all available edges for a cube.
 *
 * @return A pointer to a static table. The last element is the "invalid" edge OO.
 */
inline Edge *getEdges() {
  static Edge table[] = {LD, LU, LB, LF, RD, RU, RB, RF, DB, DF, UB, UF, OO};
  return table;
}

/**
 * Get all available vertices for a cube.
 *
 * @return A pointer to a static table. The last element is the "invalid" vertex OOO.
 */
inline Vertex *VERTICES() {
  static Vertex table[] = {LDB, LDF, LUB, LUF, RDB, RDF, RUB, RUF, OOO};
  return table;
}

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
    if (vertex == VERTICES()[i]) {
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
 * Check if a given list of elements contains the element to test.
 *
 * @tparam T The base type of the list elements
 * @param all A pointer to the start element in the list
 * @param stop A pointer to the element which is one past the last element in the list
 * @param test The element to check
 * @return true if the list contains the test element, false otherwise
 */
template <typename T>
inline bool contains(T *all, T stop, Any test) {
  bool result = false;
  for (T *t = all; *t != stop; ++t) {
    if (test == *t) {
      result = true;
      break;
    }
  }
  return result;
}

/**
 * Check if f is a face.
 *
 * @param f The parameter to check
 * @return true iff f is in the list returned from getFaces()
 */
template <typename T>
inline bool isFace(T f) {
  return std::find(Faces::table.begin(), Faces::table.end(), f) != Faces::table.end();
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
  return contains(getEdges(), OO, e);
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
  return contains(VERTICES(), OOO, v);
}

/**
 * This function implements a LUT obtained from the Samet paper:
 * "ADJ(I,O) is true if and only if octant O is adjacent to the Ith face, edge, or vertex of O's containing block"
 *
 * @param direction The direction I to search the adjacent neighbor
 * @param octant The octant O from the paper
 * @return true if the neighbor is adjacent, false otherwise
 */
inline bool ADJ(Any direction, Vertex octant) {
  static std::array<std::array<bool, 8>, 1 << 9> table;

  // Check the argument preconditions
  if (!isFace(direction) && !contains(getEdges(), OO, direction) && !contains(VERTICES(), OOO, direction)) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid direction");
  }

  if (!contains(VERTICES(), OOO, octant)) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid octant");
  }

  // Initialize if the first element is not present
  if (!table[L][0]) {
#pragma omp critical
    {
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
    }
  }

  int flatOctant = vertexToIndex(octant);
  bool result = table[direction][flatOctant];
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
  static std::array<std::array<Octant, 8>, 1 << 9> table;

  // Check the argument preconditions
  if (!isFace(direction) && !contains(getEdges(), OO, direction) && !contains(VERTICES(), OOO, direction)) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid direction");
  }

  if (!contains(VERTICES(), OOO, octant)) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid octant");
  }

  // Initialize if the first element is not present
  if (table[L][0] != RDB) {
#pragma omp critical
    {
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
    }
  }

  int flatOctant = vertexToIndex(octant);
  Octant result = table[direction][flatOctant];
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
  static std::array<std::array<Face, 8>, 1 << 9> table;

  // Check the argument preconditions
  if (!contains(getEdges(), OO, direction) && !contains(VERTICES(), OOO, direction)) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid direction");
  }

  if (!contains(VERTICES(), OOO, octant)) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid octant");
  }

  if (table[LD][2] != L) {
#pragma omp critical
    {
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
    }
  }

  int flatOctant = vertexToIndex(octant);
  Face result = table[direction][flatOctant];
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
  static std::array<std::array<Edge, 8>, 1 << 9> table;

  // Check the argument preconditions
  if (!contains(VERTICES(), OOO, direction)) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid direction");
  }

  if (!contains(VERTICES(), OOO, octant)) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid octant");
  }

  if (table[LDB][2] != LB) {
#pragma omp critical
    {
      table[LDB] = {OO, LD, LB, OO, DB, OO, OO, OO};
      table[LDF] = {LD, OO, OO, LF, OO, DF, OO, OO};
      table[LUB] = {LB, OO, OO, LU, OO, OO, UB, OO};
      table[LUF] = {OO, LF, LU, OO, OO, OO, OO, UF};
      table[RDB] = {DB, OO, OO, OO, OO, RD, RB, OO};
      table[RDF] = {OO, DF, OO, OO, RD, OO, OO, RF};
      table[RUB] = {OO, OO, UB, OO, RB, OO, OO, RU};
      table[RUF] = {OO, OO, OO, UF, OO, RF, RU, OO};
    }
  }

  int flatOctant = vertexToIndex(octant);
  Edge result = table[direction][flatOctant];

  if (not(contains(getEdges(), OO, result) or result == OO)) {
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
inline autopas::Any getOppositeDirection(autopas::Any direction) {
  using namespace autopas;

  if (!isFace(direction) && !isEdge(direction) && !isVertex(direction)) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid direction");
  }

  static std::array<Any, 1 << 9> table = {};
  if (table[L] != R) {
#pragma omp critical
    {
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
    }
  }

  Any result = table[direction];

  if (not(isFace(result) or contains(getEdges(), OO, result) or contains(VERTICES(), OOO, result))) {
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
inline std::vector<autopas::Octant> getAllowedDirections(autopas::Any along) {
  using namespace autopas;

  if (!isFace(along) && !contains(getEdges(), OO, along) && !contains(VERTICES(), OOO, along)) {
    throw std::runtime_error("[OctreeDirection.h] Received invalid direction");
  }

  static std::array<std::vector<Octant>, 1 << 9> table = {};
  if (table[L].empty()) {
#pragma omp critical
    {
      table[L] = {LDB, LDF, LUB, LUF};
      table[R] = {RDB, RDF, RUB, RUF};
      table[D] = {LDB, LDF, RDB, RDF};
      table[U] = {LUB, LUF, RUB, RUF};
      table[B] = {LDB, LUB, RDB, RUB};
      table[F] = {LDF, LUF, RDF, RUF};
      table[LD] = {LDB, LDF};
      table[LU] = {LUB, LUF};
      table[LB] = {LDB, LUB};
      table[LF] = {LDF, LUF};
      table[RD] = {RDB, RDF};
      table[RU] = {RUB, RUF};
      table[RB] = {RDB, RUB};
      table[RF] = {RDF, RUF};
      table[DB] = {LDB, RDB};
      table[DF] = {LDF, RDF};
      table[UB] = {LUB, RUB};
      table[UF] = {LUF, RUF};
      table[LDB] = {LDB};
      table[LDF] = {LDF};
      table[LUB] = {LUB};
      table[LUF] = {LUF};
      table[RDB] = {RDB};
      table[RDF] = {RDF};
      table[RUB] = {RUB};
      table[RUF] = {RUF};
    }
  }

  auto result = table[along];

  // Check post-conditions
  for (auto v : result) {
    if (!contains(VERTICES(), OOO, v)) {
      throw std::runtime_error("[OctreeDirection.h] Result contains an illegal vertex.");
    }
  }

  return result;
}

}  // namespace autopas
