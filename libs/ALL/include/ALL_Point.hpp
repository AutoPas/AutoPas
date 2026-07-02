/*
Copyright 2018-2020 Rene Halver, Forschungszentrum Juelich GmbH, Germany
Copyright 2018-2020 Godehard Sutmann, Forschungszentrum Juelich GmbH, Germany

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
   other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef ALL_POINT_HEADER_INC
#define ALL_POINT_HEADER_INC

#include "ALL_CustomExceptions.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace ALL {

template <class T> class Point;
template <typename T>
std::ostream &operator<<(std::ostream &, const Point<T> &);

/// @tparam T floating point type used, has to be identical to the type used for
/// the vertices and borders
template <class T> class Point {
public:
  /// default constructor
  Point() : dimension(0) {}
  /// constructor
  /// @param d dimension of the point object
  Point(const int d) : dimension(d) {
    coordinates.resize(d);
    weight = (T)1;
  }
  /// constructor with initialization
  /// @param d dimension of the point object
  /// @param data positions
  Point(const int d, const T *data) : Point<T>(d) {
    // copy each element of the array into the vector (insecure, no boundary
    // checks for data!)
    for (auto i = 0; i < d; ++i)
      coordinates.at(i) = data[i];
  }

  /// constructor with initialization
  /// @param data positions, dimension is the size of the
  /// vector
  Point(const std::vector<T> &data) {
    // initialize the coordinates vector with the data vector
    coordinates.insert(coordinates.begin(), data.begin(), data.end());
    // update the dimension with the size of the data vector
    dimension = data.size();
  }

  /// constructor with initialization
  /// @param d dimension of the point
  /// @param data positions
  /// @param w weight of the point
  Point(const int d, const T *data, const T w) : Point<T>(d, data) {
    weight = w;
  }

  /// constructor with initialization
  /// @param d dimension of the point
  /// @param data positions
  /// @param w weight of the point
  /// @param i index of the point
  Point(const int d, const T *data, const T w, const long i)
      : Point<T>(d, data, w) {
    id = i;
  }

  /// constructor with initialization
  /// @param data positions, dimension is the size of the
  /// vector
  /// @param w weight of the point
  Point(const std::vector<T> &data, const T w) : Point<T>(data) {
    weight = w;
  }

  /// constructor with initialization
  /// @param data positions, dimension is the size of the
  /// vector
  /// @param w weight of the point
  /// @param i index of the point
  Point(const std::vector<T> &data, const T w, const long i)
      : Point<T>(data, w) {
    id = i;
  }

  /// destructor
  ~Point(){};

  // TODO: return values to check if operation was successful?
  /// method to change the dimension of the point object
  /// @param d new dimension
  void setDimension(const int d) {
    dimension = d;
    coordinates.resize(d);
  }

  /// method to change the weight of the point object
  /// @param w new weight
  void set_weight(const T w) { weight = w; }

  /// method to change the index of the particle
  /// @param i new index
  void set_id(const long i) { id = i; }

  /// access operator to access an element of the point object
  /// @param idx the index of the element to be accessed
  /// @result reference to the indexed element
  T &operator[](const std::size_t idx) { return coordinates.at(idx); }

  /// access operator to access an element of the constant point object
  /// @param idx the index of the element to be accessed
  /// @result const reference to the indexed element
  const T &operator[](const std::size_t idx) const {
    return coordinates.at(idx);
  }

  /// method to get the weight of the point object
  /// @return the weight of the point object
  T getWeight() const { return weight; }

  /// method to get the index of the point object
  /// @return the index of the point object
  long get_id() const { return id; }

  /// method to get the dimension of the point object
  /// @return the dimension of the point object
  int getDimension() const { return dimension; }

  /// method to compute the norm of the vector described by the
  /// point object
  /// @param nd type of the norm (default: 2, i.e. the Euclidean norm)
  /// @return norm of the vector
  T norm(T nd = 2) {
    T res = 0;
    for (int d = 0; d < dimension; ++d)
      res += std::pow(coordinates.at(d), nd);
    return std::pow(res, 1.0 / nd);
  }

  /// method to compute the Euclidean distance (two-norm) between the local and
  /// the provided point object
  /// @param p the point object for which the distence to the local point object
  /// is computed
  /// @return the distance between the two points
  T d(Point<T> p) {
    int d_p = p.getDimension();
    if (d_p != dimension)
      throw PointDimensionMissmatchException(__FILE__, __func__, __LINE__);
    return dist(p).norm();
  }

  /// method to compute the Manhatten / city-block distance (one-norm) between
  /// the local and the provided point object
  /// @param p the point object for which the distance to the local point object
  /// is computed
  /// @return the distance between the two points
  T d_1(Point<T> p) {
    int d_p = p.getDimension();
    if (d_p != dimension)
      throw PointDimensionMissmatchException(__FILE__, __func__, __LINE__);
    return dist(p, 1).norm();
  }

  /// method to compute the distance vector between the local point object and
  /// the provided point object
  /// @param p the point object for which the distance vector to the local point
  /// object is computed
  /// @return the distance vector between the two points
  Point<T> dist(Point<T> &p) {
    int d_p = p.getDimension();
    if (d_p != dimension)
      throw PointDimensionMissmatchException(__FILE__, __func__, __LINE__);
    std::vector<T> d_v;
    for (int d = 0; d < dimension; ++d)
      d_v.push_back(p[d] - coordinates.at(d));

    return Point<T>(d_v);
  }

  /// method to compute the distance of the local point object from a plane
  /// spanned by provided points
  /// @param A anchor point for the plane
  /// @param B anchor point for the plane
  /// @param C anchor point for the plane
  /// @return distance between local point object and plane
  T dist_plane(const Point<T> &A, const Point<T> &B,
               const Point<T> &C) {

    if (A.getDimension() != dimension || B.getDimension() != dimension ||
        C.getDimension() != dimension)
      throw PointDimensionMissmatchException(__FILE__, __func__, __LINE__);

    // vectors spanning the plane from vertex 'a'
    Point<T> vb = B - A;
    Point<T> vc = C - A;

    // normal vector of plane
    Point<T> n = vb.cross(vc);
    n = n/n.norm();

    // return r.n - A.n = (r-A).n as distance from plane to r
    return std::abs( (coordinates.at(0) - A[0]) * n[0] +
		     (coordinates.at(1) - A[1]) * n[1] +
		     (coordinates.at(2) - A[2]) * n[2] );
  }

  /// operator for the addition of two point objects
  /// @param rhs point to add the local point object to
  /// @return sum of local point object and provided point object
  Point<T> operator+(const Point<T> &rhs) const {
    if (rhs.getDimension() != dimension) {
      throw PointDimensionMissmatchException(__FILE__, __func__, __LINE__);
    }
    Point<T> result(dimension);
    for (int d = 0; d < dimension; ++d) {
      result[d] = coordinates.at(d) + rhs[d];
    }
    return result;
  }

  /// operator for the addition of two point objects
  /// @param rhs point to subtract from the local point object
  /// @return difference vector between local point object and provided point
  /// object
  Point<T> operator-(const Point<T> &rhs) const {
    if (rhs.getDimension() != dimension) {
      throw PointDimensionMissmatchException(__FILE__, __func__, __LINE__);
    }
    Point<T> result(dimension);
    for (int d = 0; d < dimension; ++d) {
      result[d] = coordinates.at(d) - rhs[d];
    }
    return result;
  }

  /// operator to compute the dot product between two point objects
  /// @param rhs point object to compute the dot product with
  /// @return dot product between the two point objects
  T operator*(const Point<T> &rhs) const {
    if (rhs.getDimension() != dimension) {
      throw PointDimensionMissmatchException(__FILE__, __func__, __LINE__);
    }
    T result = (T)0.0;
    for (int d = 0; d < dimension; ++d) {
      result += coordinates.at(d) * rhs[d];
    }
    return result;
  }

  /// operator to scale the local point object by a provided factor
  /// @param rhs scaling factor
  /// @return scaled point object
  Point<T> operator*(const T &rhs) const {
    Point<T> result(dimension);
    for (int d = 0; d < dimension; ++d) {
      result[d] = coordinates.at(d) * rhs;
    }
    return result;
  }

  /// operator to scale the local point object by a provided factor
  /// @param rhs scaling factor
  /// @return scaled point object
  Point<T> operator/(const T &rhs) const {
    Point<T> result(dimension);
    for (int d = 0; d < dimension; ++d) {
      result[d] = coordinates.at(d) / rhs;
    }
    return result;
  }

  /// operator to compute the cross product between two point objects
  /// @param rhs point object to compute the cross product with
  /// @return cross product between the two point objects
  /// @attention only works for 3D points / vectors
  Point<T> cross(const Point<T> &rhs) const {
    if (rhs.getDimension() != dimension && dimension != 3) {
      throw PointDimensionMissmatchException(__FILE__, __func__, __LINE__);
    }
    Point<T> result(dimension);
    for (int d = 0; d < dimension; ++d) {
      result[d] = coordinates.at((d + 1) % 3) * rhs[(d + 2) % 3] -
                  coordinates.at((d + 2) % 3) * rhs[(d + 1) % 3];
    }
    return result;
  }

  /// method to determine if the local point object has the same orientation to
  /// a plane spanned by point objects A,B,C as the provided point P
  /// @param A anchor point A
  /// @param B anchor point B
  /// @param C anchor point C
  /// @param P reference point P
  /// @return the local point object has the same orientation to the plane as
  /// @attention if the reference point is located within the plane, it will not
  /// have the same orientation as the local point! the reference point
  bool same_side_plane(const Point<T> &A, const Point<T> &B,
                       const Point<T> &C, const Point<T> &P) {
    // compute difference vectors:
    Point<T> BA = B - A;
    Point<T> CA = C - A;
    Point<T> PA = P - A;
    Point<T> tA = *(this) - A;

    // compute normal vector of plane
    Point<T> n = CA.cross(BA);

    // compute scalar product of distance
    // vectors with normal vector
    T PAn = PA * n;
    T tAn = tA * n;

    return ((PAn * tAn) > 0);
  }

  /// method to check if the local point object is inside a tetrahedron
  /// described by the vertices A,B,C and D
  /// @param A vertex A
  /// @param B vertex B
  /// @param C vertex C
  /// @param D vertex D
  /// @return the point is within the tetrahedron
  /// @attention a point that is on the surface of the tetrahedron is not inside
  /// it
  bool inTetrahedron(const Point<T> &A, const Point<T> &B,
                     const Point<T> &C, const Point<T> &D) {
    /// for all surfaces of the tetrahedron check if the local point has the
    /// same orientation as the remaining vertex of the tetrahedron, to be
    /// inside the tetrahedron that must be fulfilled
    return (same_side_plane(A, B, C, D) && same_side_plane(B, C, D, A) &&
            same_side_plane(C, D, A, B) && same_side_plane(A, D, B, C));
  };

private:
  /// dimension of the point object
  int dimension;
  /// index for the point
  long id;
  // array containg the coordinates
  std::vector<T> coordinates;
  /// weight of the point to be used if points are not to be equally considered
  /// in computations, e.g. if number of interactions are considered
  T weight;
};

/// output operator for a point object
/// @param os output stream the point object should be printed to
/// @param p the point object to be printed
/// @return the modified output stream
template <class T>
std::ostream &operator<<(std::ostream &os, const Point<T> &p) {
  for (int i = 0; i < p.getDimension(); ++i)
    os << p[i] << " ";
  os << p.getWeight();
  return os;
}

/// operator to scale a point by a factor (reversed way of notation T *
/// Point<T>)
/// @param lhs scaling factor
/// @param rhs the point object to be scaled
/// @param scaled point object
template <class T>
Point<T> operator*(const T &lhs, const Point<T> &rhs) {
  return rhs * lhs;
}

// template <class T>
// Point<T> operator/(const T &lhs, const Point<T> &rhs) {
//   return rhs * ((T)1 / lhs);
// }

}//namespace ALL

#endif
