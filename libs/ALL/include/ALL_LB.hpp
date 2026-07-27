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

#ifndef ALL_LB_HEADER_INCLUDED
#define ALL_LB_HEADER_INCLUDED

#include "ALL_Point.hpp"
#include <mpi.h>
#include <vector>
#include <iostream>

namespace ALL {

/// @tparam T data for vertices and related data
/// @tparam W data for work and related data
template <class T, class W> class LB {
public:
  /// constructor for the basic load-balancing class, sets up parameters
  /// required by all balancing methods
  /// @param[in] dim dimension of the vertices used
  /// @param[in] g correction factor for the computation of shifts
  LB(const int dim, const T g) : gamma(g) {
    // set allowed minimum size to zero
    setDimension(dim);
    minSize = std::vector<T>(dim, 0.0);
    periodicity.resize(dim);
    local_coords.resize(dim);
    global_dims.resize(dim);
    neighborVertices.resize(0);
  }

  /// destructor
  virtual ~LB() = default;

  /// abstract definition of the setup method
  virtual void setup() = 0;

  /// abstract definition of the balancing method
  virtual void balance(const int step) = 0;

  /// abstract definition of the method to get the neighbors of the local domain
  /// @result std::vector<int> to store the MPI ranks of the neighbors
  /// to
  virtual std::vector<int> &getNeighbors() = 0;

  /// method to update the vertices used for the balancing step, overwrites old
  /// set of vertices
  /// @param[in] vertices_in vector containg the new vertices to be used
  virtual void setVertices(const std::vector<Point<T>> &vertices_in) {
    // as a basic assumption the number of resulting vertices
    // is equal to the number of input vertices:
    // exceptions:
    //      VORONOI
    vertices = vertices_in;
    prevVertices = vertices_in;
  }

  /// method to set the dimension of the vertices
  /// @param[in] d the dimension to be used
  /// @attention most methods currently support 3D vertices only
  virtual void setDimension(const int d) { 
    if (d != 3)
      throw InvalidArgumentException(
        __FILE__, __func__, __LINE__,
        "Currently only three dimensional vertices are supported.");
    else
        dimension = d; 
  }

  /// method to get the dimension of the vertices
  /// @return the dimension of the vertices
  virtual int getDimension() { return dimension; }

  /// method to set the correction value gamma
  /// @param[in] g the correction value to use
  virtual void setGamma(const T g) { gamma = g; }

  /// method to get the correction value currently used
  /// @return the current correction value
  virtual const T getGamma() { return gamma; }

  /// method to set a multi-dimensional work for the local domain
  /// @param[in] w vector containing all the dimensions of the work to be used
  /// for the local domain
  virtual void setWork(const std::vector<W> &w) { work = w; }

  /// method to set a scalar work for the local domain
  /// @param[in] w value containing the work to be used for the local domain
  virtual void setWork(const W w) {
    work.resize(1);
    work.at(0) = w;
  }

  /// method to set the minimum domain size in each dimension
  /// @param[in] minSize the minimum size of a domain in all dimensions
  virtual void setMinDomainSize(const std::vector<T> &minSize) {
    this->minSize = minSize;
  }

  /// method to get the minimum domain size the balancing methods have to obey
  /// @return the minimum domain size allowed
  virtual const std::vector<T> &getMinDomainSize() { return minSize; }

  /// method to set the (orthogonal) size of the system
  /// @param[in] sysSize system size in all dimensions
  virtual void setSysSize(const std::vector<T> &newSysSize) {
    sysSize = newSysSize;
  }

  /// method to get the currently stored system size
  /// @return the currently stored system size
  virtual const std::vector<T> &getSysSize() { return sysSize; }

  /// method to get the currently stored work for the local process
  /// @return the currently stored work
  /// @attention always returns a std::vector, for scalar work, it has length 1
  virtual std::vector<W> &getWork() { return work; }

  /// method to get result vertices
  /// @return the resulting vertices after a balancing step
  virtual std::vector<Point<T>> &getVertices() { return vertices; }

  /// method to get the original vertices before the last balancing step
  /// @return the original unmodified vertices before the last balancing step
  /// @attention since the original and resulting vertices are set up in the
  /// balancing method, the result is undefined before calling the balancing
  /// method at least once
  virtual std::vector<Point<T>> &getPrevVertices() { return prevVertices; };

  /// method to set the MPI communicator to be used by the balancing method
  /// @param[in] comm the MPI communicator to be used
  virtual void setCommunicator(const MPI_Comm comm) {
    // store communicator
    globalComm = comm;
    // get local rank within communicator
    MPI_Comm_rank(comm, &localRank);
  }

  /// method the get the number of vertices stored for the local domain
  /// @return number of local vertices
  int getNVertices() { return vertices.size(); }

  /// method to set undefined method specific data, which is not shared between
  /// different methods but needs to be set by a unified interface
  /// @param[in] data the data be passed to the balancing method
  virtual void setAdditionalData(const void *data) = 0;

  /// method to get an estimated LB efficiency after the balance step
  /// (currently only implemented in ALL::HISTOGRAM!)
  /// @return double providing the estimated LB after the balance step
  virtual W getEstimatedEfficiency() = 0;

  /// method to get the current LB efficiency with the given work distribution
  /// @return current LB efficiency
  double getEfficiency()
  {
    double localSum = 0.0;
    for (auto i = this->work.begin(); i < this->work.end(); ++i)
    {        
        localSum += *i;
    }        
    double globalMin;
    double globalMax;
    MPI_Allreduce(&localSum, &globalMin, 1, MPI_DOUBLE, MPI_MIN, this->globalComm);
    MPI_Allreduce(&localSum, &globalMax, 1, MPI_DOUBLE, MPI_MAX, this->globalComm);
    return (1.0 - (double)(globalMax - globalMin) / 
                  (double)(globalMax + globalMin));
  }


  /// method to provide a list of vertices describing the neighboring domains
  /// currently only implemented for VORONOI, as a means to get the anchor points
  /// of the surrounding domains
  std::vector<T> &getNeighborVertices() {
    return neighborVertices; }

protected:
  /// correction factor
  T gamma;
  /// dimension of the used vertices
  int dimension;
  /// used MPI communicator
  MPI_Comm globalComm;
  /// local rank in the used MPI communicator
  int localRank;
  /// minimum domain size
  /// @attention not all balancing method support this
  std::vector<T> minSize;
  /// (orthogonal) system size
  std::vector<T> sysSize;
  /// local work
  std::vector<W> work;
  /// local vertices after previous balancing step
  std::vector<Point<T>> vertices;
  /// original vertices before previous balancing step
  std::vector<Point<T>> prevVertices;
  /// dimensions of the global process grid
  std::vector<int> global_dims;
  /// cartesian coordinates of the local domain in the process grid
  std::vector<int> local_coords;
  /// periodicity of the MPI communicator / system
  std::vector<int> periodicity;
  /// vertices describing neighboring domains
  std::vector<T> neighborVertices;
  /// method the resize the vertex list
  /// @param [in] new_size the new size of the vertex list
  void resizeVertices(const int newSize) { vertices.resize(newSize); }

};

}//namespace ALL

#endif // ALL_LB_HEADER_INCLUDED
