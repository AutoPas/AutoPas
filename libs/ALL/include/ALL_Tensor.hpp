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

#ifndef ALL_TENSOR_HEADER_INCLUDED
#define ALL_TENSOR_HEADER_INCLUDED

#include "ALL_CustomExceptions.hpp"
#include "ALL_Defines.h"
#include "ALL_Functions.hpp"
#include "ALL_LB.hpp"
#include <exception>
#include <mpi.h>

namespace ALL {

/// Load-balancing scheme based on a independent local scheme, where the
/// workloads of domains are reduced over their cartesian position, on reduction
/// for each of the dimensions. These reduced values are then compared against
/// the values gathered for the neighbor domains. Then the boundaries between
/// the process groups are adjusted according to the ratio of cummulated work
/// loads on each of them. The border is shifted in the direction of the domain
/// with more cummulated work. This is independently done for each of the
/// dimensions.
/// @tparam T data for vertices and related data
/// @tparam W data for work and related data
template <class T, class W> class Tensor_LB : public LB<T, W> {
public:
  /// default constuctor
  Tensor_LB() {}
  /// constructor to initialize values
  /// @param[in] d dimension of the vertices used
  /// @param[in] w the local (scalar) work load
  /// @param[in] g the correction factor
  Tensor_LB(int d, W w, T g) : LB<T, W>(d, g) {
    this->setWork(w);
    // need the lower and upper bounds
    // of the domain for the tensor-based
    // load-balancing scheme

    // array of MPI communicators for each direction (to collect work
    // on each plane)
    communicators.resize(d);
    nNeighbors.resize(2 * d);
  }

  /// default destructor
  ~Tensor_LB() = default;

  /// setup internal data structures and parameters
  void setup() override;

  /// method to execute a load-balancing step
  /// @param[in] step number of the load-balancing step
  virtual void balance(int step) override {balance(step, MPI_SUM);}

  // getter for variables (passed by reference to avoid
  // direct access to private members of the object)

  /// method to provide a list of the neighbors of the local domain
  /// @param[out] list reference to a std::vector of integers where the list of
  /// neighbors will be assigned to
  virtual std::vector<int> &getNeighbors() override;

  /// method to set specific data structures (unused for tensor grid method)
  /// @param[in] data pointer to the data structure
  virtual void setAdditionalData(known_unused const void *data) override {}

  /// method to execute a load-balancing step
  /// @param[in] step number of the load-balancing step
  virtual void balance(int step, MPI_Op reductionMode);

  /// method to get an estimated work distribution after the balance step
  /// (currently only implemented in ALL::HISTOGRAM!)
  /// @param [out] double providing the estimated LB after the balance step
  virtual W getEstimatedEfficiency() override {return (W)-1;};

private:
  // type for MPI communication
  MPI_Datatype MPIDataTypeT;
  MPI_Datatype MPIDataTypeW;

  // array of MPI communicators for each direction (to collect work
  // on each plane)
  std::vector<MPI_Comm> communicators;

  // list of neighbors
  std::vector<int> neighbors;
  std::vector<int> nNeighbors;
  
protected:


};

// setup routine for the tensor-based load-balancing scheme
// requires:
//              this->globalComm (int): cartesian MPI communicator, from
//                                 which separate sub communicators
//                                 are derived in order to represent
//                                 each plane of domains in the system
template <class T, class W> void Tensor_LB<T, W>::setup() {
  int status;

  // check if Communicator is cartesian
  MPI_Topo_test(this->globalComm, &status);
  if (status != MPI_CART) {
    throw InvalidCommTypeException(
        __FILE__, __func__, __LINE__,
        "Cartesian MPI communicator required, passed communicator not "
        "cartesian");
  }

  int dim = this->getDimension();

  // get the local coordinates, periodicity and global size from the MPI
  // communicator
  MPI_Cart_get(this->globalComm, dim, this->global_dims.data(),
               this->periodicity.data(), this->local_coords.data());

  // get the local rank from the MPI communicator
  MPI_Cart_rank(this->globalComm, this->local_coords.data(), &this->localRank);

  // create sub-communicators
  for (int i = 0; i < dim; ++i) {
    MPI_Comm_split(this->globalComm, this->local_coords.at(i), 0,
                   &communicators[i]);
  }

  // determine correct MPI data type for template T
  if (std::is_same<T, double>::value)
    MPIDataTypeT = MPI_DOUBLE;
  else if (std::is_same<T, float>::value)
    MPIDataTypeT = MPI_FLOAT;
  else if (std::is_same<T, int>::value)
    MPIDataTypeT = MPI_INT;
  else if (std::is_same<T, long>::value)
    MPIDataTypeT = MPI_LONG;

  // determine correct MPI data type for template W
  if (std::is_same<W, double>::value)
    MPIDataTypeW = MPI_DOUBLE;
  else if (std::is_same<W, float>::value)
    MPIDataTypeW = MPI_FLOAT;
  else if (std::is_same<W, int>::value)
    MPIDataTypeW = MPI_INT;
  else if (std::is_same<W, long>::value)
    MPIDataTypeW = MPI_LONG;

  // calculate neighbors
  int rank_left, rank_right;

  neighbors.clear();
  for (int i = 0; i < dim; ++i) {
    MPI_Cart_shift(this->globalComm, i, 1, &rank_left, &rank_right);
    neighbors.push_back(rank_left);
    neighbors.push_back(rank_right);
    nNeighbors[2 * i] = 1;
    nNeighbors[2 * i + 1] = 1;
  }
}

template <class T, class W> void Tensor_LB<T, W>::balance(int step, MPI_Op reductionMode) {
  int dim = this->getDimension();
  std::vector<Point<T>> newVertices = this->vertices;
  this->prevVertices = this->vertices;

  bool localIsSum = reductionMode == MPI_SUM;
  bool localIsMax = reductionMode == MPI_MAX;

  // loop over all available dimensions
  for (int i = 0; i < dim; ++i) {
    W work_local_plane;
    // collect work from all processes in the same plane
    MPI_Allreduce(this->getWork().data(), &work_local_plane, 1, MPIDataTypeW,
                  reductionMode, communicators.at(i));

    // correct right border:

    W remote_work;
    T local_size;
    T remote_size;
    // determine neighbors
    int rank_left, rank_right;

    MPI_Cart_shift(this->globalComm, i, 1, &rank_left, &rank_right);

    // collect work from right neighbor plane
    MPI_Request sreq, rreq;
    MPI_Status sstat, rstat;

    MPI_Irecv(&remote_work, 1, MPIDataTypeW, rank_right, 0, this->globalComm,
              &rreq);
    MPI_Isend(&work_local_plane, 1, MPIDataTypeW, rank_left, 0,
              this->globalComm, &sreq);
    MPI_Wait(&sreq, &sstat);
    MPI_Wait(&rreq, &rstat);

    // collect size in dimension from right neighbor plane

    local_size = this->prevVertices.at(1)[i] - this->prevVertices.at(0)[i];

    MPI_Irecv(&remote_size, 1, MPIDataTypeT, rank_right, 0, this->globalComm,
              &rreq);
    MPI_Isend(&local_size, 1, MPIDataTypeT, rank_left, 0, this->globalComm,
              &sreq);
    MPI_Wait(&sreq, &sstat);
    MPI_Wait(&rreq, &rstat);

    // automatic selection of gamma to guarantee stability of method (choice of
    // user gamma ignored)
    this->gamma =
        std::max(4.1, 2.0 * (1.0 + std::max(local_size, remote_size) /
                                       std::min(local_size, remote_size)));

    T shift = Functions::borderShift1d(
        rank_right, this->local_coords.at(i), this->global_dims.at(i),
        work_local_plane, remote_work, local_size, remote_size, this->gamma,
        this->minSize[i]);

    // send shift to right neighbors
    T remote_shift = (T)0;

    MPI_Irecv(&remote_shift, 1, MPIDataTypeT, rank_left, 0, this->globalComm,
              &rreq);
    MPI_Isend(&shift, 1, MPIDataTypeT, rank_right, 0, this->globalComm, &sreq);
    MPI_Wait(&sreq, &sstat);
    MPI_Wait(&rreq, &rstat);

    // if a left neighbor exists: shift left border
    if (rank_left != MPI_PROC_NULL && this->local_coords[i] != 0)
      newVertices.at(0)[i] = this->prevVertices.at(0)[i] + remote_shift;
    else
      newVertices.at(0)[i] = this->prevVertices.at(0)[i];

    // if a right neighbor exists: shift right border
    if (rank_right != MPI_PROC_NULL &&
        this->local_coords[i] != this->global_dims[i] - 1)
      newVertices.at(1)[i] = this->prevVertices.at(1)[i] + shift;
    else
      newVertices.at(1)[i] = this->prevVertices.at(1)[i];
  }

  this->setVertices(newVertices);
}

// provide list of neighbors
template <class T, class W> std::vector<int> &Tensor_LB<T, W>::getNeighbors() {
  return neighbors;
}

/// Load-balancing scheme based on a independent local scheme, where the
/// workloads of domains are reduced over their cartesian position, on reduction
/// for each of the dimensions. These reduced values are then compared against
/// the values gathered for the neighbor domains. Then the boundaries between
/// the process groups are adjusted according to the ratio of cummulated work
/// loads on each of them. The border is shifted in the direction of the domain
/// with more cummulated work. This is independently done for each of the
/// dimensions. For this method, the reduction is a maximum over the respective
/// dimension, instead of a sum.
/// @tparam T data for vertices and related data
/// @tparam W data for work and related data
template <class T, class W> class TensorMax_LB : public Tensor_LB<T, W> {
public:
  TensorMax_LB() {}
  TensorMax_LB(int d, W w, T g) : Tensor_LB<T, W>(d, w, g) {}
  virtual void balance(int step) override {Tensor_LB<T,W>::balance(step, MPI_MAX);}
};

} // namespace ALL

#endif
