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

#ifndef ALL_STAGGERED_HEADER_INCLUDED
#define ALL_STAGGERED_HEADER_INCLUDED

#include "ALL_CustomExceptions.hpp"
#include "ALL_LB.hpp"
#include "ALL_Defines.h"
#include <exception>
#include <mpi.h>
#include <vector>

namespace ALL {

/// Load-balancing scheme based on a hierarchical scheme, which is based on a
/// local one-dimensional equilibration between neighboring domains or
/// compositions of domains. In the first step the domains in a collective plane
/// are compared to the composite plane in the chosen direction. Depending on
/// the ratio of work in both planes the border between them is shifted towards
/// the plane with more work. In the following two steps this procedure is
/// repeated for composite columns and then single domains.
/// @tparam T data for vertices and related data
/// @tparam W data for work and related data
template <class T, class W> class Staggered_LB : public LB<T, W> {
public:
  /// default constructor
  Staggered_LB() {}
  /// constructor initializing basic parameters
  /// @param[in] d the dimension of the used vertices
  /// @param[in] w the scalar work assigned to the local domain
  /// @param[in] g the correction factor gamma
  Staggered_LB(int d, W w, T g) : LB<T, W>(d, g) {
    this->setWork(w);

    // array of MPI communicators for each direction (to collect work
    // on each plane)
    communicators.resize(d);
    for (int _d = 0; _d < d; ++_d)
      communicators.at(_d) = MPI_COMM_NULL;

    nNeighbors.resize(2 * d);
  }

  /// default destructor
  ~Staggered_LB();
  ;

  /// method to setup the method specific internal parameters
  void setup() override;

  /// method to execute a load-balancing step, which updates the resulting
  /// vertices based on the stored work assigned
  /// @param[in] step the number of the load-balancing step (used in debugging
  /// output)
  void balance(int step) override;

  /// method to provide a list of the neighbors of the local domain
  /// @param[out] list reference to a std::vector of integers where the list of
  /// neighbors will be assigned to
  virtual std::vector<int> &getNeighbors() override;

  /// method to set specific data structures (unused for staggered grid method)
  /// @param[in] data pointer to the data structure
  virtual void setAdditionalData(known_unused const void *data) override {}

  /// method to get an estimated work distribution after the balance step
  /// (currently only implemented in ALL::HISTOGRAM!)
  /// @param [out] double providing the estimated LB after the balance step
  virtual W getEstimatedEfficiency() override {return (W)-1;};

private:
  /// data type for the communication of vertex-related data
  MPI_Datatype MPIDataTypeT;
  /// data type for the communication of work-related data
  MPI_Datatype MPIDataTypeW;

  /// list of internally used MPI communicators
  std::vector<MPI_Comm> communicators;

  /// list of ranks of neighbors of the local domain
  std::vector<int> neighbors;
  /// number of neighbors in each dimension
  std::vector<int> nNeighbors;

  /// method to update the lists of neighbors depending on the current vertices
  void find_neighbors();
};

template <class T, class W> Staggered_LB<T, W>::~Staggered_LB() {
  int dimension = this->getDimension();
  for (int i = 1; i < dimension; ++i) {
    if (communicators.at(i) != MPI_COMM_NULL)
      MPI_Comm_free(&communicators.at(i));
  }
}

// setup routine for the tensor-based load-balancing scheme
// requires:
//              this->globalComm (int): cartesian MPI communicator, from
//                                 which separate sub communicators
//                                 are derived in order to represent
//                                 each plane of domains in the system
template <class T, class W> void Staggered_LB<T, W>::setup() {
  int status;
  int dimension = this->getDimension();

  // check if Communicator is cartesian
  MPI_Topo_test(this->globalComm, &status);
  if (status != MPI_CART) {
    throw InvalidCommTypeException(
        __FILE__, __func__, __LINE__,
        "Cartesian MPI communicator required, passed communicator is not "
        "cartesian");
  }

  // get the local coordinates, periodicity and global size from the MPI
  // communicator
  MPI_Cart_get(this->globalComm, dimension, this->global_dims.data(),
               this->periodicity.data(), this->local_coords.data());

  // get the local rank from the MPI communicator
  MPI_Cart_rank(this->globalComm, this->local_coords.data(), &this->localRank);

  // create sub-communicators

  if (communicators.at(1) != MPI_COMM_NULL)
    MPI_Comm_free(&communicators.at(1));
  if (communicators.at(2) != MPI_COMM_NULL)
    MPI_Comm_free(&communicators.at(2));

  // z-plane
  MPI_Comm_split(this->globalComm, this->local_coords.at(2),
                 this->local_coords.at(0) +
                     this->local_coords.at(1) * this->global_dims.at(0),
                 &communicators.at(2));

  // y-column
  MPI_Comm_split(this->globalComm,
                 this->local_coords.at(2) * this->global_dims.at(1) +
                     this->local_coords.at(1),
                 this->local_coords.at(0), &communicators.at(1));

  // only cell itself
  communicators.at(0) = MPI_COMM_SELF;

  // determine correct MPI data type for template T
  if (std::is_same<T, double>::value)
    MPIDataTypeT = MPI_DOUBLE;
  else if (std::is_same<T, float>::value)
    MPIDataTypeT = MPI_FLOAT;
  else if (std::is_same<T, int>::value)
    MPIDataTypeT = MPI_INT;
  else if (std::is_same<T, long>::value)
    MPIDataTypeT = MPI_LONG;
  else {
    throw InvalidCommTypeException(
        __FILE__, __func__, __LINE__,
        "Invalid data type for boundaries given (T)");
  }

  // determine correct MPI data type for template W
  if (std::is_same<W, double>::value)
    MPIDataTypeW = MPI_DOUBLE;
  else if (std::is_same<W, float>::value)
    MPIDataTypeW = MPI_FLOAT;
  else if (std::is_same<W, int>::value)
    MPIDataTypeW = MPI_INT;
  else if (std::is_same<W, long>::value)
    MPIDataTypeW = MPI_LONG;
  else {
    throw InvalidCommTypeException(__FILE__, __func__, __LINE__,
                                       "Invalid data type for work given (W)");
  }

  // calculate neighbors
  int rank_left, rank_right;

  neighbors.clear();
  for (int i = 0; i < dimension; ++i) {
    MPI_Cart_shift(this->globalComm, i, 1, &rank_left, &rank_right);
    neighbors.push_back(rank_left);
    neighbors.push_back(rank_right);
  }
}

template <class T, class W> void Staggered_LB<T, W>::balance(int) {
  std::vector<Point<T>> newVertices = this->vertices;
  int dimension = this->getDimension();

  // store original vertices
  this->prevVertices = this->vertices;

  // loop over all available dimensions
  for (int i = 0; i < dimension; ++i) {
    W work_local_plane;
#ifdef ALL_DEBUG_ENABLED
    MPI_Barrier(this->globalComm);
    if (this->localRank == 0)
      std::cout << "ALL::Staggered_LB::balance(): before work computation..."
                << std::endl;
#endif
    // collect work from all processes in the same plane
    MPI_Allreduce(this->work.data(), &work_local_plane, 1, MPIDataTypeW,
                  MPI_SUM, communicators[i]);

#ifdef ALL_DEBUG_ENABLED
    MPI_Barrier(this->globalComm);
    if (this->localRank == 0)
      std::cout << "ALL::Staggered_LB::balance(): before work distribution..."
                << std::endl;
#endif
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

#ifdef ALL_DEBUG_ENABLED
    MPI_Barrier(this->globalComm);
    if (this->localRank == 0)
      std::cout << "ALL::Staggered_LB::balance(): before shift calculation..."
                << std::endl;
#endif

    T shift = Functions::borderShift1d(
        rank_right, this->local_coords.at(i), this->global_dims.at(i),
        work_local_plane, remote_work, local_size, remote_size, this->gamma,
        this->minSize[i]);

#ifdef ALL_DEBUG_ENABLED
    MPI_Barrier(this->globalComm);
    if (this->localRank == 0)
      std::cout << "ALL::Staggered_LB::balance(): before shift distibution..."
                << std::endl;
#endif
    // send shift to right neighbors
    T remote_shift = (T)0;

    MPI_Irecv(&remote_shift, 1, MPIDataTypeT, rank_left, 0, this->globalComm,
              &rreq);
    MPI_Isend(&shift, 1, MPIDataTypeT, rank_right, 0, this->globalComm, &sreq);
    MPI_Wait(&sreq, &sstat);
    MPI_Wait(&rreq, &rstat);

#ifdef ALL_DEBUG_ENABLED
    MPI_Barrier(this->globalComm);
    std::cout << this->localRank << ": shift = " << shift
              << " remoteShift = " << remote_shift
              << " vertices: " << this->vertices.at(0)[i] << ", "
              << this->vertices.at(1)[i] << std::endl;
#endif

    // for now: test case for simple program

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

    // check if vertices are crossed and throw exception if something went wrong
    if (newVertices.at(1)[i] < newVertices.at(0)[i]) {
      std::cout << "ERROR on process: " << this->localRank << std::endl;
      throw InternalErrorException(
          __FILE__, __func__, __LINE__,
          "Lower border of process larger than upper border of process!");
    }

    this->setVertices(newVertices);

#ifdef ALL_DEBUG_ENABLED
    MPI_Barrier(this->globalComm);
    if (this->localRank == 0)
      std::cout << "ALL::Staggered_LB::balance(): before neighbor search..."
                << std::endl;
#endif
    find_neighbors();
  }
}

template <class T, class W> void Staggered_LB<T, W>::find_neighbors() {
  auto work = this->getWork();
  auto vertices = this->prevVertices;
  auto shifted_vertices = this->vertices;

  neighbors.clear();
  // collect work from right neighbor plane
  MPI_Request sreq, rreq;
  MPI_Status sstat, rstat;
  // array to store neighbor vertices in Y/Z direction (reused)
  T *vertices_loc = new T[4];
  T *vertices_rem = new T[8 * this->global_dims[0] * this->global_dims[1]];

  int rem_rank;
  int rem_coords[3];

  // determine neighbors
  int rank_left, rank_right;

  // offset to get the correct rank
  int rank_offset;
  int offset_coords[3];

  // X-neighbors are static
  nNeighbors.at(0) = nNeighbors.at(1) = 1;

  // find X-neighbors
  MPI_Cart_shift(this->globalComm, 0, 1, &rank_left, &rank_right);

  // store X-neighbors
  neighbors.push_back(rank_left);
  neighbors.push_back(rank_right);

  // find Y-neighbors to get border information from
  MPI_Cart_shift(this->globalComm, 1, 1, &rank_left, &rank_right);

  // collect border information from local column
  vertices_loc[0] = shifted_vertices.at(0)[0];
  vertices_loc[1] = shifted_vertices.at(1)[0];
  MPI_Allgather(vertices_loc, 2, MPIDataTypeT,
                vertices_rem + 2 * this->global_dims[0], 2, MPIDataTypeT,
                communicators[1]);

  // exchange local column information with upper neighbor in Y direction (cart
  // grid)
  MPI_Irecv(vertices_rem, 2 * this->global_dims[0], MPIDataTypeT, rank_left, 0,
            this->globalComm, &rreq);
  MPI_Isend(vertices_rem + 2 * this->global_dims[0], 2 * this->global_dims[0],
            MPIDataTypeT, rank_right, 0, this->globalComm, &sreq);

  // determine the offset in ranks
  offset_coords[0] = 0;
  offset_coords[1] = this->local_coords[1] - 1;
  offset_coords[2] = this->local_coords[2];

  rem_coords[1] = offset_coords[1];
  rem_coords[2] = offset_coords[2];

  MPI_Cart_rank(this->globalComm, offset_coords, &rank_offset);

  // wait for communication
  MPI_Wait(&sreq, &sstat);
  MPI_Wait(&rreq, &rstat);

  // iterate about neighbor borders to determine the neighborship relation
  nNeighbors.at(2) = 0;
  for (int x = 0; x < this->global_dims[0]; ++x) {
    if ((vertices_rem[2 * x] <= vertices_loc[0] &&
         vertices_loc[0] < vertices_rem[2 * x + 1]) ||
        (vertices_rem[2 * x] < vertices_loc[1] &&
         vertices_loc[1] <= vertices_rem[2 * x + 1]) ||
        (vertices_rem[2 * x] >= vertices_loc[0] &&
         vertices_loc[0] < vertices_rem[2 * x + 1] &&
         vertices_loc[1] >= vertices_rem[2 * x + 1])) {
      nNeighbors.at(2)++;
      rem_coords[0] = x;
      MPI_Cart_rank(this->globalComm, rem_coords, &rem_rank);
      neighbors.push_back(rem_rank);
    }
  }

  // barrier to ensure every process concluded the calculations before
  // overwriting remote borders!
  MPI_Barrier(this->globalComm);

  // exchange local column information with lower neighbor in Y direction (cart
  // grid)
  MPI_Irecv(vertices_rem, 2 * this->global_dims[0], MPIDataTypeT, rank_right, 0,
            this->globalComm, &rreq);
  MPI_Isend(vertices_rem + 2 * this->global_dims[0], 2 * this->global_dims[0],
            MPIDataTypeT, rank_left, 0, this->globalComm, &sreq);

  // determine the offset in ranks
  offset_coords[0] = 0;
  offset_coords[1] = this->local_coords[1] + 1;
  offset_coords[2] = this->local_coords[2];

  rem_coords[1] = offset_coords[1];
  rem_coords[2] = offset_coords[2];

  MPI_Cart_rank(this->globalComm, offset_coords, &rank_offset);

  // wait for communication
  MPI_Wait(&sreq, &sstat);
  MPI_Wait(&rreq, &rstat);

  // iterate about neighbor borders to determine the neighborship relation
  nNeighbors.at(3) = 0;
  for (int x = 0; x < this->global_dims[0]; ++x) {
    if ((vertices_rem[2 * x] <= vertices_loc[0] &&
         vertices_loc[0] < vertices_rem[2 * x + 1]) ||
        (vertices_rem[2 * x] < vertices_loc[1] &&
         vertices_loc[1] <= vertices_rem[2 * x + 1]) ||
        (vertices_rem[2 * x] >= vertices_loc[0] &&
         vertices_loc[0] < vertices_rem[2 * x + 1] &&
         vertices_loc[1] >= vertices_rem[2 * x + 1])) {
      nNeighbors.at(3)++;
      rem_coords[0] = x;
      MPI_Cart_rank(this->globalComm, rem_coords, &rem_rank);
      neighbors.push_back(rem_rank);
    }
  }

  // barrier to ensure every process concluded the calculations before
  // overwriting remote borders!
  MPI_Barrier(this->globalComm);

  // find Z-neighbors to get border information from
  MPI_Cart_shift(this->globalComm, 2, 1, &rank_left, &rank_right);

  // collect border information from local column
  vertices_loc[0] = shifted_vertices.at(0)[0];
  vertices_loc[1] = shifted_vertices.at(1)[0];
  vertices_loc[2] = shifted_vertices.at(0)[1];
  vertices_loc[3] = shifted_vertices.at(1)[1];

  MPI_Barrier(this->globalComm);

  MPI_Allgather(vertices_loc, 4, MPIDataTypeT,
                vertices_rem + 4 * this->global_dims[0] * this->global_dims[1],
                4, MPIDataTypeT, communicators[2]);

  // exchange local column information with upper neighbor in Z direction (cart
  // grid)
  MPI_Irecv(vertices_rem, 4 * this->global_dims[0] * this->global_dims[1],
            MPIDataTypeT, rank_left, 0, this->globalComm, &rreq);
  MPI_Isend(vertices_rem + 4 * this->global_dims[0] * this->global_dims[1],
            4 * this->global_dims[0] * this->global_dims[1], MPIDataTypeT,
            rank_right, 0, this->globalComm, &sreq);

  // determine the offset in ranks
  offset_coords[0] = 0;
  offset_coords[1] = 0;
  offset_coords[2] = this->local_coords[2] - 1;

  rem_coords[2] = offset_coords[2];

  MPI_Cart_rank(this->globalComm, offset_coords, &rank_offset);

  // wait for communication
  MPI_Wait(&sreq, &sstat);
  MPI_Wait(&rreq, &rstat);

  // iterate about neighbor borders to determine the neighborship relation
  nNeighbors.at(4) = 0;
  for (int y = 0; y < this->global_dims[1]; ++y) {
    for (int x = 0; x < this->global_dims[0]; ++x) {
      if (((vertices_rem[4 * (x + y * this->global_dims[0]) + 2] <=
                vertices_loc[2] &&
            vertices_loc[2] <
                vertices_rem[4 * (x + y * this->global_dims[0]) + 3]) ||
           (vertices_rem[4 * (x + y * this->global_dims[0]) + 2] <
                vertices_loc[3] &&
            vertices_loc[3] <=
                vertices_rem[4 * (x + y * this->global_dims[0]) + 3]) ||
           (vertices_rem[4 * (x + y * this->global_dims[0]) + 2] >=
                vertices_loc[2] &&
            vertices_loc[2] <
                vertices_rem[4 * (x + y * this->global_dims[0]) + 3] &&
            vertices_loc[3] >=
                vertices_rem[4 * (x + y * this->global_dims[0]) + 3])))
        if (((vertices_rem[4 * (x + y * this->global_dims[0])] <=
                  vertices_loc[0] &&
              vertices_loc[0] <
                  vertices_rem[4 * (x + y * this->global_dims[0]) + 1]) ||
             (vertices_rem[4 * (x + y * this->global_dims[0])] <
                  vertices_loc[1] &&
              vertices_loc[1] <=
                  vertices_rem[4 * (x + y * this->global_dims[0]) + 1]) ||
             (vertices_rem[4 * (x + y * this->global_dims[0])] >=
                  vertices_loc[0] &&
              vertices_loc[0] <
                  vertices_rem[4 * (x + y * this->global_dims[0]) + 1] &&
              vertices_loc[1] >=
                  vertices_rem[4 * (x + y * this->global_dims[0]) + 1]))) {
          nNeighbors.at(4)++;
          rem_coords[1] = y;
          rem_coords[0] = x;
          MPI_Cart_rank(this->globalComm, rem_coords, &rem_rank);
          neighbors.push_back(rem_rank);
        }
    }
  }

  // barrier to ensure every process concluded the calculations before
  // overwriting remote borders!
  MPI_Barrier(this->globalComm);

  // exchange local column information with upper neighbor in Y direction (cart
  // grid)
  MPI_Irecv(vertices_rem, 4 * this->global_dims[0] * this->global_dims[1],
            MPIDataTypeT, rank_right, 0, this->globalComm, &rreq);
  MPI_Isend(vertices_rem + 4 * this->global_dims[0] * this->global_dims[1],
            4 * this->global_dims[0] * this->global_dims[1], MPIDataTypeT,
            rank_left, 0, this->globalComm, &sreq);

  // determine the offset in ranks
  offset_coords[0] = 0;
  offset_coords[1] = 0;
  offset_coords[2] = this->local_coords[2] + 1;

  rem_coords[2] = offset_coords[2];

  MPI_Cart_rank(this->globalComm, offset_coords, &rank_offset);

  // wait for communication
  MPI_Wait(&sreq, &sstat);
  MPI_Wait(&rreq, &rstat);

  // iterate about neighbor borders to determine the neighborship relation
  nNeighbors.at(5) = 0;
  for (int y = 0; y < this->global_dims[1]; ++y) {
    for (int x = 0; x < this->global_dims[0]; ++x) {
      if (((vertices_rem[4 * (x + y * this->global_dims[0]) + 2] <=
                vertices_loc[2] &&
            vertices_loc[2] <
                vertices_rem[4 * (x + y * this->global_dims[0]) + 3]) ||
           (vertices_rem[4 * (x + y * this->global_dims[0]) + 2] <
                vertices_loc[3] &&
            vertices_loc[3] <=
                vertices_rem[4 * (x + y * this->global_dims[0]) + 3]) ||
           (vertices_rem[4 * (x + y * this->global_dims[0]) + 2] >=
                vertices_loc[2] &&
            vertices_loc[2] <
                vertices_rem[4 * (x + y * this->global_dims[0]) + 3] &&
            vertices_loc[3] >=
                vertices_rem[4 * (x + y * this->global_dims[0]) + 3])))
        if (((vertices_rem[4 * (x + y * this->global_dims[0])] <=
                  vertices_loc[0] &&
              vertices_loc[0] <
                  vertices_rem[4 * (x + y * this->global_dims[0]) + 1]) ||
             (vertices_rem[4 * (x + y * this->global_dims[0])] <
                  vertices_loc[1] &&
              vertices_loc[1] <=
                  vertices_rem[4 * (x + y * this->global_dims[0]) + 1]) ||
             (vertices_rem[4 * (x + y * this->global_dims[0])] >=
                  vertices_loc[0] &&
              vertices_loc[0] <
                  vertices_rem[4 * (x + y * this->global_dims[0]) + 1] &&
              vertices_loc[1] >=
                  vertices_rem[4 * (x + y * this->global_dims[0]) + 1]))) {
          nNeighbors.at(5)++;
          rem_coords[1] = y;
          rem_coords[0] = x;
          MPI_Cart_rank(this->globalComm, rem_coords, &rem_rank);
          neighbors.push_back(rem_rank);
        }
    }
  }

  // barrier to ensure every process concluded the calculations before
  // overwriting remote borders!
  MPI_Barrier(this->globalComm);

  // clear up vertices array
  delete[] vertices_loc;
  delete[] vertices_rem;
}

// provide list of neighbors
template <class T, class W>
std::vector<int> &Staggered_LB<T, W>::getNeighbors() {
  return neighbors;
}

}//namespace ALL

#endif
// vim: sw=2 et ts=2
