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

#ifndef ALL_FORCEBASED_HEADER_INCLUDED
#define ALL_FORCEBASED_HEADER_INCLUDED

#define ALL_FORCE_OLD

/*

    ForceBased load-balancing scheme:

    Requirements:
        a) cartesian communicator
        b) 2^(dim) vertices describing the domain (in theory could be more,
           we assume a orthogonal grid as start point)
        c) MPI >= 3.0 due to necessary non-blocking collectives


    Method:

        For each of the vertices of a domain a communicator is created, which
   contains the surrounding processes. For one of the vertices, a domain
   collects the centres of gravitiy and the work of the neighboring domains. In
   a next step a force acting on the vertex is computed as: F = 1.0 / ( gamma *
   n * sum(W_i) ) * sum( W_i * x_i ) with:

                n:      number of neighboring domains
                W_i:    work on domain i
                x_i:    vector pointing from the vertex to
                        the center of gravity of domain i
                gamma:  correction factor to control the
                        speed of the vertex shift

        */

#include "ALL_CustomExceptions.hpp"
#include "ALL_LB.hpp"
#include "ALL_Point.hpp"
#include "ALL_Defines.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <mpi.h>
#include <vector>

namespace ALL {

/// Load-balancing scheme based on the shift of vertices according to the ratio
/// of work from all the neighboring domains. The vertex is drawn in the
/// direction of the geometric center of the domain with the most work. The
/// direction is influenced by the other domains, their work and location of
/// geometric centers. Some precautions are taken to avoid domains with too
/// small angles at their vertices. The neighboring relation of domains is
/// conserved throughout the load-balancing process.
/// @tparam T data for vertices and related data
/// @tparam W data for work and related data
template <class T, class W> class ForceBased_LB : public LB<T, W> {
public:
  /// default constructor
  ForceBased_LB() {}
  /// constructor initializing basic parameters
  /// @param[in] d the dimension of the used vertices
  /// @param[in] w the scalar work assigned to the local domain
  /// @param[in] g the correction factor gamma
  ForceBased_LB(int d, W w, T g) : LB<T, W>(d, g) {
    this->setWork(w);
  }

  /// default destructor
  ~ForceBased_LB() override;

  /// setup internal data structures and parameters
  void setup() override;

  /// method to execute a load-balancing step
  /// @param[in] step number of the load-balancing step
  virtual void balance(int step) override;

  // getter for variables (passed by reference to avoid
  // direct access to private members of the object)

  /// method to provide a list of the neighbors of the local domain
  /// @param[out] list reference to a std::vector of integers where the list of
  /// neighbors will be assigned to
  virtual std::vector<int> &getNeighbors() override;

  /// method to set specific data structures (unused for tensor grid method)
  /// @param[in] data pointer to the data structure
  virtual void setAdditionalData(known_unused const void *data) override {}

  /// method to get an estimated work distribution after the balance step
  /// (currently only implemented in ALL::HISTOGRAM!)
  /// @param [out] double providing the estimated LB after the balance step
  virtual W getEstimatedEfficiency() override {return (W)-1;}

private:
  // type for MPI communication
  MPI_Datatype MPIDataTypeT;
  MPI_Datatype MPIDataTypeW;

  // array of MPI communicators for each vertex
  std::vector<T> vertex_neighbors;

  // MPI communicator to collect the work in the
  // main dimension
  MPI_Comm main_communicator;

  // list of neighbors
  std::vector<int> neighbors;

  // find neighbors (more sophisticated than for tensor-product
  // variant of LB)
  void find_neighbors();

  // number of vertices (since it is not equal for all domains)
  int n_vertices;

  // main dimension (correction in staggered style)
  int main_dim;

  // secondary dimensions (2D force shift)
  int sec_dim[2];
};

template <class T, class W> ForceBased_LB<T, W>::~ForceBased_LB() {}

// provide list of neighbors
template <class T, class W>
std::vector<int> &ForceBased_LB<T, W>::getNeighbors() {
  return neighbors;
}

template <class T, class W> void ForceBased_LB<T, W>::setup() {
  n_vertices = this->vertices.size();
  vertex_neighbors.resize(n_vertices * 8);

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

  MPI_Group all;
  MPI_Comm_group(this->globalComm, &all);

  // check if communicator is cartesian
  int status;
  MPI_Topo_test(this->globalComm, &status);

  if (status != MPI_CART) {
    throw InvalidCommTypeException(
        __FILE__, __func__, __LINE__,
        "Cartesian MPI communicator required, passed communicator is not "
        "cartesian");
  }

  // get the local coordinates, this->periodicity and global size from the MPI
  // communicator
  MPI_Cart_get(this->globalComm, this->dimension, this->global_dims.data(),
               this->periodicity.data(), this->local_coords.data());

  // get the local rank from the MPI communicator
  MPI_Cart_rank(this->globalComm, this->local_coords.data(), &this->localRank);

  // groups required for new communicators
  MPI_Group known_unused  groups[n_vertices];
  // arrays of processes belonging to group
  int known_unused processes[n_vertices][n_vertices];

  // shifted local coordinates to find neighboring processes
  int known_unused shifted_coords[this->dimension];

  // get main and secondary dimensions
  if (this->global_dims.at(2) >= this->global_dims.at(1)) {
    if (this->global_dims.at(2) >= this->global_dims.at(0)) {
      main_dim = 2;
      sec_dim[0] = 0;
      sec_dim[1] = 1;
    } else {
      main_dim = 0;
      sec_dim[0] = 1;
      sec_dim[1] = 2;
    }
  } else {
    if (this->global_dims.at(1) >= this->global_dims.at(0)) {
      main_dim = 1;
      sec_dim[0] = 0;
      sec_dim[1] = 2;
    } else {
      main_dim = 0;
      sec_dim[0] = 1;
      sec_dim[1] = 2;
    }
  }

  if (this->localRank == 0)
    std::cout << "DEBUG: main_dim: " << main_dim << std::endl;

  // create main communicator
  MPI_Comm_split(this->globalComm, this->local_coords.at(main_dim),
                 this->local_coords.at(sec_dim[0]) +
                     this->local_coords.at(sec_dim[1]) *
                         this->global_dims.at(sec_dim[0]),
                 &main_communicator);

  /*

      // sequence of vertices in 1d
      //
      // 0 - - - 1

      // sequence of vertices in 2d
      //
      // 2 - - - 3
      // |       |
      // |       |
      // 0 - - - 1

      // sequence of vertices in 3d
      //
      //    6 - - - 7
      //   /|      /|
      //  / |     / |
      // 4 -2- - 5  3
      // | /     | /
      // |/      |/
      // 0 - - - 1

  */

#ifdef ALL_DEBUG_ENABLED
  MPI_Barrier(this->globalComm);
  if (this->localRank == 0)
    std::cout << "ALL::ForceBased_LB<T,W>::setup() preparing communicators..."
              << std::endl;
#endif
  std::vector<int> dim_vert(this->global_dims);

  MPI_Comm known_unused tmp_comm;
  int known_unused own_vertex;

#ifdef ALL_DEBUG_ENABLED
  MPI_Barrier(this->globalComm);
  if (this->localRank == 0)
    std::cout << "ALL::ForceBased_LB<T,W>::setup() computing communicators..."
              << std::endl;
  std::cout << "DEBUG: "
            << " rank: " << this->localRank << " dim_vert: " << dim_vert.at(0)
            << " " << dim_vert.at(1) << " " << dim_vert.at(2)
            << " size(local_coords): " << this->local_coords.size() << " "
            << " size(global_dims):  " << this->global_dims.size() << std::endl;
#endif
  for (int iz = 0; iz < dim_vert.at(2); ++iz) {
    for (int iy = 0; iy < dim_vert.at(1); ++iy) {
      for (int ix = 0; ix < dim_vert.at(0); ++ix) {
        bool affected[8];
        for (auto &a : affected)
          a = false;
        int v_neighbors[8];
        for (auto &vn : vertex_neighbors)
          vn = -1;
        if (ix == ((this->local_coords.at(0) + 0) % dim_vert.at(0)) &&
            iy == ((this->local_coords.at(1) + 0) % dim_vert.at(1)) &&
            iz == ((this->local_coords.at(2) + 0) % dim_vert.at(2))) {
          affected[0] = true;
          v_neighbors[0] = this->localRank;
        }
        if (ix == ((this->local_coords.at(0) + 1) % dim_vert.at(0)) &&
            iy == ((this->local_coords.at(1) + 0) % dim_vert.at(1)) &&
            iz == ((this->local_coords.at(2) + 0) % dim_vert.at(2))) {
          affected[1] = true;
          v_neighbors[1] = this->localRank;
        }
        if (ix == ((this->local_coords.at(0) + 0) % dim_vert.at(0)) &&
            iy == ((this->local_coords.at(1) + 1) % dim_vert.at(1)) &&
            iz == ((this->local_coords.at(2) + 0) % dim_vert.at(2))) {
          affected[2] = true;
          v_neighbors[2] = this->localRank;
        }
        if (ix == ((this->local_coords.at(0) + 1) % dim_vert.at(0)) &&
            iy == ((this->local_coords.at(1) + 1) % dim_vert.at(1)) &&
            iz == ((this->local_coords.at(2) + 0) % dim_vert.at(2))) {
          affected[3] = true;
          v_neighbors[3] = this->localRank;
        }
        if (ix == ((this->local_coords.at(0) + 0) % dim_vert.at(0)) &&
            iy == ((this->local_coords.at(1) + 0) % dim_vert.at(1)) &&
            iz == ((this->local_coords.at(2) + 1) % dim_vert.at(2))) {
          affected[4] = true;
          v_neighbors[4] = this->localRank;
        }
        if (ix == ((this->local_coords.at(0) + 1) % dim_vert.at(0)) &&
            iy == ((this->local_coords.at(1) + 0) % dim_vert.at(1)) &&
            iz == ((this->local_coords.at(2) + 1) % dim_vert.at(2))) {
          affected[5] = true;
          v_neighbors[5] = this->localRank;
        }
        if (ix == ((this->local_coords.at(0) + 0) % dim_vert.at(0)) &&
            iy == ((this->local_coords.at(1) + 1) % dim_vert.at(1)) &&
            iz == ((this->local_coords.at(2) + 1) % dim_vert.at(2))) {
          affected[6] = true;
          v_neighbors[6] = this->localRank;
        }
        if (ix == ((this->local_coords.at(0) + 1) % dim_vert.at(0)) &&
            iy == ((this->local_coords.at(1) + 1) % dim_vert.at(1)) &&
            iz == ((this->local_coords.at(2) + 1) % dim_vert.at(2))) {
          affected[7] = true;
          v_neighbors[7] = this->localRank;
        }

        MPI_Allreduce(MPI_IN_PLACE, v_neighbors, 8, MPI_INT, MPI_MAX,
                      this->globalComm);

        for (int v = 0; v < n_vertices; ++v) {
          if (affected[v]) {
            for (int n = 0; n < 8; ++n) {
              vertex_neighbors.at(8 * v + n) = v_neighbors[n];
            }
          }
        }
      }
    }
// ToDo: check if vertices correct
#ifdef ALL_DEBUG_ENABLED
    MPI_Barrier(this->globalComm);
    if (this->localRank == 0)
      std::cout << "ALL::ForceBased_LB<T,W>::setup() finished computing "
                   "communicators..."
                << std::endl;
#endif
  }
}

// TODO: periodic boundary conditions (would require size of the system)
//       for now: do not change the vertices along the edge of the system

template <class T, class W>
void ForceBased_LB<T, W>::balance(int /*step*/) {

  this->prevVertices = this->vertices;

  // geometrical center and work of each neighbor for each vertex
  // work is cast to vertex data type, since in the end a shift
  // of the vertex is computed
  T vertex_info[n_vertices * n_vertices * (this->dimension + 2)];
  T center[this->dimension];

  // local geometric center
  T known_unused local_info[(this->dimension + 2) * n_vertices];

  for (int i = 0; i < n_vertices * n_vertices * (this->dimension + 2); ++i)
    vertex_info[i] = -1;

  for (int d = 0; d < (this->dimension + 2) * n_vertices; ++d)
    local_info[d] = (T)0;

  for (int d = 0; d < this->dimension; ++d)
    center[d] = (T)0;

  // compute geometric center
  for (int v = 0; v < n_vertices; ++v) {
    for (int d = 0; d < this->dimension; ++d) {
      center[d] += ((this->prevVertices.at(v))[d] / (T)n_vertices);
    }
    local_info[v * (this->dimension + 2) + this->dimension] =
        (T)this->work.at(0);
  }
  for (int v = 0; v < n_vertices; ++v) {
    for (int d = 0; d < this->dimension; ++d) {
      // vector pointing to center
      local_info[v * (this->dimension + 2) + d] =
          center[d] - (this->prevVertices.at(v))[d];
    }
  }

  // compute for each vertex the maximum movement
  switch (main_dim) {
  case 0:
    local_info[0 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(0).d(this->prevVertices.at(2)),
                       this->prevVertices.at(0).d(this->prevVertices.at(4)));
    local_info[1 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(1).d(this->prevVertices.at(3)),
                       this->prevVertices.at(1).d(this->prevVertices.at(5)));
    local_info[2 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(2).d(this->prevVertices.at(0)),
                       this->prevVertices.at(2).d(this->prevVertices.at(6)));
    local_info[3 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(3).d(this->prevVertices.at(1)),
                       this->prevVertices.at(3).d(this->prevVertices.at(7)));
    local_info[4 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(4).d(this->prevVertices.at(0)),
                       this->prevVertices.at(4).d(this->prevVertices.at(6)));
    local_info[5 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(5).d(this->prevVertices.at(1)),
                       this->prevVertices.at(5).d(this->prevVertices.at(7)));
    local_info[6 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(6).d(this->prevVertices.at(2)),
                       this->prevVertices.at(6).d(this->prevVertices.at(4)));
    local_info[7 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(7).d(this->prevVertices.at(3)),
                       this->prevVertices.at(7).d(this->prevVertices.at(5)));
    break;
  case 1:
    local_info[0 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(0).d(this->prevVertices.at(1)),
                       this->prevVertices.at(0).d(this->prevVertices.at(4)));
    local_info[1 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(1).d(this->prevVertices.at(0)),
                       this->prevVertices.at(1).d(this->prevVertices.at(5)));
    local_info[2 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(2).d(this->prevVertices.at(3)),
                       this->prevVertices.at(2).d(this->prevVertices.at(6)));
    local_info[3 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(3).d(this->prevVertices.at(2)),
                       this->prevVertices.at(3).d(this->prevVertices.at(7)));
    local_info[4 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(4).d(this->prevVertices.at(0)),
                       this->prevVertices.at(4).d(this->prevVertices.at(5)));
    local_info[5 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(5).d(this->prevVertices.at(1)),
                       this->prevVertices.at(5).d(this->prevVertices.at(4)));
    local_info[6 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(6).d(this->prevVertices.at(2)),
                       this->prevVertices.at(6).d(this->prevVertices.at(7)));
    local_info[7 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(7).d(this->prevVertices.at(3)),
                       this->prevVertices.at(7).d(this->prevVertices.at(6)));
    break;
  case 2:
    local_info[0 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(0).d(this->prevVertices.at(1)),
                       this->prevVertices.at(0).d(this->prevVertices.at(2)));
    local_info[1 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(1).d(this->prevVertices.at(0)),
                       this->prevVertices.at(1).d(this->prevVertices.at(3)));
    local_info[2 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(2).d(this->prevVertices.at(0)),
                       this->prevVertices.at(2).d(this->prevVertices.at(3)));
    local_info[3 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(3).d(this->prevVertices.at(1)),
                       this->prevVertices.at(3).d(this->prevVertices.at(2)));
    local_info[4 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(4).d(this->prevVertices.at(5)),
                       this->prevVertices.at(4).d(this->prevVertices.at(6)));
    local_info[5 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(5).d(this->prevVertices.at(4)),
                       this->prevVertices.at(5).d(this->prevVertices.at(7)));
    local_info[6 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(6).d(this->prevVertices.at(4)),
                       this->prevVertices.at(6).d(this->prevVertices.at(7)));
    local_info[7 * (this->dimension + 2) + this->dimension + 1] =
        0.5 * std::min(this->prevVertices.at(7).d(this->prevVertices.at(5)),
                       this->prevVertices.at(7).d(this->prevVertices.at(6)));
    break;
  default:
    throw InternalErrorException(
        __FILE__, __func__, __LINE__,
        "Invalid main dimension provided (numerical value not 0, 1 or 2).");
    break;
  }
  // exchange information with all vertex neighbors
  MPI_Request known_unused request[n_vertices];
  MPI_Status known_unused status[n_vertices];

  // compute new position for vertex 7 (if not periodic)
  T known_unused total_work = (T)0;
  T known_unused shift_vectors[this->dimension * n_vertices];

  for (int v = 0; v < n_vertices; ++v) {
    for (int d = 0; d < this->dimension; ++d) {
      shift_vectors[v * this->dimension + d] = (T)0;
    }
  }

  // TODO:
  // 1.) collect work in main_dim communicators
  // 2.) exchange work and extension in main_dim to neighbors in main_dim
  // 3.) correct extension in main_dim
  // 4.) do force-based correction in secondary dimensions
  // 5.) find new neighbors in main_dim (hardest part, probably cell based)

  int main_up;
  int main_down;

  // as each process computes the correction of the upper layer of vertices
  // the source is 'main_up' and the target 'main_down'
  MPI_Cart_shift(this->globalComm, main_dim, 1, &main_down, &main_up);

  W local_work = this->work.at(0);
  W remote_work;

  // reduce work from local layer
  MPI_Allreduce(MPI_IN_PLACE, &local_work, 1, MPIDataTypeW, MPI_SUM,
                main_communicator);

  // exchange local work with neighbor in main direction
  MPI_Status state;
  MPI_Sendrecv(&local_work, 1, MPIDataTypeW, main_down, 1010, &remote_work, 1,
               MPIDataTypeW, main_up, 1010, this->globalComm, &state);

  // transfer width of domains as well
  T remote_width;
  T local_width;
  switch (main_dim) {
  case 0:
    local_width = this->prevVertices.at(1)[0] - this->prevVertices.at(0)[0];
    break;
  case 1:
    local_width = this->prevVertices.at(2)[1] - this->prevVertices.at(0)[1];
    break;
  case 2:
    local_width = this->prevVertices.at(4)[2] - this->prevVertices.at(0)[2];
    break;
  default:
    throw InternalErrorException(
        __FILE__, __func__, __LINE__,
        "Invalid main dimension provided (numerical value not 0, 1 or 2).");
    break;
  }
  MPI_Sendrecv(&local_width, 1, MPIDataTypeT, main_down, 1010, &remote_width, 1,
               MPIDataTypeT, main_up, 1010, this->globalComm, &state);
  T total_width = local_width + remote_width;
  T max_main_shift = std::min(0.49 * std::min(local_width, remote_width),
                              std::min(local_width, remote_width) - 1.0);

  // compute shift
  T local_shift = (remote_work - local_work) / (remote_work + local_work) /
                  this->gamma / 2.0 * total_width;
  local_shift = (std::abs(local_shift) > max_main_shift)
                    ? std::abs(local_shift) * max_main_shift / local_shift
                    : local_shift;
  T remote_shift = 0;

  // TODO: fix -> send to upper neighbor, receive from lower neighbor!!
  //              sendrecv not correct!

  // transfer shift back
  MPI_Sendrecv(&local_shift, 1, MPIDataTypeT, main_up, 2020, &remote_shift, 1,
               MPIDataTypeT, main_down, 2020, this->globalComm, &state);

  // apply shift in main direction
  switch (main_dim) {
  case 0:
    if (this->local_coords.at(0) > 0) {
      this->vertices.at(0)[0] = this->prevVertices.at(0)[0] + remote_shift;
      this->vertices.at(2)[0] = this->prevVertices.at(2)[0] + remote_shift;
      this->vertices.at(4)[0] = this->prevVertices.at(4)[0] + remote_shift;
      this->vertices.at(6)[0] = this->prevVertices.at(6)[0] + remote_shift;
    }
    if (this->local_coords.at(0) < this->global_dims.at(0) - 1) {
      this->vertices.at(1)[0] = this->prevVertices.at(1)[0] + local_shift;
      this->vertices.at(3)[0] = this->prevVertices.at(3)[0] + local_shift;
      this->vertices.at(5)[0] = this->prevVertices.at(5)[0] + local_shift;
      this->vertices.at(7)[0] = this->prevVertices.at(7)[0] + local_shift;
    }
    break;
  case 1:
    if (this->local_coords.at(1) > 0) {
      this->vertices.at(0)[1] = this->prevVertices.at(0)[1] + remote_shift;
      this->vertices.at(1)[1] = this->prevVertices.at(1)[1] + remote_shift;
      this->vertices.at(4)[1] = this->prevVertices.at(4)[1] + remote_shift;
      this->vertices.at(5)[1] = this->prevVertices.at(5)[1] + remote_shift;
    }
    if (this->local_coords.at(1) < this->global_dims.at(1) - 1) {
      this->vertices.at(2)[1] = this->prevVertices.at(2)[1] + local_shift;
      this->vertices.at(3)[1] = this->prevVertices.at(3)[1] + local_shift;
      this->vertices.at(6)[1] = this->prevVertices.at(6)[1] + local_shift;
      this->vertices.at(7)[1] = this->prevVertices.at(7)[1] + local_shift;
    }
    break;
  case 2:
    if (this->local_coords.at(2) > 0) {
      this->vertices.at(0)[2] = this->prevVertices.at(0)[2] + remote_shift;
      this->vertices.at(1)[2] = this->prevVertices.at(1)[2] + remote_shift;
      this->vertices.at(2)[2] = this->prevVertices.at(2)[2] + remote_shift;
      this->vertices.at(3)[2] = this->prevVertices.at(3)[2] + remote_shift;
    }
    if (this->local_coords.at(2) < this->global_dims.at(2) - 1) {
      this->vertices.at(4)[2] = this->prevVertices.at(4)[2] + local_shift;
      this->vertices.at(5)[2] = this->prevVertices.at(5)[2] + local_shift;
      this->vertices.at(6)[2] = this->prevVertices.at(6)[2] + local_shift;
      this->vertices.at(7)[2] = this->prevVertices.at(7)[2] + local_shift;
    }
    break;
  default:
    throw InternalErrorException(
        __FILE__, __func__, __LINE__,
        "Invalid main dimension provided (numerical value not 0, 1 or 2).");
    break;
  }

  if (this->localRank <= 1) {
    std::cout << "DEBUG (main-dim shift): " << this->localRank << " "
              << remote_work << " " << local_work << " " << total_width << " "
              << local_shift << " " << remote_shift
              << " 0: " << this->prevVertices.at(0)[2]
              << " 0: " << this->vertices.at(0)[2]
              << " 4: " << this->prevVertices.at(4)[2]
              << " 4: " << this->vertices.at(4)[2] << " " << std::endl;
  }

  int last_vertex = (n_vertices - 1) * n_vertices * (this->dimension + 2);
  int vertex_offset = this->dimension + 2;

  // compute max shift
  T max_shift = std::max(
      0.49 * vertex_info[last_vertex + 0 * vertex_offset + this->dimension + 1],
      1e-6);
  for (int v = 1; v < n_vertices; ++v) {
    max_shift = std::min(
        max_shift,
        0.49 *
            vertex_info[last_vertex + v * vertex_offset + this->dimension + 1]);
  }

  // average work of neighboring processes for last vertex
  T avg_work = 0.0;
  for (int v = 0; v < n_vertices; ++v) {
    avg_work += vertex_info[last_vertex + v * vertex_offset + this->dimension];
  }
  avg_work /= (T)n_vertices;

  // compute shift vector
  std::vector<T> vertex_shift(n_vertices * this->dimension, (T)0.0);
  for (auto &sv : vertex_shift)
    sv = (T)0.0;
  int shift_offset = (n_vertices - 1) * this->dimension;

  for (int v = 0; v < n_vertices; ++v) {
    for (int d = 0; d < this->dimension; ++d) {
      if (this->localRank == 0)
        std::cout
            << "DEBUG (shift x): " << v << ", " << d << " "
            << (vertex_info[last_vertex + v * vertex_offset + this->dimension] -
                avg_work) *
                   ((vertex_info[last_vertex + v * vertex_offset +
                                 this->dimension] > avg_work)
                        ? 1.0
                        : -1.0) /
                   this->gamma *
                   (vertex_info[last_vertex + v * vertex_offset + d] -
                    this->prevVertices.at(v)[d])
            << " => " << vertex_shift.at(shift_offset + d) <<

            " max: " << max_shift
            << " "
               " avg_work: "
            << avg_work
            << " "
               " work: "
            << vertex_info[last_vertex + v * vertex_offset + this->dimension]
            << " "
               " 1/0: "
            << ((vertex_info[last_vertex + v * vertex_offset +
                             this->dimension] > avg_work)
                    ? 1.0
                    : -1.0)

            << std::endl;
      vertex_shift.at(shift_offset + d) +=
          (vertex_info[last_vertex + v * vertex_offset + this->dimension] -
           avg_work) *
          ((vertex_info[last_vertex + v * vertex_offset + this->dimension] >
            avg_work)
               ? 1.0
               : -1.0) /
          this->gamma *
          (vertex_info[last_vertex + v * vertex_offset + d] -
           this->prevVertices.at(v)[d]);
    }
  }

  Point<T> shift_vector(3);
  shift_vector[0] = vertex_shift.at(shift_offset);
  shift_vector[1] = vertex_shift.at(shift_offset + 1);
  shift_vector[2] = vertex_shift.at(shift_offset + 2);
  T shift_length = shift_vector.norm();
  if (this->localRank == 0)
    std::cout << "DEBUG (shift length): " << shift_length
              << " shift vector: " << shift_vector[0] << " " << shift_vector[1]
              << " " << shift_vector[2] << " " << std::endl;

  // apply correct length, if the shift vector is too large
  if (shift_length > max_shift) {
    for (int d = 0; d < this->dimension; ++d) {
      vertex_shift.at(shift_offset + d) *= max_shift / shift_length;
    }
  }

  if (this->localRank == 0) {
    for (int v = 0; v < n_vertices; ++v) {
      std::cout << "DEBUG shift vector vertex (before): " << v << ": ";
      for (int d = 0; d < this->dimension; ++d) {
        std::cout << vertex_shift.at(v * this->dimension + d) << " ";
      }
      std::cout << std::endl;
    }
  }

  if (this->localRank == 0 || this->localRank == 1) {
    for (int v = 0; v < n_vertices; ++v) {
      std::cout << "DEBUG shift vector vertex " << v << ": ";
      for (int d = 0; d < this->dimension; ++d) {
        std::cout << vertex_shift.at(v * this->dimension + d) << " ";
      }
      std::cout << std::endl;
    }
  }

  bool shift[3];
  for (int z = 0; z < 2; ++z) {
    if ((z == 0 && this->local_coords.at(2) == 0) ||
        (z == 1 && this->local_coords.at(2) == this->global_dims.at(2) - 1))
      shift[2] = false;
    else
      shift[2] = true;
    for (int y = 0; y < 2; ++y) {
      if ((y == 0 && this->local_coords.at(1) == 0) ||
          (y == 1 && this->local_coords.at(1) == this->global_dims.at(1) - 1))
        shift[1] = false;
      else
        shift[1] = true;
      for (int x = 0; x < 2; ++x) {
        if ((x == 0 && this->local_coords.at(0) == 0) ||
            (x == 1 && this->local_coords.at(0) == this->global_dims.at(0) - 1))
          shift[0] = false;
        else
          shift[0] = true;
        for (int d = 0; d < this->dimension; ++d) {
          if (d == main_dim)
            continue;
          int v = 4 * z + 2 * y + x;
          if (shift[d])
            this->vertices.at(v)[d] = this->prevVertices.at(v)[d] +
                                      vertex_shift.at(v * this->dimension + d);
          else
            this->vertices.at(v)[d] = this->prevVertices.at(v)[d];
        }
      }
    }
  }
}

}//namespace ALL

#endif
