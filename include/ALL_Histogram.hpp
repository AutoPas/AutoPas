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

#ifndef ALL_HISTOGRAM_HEADER_INCLUDED
#define ALL_HISTOGRAM_HEADER_INCLUDED

#include "ALL_CustomExceptions.hpp"
#include "ALL_LB.hpp"
#include <exception>
#include <mpi.h>
#include <vector>

namespace ALL {

/// Load-balancing scheme based on the usage of global histograms in order to
/// provide a most optimal decomposition of data. The decomposition is created
/// by calling the balancing method three times. For each of the different
/// dimensions the data is sorted into histograms (the correction factor gamma
/// determines the width of each data bin) and using this as a cumulative
/// distribution function a decomposition is computed where each chunk contains
/// roughly the same amount of data. For the next dimensions this operation is
/// repeated within the chunks provided by the previous call, therefore
/// resulting after three call in a full three-dimensional decomposition.
/// Between each call to the balancing method the data needs to be transferred
/// according to the newly created domain boundaries.
/// @tparam T data for vertices and related data
/// @tparam W data for work and related data
template <class T, class W> class Histogram_LB : public LB<T, W> {
public:
  /// default constructor
  Histogram_LB() {}
  /// constructor to initialize the main parameters
  /// @param d the dimension of the used vertices
  /// @param w the multi-dimensional work of the local domain, the number of
  /// entries has to be equal to the number of bins overlapping the local domain
  /// in the current direction (bin width = correction factor gamma)
  /// @param g the correction factor, i.e. the width of a single bin
  Histogram_LB(int d, std::vector<W> w, T g) : LB<T, W>(d, g) {
    this->setWork(w);
    // need the lower and upper bounds
    // of the domain for the tensor-based
    // load-balancing scheme

    // array of MPI communicators for each direction (to collect work
    // on each plane)
    communicators.resize(6);
    for (int i=0; i<6; i++) communicators.at(i) = MPI_COMM_NULL;

    nNeighbors.resize(2 * this->dimension);
    nBins.resize(this->dimension);
  }

  /// default destructor
  ~Histogram_LB() {
    for (int i=0; i<communicators.size(); i++) {
      if (communicators.at(i) != MPI_COMM_NULL &&
          communicators.at(i) != MPI_COMM_SELF &&
          communicators.at(i) != MPI_COMM_WORLD)
        MPI_Comm_free(&(communicators.at(i)));
    }
  }

  /// method to setup the loac-balancing method
  void setup() override;

  /// method to perform a single load-balancing step
  /// @param step the number the actual load-balancing step should get (e.g. if
  /// counting the number of loadbalancing steps)
  void balance(int step) override;

  // getter for variables (passed by reference to avoid
  // direct access to private members of the object)

  // neighbors
  /// method to provide a list of neighbors by MPI rank
  /// @param [out] list the std::vector the list of neighbors is stored to
  std::vector<int> &getNeighbors() override;

  /// method to set method specific data
  /// @param data pointer to an array of integers (int) containing the number of
  /// bins in the system for all three dimensions (int, int, int)
  virtual void setAdditionalData(const void *data) override;


  /// method to get an estimated work distribution after the balance step
  /// (currently only implemented in ALL::HISTOGRAM!)
  /// @param [out] double providing the estimated LB after the balance step
  virtual W getEstimatedEfficiency() override;

private:
  /// number of bins in each dimension
  std::vector<int> nBins;

  /// MPI data type to transfer data related to vertices
  MPI_Datatype MPIDataTypeT;
  /// MPI data type to transfer data related to work
  MPI_Datatype MPIDataTypeW;

  /// set of communicators used to collect data over the different dimensions,
  /// e.g. to compute the total amount of work in it and to compute the amount
  /// of work each chunk in that dimension should receive
  std::vector<MPI_Comm> communicators;

  /// list of MPI ranks of neighboring domains
  std::vector<int> neighbors;
  /// number of neighbors in each dimension
  std::vector<int> nNeighbors;

  /// method to determine the MPI ranks of the neighbors of the local process
  void find_neighbors();

  /// estimated efficiency after a load-balancing step
  W estimatedEfficiency;
};

template <class T, class W>
void Histogram_LB<T, W>::setAdditionalData(const void *data) {
  if (nBins.size() < 3)
    nBins.resize(3);

  for (int i = 0; i < 3; ++i)
    nBins.at(i) = *((int *)data + i);
}

// setup routine for the tensor-based load-balancing scheme
// requires:
//              this->globalComm (int): cartesian MPI communicator, from
//                                 which separate sub communicators
//                                 are derived in order to represent
//                                 each plane of domains in the system
template <class T, class W> void Histogram_LB<T, W>::setup() {
  int status;

  // check if Communicator is cartesian
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

  // free sub-communicators
  for (int i=0; i<6; i++) {
    if (communicators.at(i) != MPI_COMM_SELF &&
        communicators.at(i) != MPI_COMM_WORLD &&
        communicators.at(i) != MPI_COMM_NULL)
      MPI_Comm_free(&(communicators.at(i)));
  }

  // create sub-communicators

  // communicators 0 - 2: reduction of partial histograms
  // communicators 3 - 5: gathering of complete histograms

  // reduction sub-communicators (equal to staggered grid)

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
  communicators[0] = MPI_COMM_SELF;

  // gathering sub-communicators

  // x-gathering (same y/z coordinates)
  MPI_Comm_split(this->globalComm,
                 this->local_coords.at(1) +
                     this->local_coords.at(2) * this->global_dims.at(1),
                 this->local_coords.at(0), &communicators.at(3));

  // y-gathering
  MPI_Comm_split(this->globalComm,
                 this->local_coords.at(0) +
                     this->local_coords.at(2) * this->global_dims.at(0),
                 this->local_coords.at(1), &communicators.at(4));

  // z-gathering
  MPI_Comm_split(this->globalComm,
                 this->local_coords.at(0) +
                     this->local_coords.at(1) * this->global_dims.at(0),
                 this->local_coords.at(2), &communicators.at(5));

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
  for (int i = 0; i < this->dimension; ++i) {
    MPI_Cart_shift(this->globalComm, i, 1, &rank_left, &rank_right);
    neighbors.push_back(rank_left);
    neighbors.push_back(rank_right);
  }
}

template <class T, class W> void Histogram_LB<T, W>::balance(int step) {

  int i = 2 - step % 3;

  // required work values for the scheme
  W workSumLocal = (W)0;
  W workSumDimension;
  W workTotalDimension;
  W workAvgDimension;
  W workMinDimension;
  W workMaxDimension;

  // vector to store all values of the histogram for the current dimension
  std::vector<W> work_dimension(nBins.at(i));

  // compute total number of bins in dimension
  int nBinsDimension;
  MPI_Allreduce(&(nBins.at(i)), &nBinsDimension, 1, MPI_INT, MPI_SUM,
                communicators.at(i + 3));

  std::vector<W> work_collection(nBinsDimension);
  std::vector<int> histogramSlices(this->global_dims.at(i));
  std::vector<W> work_new_dimension(this->global_dims.at(i), 0.0);

  // collect how many bins from each process will be received
  MPI_Allgather(&(nBins.at(i)), 1, MPI_INT, histogramSlices.data(), 1, MPI_INT,
                communicators.at(i + 3));

  // add up work to get total work on domain
  for (int n = 0; n < nBins.at(i); ++n) {
    //        std::cout << "DEBUG: " << this->work.size() << " " << n << " " <<
    //        i << " on " << this->localRank << std::endl;
    workSumLocal += this->work.at(n);
  }

  // compute total work in current dimension
  MPI_Allreduce(&workSumLocal, &workSumDimension, 1, MPIDataTypeW, MPI_SUM,
                communicators.at(i));

  // compute total work for current dimension
  MPI_Allreduce(&workSumDimension, &workTotalDimension, 1, MPIDataTypeW,
                MPI_SUM, communicators.at(i + 3));

  int avg_num = this->global_dims.at(i);
  workAvgDimension = workTotalDimension / (W)avg_num;

  // compute local slice of the histogram
  MPI_Allreduce(this->work.data(), work_dimension.data(), nBins.at(i),
                MPIDataTypeW, MPI_SUM, communicators.at(i));
  // displacement array
  std::vector<int> displs(this->global_dims.at(i), 0);
  int tmp = 0;

  for (int n = 0; n < this->global_dims.at(i); ++n) {
    displs[n] = tmp;
    tmp += histogramSlices.at(n);
  }

  // gather complete histogram in current dimension
  MPI_Allgatherv(work_dimension.data(), nBins.at(i), MPIDataTypeW,
                 work_collection.data(), histogramSlices.data(), displs.data(),
                 MPIDataTypeW, communicators[i + 3]);

  int current_slice = 0;

  // size of one bin
  T size = (this->sysSize[2 * i + 1] - this->sysSize[2 * i]) / (T)work_collection.size();

  // minimum size of one slice due to minimum domain width
  int minBinsPerSlice = ceil(this->minSize[i]/size);
  if (minBinsPerSlice < 1) minBinsPerSlice = 1;

  if (minBinsPerSlice * this->global_dims.at(i) > (int)work_collection.size()) {
    std::cout << "ERROR on process: " << this->localRank << std::endl;
    throw InternalErrorException(
        __FILE__, __func__, __LINE__,
        "Less histogram bins than required due to minimum domain size!");
  }

  // TODO: compute cumulative function - up to workAvgDimension work in each
  // box

  // consider minimum domain size
  for (int idx = 0; idx < minBinsPerSlice - 1; ++idx)
    work_new_dimension.at(current_slice) += work_collection.at(idx);

  for (int idx = minBinsPerSlice; idx < (int)work_collection.size(); ++idx) {
    work_new_dimension.at(current_slice) += work_collection.at(idx);
    if ((idx < (int)work_collection.size() - 1 &&
        work_new_dimension.at(current_slice) + work_collection.at(idx+1)
           > workAvgDimension)
        || idx + (this->global_dims.at(i) - current_slice) * minBinsPerSlice >= work_collection.size()
       ) {
      histogramSlices.at(current_slice) = idx;

      // calculate last slice with histogram size
      if (current_slice == (this->global_dims.at(i) - 1)) {
        histogramSlices.at(current_slice) = work_collection.size();
        W tmp_work = (W)0;
        for (int j = 0; j < this->global_dims.at(i) - 1; ++j)
          tmp_work += work_new_dimension.at(j);
        work_new_dimension.at(current_slice) = workTotalDimension - tmp_work;
        break;
      }

      current_slice++;

      // consider minimum domain size
      for (int i_bin=0; i_bin<minBinsPerSlice-1; ++i_bin) {
        idx++;
        work_new_dimension.at(current_slice) += work_collection.at(idx);
      }
    }
  }

  // TODO: compute width of domain
  T up = (this->local_coords.at(i) == this->global_dims.at(i) - 1)
             ? (T)work_collection.size()
             : (T)histogramSlices.at(this->local_coords.at(i));
  T down = (this->local_coords.at(i) == 0)
               ? (T)0
               : histogramSlices.at(this->local_coords.at(i) - 1);

  this->prevVertices = this->vertices;

  this->vertices.at(0)[i] = (T)down * size + this->sysSize[2 * i];
  this->vertices.at(1)[i] = (T)up * size + this->sysSize[2 * i];

  // compute min / max work in current dimension
  MPI_Allreduce(work_new_dimension.data() + this->local_coords.at(i),
                &workMinDimension, 1, MPIDataTypeW, MPI_MIN, this->globalComm);
  MPI_Allreduce(work_new_dimension.data() + this->local_coords.at(i),
                &workMaxDimension, 1, MPIDataTypeW, MPI_MAX, this->globalComm);

  estimatedEfficiency = (W)1.0 - 
                        (workMaxDimension - workMinDimension) /
                        (workMaxDimension + workMinDimension);
  /*
  if (this->localRank == 0)
    std::cout << "HISTOGRAM: " << i << " " << workMinDimension << " "
              << workMaxDimension << " "
              << (workMaxDimension - workMinDimension) /
                     (workMaxDimension + workMinDimension)
              << std::endl;
  */
  find_neighbors();
}

template <class T, class W> void Histogram_LB<T, W>::find_neighbors() {
  neighbors.clear();
  // collect work from right neighbor plane
  MPI_Request sreq, rreq;
  MPI_Status sstat, rstat;
  // array to store neighbor vertices in Y/Z direction (reused)
  T *vertices_loc = new T[4];
  T *vertices_rem =
      new T[8 * this->global_dims.at(0) * this->global_dims.at(1)];

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
  vertices_loc[0] = this->vertices.at(0)[0];
  vertices_loc[1] = this->vertices.at(1)[0];
  MPI_Allgather(vertices_loc, 2, MPIDataTypeT,
                vertices_rem + 2 * this->global_dims.at(0), 2, MPIDataTypeT,
                communicators.at(1));

  // exchange local column information with upper neighbor in Y direction (cart
  // grid)
  MPI_Irecv(vertices_rem, 2 * this->global_dims.at(0), MPIDataTypeT, rank_left,
            0, this->globalComm, &rreq);
  MPI_Isend(vertices_rem + 2 * this->global_dims.at(0),
            2 * this->global_dims.at(0), MPIDataTypeT, rank_right, 0,
            this->globalComm, &sreq);

  // determine the offset in ranks
  offset_coords[0] = 0;
  offset_coords[1] = this->local_coords.at(1) - 1;
  offset_coords[2] = this->local_coords.at(2);

  rem_coords[1] = offset_coords[1];
  rem_coords[2] = offset_coords[2];

  MPI_Cart_rank(this->globalComm, offset_coords, &rank_offset);

#ifdef ALL_DEBUG_ENABLED  
  MPI_Barrier(this->globalComm);
  if (this->localRank == 0)
    std::cout << "HISTOGRAM: neighbor_find y-communication" << std::endl;
#endif  
  // wait for communication
  MPI_Wait(&sreq, &sstat);
  MPI_Wait(&rreq, &rstat);

  // iterate about neighbor borders to determine the neighborship relation
  nNeighbors.at(2) = 0;
  for (int x = 0; x < this->global_dims.at(0); ++x) {
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
  MPI_Irecv(vertices_rem, 2 * this->global_dims.at(0), MPIDataTypeT, rank_right,
            0, this->globalComm, &rreq);
  MPI_Isend(vertices_rem + 2 * this->global_dims.at(0),
            2 * this->global_dims.at(0), MPIDataTypeT, rank_left, 0,
            this->globalComm, &sreq);

  // determine the offset in ranks
  offset_coords[0] = 0;
  offset_coords[1] = this->local_coords.at(1) + 1;
  offset_coords[2] = this->local_coords.at(2);

  rem_coords[1] = offset_coords[1];
  rem_coords[2] = offset_coords[2];

  MPI_Cart_rank(this->globalComm, offset_coords, &rank_offset);

#ifdef ALL_DEBUG_ENABLED  
  MPI_Barrier(this->globalComm);
  if (this->localRank == 0)
    std::cout << "HISTOGRAM: neighbor_find y-communication 2" << std::endl;
#endif
  // wait for communication
  MPI_Wait(&sreq, &sstat);
  MPI_Wait(&rreq, &rstat);

  // iterate about neighbor borders to determine the neighborship relation
  nNeighbors.at(3) = 0;
  for (int x = 0; x < this->global_dims.at(0); ++x) {
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
  vertices_loc[0] = this->vertices.at(0)[0];
  vertices_loc[1] = this->vertices.at(1)[0];
  vertices_loc[2] = this->vertices.at(0)[1];
  vertices_loc[3] = this->vertices.at(1)[1];

  MPI_Barrier(this->globalComm);

  MPI_Allgather(vertices_loc, 4, MPIDataTypeT,
                vertices_rem +
                    4 * this->global_dims.at(0) * this->global_dims.at(1),
                4, MPIDataTypeT, communicators.at(2));

  // exchange local column information with upper neighbor in Z direction (cart
  // grid)
  MPI_Irecv(vertices_rem, 4 * this->global_dims.at(0) * this->global_dims.at(1),
            MPIDataTypeT, rank_left, 0, this->globalComm, &rreq);
  MPI_Isend(vertices_rem +
                4 * this->global_dims.at(0) * this->global_dims.at(1),
            4 * this->global_dims.at(0) * this->global_dims.at(1), MPIDataTypeT,
            rank_right, 0, this->globalComm, &sreq);

  // determine the offset in ranks
  offset_coords[0] = 0;
  offset_coords[1] = 0;
  offset_coords[2] = this->local_coords.at(2) - 1;

  rem_coords[2] = offset_coords[2];

  MPI_Cart_rank(this->globalComm, offset_coords, &rank_offset);

#ifdef ALL_DEBUG_ENABLED  
  MPI_Barrier(this->globalComm);
  if (this->localRank == 0)
    std::cout << "HISTOGRAM: neighbor_find z-communication" << std::endl;
#endif
  // wait for communication
  MPI_Wait(&sreq, &sstat);
  MPI_Wait(&rreq, &rstat);

  // iterate about neighbor borders to determine the neighborship relation
  nNeighbors.at(4) = 0;
  for (int y = 0; y < this->global_dims.at(1); ++y) {
    for (int x = 0; x < this->global_dims.at(0); ++x) {
      if (((vertices_rem[4 * (x + y * this->global_dims.at(0)) + 2] <=
                vertices_loc[2] &&
            vertices_loc[2] <
                vertices_rem[4 * (x + y * this->global_dims.at(0)) + 3]) ||
           (vertices_rem[4 * (x + y * this->global_dims.at(0)) + 2] <
                vertices_loc[3] &&
            vertices_loc[3] <=
                vertices_rem[4 * (x + y * this->global_dims.at(0)) + 3]) ||
           (vertices_rem[4 * (x + y * this->global_dims.at(0)) + 2] >=
                vertices_loc[2] &&
            vertices_loc[2] <
                vertices_rem[4 * (x + y * this->global_dims.at(0)) + 3] &&
            vertices_loc[3] >=
                vertices_rem[4 * (x + y * this->global_dims.at(0)) + 3])))
        if (((vertices_rem[4 * (x + y * this->global_dims.at(0))] <=
                  vertices_loc[0] &&
              vertices_loc[0] <
                  vertices_rem[4 * (x + y * this->global_dims.at(0)) + 1]) ||
             (vertices_rem[4 * (x + y * this->global_dims.at(0))] <
                  vertices_loc[1] &&
              vertices_loc[1] <=
                  vertices_rem[4 * (x + y * this->global_dims.at(0)) + 1]) ||
             (vertices_rem[4 * (x + y * this->global_dims.at(0))] >=
                  vertices_loc[0] &&
              vertices_loc[0] <
                  vertices_rem[4 * (x + y * this->global_dims.at(0)) + 1] &&
              vertices_loc[1] >=
                  vertices_rem[4 * (x + y * this->global_dims.at(0)) + 1]))) {
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
  MPI_Irecv(vertices_rem, 4 * this->global_dims.at(0) * this->global_dims.at(1),
            MPIDataTypeT, rank_right, 0, this->globalComm, &rreq);
  MPI_Isend(vertices_rem +
                4 * this->global_dims.at(0) * this->global_dims.at(1),
            4 * this->global_dims.at(0) * this->global_dims.at(1), MPIDataTypeT,
            rank_left, 0, this->globalComm, &sreq);

  // determine the offset in ranks
  offset_coords[0] = 0;
  offset_coords[1] = 0;
  offset_coords[2] = this->local_coords.at(2) + 1;

  rem_coords[2] = offset_coords[2];

  MPI_Cart_rank(this->globalComm, offset_coords, &rank_offset);

#ifdef ALL_DEBUG_ENABLED  
  MPI_Barrier(this->globalComm);
  if (this->localRank == 0)
    std::cout << "HISTOGRAM: neighbor_find z-communication 2" << std::endl;
#endif

  // wait for communication
  MPI_Wait(&sreq, &sstat);
  MPI_Wait(&rreq, &rstat);
#ifdef ALL_DEBUG_ENABLED  
  if (this->localRank == 0)
    std::cout << "HISTOGRAM: neighbor_find z-communication 2" << std::endl;
#endif
  // iterate about neighbor borders to determine the neighborship relation
  nNeighbors.at(5) = 0;
  for (int y = 0; y < this->global_dims.at(1); ++y) {
    for (int x = 0; x < this->global_dims.at(0); ++x) {
      if (((vertices_rem[4 * (x + y * this->global_dims.at(0)) + 2] <=
                vertices_loc[2] &&
            vertices_loc[2] <
                vertices_rem[4 * (x + y * this->global_dims.at(0)) + 3]) ||
           (vertices_rem[4 * (x + y * this->global_dims.at(0)) + 2] <
                vertices_loc[3] &&
            vertices_loc[3] <=
                vertices_rem[4 * (x + y * this->global_dims.at(0)) + 3]) ||
           (vertices_rem[4 * (x + y * this->global_dims.at(0)) + 2] >=
                vertices_loc[2] &&
            vertices_loc[2] <
                vertices_rem[4 * (x + y * this->global_dims.at(0)) + 3] &&
            vertices_loc[3] >=
                vertices_rem[4 * (x + y * this->global_dims.at(0)) + 3])))
        if (((vertices_rem[4 * (x + y * this->global_dims.at(0))] <=
                  vertices_loc[0] &&
              vertices_loc[0] <
                  vertices_rem[4 * (x + y * this->global_dims.at(0)) + 1]) ||
             (vertices_rem[4 * (x + y * this->global_dims.at(0))] <
                  vertices_loc[1] &&
              vertices_loc[1] <=
                  vertices_rem[4 * (x + y * this->global_dims.at(0)) + 1]) ||
             (vertices_rem[4 * (x + y * this->global_dims.at(0))] >=
                  vertices_loc[0] &&
              vertices_loc[0] <
                  vertices_rem[4 * (x + y * this->global_dims.at(0)) + 1] &&
              vertices_loc[1] >=
                  vertices_rem[4 * (x + y * this->global_dims.at(0)) + 1]))) {
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
std::vector<int> &Histogram_LB<T, W>::getNeighbors() {
  return neighbors;
}

template <class T, class W>
W Histogram_LB<T, W>::getEstimatedEfficiency()
{
    return estimatedEfficiency;
}

}//namespace ALL

#endif
