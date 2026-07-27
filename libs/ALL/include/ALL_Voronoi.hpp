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

#ifndef ALL_VORONOI_HEADER_INCLUDED
#define ALL_VORONOI_HEADER_INCLUDED

// number of maximum neighboring Voronoi cells
#define ALL_VORONOI_MAX_NEIGHBORS 32
// number of subblock division for Voronoi cell creation
#define ALL_VORONOI_SUBBLOCKS 20

// depth used to find neighbor cells
#define ALL_VORONOI_NEIGHBOR_SEARCH_DEPTH 2

#ifdef ALL_VORONOI_ACTIVE
#include "ALL_CustomExceptions.hpp"
#include "ALL_LB.hpp"
#include "ALL_Point.hpp"
#include <algorithm>
#include <exception>
#include <iomanip>
#include <map>
#include <mpi.h>
#include <sstream>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <voro++.hh>
#pragma GCC diagnostic pop

namespace ALL {

/// Load-balancing scheme based on Voronoi cells and the shift of their anchor
/// points due to the ratio of work assigned to the local and neighboring
/// domains. For each domain a center of work, i.e. the average of the particle
/// locations, weighted by their indivual work share is computed (vertex 2).
/// This and the work load of the local and the neighboring domains are used to
/// compute a shift of the anchor point in the direction of the largest work
/// load.
/// @tparam T data for vertices and related data
/// @tparam W data for work and related data
template <class T, class W> class Voronoi_LB : public LB<T, W> {
public:
  /// default constructor
  Voronoi_LB() {}
  /// constructor initializing basic parameters
  /// @param[in] d the dimension of the used vertices
  /// @param[in] w the scalar work assigned to the local domain
  /// @param[in] g the correction factor gamma
  Voronoi_LB(int d, W w, T g) : LB<T, W>(d, g) { this->setWork(w); }

  /// default destructor
  ~Voronoi_LB() override;

  voro::voronoicell vc;

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

  /// method to get a list of the vertices of the neighboring domains
  /// @param[out] ret list of vertices of neighboring domains
  std::vector<T> &getNeighborVertices();

  /// method to get an estimated work distribution after the balance step
  /// (currently only implemented in ALL::HISTOGRAM!)
  /// @param [out] double providing the estimated LB after the balance step
  virtual W getEstimatedEfficiency() override {return (W)-1;};

private:
  // MPI values

  // type for MPI communication
  MPI_Datatype MPIDataTypeT;
  MPI_Datatype MPIDataTypeW;

  // list of neighbors
  std::vector<int> neighbors;
  int nNeighbors[1];

  // collection of generator points
  std::vector<T> generator_points;

  // number of ranks in global communicator
  int nDomains;
};

template <class T, class W> Voronoi_LB<T, W>::~Voronoi_LB() {}

// setup routine for the tensor-based load-balancing scheme
// requires:
//              this->globalComm (int): cartesian MPI communicator, from
//                                 which separate sub communicators
//                                 are derived in order to represent
//                                 each plane of domains in the system
template <class T, class W> void Voronoi_LB<T, W>::setup() {
  // store global communicator
  int dimension = this->getDimension();

  // no special communicator required

  // create array to store information about generator points
  // TODO: hierarchical scheme or better preselection which
  //       generator points are required for correct creation
  //       of the local cell -> large-scale Voronoi grid creation
  //                            is too expensive to be done
  //                            every step
  MPI_Comm_size(this->globalComm, &nDomains);
  MPI_Comm_rank(this->globalComm, &this->localRank);
  generator_points.resize(nDomains * (2 * dimension + 1));

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
}

template <class T, class W>
void Voronoi_LB<T, W>::balance(int known_unused loadbalancing_step) {
  int dimension = this->getDimension();
  // collect local information in array
  int info_length = 2 * dimension + 1;
  std::vector<T> local_info(info_length);
  // Todo(s.schulz): prevVertices is unused, might be a bug?
  std::vector<Point<T>> &prevVertices = this->getVertices();
  std::vector<Point<T>> &vertices = this->getVertices();
  std::vector<W> &work = this->getWork();

  for (auto i = 0; i < 2; ++i)
    for (auto d = 0; d < dimension; ++d) {
      local_info.at(d) = vertices.at(i)[d];
    }
  local_info.at(2 * dimension) = (T)work.at(0);

  MPI_Allgather(local_info.data(), info_length, MPIDataTypeT,
                generator_points.data(), info_length, MPIDataTypeT,
                this->globalComm);

  // create Voronoi-cells
  // TODO: check if system is periodic or not -> true / false!
  voro::container con_old(
      this->sysSize.at(0), this->sysSize.at(1), this->sysSize.at(2),
      this->sysSize.at(3), this->sysSize.at(4), this->sysSize.at(5),
      ALL_VORONOI_SUBBLOCKS, ALL_VORONOI_SUBBLOCKS, ALL_VORONOI_SUBBLOCKS,
      false, false, false, nDomains + 10);

  // add generator points to container
  for (auto d = 0; d < nDomains; ++d) {
    con_old.put(d, generator_points.at(info_length * d),
                generator_points.at(info_length * d + 1),
                generator_points.at(info_length * d + 2));
  }

#ifdef ALL_DEBUG_ENABLED
  // print old voronoi cell structure
  std::ostringstream ss_local_gp_pov;
  ss_local_gp_pov << "voronoi/generator_points_" << std::setw(7)
                  << std::setfill('0') << loadbalancing_step << ".pov";
  std::ostringstream ss_local_gp_gnu;
  ss_local_gp_gnu << "voronoi/generator_points_" << std::setw(7)
                  << std::setfill('0') << loadbalancing_step << ".gnu";
  std::ostringstream ss_local_vc_pov;
  ss_local_vc_pov << "voronoi/voronoi_cells_" << std::setw(7)
                  << std::setfill('0') << loadbalancing_step << ".pov";
  std::ostringstream ss_local_vc_gnu;
  ss_local_vc_gnu << "voronoi/voronoi_cells_" << std::setw(7)
                  << std::setfill('0') << loadbalancing_step << ".gnu";

  loadbalancing_step++;

  if (this->localRank == 0) {
    con_old.draw_particles(ss_local_gp_gnu.str().c_str());
    con_old.draw_cells_gnuplot(ss_local_vc_gnu.str().c_str());
    con_old.draw_particles_pov(ss_local_gp_pov.str().c_str());
    con_old.draw_cells_pov(ss_local_vc_pov.str().c_str());
  }
#endif

  // vector to store neighbor information
  std::vector<double> cell_vertices;
  voro::voronoicell_neighbor c;

  // generate loop over all old generator points
  voro::c_loop_all cl_old(con_old);

  int pid;
  double x, y, z, r;
  // reset list of neighbors
  neighbors.clear();

  std::map<int, double> neig_area;

  for (auto i = 0; i < 1; ++i) {
    // find next neighbors
    std::vector<std::vector<int>> next_neighbors(
        std::max((int)neighbors.size(), 1));
    std::vector<std::vector<double>> neighbors_area(
        std::max((int)neighbors.size(), 1));
    int idx = 0;
    if (cl_old.start()) {
      do {
        // compute only local cell
        cl_old.pos(pid, x, y, z, r);
        if (i == 0) {
          if (pid == this->localRank) {
            con_old.compute_cell(c, cl_old);
            c.neighbors(next_neighbors.at(idx));
            c.face_areas(neighbors_area.at(idx));
            for (int j = 0; j < (int)next_neighbors.at(idx).size(); ++j)
              neig_area.insert(std::pair<int, double>(
                  next_neighbors.at(idx).at(j), neighbors_area.at(idx).at(j)));
            idx++;
            if (idx == (int)neighbors.size())
              break;
          }
        } else {
          if (std::count(neighbors.begin(), neighbors.end(), pid) == 1) {
            con_old.compute_cell(c, cl_old);
            c.neighbors(next_neighbors.at(idx));
            idx++;
            if (idx == (int)neighbors.size())
              break;
          }
        }
      } while (cl_old.inc());
    }
    for (auto it = next_neighbors.begin(); it != next_neighbors.end(); ++it) {
      neighbors.insert(neighbors.begin(), it->begin(), it->end());
    }

    std::sort(neighbors.begin(), neighbors.end());
    auto uniq_it = std::unique(neighbors.begin(), neighbors.end());
    neighbors.resize(std::distance(neighbors.begin(), uniq_it));
    neighbors.erase(
        std::remove(neighbors.begin(), neighbors.end(), this->localRank),
        neighbors.end());
    neighbors.erase(std::remove_if(neighbors.begin(), neighbors.end(),
                                   [](int x) { return x < 0; }),
                    neighbors.end());
  }

  std::vector<double> volumes(neighbors.size());
  std::vector<double> surfaces(neighbors.size());

  for (int i = 0; i < (int)neighbors.size(); ++i) {
    auto it = neig_area.find(neighbors.at(i));
    surfaces.at(i) = it->second;
  }

  // Todo(s.schulz): is local_volume needed for future work or can we remove it?
  double local_volume;
  // get volumes of each neighbor cell
  if (cl_old.start()) {
    do {
      cl_old.pos(pid, x, y, z, r);
      auto it = std::find(neighbors.begin(), neighbors.end(), pid);

      if (it != neighbors.end()) {
        int idx = std::distance(neighbors.begin(), it);
        con_old.compute_cell(c, cl_old);
        volumes.at(idx) = c.volume();
      }
      if (pid == this->localRank) {
        con_old.compute_cell(c, cl_old);
        local_volume = c.volume();
      }
    } while (cl_old.inc());
  }
#ifdef ALL_DEBUG_ENABLED
  MPI_Barrier(this->globalComm);
  if (this->localRank == 0)
    std::cout << "ALL::Voronoi_LB<T,W>::balance begin checking..." << std::endl;
#endif

  // compute shift
  std::vector<T> shift(dimension, 0.0);
  T norm;

  T work_load = (T)local_info.at(info_length - 1);
  for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
    int neighbor = *it;
    work_load += generator_points.at(info_length * neighbor + info_length - 1);
  }

#ifdef ALL_DEBUG_ENABLED
  MPI_Barrier(this->globalComm);
  if (this->localRank == 0)
    std::cout << "ALL::Voronoi_LB<T,W>::balance find neighbors..." << std::endl;
#endif

  T max_diff = 20.0;

  for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
    int neighbor = *it;
    std::vector<T> diff(dimension);
    bool correct = true;
    for (int d = 0; d < dimension && correct; ++d) {
      diff.at(d) = generator_points.at(info_length * neighbor + dimension + d) -
                   local_info.at(d);

      if (diff.at(d) >
          0.5 * (this->sysSize[2 * d + 1] - this->sysSize[2 * d])) {
        diff.at(d) -= (this->sysSize[2 * d + 1] - this->sysSize[2 * d]);
        correct = false;
      } else if (diff.at(d) <
                 -0.5 * (this->sysSize[2 * d + 1] - this->sysSize[2 * d])) {
        diff.at(d) += (this->sysSize[2 * d + 1] - this->sysSize[2 * d]);
        correct = false;
      }
    }
    max_diff = std::min(max_diff,
                        sqrt(diff.at(0) * diff.at(0) + diff.at(1) * diff.at(1) +
                             diff.at(2) * diff.at(2)));
    if (correct) {
      // compute difference in work load
      T work_diff = (T)0.0;
      work_diff =
          (generator_points.at(info_length * neighbor + info_length - 1) -
           local_info.at(info_length - 1));
      T work_density =
          (generator_points.at(info_length * neighbor + info_length - 1) +
           local_info.at(info_length - 1));

      if (work_density < 1e-6)
        work_diff = (T)0;
      else
        work_diff /= work_density;
      for (int d = 0; d < dimension; ++d) {
        shift.at(d) += 0.5 * work_diff * diff.at(d) / this->gamma;
      }
    }
  }

#ifdef ALL_DEBUG_ENABLED
  MPI_Barrier(this->globalComm);
  if (this->localRank == 0)
    std::cout << "ALL::Voronoi_LB<T,W>::balance after neighbors found..."
              << std::endl;
#endif

  norm = sqrt(shift.at(0) * shift.at(0) + shift.at(1) * shift.at(1) +
              shift.at(2) * shift.at(2));

  T scale = 1.0;

  if (norm > 0.45 * max_diff)
    scale = 0.45 * max_diff / norm;

#ifdef ALL_DEBUG_ENABLED
  MPI_Barrier(this->globalComm);
  if (this->localRank == 0)
    std::cout << "ALL::Voronoi_LB<T,W>::balance find new neighbors..." << std::endl;
#endif

  // to find new neighbors
  for (int d = 0; d < dimension; ++d) {
    local_info.at(d) += scale * shift.at(d);

    // periodic correction of points
    local_info.at(d) = (local_info.at(d) < this->sysSize.at(2 * d))
                           ? this->sysSize.at(2 * d) + 1.0
                           : ((local_info.at(d) >= this->sysSize.at(2 * d + 1))
                                  ? this->sysSize.at(2 * d + 1) - 1.0
                                  : local_info.at(d));
  }

  MPI_Allgather(local_info.data(), info_length, MPIDataTypeT,
                generator_points.data(), info_length, MPIDataTypeT,
                this->globalComm);

#ifdef ALL_DEBUG_ENABLED
  MPI_Barrier(this->globalComm);
  if (this->localRank == 0)
    std::cout << "ALL::Voronoi_LB<T,W>::balance create new containers..."
              << std::endl;
#endif

  // create Voronoi-cells
  // TODO: check if system is periodic or not -> true / false!
  voro::container con_new(
      this->sysSize.at(0), this->sysSize.at(1), this->sysSize.at(2),
      this->sysSize.at(3), this->sysSize.at(4), this->sysSize.at(5),
      ALL_VORONOI_SUBBLOCKS, ALL_VORONOI_SUBBLOCKS, ALL_VORONOI_SUBBLOCKS,
      false, false, false, nDomains + 10);

  // add generator points to container
  for (int d = 0; d < nDomains; ++d) {
    con_new.put(d, generator_points.at(info_length * d),
                generator_points.at(info_length * d + 1),
                generator_points.at(info_length * d + 2));
  }

#ifdef ALL_DEBUG_ENABLED
  std::ostringstream ss_local_gp_pov2;
  ss_local_gp_pov2 << "voronoi/generator_points_" << std::setw(7)
                   << std::setfill('0') << loadbalancing_step << ".pov";
  std::ostringstream ss_local_gp_gnu2;
  ss_local_gp_gnu2 << "voronoi/generator_points_" << std::setw(7)
                   << std::setfill('0') << loadbalancing_step << ".gnu";
  std::ostringstream ss_local_vc_pov2;
  ss_local_vc_pov2 << "voronoi/voronoi_cells_" << std::setw(7)
                   << std::setfill('0') << loadbalancing_step << ".pov";
  std::ostringstream ss_local_vc_gnu2;
  ss_local_vc_gnu2 << "voronoi/voronoi_cells_" << std::setw(7)
                   << std::setfill('0') << loadbalancing_step << ".gnu";

  loadbalancing_step++;

  if (this->localRank == 0) {
    con_new.draw_particles(ss_local_gp_gnu2.str().c_str());
    con_new.draw_cells_gnuplot(ss_local_vc_gnu2.str().c_str());
    con_new.draw_particles_pov(ss_local_gp_pov2.str().c_str());
    con_new.draw_cells_pov(ss_local_vc_pov2.str().c_str());
  }
#endif

#ifdef ALL_DEBUG_ENABLED
  MPI_Barrier(this->globalComm);
  if (this->localRank == 0)
    std::cout << "ALL::Voronoi_LB<T,W>::balance storing shifted vertices..."
              << " " << vertices.size() << " " << dimension << std::endl;
#endif

  Point<T> tmp_pnt(dimension, local_info.data());

  // compute new neighboring cells and generator points

  // generate loop over all new generator points
  voro::c_loop_all cl_new(con_new);

  // reset list of neighbors
  neighbors.clear();

#ifdef ALL_DEBUG_ENABLED
  MPI_Barrier(this->globalComm);
  if (this->localRank == 0)
    std::cout << "ALL::Voronoi_LB<T,W>::balance searching neighboring vertices..."
              << std::endl;
#endif

  for (auto i = 0; i < ALL_VORONOI_NEIGHBOR_SEARCH_DEPTH; ++i) {
    // find next neighbors
    std::vector<std::vector<int>> next_neighbors(
        std::max((int)neighbors.size(), 1));
    int idx = 0;
    if (cl_new.start()) {
      do {
        // compute next voronoi cell
        cl_new.pos(pid, x, y, z, r);
        if (i == 0) {
          if (pid == this->localRank) {
            con_new.compute_cell(c, cl_new);
            c.neighbors(next_neighbors.at(idx));
            idx++;
            if (idx == (int)neighbors.size())
              break;
          }
        } else {
          if (std::count(neighbors.begin(), neighbors.end(), pid) == 1) {
            con_new.compute_cell(c, cl_new);
            c.neighbors(next_neighbors.at(idx));
            idx++;
            if (idx == (int)neighbors.size())
              break;
          }
        }
      } while (cl_new.inc());
    }
    for (auto it = next_neighbors.begin(); it != next_neighbors.end(); ++it) {
      neighbors.insert(neighbors.begin(), it->begin(), it->end());
    }

    std::sort(neighbors.begin(), neighbors.end());
    auto uniq_it = std::unique(neighbors.begin(), neighbors.end());
    neighbors.resize(std::distance(neighbors.begin(), uniq_it));
    neighbors.erase(
        std::remove(neighbors.begin(), neighbors.end(), this->localRank),
        neighbors.end());
    neighbors.erase(std::remove_if(neighbors.begin(), neighbors.end(),
                                   [](int x) { return x < 0; }),
                    neighbors.end());
  }

  // determine number of neighbors
  nNeighbors[0] = neighbors.size();

  // clear old neighbor vertices
  this->neighborVertices.clear();

  // find vertices of neighbors

  for (auto n : neighbors) {
    for (int i = 0; i < dimension; ++i)
      this->neighborVertices.push_back(generator_points.at(info_length * n + i));
  }
}

// provide list of neighbors
template <class T, class W>
std::vector<int> &Voronoi_LB<T, W>::getNeighbors() {
  return neighbors;
}


}//namespace ALL

#endif
#endif
