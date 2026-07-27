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

#ifndef ALL_MAIN_HEADER_INC
#define ALL_MAIN_HEADER_INC

#include "ALL_CustomExceptions.hpp"
#include "ALL_LB.hpp"
#include "ALL_Point.hpp"
#include "ALL_Tensor.hpp"
#ifdef ALL_VORONOI_ACTIVE
#include "ALL_Voronoi.hpp"
#endif
#include "ALL_ForceBased.hpp"
#include "ALL_Histogram.hpp"
#include "ALL_Staggered.hpp"
#include <algorithm>
#include <iomanip>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <cerrno>
#include <sys/stat.h>

#ifdef ALL_VTK_OUTPUT
#include <limits>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkMPIController.h>
#include <vtkPolyhedron.h>
#include <vtkProgrammableFilter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVoxel.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#ifdef VTK_CELL_ARRAY_V2
#include <vtkNew.h>
#endif
#endif

namespace ALL {

#define ALL_ESTIMATION_LIMIT_DEFAULT 0

/// enum type to describe the different balancing methods
enum LB_t : int {
  /// staggered grid load balancing
  STAGGERED = 0,
  /// tensor based load balancing
  TENSOR = 1,
  /// unstructured-mesh load balancing
  FORCEBASED = 2,
#ifdef ALL_VORONOI_ACTIVE
  /// voronoi cell based load balancing
  VORONOI = 3,
#endif
  /// histogram based load balancing
  HISTOGRAM = 4,
  /// tensor based load balancing using maximum values
  TENSOR_MAX = 5,
  /// Unimplemented load balancing NEVER SUPPORTED
  UNIMPLEMENTED = 128
};

/// @tparam T data type for vertices and related data
/// @tparam W data type for work and related data

template <class T, class W> class ALL {
public:
  /// default constructor, that sets up internal data structures and initializes
  /// internal values shared between all balancing methods
  ALL()
      : loadbalancing_step(0)
#ifdef ALL_VTK_OUTPUT
        ,
        vtk_init(false)
#endif
  {
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

  /// constructor providing the balancing method to use, the dimensions of the
  /// vertices and the gamma correction value and setting the chosen balancing
  /// methods up
  /// @param[in] m the balancing method to be used
  /// @param[in] d the dimension of the vertices to be used (most methods
  /// currently only support 3D vertices)
  /// @param[in] g the value used for the correction value, if required by the
  /// chosen balancing method
  ALL(const LB_t m, const int d, const T g) : ALL() {
    method = m;
    switch (method) {
    case LB_t::TENSOR:
      balancer.reset(new Tensor_LB<T, W>(d, (W)0, g));
      break;
    case LB_t::STAGGERED:
      balancer.reset(new Staggered_LB<T, W>(d, (W)0, g));
      break;
    case LB_t::FORCEBASED:
      balancer.reset(new ForceBased_LB<T, W>(d, (W)0, g));
      break;
#ifdef ALL_VORONOI_ACTIVE
    case LB_t::VORONOI:
      balancer.reset(new Voronoi_LB<T, W>(d, (W)0, g));
      break;
#endif
    case LB_t::HISTOGRAM:
      balancer.reset(new Histogram_LB<T, W>(d, std::vector<W>(10), g));
      break;
    case LB_t::TENSOR_MAX:
      balancer.reset(new TensorMax_LB<T, W>(d, (W)0, g));
      break;
    default:
      throw InvalidArgumentException(__FILE__, __func__, __LINE__,
                                     "Unknown type of loadbalancing passed.");
    }
    balancer->setDimension(d);
    balancer->setGamma(g);
    estimation_limit = ALL_ESTIMATION_LIMIT_DEFAULT;
  }

  /// constructor to setup the method, the dimension of the used vertices, the
  /// correction value and also providing a first set of vertices to start from
  /// @param[in] m the balancing method to be used
  /// @param[in] d the dimension of the vertices to be used (most methods
  /// currently only support 3D vertices)
  /// @param[in] inp the set of vertices to be used in the
  /// balancing step
  /// @param[in] g the value used for the correction value, if required by the
  /// chosen balancing method
  ALL(const LB_t m, const int d, const std::vector<Point<T>> &inp, const T g)
      : ALL(m, d, g) {
    balancer->setVertices(inp);
    calculate_outline();
  }

  /// destructor
  ~ALL() = default;

  /// method to provide a new set of vertices
  /// @param[in] inp reference to a vector of ALL::Point objects describing
  /// the vertices to be used in the balancing step
  void setVertices(const std::vector<Point<T>> &inp) {
    int dim = inp.at(0).getDimension();
    if (dim != balancer->getDimension())
      throw PointDimensionMissmatchException(
          __FILE__, __func__, __LINE__,
          "Dimension of ALL::Points in input vector do not match dimension of "
          "ALL "
          "object.");
    balancer->setVertices(inp);
    calculate_outline();
  }

  /// method to set the communicator to be used in the balancing step, if a non
  /// cartesian communicator is provided, a cartesian communicator is created as
  /// a copy and used internally
  /// @param[in] comm the communicator to be used
  void setCommunicator(const MPI_Comm comm) {
    int comm_type;
    MPI_Topo_test(comm, &comm_type);
    if (comm_type == MPI_CART) {
      communicator = comm;
      balancer->setCommunicator(communicator);
      // set procGridLoc and procGridDim
      const int dimension = balancer->getDimension();
      int location[dimension];
      int size[dimension];
      int periodicity[dimension];
      MPI_Cart_get(communicator, dimension, size, periodicity, location);
      procGridLoc.assign(location, location + dimension);
      procGridDim.assign(size, size + dimension);
      procGridSet = true;
    } else {
      int rank;
      MPI_Comm_rank(comm, &rank);
      if (procGridSet) {
        // compute 1D index from the location in the process grid (using z-y-x
        // ordering, as MPI_Dims_create)
        int idx;
        switch (balancer->getDimension()) {
        case 3:
          idx = procGridLoc.at(2) + procGridLoc.at(1) * procGridDim.at(2) +
                procGridLoc.at(0) * procGridDim.at(2) * procGridDim.at(1);
          break;
        case 2:
          idx = procGridLoc.at(1) + procGridLoc.at(0) * procGridDim.at(1);
          break;
        case 1:
          idx = procGridLoc.at(0);
          break;
        default:
          throw InternalErrorException(
              __FILE__, __func__, __LINE__,
              "ALL cannot deal with more then three dimensions (or less then "
              "one).");
          break;
        }

        // create new temporary communicator with correct ranks
        MPI_Comm temp_comm;
        MPI_Comm_split(comm, 0, idx, &temp_comm);

        std::vector<int> periods(3, 1);

        // transform temporary communicator to cartesian communicator
        MPI_Cart_create(temp_comm, balancer->getDimension(), procGridDim.data(),
                        periods.data(), 0, &communicator);
        balancer->setCommunicator(communicator);
      } else {
        throw InternalErrorException(
            __FILE__, __func__, __LINE__,
            "When using a non-cartesian communicator process grid parameters "
            "required (not set).");
      }
    }
  }

  /// method the get the correction value used in the balancing method
  ///@result T the correction value for the balancing method
  T getGamma() { return balancer->getGamma(); };

  /// method to set the correction value to be used in the balancing method
  /// @param[in] g the correctin value to be used
  void setGamma(const T g) { balancer->setGamma(g); };

  /// method to set a scalar work for the local domain
  /// @param[in] work the scalar work for the local domain
  void setWork(const W work) {
    // set new value for work (single value for whole domain)
    balancer->setWork(work);
  }

  // method to set a multi-dimensional work for the local domain
  ///@param[in] work std::vector<W> containing the work for the local domain
  void setWork(const std::vector<W> &work) { balancer->setWork(work); }

  /// method to call the setup of the chosen balancing method (not all methods
  /// require a setup)
  void setup() { balancer->setup(); }

  /// method the trigger the balancing step, that updates the vertices according
  /// to the previously provided work and chosen method
  /// @param[in] internal toggles if internal steps are performed, needed
  /// for some recursive calls of the routine by some methods
  /// @result std::vector<ALL::Point<T>> containing the shifted set of vertices
  std::vector<Point<T>> &balance() {
    nVertices = balancer->getNVertices();
    switch (method) {
    case LB_t::TENSOR:
      balancer->balance(loadbalancing_step);
      calculate_outline();
      break;
    case LB_t::STAGGERED:

      /* ToDo: Check if repetition is useful at all and
       *       change it to be included for all methods,
       *       not only STAGGERED */
      /*
      // estimation to improve speed of convergence
      // changing vertices during estimation
      std::vector<Point<T>> tmp_vertices = balancer->getVertices();
      // saved original vertices
      std::vector<Point<T>> old_vertices = balancer->getVertices();
      // old temporary vertices
      std::vector<Point<T>> old_tmp_vertices = balancer->getVertices();
      // result temporary vertices
      std::vector<Point<T>> result_vertices = balancer->getVertices();
      // temporary work
      W tmp_work = balancer->getWork().at(0);
      // old temporary work
      W old_tmp_work = balancer->getWork().at(0);

      W max_work;
      W sum_work;
      double d_max;
      int n_ranks;

      // compute work density on local domain
      double V = 1.0;
      for (int i = 0; i < balancer->getDimension(); ++i)
          V *= ( outline.at(1)[i] - outline.at(0)[i] );
      double rho = workArray.at(0) / V;

      // collect maximum work in system
      MPI_Allreduce(workArray.data(), &max_work, 1, MPIDataTypeW, MPI_MAX,
      communicator); MPI_Allreduce(workArray.data(), &sum_work, 1,
      MPIDataTypeW, MPI_SUM, communicator);
      MPI_Comm_size(communicator,&n_ranks);
      d_max = (double)(max_work) * (double)n_ranks / (double)sum_work - 1.0;


      // internal needs to be readded to argument list
      const bool internal = false
      for (int i = 0; i < estimation_limit && d_max > 0.1 && !internal; ++i)
      {
          balancer->setWork(tmp_work);
          balancer->setup();
          balancer->balance(true);
          bool sane = true;
          // check if the estimated boundaries are not too deformed
          for (int j = 0; j < balancer->getDimension(); ++j)
              sane = sane && (result_vertices.at(0)[j] <=
      old_vertices.at(1)[j]);
          MPI_Allreduce(MPI_IN_PLACE,&sane,1,MPI_CXX_BOOL,MPI_LAND,communicator);
          if (sane)
          {
              old_tmp_vertices = tmp_vertices;
              tmp_vertices = result_vertices;
              V = 1.0;
              for (int i = 0; i < balancer->getDimension(); ++i)
                  V *= ( tmp_vertices.at(1)[i] - tmp_vertices.at(0)[i] );
              old_tmp_work = tmp_work;
              tmp_work = rho * V;
          }
          else if (!sane || i == estimation_limit - 1)
          {
              balancer->getVertices() = old_tmp_vertices;
              workArray->at(0) = old_tmp_work;
              calculate_outline();
              i = estimation_limit;
          }
      }

      ((Staggered_LB<T,W>*)balancer.get())->setVertices(outline->data());
      */
      balancer->balance(loadbalancing_step);
      calculate_outline();
      break;
    case LB_t::FORCEBASED:
      balancer->balance(loadbalancing_step);
      calculate_outline();
      break;
#ifdef ALL_VORONOI_ACTIVE
    case LB_t::VORONOI:
      balancer->balance(loadbalancing_step);
      break;
#endif
    case LB_t::HISTOGRAM:
      balancer->balance(loadbalancing_step);
      calculate_outline();
      break;
    case LB_t::TENSOR_MAX:
      balancer->balance(loadbalancing_step);
      calculate_outline();
      break;

    default:
      throw InvalidArgumentException(__FILE__, __func__, __LINE__,
                                     "Unknown type of loadbalancing passed.");
    }
    loadbalancing_step++;
    return balancer->getVertices();
  }

  /**********************/
  /*  getter functions  */
  /**********************/

  /// get the vertices before performing the load-balancing step
  /// @result std::vector<ALL::Point<T>> containing the vertices before the
  /// balancing step
  std::vector<Point<T>> &getPrevVertices() {
    return balancer->getPrevVertices();
  };

  /// get the resulting vertices after the load-balancing step, if it has been
  /// performed
  /// @result std::vector<ALL::Point<T>> containing the resulting vertices after
  /// the balancing step
  std::vector<Point<T>> &getVertices() { return balancer->getVertices(); };

  /// get the dimension of the ALL::Point<T> objects used to describe the
  /// vertices of the domains
  /// @result int containing the dimension of the vertices
  int getDimension() { return balancer->getDimension(); }

  /// method to get the work provided to the method
  /// @param[out] result reference to std::vector<W> to store the vector of work
  /// if an array of work was provided, e.g. for the histogram method
  void getWork(std::vector<W> &result) { result = balancer->getWork(); }

  /// method to get the work provided to the method
  /// @result scalar work or first value of vector work
  W getWork() { return balancer->getWork().at(0); }

  /// method to provide a list of the ranks of the neighbors the local domain
  /// has in all directions
  /// @result vector if neighboring ranks
  std::vector<int> &getNeighbors() { return balancer->getNeighbors(); }

  /// method to provide a list of neighboring vertices, e.g. required for
  /// VORONOI
  /// @result vector of neighboring vertices
  /// neighboring vertices are stored in
  std::vector<T> &getNeighborVertices() {
    return balancer->getNeighborVertices();
  }

  /// method to provide the current load-balancing efficiency
  /// @result the current LB efficiency (only valid before balance() was called)
  W getEfficiency(){
    return balancer->getEfficiency();
  }

  /// method to provide an estimated work efficiency after the balancing
  /// @result the estimated LB efficieny (only valid after balance() was called)
  W getEstimatedEfficiency() {
    // currently only implemented for HISTOGRAM
    return balancer->getEstimatedEfficiency();
  }

  /**********************/
  /*  setter functions  */
  /**********************/

  /// method to set the size of the system, e.g. required for the HISTOGRAM
  /// balancing method
  /// @param[in] sysSize reference to a vector of T containing the size of
  /// the orthogonal system in the following format: (x_min, x_max, y_min,
  /// y_max, z_min, z_max)
  void setSysSize(const std::vector<T> &sysSize) {
    balancer.get()->setSysSize(sysSize);
  }

  /// method to set method specific data, that is not required by all different
  /// balancing methods
  /// @param[in] data pointer to a struct or other data object in the format the
  /// methods requires
  ///
  /// For the histogram method the number of bins per dimensions can be given
  /// as a C array of `int`s. The length of the array must be the number of
  /// dimensions given to the load balancer. BE CAREFUL this is not enforced!
  void setMethodData(const void *data) {
    balancer.get()->setAdditionalData(data);
  }

  /// method to set the parameters of the used cartesian communicator
  /// @param[in] loc the cartesian coordinates of the local domain in the
  /// process grid
  /// @param[in] size the size of the cartesian process grid
  void setProcGridParams(const std::vector<int> &loc,
                         const std::vector<int> &size) {
    // todo(s.schulz): We should probably throw if the dimension does not match
    procGridLoc.insert(procGridLoc.begin(), loc.begin(), loc.end());
    procGridDim.insert(procGridDim.begin(), size.begin(), size.end());
    procGridSet = true;
  }

  /// method to set the process tag output in VTK output
  /// @param[in] tag
  void setProcTag(int tag) { procTag = tag; }

  /// method the set the minimum domain sizes in all directions, can be required
  /// if using linked cell algorithms and wanting to prevent data exchange with
  /// next-near neighbors instead of only with next neighbors
  /// @param[in] minSize the minimum size of a domain in each dimension
  void setMinDomainSize(const std::vector<T> &minSize) {
    (balancer.get())->setMinDomainSize(minSize);
  }

#ifdef ALL_VTK_OUTPUT
  /// method to create VTK based output of the domains used in the
  /// load-balancing process (for now only orthogonal domains are supported)
  /// @param[in] step the number of the loadbalancing step used for numbering
  /// the output files
  void printVTKoutlines(const int step) {
    // define grid points, i.e. vertices of local domain
    auto points = vtkSmartPointer<vtkPoints>::New();
    // allocate space for eight points (describing the cuboid around the domain)
    points->Allocate(8 * balancer->getDimension());

    int n_ranks;
    int local_rank;

    if (!vtk_init) {
      controller = vtkMPIController::New();
      controller->Initialize();
      controller->SetGlobalController(controller);
      vtk_init = true;
    }

    MPI_Comm_rank(communicator, &local_rank);
    MPI_Comm_size(communicator, &n_ranks);

    std::vector<Point<W>> tmp_outline(2);
    std::vector<W> tmp_0(
        {outline.at(0)[0], outline.at(0)[1], outline.at(0)[2]});
    std::vector<W> tmp_1(
        {outline.at(1)[0], outline.at(1)[1], outline.at(1)[2]});
    tmp_outline.at(0) = Point<W>(tmp_0);
    tmp_outline.at(1) = Point<W>(tmp_1);

    for (auto z = 0; z <= 1; ++z)
      for (auto y = 0; y <= 1; ++y)
        for (auto x = 0; x <= 1; ++x) {
          points->InsertPoint(x + 2 * y + 4 * z, tmp_outline.at(x)[0],
                              tmp_outline.at(y)[1], tmp_outline.at(z)[2]);
        }

    auto hexa = vtkSmartPointer<vtkVoxel>::New();
    for (int i = 0; i < 8; ++i)
      hexa->GetPointIds()->SetId(i, i);

    auto cellArray = vtkSmartPointer<vtkCellArray>::New();
    cellArray->InsertNextCell(hexa);

    // define work array (length = 1)
    auto work = vtkSmartPointer<vtkFloatArray>::New();
    work->SetNumberOfComponents(1);
    work->SetNumberOfTuples(1);
    work->SetName("work");
    W total_work = std::accumulate(balancer->getWork().begin(),
                                   balancer->getWork().end(), (W)0);
    work->SetValue(0, total_work);

    // define rank array (length = 4, x,y,z, rank)
    int rank = 0;
    MPI_Comm_rank(communicator, &rank);
    auto coords = vtkSmartPointer<vtkFloatArray>::New();
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(1);
    coords->SetName("coords");
    coords->SetValue(0, procGridLoc.at(0));
    coords->SetValue(1, procGridLoc.at(1));
    coords->SetValue(2, procGridLoc.at(2));

    auto rnk = vtkSmartPointer<vtkFloatArray>::New();
    rnk->SetNumberOfComponents(1);
    rnk->SetNumberOfTuples(1);
    rnk->SetName("MPI rank");
    rnk->SetValue(0, rank);

    // define tag array (length = 1)
    auto tag = vtkSmartPointer<vtkIntArray>::New();
    tag->SetNumberOfComponents(1);
    tag->SetNumberOfTuples(1);
    tag->SetName("tag");
    tag->SetValue(0, procTag);

    // determine extent of system
    W known_unused global_extent[6];
    W local_min[3];
    W local_max[3];
    W global_min[3];
    W global_max[3];
    W width[3];
    W volume = (W)1.0;

    for (int i = 0; i < 3; ++i) {
      local_min[i] = outline.at(0)[i];
      local_max[i] = outline.at(1)[i];
      width[i] = local_max[i] - local_min[i];
      volume *= width[i];
    }

    W surface = (W)2.0 * width[0] * width[1] + (W)2.0 * width[1] * width[2] +
                (W)2.0 * width[0] * width[2];

    auto expanse = vtkSmartPointer<vtkFloatArray>::New();
    expanse->SetNumberOfComponents(6);
    expanse->SetNumberOfTuples(1);
    expanse->SetName("expanse");
    expanse->SetValue(0, width[0]);
    expanse->SetValue(1, width[0]);
    expanse->SetValue(2, width[0]);
    expanse->SetValue(3, volume);
    expanse->SetValue(4, surface);
    expanse->SetValue(5, surface / volume);

    MPI_Allreduce(&local_min, &global_min, 3, MPIDataTypeW, MPI_MIN,
                  communicator);
    MPI_Allreduce(&local_max, &global_max, 3, MPIDataTypeW, MPI_MAX,
                  communicator);

    for (int i = 0; i < 3; ++i) {
      global_extent[2 * i] = global_min[i];
      global_extent[2 * i + 1] = global_max[i];
    }

    // create a structured grid and assign data to it
    auto unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_VOXEL, cellArray);
    unstructuredGrid->GetCellData()->AddArray(work);
    unstructuredGrid->GetCellData()->AddArray(expanse);
    unstructuredGrid->GetCellData()->AddArray(rnk);
    unstructuredGrid->GetCellData()->AddArray(coords);
    unstructuredGrid->GetCellData()->AddArray(tag);

    createDirectory("vtk_outline");
    std::ostringstream ss_local;
    ss_local << "vtk_outline/ALL_vtk_outline_" << std::setw(7)
             << std::setfill('0') << step << "_" << local_rank << ".vtu";

    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetInputData(unstructuredGrid);
    writer->SetFileName(ss_local.str().c_str());
    writer->SetDataModeToAscii();
    // writer->SetDataModeToBinary();
    writer->Write();

    // if (local_rank == 0)
    //{
    std::ostringstream ss_para;
    ss_para << "vtk_outline/ALL_vtk_outline_" << std::setw(7)
            << std::setfill('0') << step << ".pvtu";
    // create the parallel writer
    auto parallel_writer =
        vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
    parallel_writer->SetFileName(ss_para.str().c_str());
    parallel_writer->SetNumberOfPieces(n_ranks);
    parallel_writer->SetStartPiece(local_rank);
    parallel_writer->SetEndPiece(local_rank);
    parallel_writer->SetInputData(unstructuredGrid);
    parallel_writer->SetDataModeToAscii();
    // parallel_writer->SetDataModeToBinary();
    parallel_writer->Write();
    //}
  }

  /// method to create VTK based output of the vertices used in the
  /// load-balancing process
  /// @param[in] step the number of the loadbalancing step used for numbering
  /// the output files
  void printVTKvertices(const int step) {
    int n_ranks;
    int local_rank;

    MPI_Comm_rank(communicator, &local_rank);
    MPI_Comm_size(communicator, &n_ranks);

    // local vertices
    // (vertices + work)
    T local_vertices[nVertices * balancer->getDimension() + 1];

    for (int v = 0; v < nVertices; ++v) {
      for (int d = 0; d < balancer->getDimension(); ++d) {
        local_vertices[v * balancer->getDimension() + d] =
            balancer->getVertices().at(v)[d];
      }
    }
    local_vertices[nVertices * balancer->getDimension()] =
        (T)balancer->getWork().at(0);

    /*
    T *global_vertices;
    if (local_rank == 0) {
      global_vertices =
          new T[n_ranks * (nVertices * balancer->getDimension() + 1)];
    }
    */
    T global_vertices[n_ranks * (nVertices * balancer->getDimension() + 1)];

    // collect all works and vertices on a single process
    MPI_Gather(local_vertices, nVertices * balancer->getDimension() + 1,
               MPIDataTypeT, global_vertices,
               nVertices * balancer->getDimension() + 1, MPIDataTypeT, 0,
               communicator);

    if (local_rank == 0) {
      auto vtkpoints = vtkSmartPointer<vtkPoints>::New();
#ifdef VTK_CELL_ARRAY_V2
      vtkNew<vtkUnstructuredGrid> unstructuredGrid;
      unstructuredGrid->Allocate(n_ranks + 10);
#else
      auto unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
#endif

      // enter vertices into unstructured grid
      for (int i = 0; i < n_ranks; ++i) {
        for (int v = 0; v < nVertices; ++v) {
          vtkpoints->InsertNextPoint(
              global_vertices[i * (nVertices * balancer->getDimension() + 1) +
                              v * balancer->getDimension() + 0],
              global_vertices[i * (nVertices * balancer->getDimension() + 1) +
                              v * balancer->getDimension() + 1],
              global_vertices[i * (nVertices * balancer->getDimension() + 1) +
                              v * balancer->getDimension() + 2]);
        }
      }
      unstructuredGrid->SetPoints(vtkpoints);

      // data arrays for work and cell id
      auto work = vtkSmartPointer<vtkFloatArray>::New();
      auto cell = vtkSmartPointer<vtkFloatArray>::New();
      work->SetNumberOfComponents(1);
      work->SetNumberOfTuples(n_ranks);
      work->SetName("work");
      cell->SetNumberOfComponents(1);
      cell->SetNumberOfTuples(n_ranks);
      cell->SetName("cell id");

      for (int n = 0; n < n_ranks; ++n) {

#ifdef VTK_CELL_ARRAY_V2
        // define grid points, i.e. vertices of local domain
        vtkIdType pointIds[8] = {8 * n + 0, 8 * n + 1, 8 * n + 2, 8 * n + 3,
                                 8 * n + 4, 8 * n + 5, 8 * n + 6, 8 * n + 7};

        vtkIdType faces[48] = {
            3,         8 * n + 0, 8 * n + 2, 8 * n + 1, 
            3,         8 * n + 1, 8 * n + 2, 8 * n + 3, 
            3,         8 * n + 0, 8 * n + 4, 8 * n + 2, 
            3,         8 * n + 2, 8 * n + 4, 8 * n + 6, 
            3,         8 * n + 2, 8 * n + 6, 8 * n + 3, 
            3,         8 * n + 3, 8 * n + 6, 8 * n + 7, 
            3,         8 * n + 1, 8 * n + 5, 8 * n + 3, 
            3,         8 * n + 3, 8 * n + 5, 8 * n + 7, 
            3,         8 * n + 0, 8 * n + 4, 8 * n + 1, 
            3,         8 * n + 1, 8 * n + 4, 8 * n + 5, 
            3,         8 * n + 4, 8 * n + 6, 8 * n + 5, 
            3,         8 * n + 5, 8 * n + 6, 8 * n + 7};

        unstructuredGrid->InsertNextCell(VTK_POLYHEDRON, 8, pointIds, 12,
                                         faces);
#else
        // define grid points, i.e. vertices of local domain
        vtkIdType pointIds[8] = {8 * n + 0, 8 * n + 1, 8 * n + 2, 8 * n + 3,
                                 8 * n + 4, 8 * n + 5, 8 * n + 6, 8 * n + 7};

        auto faces = vtkSmartPointer<vtkCellArray>::New();
        // setup faces of polyhedron
        vtkIdType f0[3] = {8 * n + 0, 8 * n + 2, 8 * n + 1};
        vtkIdType f1[3] = {8 * n + 1, 8 * n + 2, 8 * n + 3};

        vtkIdType f2[3] = {8 * n + 0, 8 * n + 4, 8 * n + 2};
        vtkIdType f3[3] = {8 * n + 2, 8 * n + 4, 8 * n + 6};

        vtkIdType f4[3] = {8 * n + 2, 8 * n + 6, 8 * n + 3};
        vtkIdType f5[3] = {8 * n + 3, 8 * n + 6, 8 * n + 7};

        vtkIdType f6[3] = {8 * n + 1, 8 * n + 5, 8 * n + 3};
        vtkIdType f7[3] = {8 * n + 3, 8 * n + 5, 8 * n + 7};

        vtkIdType f8[3] = {8 * n + 0, 8 * n + 4, 8 * n + 1};
        vtkIdType f9[3] = {8 * n + 1, 8 * n + 4, 8 * n + 5};

        vtkIdType fa[3] = {8 * n + 4, 8 * n + 6, 8 * n + 5};
        vtkIdType fb[3] = {8 * n + 5, 8 * n + 6, 8 * n + 7};

        faces->InsertNextCell(3, f0);
        faces->InsertNextCell(3, f1);
        faces->InsertNextCell(3, f2);
        faces->InsertNextCell(3, f3);
        faces->InsertNextCell(3, f4);
        faces->InsertNextCell(3, f5);
        faces->InsertNextCell(3, f6);
        faces->InsertNextCell(3, f7);
        faces->InsertNextCell(3, f8);
        faces->InsertNextCell(3, f9);
        faces->InsertNextCell(3, fa);
        faces->InsertNextCell(3, fb);

        unstructuredGrid->InsertNextCell(VTK_POLYHEDRON, 8, pointIds, 12,
                                         faces->GetPointer());
#endif
        work->SetValue(
            n, global_vertices[n * (nVertices * balancer->getDimension() + 1) +
                               8 * balancer->getDimension()]);
        cell->SetValue(n, (T)n);
      }
      unstructuredGrid->GetCellData()->AddArray(work);
      unstructuredGrid->GetCellData()->AddArray(cell);

      createDirectory("vtk_vertices");
      std::ostringstream filename;
      filename << "vtk_vertices/ALL_vtk_vertices_" << std::setw(7)
               << std::setfill('0') << step << ".vtu";
      auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
      writer->SetInputData(unstructuredGrid);
      writer->SetFileName(filename.str().c_str());
      writer->SetDataModeToAscii();
      // writer->SetDataModeToBinary();
      writer->Write();

      // delete[] global_vertices;
    }
  }
#endif

private:
  /// the used balancing method
  LB_t method;

  /// outer cuboid encasing the domain (for cuboid domains identical to domain)
  /// defined by front left lower [0][*] and back upper right points [1][*]
  /// where * is 0,...,dim-1
  std::vector<Point<T>> outline;

  /// the interal balancer object
  std::unique_ptr<LB<T, W>> balancer;

  /// storage of local work (scalar work is stored in the first entry of the
  /// vector)
  std::vector<W> workArray;

  /// internal MPI communicator used in library
  MPI_Comm communicator;

  /// number of vertices (for non-cuboid domains)
  int nVertices;

  /// number of neighbors
  int nNeighbors;

  /// calculate the outline of the domain using the vertices, i.e. the smallest
  /// orthogonal domain containing the vertices used
  void calculate_outline() {
    // calculation only possible if there are vertices defining the domain
    if (balancer->getVertices().size() > 0) {
      outline.resize(2);
      // setup the outline with the first point
      for (int i = 0; i < balancer->getDimension(); ++i) {
        outline.at(0) = balancer->getVertices().at(0);
        outline.at(1) = balancer->getVertices().at(0);
      }
      // compare value of each outline point with all vertices to find the
      // maximum extend of the domain in each dimension
      for (int i = 1; i < (int)balancer->getVertices().size(); ++i) {
        Point<T> p = balancer->getVertices().at(i);
        for (int j = 0; j < balancer->getDimension(); ++j) {
          outline.at(0)[j] = std::min(outline.at(0)[j], p[j]);
          outline.at(1)[j] = std::max(outline.at(1)[j], p[j]);
        }
      }
    }
  }

  /// create directory if it does not already exist
  ///
  /// May throw FilesystemErrorException
  void createDirectory(const char *filename) {
    errno = 0;
    mode_t mode = S_IRWXU | S_IRWXG; // rwx for user and group
    if (mkdir(filename, mode)) {
      if (errno != EEXIST) {
        throw FilesystemErrorException(__FILE__, __func__, __LINE__,
                                       "Could not create output directory.");
      }
    }
    // check permission in case directory already existed, but had wrong
    // permissions
    struct stat attr;
    stat(filename, &attr);
    if ((attr.st_mode & S_IRWXU) != S_IRWXU) {
      throw FilesystemErrorException(
          __FILE__, __func__, __LINE__,
          "Possibly already existing directory does not have correct "
          "permissions (rwx) for user");
    }
  }

  /// limit for the estimation procedure, which is an experimental way to use
  /// non-updated work-loads to provide better results for the final vertex
  /// positions
  int estimation_limit;

  /// cartesian coordinates of the local domain in the process grid
  std::vector<int> procGridLoc;
  /// size of the process grid
  std::vector<int> procGridDim;
  /// check if the process grid data is set
  bool procGridSet;
  /// tag of the process
  int procTag;

  /// datatype for MPI to use when communicating vertex related data
  MPI_Datatype MPIDataTypeT;

  /// datatype for MPI to use when communicating work related data
  MPI_Datatype MPIDataTypeW;

  /// counter for the number of loadbalancing steps already performed
  int loadbalancing_step;

#ifdef ALL_VTK_OUTPUT
  /// check if the initialization of parallel VTK is already performed
  bool vtk_init;

  /// the parallel VTK controller for the domain output
  vtkMPIController *controller;

  /// the parallel VTK controller for the vertex output
  vtkMPIController *controller_vertices;
#endif
};

} // namespace ALL

#endif

// VIM modeline for indentation
// vim: sw=2 et
