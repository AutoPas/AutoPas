/*
   Copyright 2018 Rene Halver, Forschungszentrum Juelich GmbH, Germany
   Copyright 2018 Godehard Sutmann, Forschungszentrum Juelich GmbH, Germany

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

   3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
 */

#include "ALL.hpp"
#include "ALL_CustomExceptions.hpp"
#include "ALL_Point.hpp"
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

#ifdef ALL_VTK_OUTPUT
#include <vtkCellArray.h>
#include <vtkCellLocator.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkVersion.h>
#include <vtkXMLPPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#endif

#define BOX_SIZE 300.0
#define N_PARTICLES 2080
#define N_GENERATE 1000
#define SEED 123789456u
#define N_LOOP 500
#define OUTPUT_INTV 1
#define MAX_NEIG 1024
#define ALL_HISTOGRAM_DEFAULT_WIDTH 1.0

#ifdef ALL_VORONOI_ACTIVE
#define ALL_VORONOI_PREP_STEPS 50
#endif

void generate_points(std::vector<ALL::Point<double>> &points,
                     std::vector<double> l, std::vector<double> u, int *lcoords,
                     int *gcoords, int dimension, int &n, int rank) {

  unsigned seed = SEED + rank;
  std::default_random_engine generator(seed);
  std::vector<std::uniform_real_distribution<double>> dist(dimension);

  std::uniform_int_distribution<int> dist2(0, 2 * n);
  std::uniform_real_distribution<double> weight_dist(0.1, 1.9);
  n = dist2(generator);
  n = dist2(generator);
  int offset = 1;
  double x = (lcoords[0] < gcoords[0] / 2)
                 ? (double)lcoords[0]
                 : (double)((gcoords[0] - 1) - lcoords[0]);
  double y = (lcoords[1] < gcoords[1] / 2)
                 ? (double)lcoords[1]
                 : (double)((gcoords[1] - 1) - lcoords[1]);
  double z = (lcoords[2] < gcoords[2] / 2)
                 ? (double)lcoords[2]
                 : (double)((gcoords[2] - 1) - lcoords[2]);
  // offset = (int) ( 15.0 - 0.5 * x * y - (5.0 - x) * (5.0 - y) / 2.0 - 2.0 *
  // std::abs(2.5 - z) ); offset = std::max(offset,0); offset =
  // dist2(generator); offset = 1 + (int)x * (int)y * (int)z;
  offset = 1 + (int)x + (int)y + (int)z;
  n = offset * 64;
  double weight = weight_dist(generator);

  for (int i = 0; i < dimension; ++i) {
    dist[i] = std::uniform_real_distribution<double>(l.at(i), u.at(i));
  }

  double coords[dimension];
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < dimension; ++i) {
      coords[i] = dist.at(i)(generator);
    }
    ALL::Point<double> p(dimension, coords, weight);
    points.push_back(p);
  }
}

/****************************************************************
 *
 * input routine for example files
 *
 * format of MPI I/O files:
 *
 *   n_proc int                     : offset in point data
 *   n_proc int                     : number of points on domain
 *   long + 4*n_point double        : point data
 ****************************************************************/

void read_points(std::vector<ALL::Point<double>> &points, std::vector<double> l,
                 std::vector<double> u, int &n_points, char *filename,
                 int dimension, int rank, MPI_Comm comm) {
  MPI_File file;
  MPI_Barrier(comm);
  int err =
      MPI_File_open(comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
  MPI_Barrier(comm);

  n_points = 0;

  int n_ranks;
  long offset;
  MPI_Comm_size(comm, &n_ranks);

  if (err) {
    if (rank == 0)
      std::cout << "File could not be opened: " << filename << std::endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  // read offset from file
  MPI_File_read_at(file, (MPI_Offset)(rank * sizeof(long)), &offset, 1,
                   MPI_LONG, MPI_STATUS_IGNORE);
  // read number of points from file
  MPI_File_read_at(file,
                   (MPI_Offset)(n_ranks * sizeof(long) + rank * sizeof(int)),
                   &n_points, 1, MPI_INT, MPI_STATUS_IGNORE);

  double values[4];
  long ID;
  long block_size = sizeof(long) + 4 * sizeof(double);
  for (int i = 0; i < n_points; ++i) {
    MPI_File_read_at(file,
                     (MPI_Offset)((offset + i) * block_size +
                                  n_ranks * (sizeof(int) + sizeof(long))),
                     &ID, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_read_at(file,
                     (MPI_Offset)((offset + i) * block_size +
                                  n_ranks * (sizeof(int) + sizeof(long)) +
                                  sizeof(long)),
                     values, 4, MPI_DOUBLE, MPI_STATUS_IGNORE);
    ALL::Point<double> p(dimension, &values[0], values[3], ID);
    points.push_back(p);
  }

  MPI_File_close(&file);

  /*

     std::ifstream ifile;
     ifile.open(filename);

     double coords[3];
     double weight;
     int id;

     n_points = 0;
     char line[256];

     while (ifile.getline(line,256))
     {
     if (ifile.eof()) break;
     std::string str(line);
     std::istringstream istr(str);
     istr >> id >> coords[0] >> coords[1] >> coords[2] >> weight;

  // check if point is within the local domain
  if (
  coords[0] >= l.at(0) && coords[0] < u.at(0) &&
  coords[1] >= l.at(1) && coords[1] < u.at(1) &&
  coords[2] >= l.at(2) && coords[2] < u.at(2)
  )
  {
  ALL::Point<double> p(dimension,coords,weight);
  points.push_back(p);
  n_points++;
  }
  }

  ifile.close();
   */
}

#ifdef ALL_VTK_OUTPUT
// function to create a VTK output of the points in the system
void print_points(std::vector<ALL::Point<double>> plist, int step,
                  ALL::LB_t method, MPI_Comm comm) {
  int rank, n_ranks;
  static bool vtk_init = false;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &n_ranks);
  vtkMPIController *controller;

  // seperate init step required, since VORONOI does not
  // use the MPI initialization of VTK within the library
  if (
#ifdef ALL_VORONOI_ACTIVE
      method == ALL::LB_t::VORONOI ||
#endif
      method == ALL::LB_t::FORCEBASED) {
    if (!vtk_init) {
      controller = vtkMPIController::New();
      controller->Initialize();
      controller->SetGlobalController(controller);
      vtk_init = true;
    }
  }

  auto points = vtkSmartPointer<vtkPoints>::New();

  for (auto p : plist)
    points->InsertNextPoint(p[0], p[1], p[2]);

  auto polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);

  auto weight = vtkSmartPointer<vtkFloatArray>::New();
  weight->SetNumberOfComponents(1);
  weight->SetNumberOfTuples(plist.size());
  weight->SetName("weight");
  for (auto i = 0; i < (int)plist.size(); ++i)
    weight->SetValue(i, plist.at(i).getWeight());

  auto domain = vtkSmartPointer<vtkIntArray>::New();
  domain->SetNumberOfComponents(1);
  domain->SetNumberOfTuples(plist.size());
  domain->SetName("domain");
  for (auto i = 0; i < (int)plist.size(); ++i)
    domain->SetValue(i, rank);

  polydata->GetPointData()->AddArray(weight);
  polydata->GetPointData()->AddArray(domain);

  std::ostringstream ss_local;
  ss_local << "vtk_points/ALL_vtk_points_" << std::setw(7) << std::setfill('0')
           << step << "_" << rank << ".vtp";

  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputData(polydata);
  writer->SetFileName(ss_local.str().c_str());
  writer->SetDataModeToAscii();
  writer->Write();

  std::ostringstream ss_para;
  ss_para << "vtk_points/ALL_vtk_points_" << std::setw(7) << std::setfill('0')
          << step << ".pvtp";

  auto parallel_writer = vtkSmartPointer<vtkXMLPPolyDataWriter>::New();
  parallel_writer->SetFileName(ss_para.str().c_str());
  parallel_writer->SetNumberOfPieces(n_ranks);
  parallel_writer->SetStartPiece(rank);
  parallel_writer->SetEndPiece(rank);
  parallel_writer->SetInputData(polydata);
  parallel_writer->SetDataModeToAscii();
  parallel_writer->Update();
  MPI_Barrier(comm);
  parallel_writer->Write();
}
#endif

int main(int argc, char **argv) {
  try {
    MPI_Init(&argc, &argv);

    const int sys_dim = 3;
    int max_loop = N_LOOP;
    double box_size[sys_dim] = {BOX_SIZE, BOX_SIZE, BOX_SIZE};
    std::vector<double> sys_size(6);
    sys_size.at(0) = 0;
    sys_size.at(1) = BOX_SIZE;
    sys_size.at(2) = 0;
    sys_size.at(3) = BOX_SIZE;
    sys_size.at(4) = 0;
    sys_size.at(5) = BOX_SIZE;

    std::vector<double> minSize(3, 2.0);

    ALL::LB_t chosen_method = ALL::LB_t::STAGGERED;
    bool weighted_points = false;
    char *filename = NULL;
    double gamma = 16.0;
    int output_step = 0;
    int global_dim[sys_dim];
    for (int d = 0; d < sys_dim; ++d)
      global_dim[d] = 0;

    // check arguments
    if (argc < 2) {
      throw ALL::InvalidArgumentException(__FILE__, __func__, __LINE__,
                                         "Wrong number of arguments: ALL_test \
                    [ <int:method> [ <int:max_loop> \
                    [ <double:gamma> [ <int:weighted> \
                    [ <char*:input_file> \
                    [ <double:lx> <double:ly> <double:lz> \
                    [ <int:dx> <int:dy> <int:dz> ]]]]]]]");
    }
    if (argc < 3) {
      // read system dimension from command line
      chosen_method = (ALL::LB_t)atoi(argv[1]);
    } else if (argc < 4) {
      chosen_method = (ALL::LB_t)atoi(argv[1]);
      max_loop = atoi(argv[2]);
    } else if (argc < 5) {
      chosen_method = (ALL::LB_t)atoi(argv[1]);
      max_loop = atoi(argv[2]);
      gamma = atof(argv[3]);
    } else if (argc < 6) {
      chosen_method = (ALL::LB_t)atoi(argv[1]);
      max_loop = atoi(argv[2]);
      gamma = atof(argv[3]);
      weighted_points = atoi(argv[4]) == 1;
    } else if (argc < 7) {
      chosen_method = (ALL::LB_t)atoi(argv[1]);
      max_loop = atoi(argv[2]);
      gamma = atof(argv[3]);
      weighted_points = atoi(argv[4]) == 1;
      filename = argv[5];
    } else if (argc == 9) {
      chosen_method = (ALL::LB_t)atoi(argv[1]);
      max_loop = atoi(argv[2]);
      gamma = atof(argv[3]);
      weighted_points = atoi(argv[4]) == 1;
      filename = argv[5];
      box_size[0] = atof(argv[6]);
      box_size[1] = atof(argv[7]);
      box_size[2] = atof(argv[8]);
      sys_size.at(1) = atof(argv[6]);
      sys_size.at(3) = atof(argv[7]);
      sys_size.at(5) = atof(argv[8]);
    } else if (argc == 12) {
      chosen_method = (ALL::LB_t)atoi(argv[1]);
      max_loop = atoi(argv[2]);
      gamma = atof(argv[3]);
      weighted_points = atoi(argv[4]) == 1;
      filename = argv[5];
      box_size[0] = atof(argv[6]);
      box_size[1] = atof(argv[7]);
      box_size[2] = atof(argv[8]);
      global_dim[0] = atoi(argv[9]);
      global_dim[1] = atoi(argv[10]);
      global_dim[2] = atoi(argv[11]);
      sys_size.at(1) = atof(argv[6]);
      sys_size.at(3) = atof(argv[7]);
      sys_size.at(5) = atof(argv[8]);
    }

    // setup of vector of points on each MPI rank
    std::vector<double> dummy(sys_dim);
    std::vector<ALL::Point<double>> points;
    int max_neighbors, loc_neighbors;

    // setup of cartesian communicator
    int localRank;
    int n_ranks;
    MPI_Comm cart_comm;
    int local_coords[sys_dim];
    int periodicity[sys_dim];
    // domain sizes
    std::vector<double> ds(sys_dim);
    // domain boundaries
    std::vector<double> l(sys_dim);
    std::vector<double> u(sys_dim);

    double min_ratio = 1.01;
    int min_step = -1;

    for (int i = 0; i < sys_dim; ++i) {
      periodicity[i] = 1;
    }

    // get number of total ranks
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

    if (global_dim[0] == 0) {
      // get distribution into number of dimensions
      MPI_Dims_create(n_ranks, sys_dim, global_dim);
    }

    // create cartesian MPI communicator
    MPI_Cart_create(MPI_COMM_WORLD, sys_dim, global_dim, periodicity, 1,
                    &cart_comm);

    // get local coordinates
    MPI_Cart_get(cart_comm, sys_dim, global_dim, periodicity, local_coords);

    // get local rank
    MPI_Cart_rank(cart_comm, local_coords, &localRank);

    int coords[sys_dim];
    MPI_Cart_coords(cart_comm, localRank, sys_dim, coords);

    // print out chosen method
    if (localRank == 0) {
      switch (chosen_method) {
      case ALL::LB_t::TENSOR_MAX:
      case ALL::LB_t::TENSOR:
        std::cout << "chosen method: TENSOR" << std::endl;
        break;
      case ALL::LB_t::STAGGERED:
        std::cout << "chosen method: STAGGERED" << std::endl;
        break;
      case ALL::LB_t::FORCEBASED:
        std::cout << "chosen method: FORCEBASED" << std::endl;
        break;
#ifdef ALL_VORONOI_ACTIVE
      case ALL::LB_t::VORONOI:
        std::cout << "chosen method: VORONOI" << std::endl;
        break;
#endif
      case ALL::LB_t::HISTOGRAM:
        std::cout << "chosen method: HISTOGRAM" << std::endl;
        break;
      default:
        throw ALL::InvalidArgumentException(__FILE__, __func__, __LINE__,
                                           "Invalid method chosen.");
      }
      std::cout << "chosen gamma: " << gamma << std::endl;
    }
    // calculate domain extension
    for (int i = 0; i < sys_dim; ++i) {
      ds.at(i) = box_size[i] / (double)global_dim[i];
      l.at(i) = local_coords[i] * ds.at(i);
      u.at(i) = (1.0 + local_coords[i]) * ds.at(i);
    }

    // generate vertices (equal to outline (for tensor only))
    ALL::Point<double> *tp = new ALL::Point<double>(2);
    std::vector<ALL::Point<double>> *tv = new std::vector<ALL::Point<double>>(2);
    tv->at(0) = *tp;
    delete tv;
    delete tp;

    int nvertices = 2;

    switch (chosen_method) {
    case ALL::LB_t::TENSOR_MAX:
    case ALL::LB_t::TENSOR:
      nvertices = 2;
      break;
    case ALL::LB_t::STAGGERED:
      nvertices = 2;
      break;
    case ALL::LB_t::FORCEBASED:
      nvertices = 8;
      break;
#ifdef ALL_VORONOI_ACTIVE
    case ALL::LB_t::VORONOI:
      // the generator point is stored as 'vertex(0)'
      // the center of work is stored as 'vertex(1)'
      nvertices = 2;
      break;
#endif
    case ALL::LB_t::HISTOGRAM:
      nvertices = 2;
      break;
    default:
      nvertices = 2;
      break;
    }

    ALL::Point<double> lp(l);
    ALL::Point<double> up(u);
    std::vector<ALL::Point<double>> vertices(nvertices, lp);
    std::vector<ALL::Point<double>> new_vertices(nvertices, lp);
    std::vector<ALL::Point<double>> old_vertices(nvertices, lp);

    switch (chosen_method) {
    case ALL::LB_t::TENSOR_MAX:
    case ALL::LB_t::TENSOR:
      vertices.at(0) = lp;
      vertices.at(1) = up;
      break;
    case ALL::LB_t::STAGGERED:
      vertices.at(0) = lp;
      vertices.at(1) = up;
      break;
    case ALL::LB_t::FORCEBASED:
      for (int i = 0; i < nvertices; ++i)
        vertices.at(0) = lp;
      for (int z = 0; z < 2; ++z)
        for (int y = 0; y < 2; ++y)
          for (int x = 0; x < 2; ++x) {
            int vertex = 4 * z + 2 * y + x;
            vertices.at(vertex)[0] =
                vertices.at(vertex)[0] + (double)x * ds.at(0);
            vertices.at(vertex)[1] =
                vertices.at(vertex)[1] + (double)y * ds.at(1);
            vertices.at(vertex)[2] =
                vertices.at(vertex)[2] + (double)z * ds.at(2);
          }
      break;
#ifdef ALL_VORONOI_ACTIVE
    case ALL::LB_t::VORONOI:
      // first generator point is in the center of the orthogonal,
      // equidistant decompositon
      for (int d = 0; d < sys_dim; ++d)
        vertices.at(0)[d] = 0.5 * (lp[d] + up[d]);
      vertices.at(1) = vertices.at(0);
      break;
#endif
    case ALL::LB_t::HISTOGRAM:
      vertices.at(0) = lp;
      vertices.at(1) = up;
      break;
    default:
      throw ALL::InvalidArgumentException(__FILE__, __func__, __LINE__,
                                         "Invalid method chosen.");
      break;
    }

    // generate points on domains
    int n_points = N_GENERATE;
    int max_particles = 1;
    if (localRank == 0)
      std::cout << "creating / reading point data" << std::endl;
    if (argc <= 4) {
      generate_points(points, l, u, coords, global_dim, sys_dim, n_points,
                      localRank);
    } else {
      read_points(points, l, u, n_points, filename, sys_dim, localRank,
                  cart_comm);
    }
    double *transfer;
    double *recv;

    MPI_Allreduce(&n_points, &max_particles, 1, MPI_INT, MPI_MAX, cart_comm);

    max_particles = (int)std::ceil((double)max_particles * 1.5);

    if (localRank == 0)
      std::cout << "maximum number of points on any process: " << max_particles
                << std::endl;
    int max_neig = 27;

    recv = new double[max_neig * (sys_dim + 1) * max_particles];
    transfer = new double[max_neig * (sys_dim + 1) * max_particles];

    // find neighbors on cartesian communicator
    int l_neig[sys_dim];
    int r_neig[sys_dim];
    int self;
    for (int i = 0; i < sys_dim; ++i) {
      MPI_Cart_shift(cart_comm, i, 1, &self, &r_neig[i]);
      MPI_Cart_shift(cart_comm, i, -1, &self, &l_neig[i]);
      if (local_coords[i] == 0)
        l_neig[i] = MPI_PROC_NULL;
      if (local_coords[i] == global_dim[i] - 1)
        r_neig[i] = MPI_PROC_NULL;
    }

    double d_min, d_max, d_ratio;
    double n_local;
    if (!weighted_points) {
      n_local = (double)n_points;
    } else {
      n_local = 0.0;
      for (auto p = points.begin(); p != points.end(); ++p)
        n_local += p->getWeight();
    }
    double n_total;
    MPI_Allreduce(&n_local, &n_total, 1, MPI_DOUBLE, MPI_SUM, cart_comm);
    double avg_work = (double)n_total / (double)n_ranks;
    double n_min, n_max;
    MPI_Allreduce(&n_local, &n_min, 1, MPI_DOUBLE, MPI_MIN, cart_comm);
    MPI_Allreduce(&n_local, &n_max, 1, MPI_DOUBLE, MPI_MAX, cart_comm);
    d_min = n_min / avg_work;
    d_max = n_max / avg_work;
    d_ratio = (d_max - d_min) / (d_max + d_min);
    min_ratio = d_ratio;
    min_step = 0;

    double total_points;
    // get starting number of particles
    MPI_Allreduce(&n_local, &total_points, 1, MPI_DOUBLE, MPI_SUM, cart_comm);

    // output of borders / contents
    if (n_ranks < 216) {
      for (int i = 0; i < n_ranks; ++i) {
        if (localRank == i) {
          std::ofstream of;
          of.open("domain_data.dat", std::ios::out | std::ios::app);
          of << 0 << " " << localRank << " ";

          for (int j = 0; j < sys_dim; ++j) {
            of << " " << local_coords[j] << " " << lp[j] << " " << up[j] << " ";
          }

          of << " " << n_local << " ";

          of << std::endl;
          if (i == n_ranks - 1)
            of << std::endl;
          of.close();
          MPI_Barrier(cart_comm);
        } else
          MPI_Barrier(cart_comm);
      }
    }

#ifdef ALL_VORONOI_ACTIVE
    if (chosen_method == ALL::LB_t::VORONOI) {
      // one-time particle output to voronoi/particles.pov

      for (int i = 0; i < n_ranks; ++i) {
        if (localRank == i) {
          std::ofstream of;
          if (i != 0)
            of.open("voronoi/particles.pov", std::ios::out | std::ios::app);
          else
            of.open("voronoi/particles.pov", std::ios::out);
          int j = 0;
          for (auto it = points.begin(); it != points.end(); ++it) {
            of << "// rank " << localRank << " particle " << j << std::endl;
            of << "sphere{<" << (*it)[0] << "," << (*it)[1] << "," << (*it)[2]
               << ">, p}" << std::endl;
            ++j;
          }
        }
        MPI_Barrier(cart_comm);
      }

      int minpoints = 0;
      double cow_sys[sys_dim + 1];
      double target_point[sys_dim + 1];

      cow_sys[0] = cow_sys[1] = cow_sys[2] = cow_sys[3] = 0.0;

      for (int i = 0; i < n_points; ++i) {
        for (int d = 0; d < sys_dim; ++d)
          cow_sys[d] += points.at(i)[d] * points.at(i).getWeight();
        cow_sys[sys_dim] += points.at(i).getWeight();
      }

      MPI_Allreduce(cow_sys, target_point, sys_dim + 1, MPI_DOUBLE, MPI_SUM,
                    cart_comm);

      for (int d = 0; d < sys_dim; ++d)
        target_point[d] /= target_point[sys_dim];

      std::default_random_engine generator(SEED + localRank - 42 * 23);
      std::uniform_real_distribution<double> dist(-0.25, 0.25);

      // experimental: as a first step try to find more optimal start points
      for (int i_preloop = 0; i_preloop < ALL_VORONOI_PREP_STEPS; ++i_preloop) {
        // get neighbor information
        int nNeighbors = n_ranks;
        std::vector<double> neighbor_vertices(sys_dim * n_ranks);
        std::vector<int> neighbors(n_ranks);

        for (int n = 0; n < n_ranks; ++n)
          neighbors.at(n) = n;

        double local_vertex[sys_dim];
        if (n_points > 0) {
          for (auto &v : local_vertex)
            v = 0.0;
          double local_work = 0.0;

          for (int p = 0; p < n_points; ++p) {
            local_work += points.at(p).getWeight();
            for (int d = 0; d < sys_dim; ++d)
              local_vertex[d] += points.at(p)[d] * points.at(p).getWeight();
          }

          for (int d = 0; d < sys_dim; ++d) {
            if (local_work < 1e-6)
              local_work = 1.0;
            local_vertex[d] /= local_work;
            vertices.at(0)[d] = local_vertex[d];
          }
        } else {
          double diff[3];

          for (int d = 0; d < sys_dim; ++d) {
            // diff[d] = dist(generator);
            diff[d] = target_point[d] - local_vertex[d] + dist(generator);
          }

          double dd = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];

          dd = sqrt(dd);
          if (std::fabs(dd) < 1e-6)
            dd = 1.0;

          for (int d = 0; d < sys_dim; ++d) {
            local_vertex[d] = vertices.at(0)[d] + diff[d] / dd;
            local_vertex[d] = (local_vertex[d] < 0.0)
                                  ? local_vertex[d] + 1.0
                                  : ((local_vertex[d] > box_size[d])
                                         ? local_vertex[d] - box_size[d]
                                         : local_vertex[d]);
            vertices.at(0)[d] = local_vertex[d];
          }
        }
        MPI_Allgather(local_vertex, sys_dim, MPI_DOUBLE,
                      neighbor_vertices.data(), sys_dim, MPI_DOUBLE, cart_comm);

        // compute voronoi cells

        voro::container con_prep(0.0, box_size[0], 0.0, box_size[1], 0.0,
                                 box_size[2], ALL_VORONOI_SUBBLOCKS,
                                 ALL_VORONOI_SUBBLOCKS, ALL_VORONOI_SUBBLOCKS,
                                 false, false, false, 10);
        // add neighbor points first to maintain
        // mapping to neighbors array
        for (auto i = 0; i < nNeighbors; ++i) {
          con_prep.put(i, neighbor_vertices.at(3 * i),
                       neighbor_vertices.at(3 * i + 1),
                       neighbor_vertices.at(3 * i + 2));
        }

        /* OUTPUT: Preparation steps */

        /*
           std::ostringstream ss_local_gp_pov2;
           ss_local_gp_pov2 << "voronoi/prep_points_"
           << std::setw(7) << std::setfill('0')
           << i_preloop << ".pov";
           std::ostringstream ss_local_vc_pov2;
           ss_local_vc_pov2 << "voronoi/prep_cells_"
           << std::setw(7) << std::setfill('0')
           << i_preloop << ".pov";

           if (localRank == 0)
           {
           con_prep.draw_particles_pov(ss_local_gp_pov2.str().c_str());
           con_prep.draw_cells_pov(ss_local_vc_pov2.str().c_str());
           }
         */

        // collect particles that left the domain
        // and determine to which domain they will be transfered

        std::vector<std::vector<double>> remote_particles(nNeighbors);

        for (auto it = points.begin(); it != points.end(); ++it) {
          double x, y, z;
          int pos;

          con_prep.find_voronoi_cell((*it)[0], (*it)[1], (*it)[2], x, y, z,
                                     pos);
          remote_particles.at(pos).push_back((*it)[0]);
          remote_particles.at(pos).push_back((*it)[1]);
          remote_particles.at(pos).push_back((*it)[2]);
          remote_particles.at(pos).push_back(it->getWeight());
        }

        points.clear();
        n_points = 0;
        // exchange number of particles to be exchanged

        int remote_s[nNeighbors];
        int remote_r[nNeighbors];
        MPI_Request request_s[nNeighbors];
        MPI_Request request_r[nNeighbors];
        MPI_Status status_s[nNeighbors];
        MPI_Status status_r[nNeighbors];

        for (auto i = 0; i < nNeighbors; ++i) {
          // send number of values to be send
          remote_s[i] = remote_particles.at(i).size();
          MPI_Isend(&remote_s[i], 1, MPI_INT, neighbors.at(i), 3000, cart_comm,
                    &request_s[i]);
          MPI_Irecv(&remote_r[i], 1, MPI_INT, neighbors.at(i), 3000, cart_comm,
                    &request_r[i]);
        }

        MPI_Waitall(nNeighbors, request_s, status_s);
        MPI_Waitall(nNeighbors, request_r, status_r);

        /*
           for (int i = 0; i < 25; ++i)
           {
           if (localRank == i)
           {
           std::cout << localRank << ": ";
           for (int j = 0; j < nNeighbors; ++j)
           std::cout << " " << neighbors.at(j) << " ";
           std::cout << "| " << nNeighbors << std::endl;
           }
           MPI_Barrier(cart_comm);
           }

           for (int i = 0; i < 25; ++i)
           {
           if (localRank == i)
           {
           std::cout << localRank << ": ";
           for (int j = 0; j < nNeighbors; ++j)
           std::cout << " " << remote_s[j] << " / " << remote_r[j] << " ";
           std::cout << n_points << std::endl;
           }
           MPI_Barrier(cart_comm);
           }
         */

        std::vector<std::vector<double>> received_particles(nNeighbors);
        for (auto i = 0; i < nNeighbors; ++i) {
          if (remote_r[i] > 0) {
            received_particles.at(i).resize(remote_r[i]);
            MPI_Irecv(received_particles.at(i).data(), remote_r[i], MPI_DOUBLE,
                      neighbors.at(i), 4000, cart_comm, &request_r[i]);
          } else
            request_r[i] = MPI_REQUEST_NULL;
          if (remote_s[i] > 0)
            MPI_Isend(remote_particles.at(i).data(), remote_s[i], MPI_DOUBLE,
                      neighbors.at(i), 4000, cart_comm, &request_s[i]);
          else
            request_s[i] = MPI_REQUEST_NULL;
        }

        MPI_Waitall(nNeighbors, request_s, status_s);
        MPI_Waitall(nNeighbors, request_r, status_r);

        for (auto i = 0; i < nNeighbors; ++i) {
          for (auto j = 0; j < remote_r[i]; j += 4) {
            ALL::Point<double> tmp_point(3, received_particles.at(i).data() + j,
                                        received_particles.at(i).at(j + 3));
            points.push_back(tmp_point);
            n_points++;
          }
        }

        /*
           int check;
           MPI_Reduce(&n_points,&check,1,MPI_INT,MPI_SUM,0,cart_comm);
           MPI_Allreduce(&n_points,&minpoints,1,MPI_INT,MPI_MIN,cart_comm);
           if (localRank == 0) std::cout << "minpoints: " << minpoints <<
           std::endl;
         */

        // output of borders / contents
        for (int i = 0; i < n_ranks; ++i) {
          if (localRank == i) {
            std::ofstream of;
            if (!weighted_points)
              of.open("prep_data.dat", std::ios::out | std::ios::app);
            else
              of.open("prep_data_w.dat", std::ios::out | std::ios::app);
            of << (i_preloop + 1) << " " << localRank << " ";

            of << " " << vertices.at(0)[0] << " " << vertices.at(0)[1] << " "
               << vertices.at(0)[2] << " " << n_points << std::endl;

            if (i == n_ranks - 1)
              of << std::endl;
            of.close();
            MPI_Barrier(cart_comm);
          } else
            MPI_Barrier(cart_comm);
        }
      }
    }
#endif
    double limit_efficiency = 0.5;
    // create ALL object
    ALL::ALL<double, double> lb_obj(chosen_method, sys_dim, vertices, gamma);
    std::vector<int> local_coordinates(local_coords, local_coords+sys_dim);
    std::vector<int> global_dimensions(global_dim, global_dim+sys_dim);
    lb_obj.setProcGridParams(local_coordinates, global_dimensions);
    lb_obj.setCommunicator(MPI_COMM_WORLD);
    lb_obj.setup();
    for (int i_loop = 0; i_loop < max_loop; ++i_loop) {
      MPI_Barrier(cart_comm);
      if (d_ratio < limit_efficiency) {
        gamma *= 2.0;
        limit_efficiency /= 2.0;
      }
      if (localRank == 0)
        std::cout << "loop " << i_loop << ": " << std::endl;
      std::vector<double> work;
      std::vector<int> n_bins(3, -1);

      // double histogram_width = ALL_HISTOGRAM_DEFAULT_WIDTH / (double)(i_loop
      // + 1);
      double histogram_width = ALL_HISTOGRAM_DEFAULT_WIDTH / gamma;

      if (!weighted_points) {
        // work -> number of points on domain
#ifdef ALL_VORONOI_ACTIVE
        if (chosen_method == ALL::LB_t::VORONOI) {
          ALL::Point<double> cow(3);
          for (int i = 0; i < 3; ++i)
            cow[i] = 0.0;
          for (auto p : points) {
            cow = cow + (p * (1.0 / n_points));
          }
          if (points.size() != 0)
            vertices.at(1) = cow;
          else
            vertices.at(1) = vertices.at(0);
        }
#endif
        // compute work histogram
        if (chosen_method != ALL::LB_t::HISTOGRAM)
          work = std::vector<double>(1, (double)n_points);
        else {
          double lb(-1.0);
          double ub(-1.0);
          double overlap(0.0);
          int d = 2 - i_loop % 3;
          // compute number of bins in each direction
          lb = std::ceil(lp[d] / histogram_width) * histogram_width;
          ub = std::ceil(up[d] / histogram_width) * histogram_width;
          n_bins.at(d) = (ub - lb) / histogram_width;

          work = std::vector<double>(n_bins.at(d), 0.0);
          // compute histogram of work load
          for (auto p : points) {
            int idx = (int)std::floor(((p[d] - lb) / histogram_width));
            if (idx >= 0) {
              work.at(idx) += 1.0;
            } else
              overlap += 1.0;
          }

          // exchange overlapping workload (histograms might overlap
          // over the domain boundaries
          int rank_left, rank_right;
          MPI_Cart_shift(cart_comm, 0, 1, &rank_left, &rank_right);

          MPI_Request sreq, rreq;
          MPI_Status ssta, rsta;

          double recv_work;

          MPI_Isend(&overlap, 1, MPI_DOUBLE, rank_left, 0, cart_comm, &sreq);
          MPI_Irecv(&recv_work, 1, MPI_DOUBLE, rank_right, 0, cart_comm, &rreq);
          MPI_Wait(&sreq, &ssta);
          MPI_Wait(&rreq, &rsta);

          if (local_coords[d] != global_dim[d] - 1)
            work.at(n_bins.at(d) - 1) += recv_work;
        }
      } else {
        // work -> weighted number of points on domain
        if (chosen_method != ALL::LB_t::HISTOGRAM) {
          work = std::vector<double>(1, 0.0);
          for (auto p = points.begin(); p != points.end(); ++p)
            work.at(0) += p->getWeight();
        } else {
          double lb(-1.0);
          double ub(-1.0);
          double overlap(0.0);
          int d = 2 - i_loop % 3;
          // compute number of bins in each direction
          lb = std::ceil(lp[d] / histogram_width) * histogram_width;
          ub = std::ceil(up[d] / histogram_width) * histogram_width;
          n_bins.at(d) = (ub - lb) / histogram_width;

          work = std::vector<double>(n_bins.at(d), 0.0);
          // compute histogram of work load
          for (auto p : points) {
            int idx = (int)std::floor(((p[d] - lb) / histogram_width));
            if (idx >= 0) {
              work.at(idx) += p.getWeight();
            } else
              overlap += p.getWeight();
          }

          // exchange overlapping workload (histograms might overlap
          // over the domain boundaries
          int rank_left, rank_right;
          MPI_Cart_shift(cart_comm, 0, 1, &rank_left, &rank_right);

          MPI_Request sreq, rreq;
          MPI_Status ssta, rsta;

          double recv_work;

          MPI_Isend(&overlap, 1, MPI_DOUBLE, rank_left, 0, cart_comm, &sreq);
          MPI_Irecv(&recv_work, 1, MPI_DOUBLE, rank_right, 0, cart_comm, &rreq);
          MPI_Wait(&sreq, &ssta);
          MPI_Wait(&rreq, &rsta);

          if (local_coords[d] != global_dim[d] - 1)
            work.at(n_bins.at(d) - 1) += recv_work;
        }
#ifdef ALL_VORONOI_ACTIVE
        if (chosen_method == ALL::LB_t::VORONOI) {
          if (points.size() > 0) {
            ALL::Point<double> cow(3);
            for (int i = 0; i < 3; ++i)
              cow[i] = 0.0;
            for (auto p : points) {
              cow = cow + (p * (1.0 / (work.at(0) * p.getWeight())));
            }
            vertices.at(1) = cow;
          } else
            vertices.at(1) = vertices.at(0);
        }
#endif
      }
#ifdef ALL_DEBUG_ENABLED
      MPI_Barrier(cart_comm);
      if (localRank == 0)
        std::cout << "finished computation of work" << std::endl;
      // lb_obj.setWork((double)n_points);

      if (localRank == 0)
        std::cout << "Setting work..." << std::endl;
#endif
      lb_obj.setWork(work);
#ifdef ALL_DEBUG_ENABLED
      MPI_Barrier(cart_comm);
      if (localRank == 0)
        std::cout << "Setting process grid information..." << std::endl;
#endif
#ifdef ALL_DEBUG_ENABLED
      MPI_Barrier(cart_comm);
      if (localRank == 0)
        std::cout << "Setting communicator..." << std::endl;
#endif
      lb_obj.setCommunicator(MPI_COMM_WORLD);
#ifdef ALL_DEBUG_ENABLED
      MPI_Barrier(cart_comm);
      if (localRank == 0)
        std::cout << "Setting up method..." << std::endl;
#endif
#ifdef ALL_DEBUG_ENABLED
      MPI_Barrier(cart_comm);
      if (localRank == 0)
        std::cout << "Setting method data..." << std::endl;
#endif
      if (chosen_method == ALL::LB_t::HISTOGRAM)
        lb_obj.setMethodData(n_bins.data());

#ifdef ALL_DEBUG_ENABLED
      MPI_Barrier(cart_comm);
      if (localRank == 0)
        std::cout << "Setting system size..." << std::endl;
#endif
      lb_obj.setSysSize(sys_size);
#ifdef ALL_DEBUG_ENABLED
      MPI_Barrier(cart_comm);
      if (localRank == 0)
        std::cout << "Setting domain vertices..." << std::endl;
#endif
      lb_obj.setVertices(vertices);
#ifdef ALL_DEBUG_ENABLED
      MPI_Barrier(cart_comm);
      if (localRank == 0)
        std::cout << "Setting domain minimum size..." << std::endl;
#endif
      lb_obj.setMinDomainSize(minSize);
      // print out number of particles to check communication!
      int n_points_global = 0;
      MPI_Reduce(&n_points, &n_points_global, 1, MPI_INT, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      if (localRank == 0)
        std::cout << "number of particles in step " << i_loop << ": "
                  << n_points_global << std::endl;
#ifdef ALL_DEBUG_ENABLED
      MPI_Barrier(cart_comm);
      if(localRank == 0)
          std::cout << "Providing current LB efficiency" << std::endl;
#endif 
      double LBefficiency = lb_obj.getEfficiency();
      if (localRank == 0)
          std::cout << "current LB efficieny: " << LBefficiency << std::endl;

#ifdef ALL_DEBUG_ENABLED
      MPI_Barrier(cart_comm);
      if (localRank == 0)
        std::cout << "Balancing domains..." << std::endl;
#endif
      lb_obj.balance();
#ifdef ALL_DEBUG_ENABLED
      MPI_Barrier(cart_comm);
      if (localRank == 0)
          std::cout << "Providing estimated new LB efficiency" << std::endl;
#endif
      if (localRank == 0)
      std::cout << "Estimated LB efficency: " 
                << lb_obj.getEstimatedEfficiency()
                << std::endl;

#ifdef ALL_DEBUG_ENABLED
      MPI_Barrier(cart_comm);
      if (localRank == 0)
        std::cout << "Updating vertices..." << std::endl;
#endif
      new_vertices = lb_obj.getVertices();
      old_vertices = vertices;
      vertices = new_vertices;

#ifdef ALL_VORONOI_ACTIVE
      if (chosen_method == ALL::LB_t::VORONOI) {
#ifdef ALL_DEBUG_ENABLED
        MPI_Barrier(cart_comm);
        if (localRank == 0)
          std::cout << "Updating vertices (Voronoi)..." << std::endl;
#endif
        vertices.resize(2);
        vertices.at(1) = vertices.at(0);
      }
#endif
#ifdef ALL_DEBUG_ENABLED
      MPI_Barrier(cart_comm);
      if (localRank == 0)
        std::cout << "Updating upper / lower vertices..." << std::endl;
#endif
      lp = vertices.at(0);
      up = vertices.at(vertices.size() - 1);

#ifdef ALL_DEBUG_ENABLED
      MPI_Barrier(cart_comm);
      if (localRank == 0)
        std::cout << "beginning particles transfer..." << std::endl;
#endif
      switch (chosen_method) {
      case ALL::LB_t::TENSOR_MAX:
      case ALL::LB_t::TENSOR: {
        int n_transfer[2];
        int n_recv[2];
        for (int i = 0; i < sys_dim; ++i) {
          for (int j = 0; j < 2; ++j) {
            n_transfer[j] = 0;
            n_recv[j] = 0;
          }
          for (auto P = points.begin(); P != points.end(); ++P) {
            ALL::Point<double> p = *P;
            if (p[i] < vertices.at(0)[i]) {
              // copy particles that left the domain
              // to the transfer array
              for (int j = 0; j < sys_dim; j++) {
                if (p[j] < 0.0)
                  p[j] = p[j] + box_size[j];
                transfer[n_transfer[0] * (sys_dim + 1) + j] = p[j];
              }
              transfer[n_transfer[0] * (sys_dim + 1) + sys_dim] = p.getWeight();
              n_transfer[0]++;
            } else if (p[i] >= vertices.at(1)[i]) {
              // copy particles that left the domain to the transfer array
              for (int j = 0; j < sys_dim; j++) {
                if (p[j] > box_size[j])
                  p[j] = p[j] - box_size[j];
                transfer[(sys_dim + 1) * max_particles +
                         n_transfer[1] * (sys_dim + 1) + j] = p[j];
              }
              transfer[(sys_dim + 1) * max_particles +
                       n_transfer[1] * (sys_dim + 1) + sys_dim] = p.getWeight();
              n_transfer[1]++;
            }
          }
          for (auto P = points.begin(); P != points.end(); ++P) {
            ALL::Point<double> p = *P;
            if (p[i] < vertices.at(0)[i] || p[i] >= vertices.at(1)[i]) {
              points.erase(P);
              n_points--;
              P--;
            }
          }
          MPI_Status status;

          MPI_Request sreq_r, rreq_r;
          MPI_Request sreq_l, rreq_l;
          MPI_Status sstatus, rstatus;

          MPI_Irecv(&n_recv[1], 1, MPI_INT, r_neig[i], 10, cart_comm, &rreq_r);
          MPI_Irecv(&n_recv[0], 1, MPI_INT, l_neig[i], 20, cart_comm, &rreq_l);

          MPI_Isend(&n_transfer[0], 1, MPI_INT, l_neig[i], 10, cart_comm,
                    &sreq_l);
          MPI_Isend(&n_transfer[1], 1, MPI_INT, r_neig[i], 20, cart_comm,
                    &sreq_r);

          MPI_Wait(&sreq_l, &sstatus);
          MPI_Wait(&sreq_r, &rstatus);

          MPI_Wait(&rreq_l, &sstatus);
          MPI_Wait(&rreq_r, &rstatus);

          // send particles to corresponding neighbor
          MPI_Irecv(&recv[max_particles * (sys_dim + 1)],
                    (sys_dim + 1) * n_recv[1], MPI_DOUBLE, r_neig[i], 30,
                    cart_comm, &rreq_r);
          MPI_Irecv(&recv[0], (sys_dim + 1) * n_recv[0], MPI_DOUBLE, l_neig[i],
                    40, cart_comm, &rreq_l);

          MPI_Isend(&transfer[0], (sys_dim + 1) * n_transfer[0], MPI_DOUBLE,
                    l_neig[i], 30, cart_comm, &sreq_l);
          MPI_Isend(&transfer[max_particles * (sys_dim + 1)],
                    (sys_dim + 1) * n_transfer[1], MPI_DOUBLE, r_neig[i], 40,
                    cart_comm, &sreq_r);

          MPI_Wait(&sreq_l, &sstatus);
          MPI_Wait(&sreq_r, &rstatus);

          MPI_Wait(&rreq_l, &sstatus);
          MPI_Wait(&rreq_r, &rstatus);

          for (int j = 0; j < n_recv[0]; ++j) {
            ALL::Point<double> p(sys_dim, &recv[j * (sys_dim + 1)],
                                recv[j * (sys_dim + 1) + sys_dim]);
            points.push_back(p);
            n_points++;
          }
          for (int j = 0; j < n_recv[1]; ++j) {
            ALL::Point<double> p(
                sys_dim, &recv[(max_particles + j) * (sys_dim + 1)],
                recv[(max_particles + j) * (sys_dim + 1) + sys_dim]);
            points.push_back(p);
            n_points++;
          }
        }
      } break;
      case ALL::LB_t::STAGGERED: {
        int n_transfer[2];
        int n_recv[2 * MAX_NEIG];
        int offset_neig[2];
        std::vector<int> neighbors = lb_obj.getNeighbors();
        int *nNeighbors;
        nNeighbors = neighbors.data();

        offset_neig[0] = 0;
        offset_neig[1] = nNeighbors[0];

        loc_neighbors = 0;
        for (int n = 0; n < 6; ++n)
          loc_neighbors += nNeighbors[n];

        MPI_Allreduce(&loc_neighbors, &max_neighbors, 1, MPI_INT, MPI_MAX,
                      cart_comm);

        for (int i = 0; i < sys_dim; ++i) {

          if (nNeighbors[i] + nNeighbors[i + 1] > max_neig) {
            std::cout << ">>> resizing transfer buffers from " << max_neig
                      << " neighbors ";
            max_neig = (int)(1.5 * (double)(nNeighbors[i] + nNeighbors[i + 1]));
            std::cout << "to " << max_neig << " neighbors on process "
                      << localRank << std::endl;
            delete[] transfer;
            delete[] recv;
            recv = new double[max_neig * (sys_dim + 1) * max_particles];
            transfer = new double[max_neig * (sys_dim + 1) * max_particles];
          }
          // MPI_Barrier(cart_comm);
          // if (localRank == 0)
          //  std::cout << "filling buffers with remote points, dim "
          //            << i << std::endl;

          for (int j = 0; j < 2; ++j) {
            n_transfer[j] = 0;
          }

          for (auto P = points.begin(); P != points.end(); ++P) {
            ALL::Point<double> p = *P;
            if (p[i] < vertices.at(0)[i]) {
              // periodic correction
              if (p[i] < 0.0)
                p[i] += box_size[i];
              // copy particles that left the domain
              // to the transfer array
              for (int d = 0; d < sys_dim; ++d)
                transfer[n_transfer[0] * (sys_dim + 1) + d] = p[d];
              transfer[n_transfer[0] * (sys_dim + 1) + sys_dim] = p.getWeight();
              n_transfer[0]++;
              if (n_transfer[0] > max_particles) {
                std::stringstream ss;
                ss << "Trying to send more particles than buffer size allows! "
                   << " n_transfer: " << n_transfer[0]
                   << " max_particles: " << max_particles;
                throw ALL::InvalidArgumentException(__FILE__, __func__, __LINE__,
                                                   ss.str().c_str());
              }
            } else if (p[i] >= vertices.at(1)[i]) {
              // periodic correction
              if (p[i] >= box_size[i])
                p[i] -= box_size[i];
              // copy particles that left the domain
              // to the transfer array
              for (int d = 0; d < sys_dim; ++d) {
                transfer[(sys_dim + 1) * max_particles +
                         n_transfer[1] * (sys_dim + 1) + d] = p[d];
              }
              transfer[(sys_dim + 1) * max_particles +
                       n_transfer[1] * (sys_dim + 1) + sys_dim] = p.getWeight();
              n_transfer[1]++;
              if (n_transfer[1] > max_particles) {
                std::stringstream ss;
                ss << "Trying to send more particles than buffer size allows! "
                   << " n_transfer: " << n_transfer[0]
                   << " max_particles: " << max_particles;
                throw ALL::InvalidArgumentException(__FILE__, __func__, __LINE__,
                                                   ss.str().c_str());
              }
            }
          }

          // MPI_Barrier(cart_comm);
          // if (localRank == 0)
          //  std::cout << "cleaning point vector, dim " << i << std::endl;

          auto rem_it = std::remove_if(
              points.begin(), points.end(), [vertices, i](ALL::Point<double> p) {
                return p[i] < vertices.at(0)[i] || p[i] >= vertices.at(1)[i];
              });
          points.erase(rem_it, points.end());
          n_points = points.size();

          // MPI_Barrier(cart_comm);
          // if (localRank == 0)
          //  std::cout << "exchanging number of sent points, dim "
          //            << i << std::endl;

          // TODO -> better estimation than max 1024 neighbors ...
          MPI_Request sreq_r[MAX_NEIG], rreq_r[MAX_NEIG];
          MPI_Request sreq_l[MAX_NEIG], rreq_l[MAX_NEIG];
          MPI_Status lsstatus[MAX_NEIG], rsstatus[MAX_NEIG];
          MPI_Status lrstatus[MAX_NEIG], rrstatus[MAX_NEIG];

          for (int d = 0; d < nNeighbors[2 * i]; ++d) {
            MPI_Irecv(&n_recv[d], 1, MPI_INT, neighbors.at(offset_neig[0] + d),
                      20, cart_comm, &rreq_l[d]);
            MPI_Isend(&n_transfer[0], 1, MPI_INT,
                      neighbors.at(offset_neig[0] + d), 10, cart_comm,
                      &sreq_l[d]);
          }
          for (int d = 0; d < nNeighbors[2 * i + 1]; ++d) {
            MPI_Irecv(&n_recv[MAX_NEIG + d], 1, MPI_INT,
                      neighbors.at(offset_neig[1] + d), 10, cart_comm,
                      &rreq_r[d]);
            MPI_Isend(&n_transfer[1], 1, MPI_INT,
                      neighbors.at(offset_neig[1] + d), 20, cart_comm,
                      &sreq_r[d]);
          }

          if (nNeighbors[2 * i] > MAX_NEIG) {
            throw ALL::InternalErrorException(
                __FILE__, __func__, __LINE__,
                "Number of neighbors higher than expected \
                                        maximum number of neighbors.");
          }

          if (nNeighbors[2 * i + 1] > MAX_NEIG) {
            throw ALL::InvalidArgumentException(
                __FILE__, __func__, __LINE__,
                "Number of neighbors higher than expected \
                                        maximum number of neighbors.");
          }

          MPI_Waitall(nNeighbors[2 * i], sreq_l, lsstatus);
          MPI_Waitall(nNeighbors[2 * i + 1], sreq_r, lrstatus);

          MPI_Waitall(nNeighbors[2 * i], rreq_l, rsstatus);
          MPI_Waitall(nNeighbors[2 * i + 1], rreq_r, rrstatus);

          // MPI_Barrier(cart_comm);
          // if (localRank == 0)
          //  std::cout << "exchanging point data, dim " << i << std::endl;

          int offset_l = 0;
          // send particles to corresponding neighbor
          for (int d = 0; d < nNeighbors[2 * i]; ++d) {
            MPI_Irecv(&recv[offset_l * (sys_dim + 1)],
                      (sys_dim + 1) * n_recv[d], MPI_DOUBLE,
                      neighbors.at(offset_neig[0] + d), 40, cart_comm,
                      &rreq_l[d]);
            MPI_Isend(&transfer[0], (sys_dim + 1) * n_transfer[0], MPI_DOUBLE,
                      neighbors.at(offset_neig[0] + d), 30, cart_comm,
                      &sreq_l[d]);
            offset_l += n_recv[d];
          }

          int offset_r = 0;
          for (int d = 0; d < nNeighbors[2 * i + 1]; ++d) {
            MPI_Irecv(
                &recv[max_particles * (sys_dim + 1) + offset_r * (sys_dim + 1)],
                (sys_dim + 1) * n_recv[MAX_NEIG + d], MPI_DOUBLE,
                neighbors.at(offset_neig[1] + d), 30, cart_comm, &rreq_r[d]);
            MPI_Isend(&transfer[max_particles * (sys_dim + 1)],
                      (sys_dim + 1) * n_transfer[1], MPI_DOUBLE,
                      neighbors.at(offset_neig[1] + d), 40, cart_comm,
                      &sreq_r[d]);
            offset_r += n_recv[MAX_NEIG + d];
          }

          MPI_Waitall(nNeighbors[2 * i], sreq_l, lsstatus);
          MPI_Waitall(nNeighbors[2 * i + 1], sreq_r, lrstatus);

          MPI_Waitall(nNeighbors[2 * i], rreq_l, rsstatus);
          MPI_Waitall(nNeighbors[2 * i + 1], rreq_r, rrstatus);

          // MPI_Barrier(cart_comm);
          // if (localRank == 0)
          //  std::cout << "adding received point data, dim "
          //            << i << std::endl;
          for (int j = 0; j < offset_l; ++j) {
            ALL::Point<double> p(sys_dim, &recv[j * (sys_dim + 1)],
                                recv[j * (sys_dim + 1) + sys_dim]);
            // if the particle does not belong to the local domain skip it
            bool outside = false;
            for (int d = 0; d < i; ++d)
              outside = outside ||
                        (p[d] < vertices.at(0)[d] || p[d] >= vertices.at(1)[d]);
            if (outside)
              continue;
            points.push_back(p);
            n_points++;
          }
          for (int j = 0; j < offset_r; ++j) {
            ALL::Point<double> p(
                sys_dim, &recv[(max_particles + j) * (sys_dim + 1)],
                recv[(max_particles + j) * (sys_dim + 1) + sys_dim]);
            // if the particle does not belong to the local domain skip it
            bool outside = false;
            for (int d = 0; d < i; ++d)
              outside = outside ||
                        (p[d] < vertices.at(0)[d] || p[d] >= vertices.at(1)[d]);
            if (outside)
              continue;
            points.push_back(p);
            n_points++;
          }

          offset_neig[0] = offset_neig[1] + nNeighbors[2 * i + 1];
          if (i < sys_dim - 1)
            offset_neig[1] = offset_neig[0] + nNeighbors[2 * (i + 1)];
        }
      } break;
      case ALL::LB_t::FORCEBASED: {

        // array of all vertices of all neighbors (and self)
        double comm_vertices[27 * 8 * 3];

        std::vector<int> neighbors(27);
        int neig = 0;
        int coords[3];
        bool exists[27];

        // get list of neighbors
        for (int z = -1; z <= 1; ++z) {
          coords[2] = local_coords[2] + z;
          for (int y = -1; y <= 1; ++y) {
            coords[1] = local_coords[1] + y;
            for (int x = -1; x <= 1; ++x) {
              coords[0] = local_coords[0] + x;
              if (coords[0] >= 0 && coords[0] < global_dim[0] &&
                  coords[1] >= 0 && coords[1] < global_dim[1] &&
                  coords[2] >= 0 && coords[2] < global_dim[2]) {
                MPI_Cart_rank(cart_comm, coords, &neighbors.at(neig));
                exists[neig] = true;
              } else {
                neighbors.at(neig) = MPI_PROC_NULL;
                exists[neig] = false;
              }
              ++neig;
            }
          }
        }
        // no communication to self
        neighbors.at(13) = MPI_PROC_NULL;

        for (int i = 0; i < 27 * 8 * 3; ++i)
          comm_vertices[i] = -1.0;

        // copy local indices to position 13
        for (int i = 0; i < 8; ++i) {
          for (int d = 0; d < 3; ++d) {
            comm_vertices[13 * 8 * 3 + 3 * i + d] = vertices.at(i)[d];
          }
        }

        int offset = 0;
        int comm_size = 0;
        MPI_Request request[54];
        MPI_Status status[54];

        for (int i = 0; i < 54; ++i)
          request[i] = MPI_REQUEST_NULL;

        // send local indices to neighbors
        for (int d = 0; d < 3; ++d) {
          // offset for the buffer
          comm_size = (int)std::pow(3, d);

          // local coords
          coords[0] = local_coords[0];
          coords[1] = local_coords[1];
          coords[2] = local_coords[2];

          int low_neig;
          int up_neig;

          coords[d] -= 1;
          if (coords[d] >= 0)
            MPI_Cart_rank(cart_comm, coords, &low_neig);
          else
            low_neig = MPI_PROC_NULL;
          coords[d] += 2;
          if (coords[d] < global_dim[d])
            MPI_Cart_rank(cart_comm, coords, &up_neig);
          else
            up_neig = MPI_PROC_NULL;

          MPI_Isend(&comm_vertices[(13 - offset) * 24], comm_size * 24,
                    MPI_DOUBLE, low_neig, 1010, cart_comm, &request[0]);
          MPI_Isend(&comm_vertices[(13 - offset) * 24], comm_size * 24,
                    MPI_DOUBLE, up_neig, 2010, cart_comm, &request[1]);

          // increase offset here, as the previous offset gives
          // the start address of the whole dimension
          // data to be transfered

          offset += comm_size;

          MPI_Irecv(&comm_vertices[(13 - offset) * 24], comm_size * 24,
                    MPI_DOUBLE, low_neig, 2010, cart_comm, &request[28]);
          MPI_Irecv(&comm_vertices[(14 + offset - comm_size) * 24],
                    comm_size * 24, MPI_DOUBLE, up_neig, 1010, cart_comm,
                    &request[27]);

          MPI_Waitall(54, request, status);
        }

#ifdef ALL_VTK_OUTPUT
#ifdef ALL_VTK_FORCE_SORT
        // creating an unstructured grid for local domain and neighbors
        auto vtkpoints = vtkSmartPointer<vtkPoints>::New();
        auto unstructuredGrid = vtkSmartPointer<vtkForceBasedGrid>::New();
        for (int i = 0; i < 27; ++i) {
          for (int v = 0; v < 8; ++v) {
            vtkpoints->InsertNextPoint(comm_vertices[i * 24 + v * 3],
                                       comm_vertices[i * 24 + v * 3 + 1],
                                       comm_vertices[i * 24 + v * 3 + 2]);
          }
        }
        unstructuredGrid->SetPoints(vtkpoints);

        auto work = vtkSmartPointer<vtkFloatArray>::New();
        work->SetNumberOfComponents(1);
        work->SetNumberOfTuples(27);
        work->SetName("Cell");

        for (int n = 0; n < 27; ++n) {
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
          work->SetValue(n, (double)n);
        }
        unstructuredGrid->GetCellData()->AddArray(work);

        /* Debug output: print local cell and neighbors */
        /*
           if (localRank == 26)
           {
           auto writer = vtkSmartPointer<vtkXMLForceBasedGridWriter>::New();
           writer->SetInputData(unstructuredGrid);
           writer->SetFileName("test.vtu");
           writer->SetDataModeToAscii();
        //writer->SetDataModeToBinary();
        writer->Write();
        }
         */

        MPI_Allreduce(&n_points, &check_np, 1, MPI_INT, MPI_SUM, cart_comm);

        auto locator = vtkSmartPointer<vtkCellLocator>::New();
        locator->SetDataSet(unstructuredGrid);
        locator->BuildLocator();
#else
        auto in_triangle = [=](double nv[27 * 8 * 3], int n, int A, int B,
                               int C, double x, double y, double z) {
          ALL::Point<double> a(3, nv + n * 24 + A * 8);
          ALL::Point<double> b(3, nv + n * 24 + B * 8);
          ALL::Point<double> c(3, nv + n * 24 + C * 8 + 0);
          ALL::Point<double> p(3);
          p[0] = x;
          p[1] = y;
          p[2] = z;
          return false;
          // return p.same_side_line(a, b, c) &&
          //       p.same_side_line(b, c, a) &&
          //       p.same_side_line(c, a, b);
        };

        // check if (x,y,z) is on the same side as A in
        // regard to the plane spanned by (B,C,D)
        auto same_side_plane = [=](double nv[27 * 8 * 3], int n, int A, int B,
                                   int C, int D, double x, double y, double z) {
          // compute difference vectors required:
          // C-B
          double CBx = nv[n * 24 + C * 8 + 0] - nv[n * 24 + B * 8 + 0];
          double CBy = nv[n * 24 + C * 8 + 1] - nv[n * 24 + B * 8 + 1];
          double CBz = nv[n * 24 + C * 8 + 2] - nv[n * 24 + B * 8 + 2];
          // D-B
          double DBx = nv[n * 24 + D * 8 + 0] - nv[n * 24 + B * 8 + 0];
          double DBy = nv[n * 24 + D * 8 + 1] - nv[n * 24 + B * 8 + 1];
          double DBz = nv[n * 24 + D * 8 + 2] - nv[n * 24 + B * 8 + 2];
          // A-B
          double ABx = nv[n * 24 + A * 8 + 0] - nv[n * 24 + B * 8 + 0];
          double ABy = nv[n * 24 + A * 8 + 1] - nv[n * 24 + B * 8 + 1];
          double ABz = nv[n * 24 + A * 8 + 2] - nv[n * 24 + B * 8 + 2];
          // P-B
          double PBx = x - nv[n * 24 + B * 8 + 0];
          double PBy = y - nv[n * 24 + B * 8 + 1];
          double PBz = z - nv[n * 24 + B * 8 + 2];

          // compute normal vector of plane
          // n = (C-B) x (D-B)
          double nx = CBy * DBz - CBz * DBy;
          double ny = CBz * DBx - CBx * DBz;
          double nz = CBx * DBy - CBy * DBx;

          // compute dot product <A-B,n>
          double ABn = ABx * nx + ABy * ny + ABz * nz;

          // compute dot product <P-B,n>
          double PBn = PBx * nx + PBy * ny + PBz * nz;

          /*
          if (
                  x <= TRACK_X+0.01 && x >= TRACK_X &&
                  y <= TRACK_Y+0.01 && y >= TRACK_Y &&
                  z <= TRACK_Z+0.01 && z >= TRACK_Z
             )
          {
              std::cerr << "Tracked Particle: "
                        << localRank << " ( "
                        << x << " , "
                        << y << " , "
                        << z << " ) "
                        << n << ", "
                        << A << ", "
                        << B << ", "
                        << C << ", "
                        << D << ") "
                        << ( ( std::abs(PBn) < 1e-12 ) && in_triangle(nv, n, B,
          C, D, x, y, z) ) << " "
                        << ( ( ABn * PBn ) > 0 )
                        << std::endl;
          }
          */

          if (std::abs(PBn) < 1e-12) {
            return in_triangle(nv, n, B, C, D, x, y, z);
          } else
            // check if sign matches
            return ((ABn * PBn) > 0);
        };

        // check if (x,y,z) is in the tetrahedron described
        // by the points within the array t
        auto in_tetraeder = [=](double nv[27 * 8 * 3], int n, int t[4],
                                double x, double y, double z) {
          return same_side_plane(nv, n, t[0], t[1], t[2], t[3], x, y, z) &&
                 same_side_plane(nv, n, t[3], t[0], t[1], t[2], x, y, z) &&
                 same_side_plane(nv, n, t[2], t[3], t[0], t[1], x, y, z) &&
                 same_side_plane(nv, n, t[1], t[2], t[3], t[0], x, y, z);
        };

        auto containing_neighbor = [=](double nv[27 * 8 * 3], bool e[27],
                                       double x, double y, double z) {
          // definition of tetrahedra
          /*
          int tetrahedra[8][4] =
          {
              { 0, 1, 2, 4 },
              { 1, 0, 3, 5 },
              { 2, 3, 0, 6 },
              { 3, 2, 1, 7 },
              { 4, 5, 6, 0 },
              { 5, 4, 7, 1 },
              { 6, 7, 4, 2 },
              { 7, 6, 5, 3 }
          };
          */
          int tetrahedra[5][4] = {{1, 2, 4, 7},
                                  {0, 1, 2, 4},
                                  {1, 2, 3, 7},
                                  {1, 4, 5, 7},
                                  {2, 4, 6, 7}};

          // check if the particle is in the local
          // tetrahedron, if so keep it
          // (needed to decide what about particles
          //  on surfaces)
          for (int t = 0; t < 5; ++t) {
            if (in_tetraeder(nv, 13, tetrahedra[t], x, y, z))
              return 13;
          }

          // loop over neighbors
          for (int n = 0; n < 27; ++n) {
            if (!e[n] || n == 13)
              continue;
            // loop over tetrahedra
            for (int t = 0; t < 5; ++t) {
              if (in_tetraeder(nv, n, tetrahedra[t], x, y, z))
                return n;
            }
          }

          return -1;
        };

#endif

        int n_send[27];
        int n_recv[27];
        for (int i = 0; i < 27; ++i) {
          n_send[i] = 0;
          n_recv[i] = 0;
        }

        int check_np = 0;
        int check_np_new = 0;

        for (auto P = points.begin(); P != points.end(); ++P) {
          ALL::Point<double> p = *P;
          double pcoords[3];
          double pccoords[3];
          int subId;
          double weights[4];
          for (int d = 0; d < 3; ++d)
            pcoords[d] = p[d];
#ifdef ALL_VTK_FORCE_SORT
          vtkIdType cellId = locator->FindCell(pcoords);
#else
          int cellId = containing_neighbor(comm_vertices, exists, pcoords[0],
                                           pcoords[1], pcoords[2]);
#endif
          /*
                                              vtkIdType cellId =
             unstructuredGrid->FindCell(pcoords, NULL, 0, 1e-6, subId, pccoords,
                                              weights);
          */

          // if the particle is in a valid neighboring cell
          if (cellId >= 0 && cellId != 13 && cellId <= 26) {
            if (n_send[cellId] == max_particles) {
              throw ALL::InvalidArgumentException(
                  __FILE__, __func__, __LINE__,
                  "Trying to send more particles than \
                                            buffer size allows!");
            }
            for (int d = 0; d < 3; ++d) {
              transfer[cellId * (sys_dim + 1) * max_particles +
                       n_send[cellId] * (sys_dim + 1) + d] = p[d];
            }
            transfer[cellId * (sys_dim + 1) * max_particles +
                     n_send[cellId] * (sys_dim + 1) + sys_dim] = p.getWeight();
            n_send[cellId]++;
            points.erase(P);
            --P;
            --n_points;
          }
        }

        for (int n = 0; n < 27; ++n) {
          MPI_Isend(&n_send[n], 1, MPI_INT, neighbors.at(n), 1020, cart_comm,
                    &request[n]);
          MPI_Irecv(&n_recv[n], 1, MPI_INT, neighbors.at(n), 1020, cart_comm,
                    &request[27 + n]);
        }
        MPI_Waitall(54, request, status);

        for (int n = 0; n < 27; ++n) {
          MPI_Isend(&transfer[n * (sys_dim + 1) * max_particles],
                    (sys_dim + 1) * n_send[n], MPI_DOUBLE, neighbors.at(n),
                    1030, cart_comm, &request[n]);
          MPI_Irecv(&recv[n * (sys_dim + 1) * max_particles],
                    (sys_dim + 1) * n_recv[n], MPI_DOUBLE, neighbors.at(n),
                    1030, cart_comm, &request[27 + n]);
        }
        MPI_Waitall(54, request, status);

        for (int n = 0; n < 27; ++n) {
          if (exists[n]) {
            for (int i = 0; i < n_recv[n]; ++i) {
              ALL::Point<double> p(
                  sys_dim,
                  &recv[n * (sys_dim + 1) * max_particles + i * (sys_dim + 1)],
                  recv[n * (sys_dim + 1) * max_particles + i * (sys_dim + 1) +
                       sys_dim]);
              points.push_back(p);
              ++n_points;
            }
          }
        }

#else
        if (localRank == 0)
          std::cout << "Currently no FORCEBASED test without VTK!"
                    << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
      } break;
#ifdef ALL_VORONOI_ACTIVE
      case ALL::LB_t::VORONOI: {
        // get neighbor information
        std::vector<double> neighbor_vertices = lb_obj.getNeighborVertices();
        std::vector<int> neighbors = lb_obj.getNeighbors();
        int nNeighbors = neighbors.size();

        // compute voronoi cells

        voro::container con_copy(0.0, box_size[0], 0.0, box_size[1], 0.0,
                                 box_size[2], ALL_VORONOI_SUBBLOCKS,
                                 ALL_VORONOI_SUBBLOCKS, ALL_VORONOI_SUBBLOCKS,
                                 false, false, false, 10);
        // add neighbor points first to maintain
        // mapping to neighbors array
        for (auto i = 0; i < nNeighbors; ++i) {
          con_copy.put(i, neighbor_vertices.at(3 * i),
                       neighbor_vertices.at(3 * i + 1),
                       neighbor_vertices.at(3 * i + 2));
        }

        // add local vertex to complete map of cells
        con_copy.put(nNeighbors, vertices.at(0)[0], vertices.at(0)[1],
                     vertices.at(0)[2]);

        /*
                              std::ostringstream ss_local_gp_pov2;
                              ss_local_gp_pov2 << "voronoi/copy_points_"
                              << std::setw(7) << std::setfill('0')
                              << i_loop << ".pov";
                              std::ostringstream ss_local_vc_pov2;
                              ss_local_vc_pov2 << "voronoi/copy_cells_"
                              << std::setw(7) << std::setfill('0')
                              << i_loop << ".pov";

                              if (localRank == 0)
                              {
                              con_copy.draw_particles_pov(ss_local_gp_pov2.str().c_str());
                              con_copy.draw_cells_pov(ss_local_vc_pov2.str().c_str());
                              }
         */
        // collect particles that left the domain
        // and determine to which domain they will be transfered

        std::vector<std::vector<double>> remote_particles(nNeighbors);

        for (auto P = points.begin(); P != points.end(); ++P) {
          ALL::Point<double> p = *P;
          double x, y, z;
          int pos;

          con_copy.find_voronoi_cell(p[0], p[1], p[2], x, y, z, pos);
          if (pos < nNeighbors) {
            remote_particles.at(pos).push_back(p[0]);
            remote_particles.at(pos).push_back(p[1]);
            remote_particles.at(pos).push_back(p[2]);
            remote_particles.at(pos).push_back(p.getWeight());
            p[0] = -1000.0;
          }
        }

        for (auto P = points.begin(); P != points.end(); ++P) {
          ALL::Point<double> p = *P;
          if (p[0] < -500.0) {
            points.erase(P);
            n_points--;
            P--;
          }
        }

        // exchange number of particles to be exchanged

        int remote_s[nNeighbors];
        int remote_r[nNeighbors];
        MPI_Request request_s[nNeighbors];
        MPI_Request request_r[nNeighbors];
        MPI_Status status_s[nNeighbors];
        MPI_Status status_r[nNeighbors];

        for (auto i = 0; i < nNeighbors; ++i) {
          // send number of values to be send
          remote_s[i] = remote_particles.at(i).size();
          MPI_Isend(&remote_s[i], 1, MPI_INT, neighbors.at(i), 3000, cart_comm,
                    &request_s[i]);
          MPI_Irecv(&remote_r[i], 1, MPI_INT, neighbors.at(i), 3000, cart_comm,
                    &request_r[i]);
        }

        MPI_Waitall(nNeighbors, request_s, status_s);
        MPI_Waitall(nNeighbors, request_r, status_r);

        /*
           for (int i = 0; i < 25; ++i)
           {
           if (localRank == i)
           {
           std::cout << localRank << ": ";
           for (int j = 0; j < nNeighbors; ++j)
           std::cout << " " << neighbors.at(j) << " ";
           std::cout << "| " << nNeighbors << std::endl;
           }
           MPI_Barrier(cart_comm);
           }

           for (int i = 0; i < 25; ++i)
           {
           if (localRank == i)
           {
           std::cout << localRank << ": ";
           for (int j = 0; j < nNeighbors; ++j)
           std::cout << " " << remote_s[j]
           << " / " << remote_r[j] << " ";
           std::cout << n_points << std::endl;
           }
           MPI_Barrier(cart_comm);
           }
         */

        std::vector<std::vector<double>> received_particles(nNeighbors);
        for (auto i = 0; i < nNeighbors; ++i) {
          if (remote_r[i] > 0) {
            received_particles.at(i).resize(remote_r[i]);
            MPI_Irecv(received_particles.at(i).data(), remote_r[i], MPI_DOUBLE,
                      neighbors.at(i), 4000, cart_comm, &request_r[i]);
          } else
            request_r[i] = MPI_REQUEST_NULL;
          if (remote_s[i] > 0)
            MPI_Isend(remote_particles.at(i).data(), remote_s[i], MPI_DOUBLE,
                      neighbors.at(i), 4000, cart_comm, &request_s[i]);
          else
            request_s[i] = MPI_REQUEST_NULL;
        }

        MPI_Waitall(nNeighbors, request_s, status_s);
        MPI_Waitall(nNeighbors, request_r, status_r);

        for (auto i = 0; i < nNeighbors; ++i) {
          for (auto j = 0; j < remote_r[i]; j += 4) {
            ALL::Point<double> tmp_point(3, received_particles.at(i).data() + j,
                                        received_particles.at(i).at(j + 3));
            points.push_back(tmp_point);
            n_points++;
          }
        }

        // exchange particles
      } break;
#endif
      case ALL::LB_t::HISTOGRAM: {
        // determine current dimension
        int curr_dim = 2 - (i_loop % 3);

        MPI_Comm comm_col;
        // create temporary communicator to exchange borders
        switch (curr_dim) {
        case 0: {
          // x-plane
          MPI_Comm_split(cart_comm,
                         local_coords[1] + local_coords[2] * global_dim[1],
                         local_coords[0], &comm_col);
          break;
        }
        case 1: {
          // y-plane
          MPI_Comm_split(cart_comm,
                         local_coords[0] + local_coords[2] * global_dim[0],
                         local_coords[1], &comm_col);
          break;
        }
        case 2: {
          // z-plane
          MPI_Comm_split(cart_comm,
                         local_coords[0] + local_coords[1] * global_dim[0],
                         local_coords[2], &comm_col);
          break;
        }
        }

        // vector to collect the borders into

        std::vector<double> borders(2 * global_dim[curr_dim]);
        std::vector<double> old_borders(2 * global_dim[curr_dim]);
        std::vector<double> local_borders(2);
        std::vector<double> old_border(2);

        local_borders.at(0) = vertices.at(0)[curr_dim];
        local_borders.at(1) = vertices.at(1)[curr_dim];
        old_border.at(0) = old_vertices.at(0)[curr_dim];
        old_border.at(1) = old_vertices.at(1)[curr_dim];

        int size, rank;
        MPI_Comm_rank(comm_col, &rank);
        MPI_Comm_size(comm_col, &size);

        // collect borders
        MPI_Allgather(local_borders.data(), 2, MPI_DOUBLE, borders.data(), 2,
                      MPI_DOUBLE, comm_col);

        // collect borders
        MPI_Allgather(old_border.data(), 2, MPI_DOUBLE, old_borders.data(), 2,
                      MPI_DOUBLE, comm_col);

        // compare old domains with new domains
        std::vector<int> send_neig;
        std::vector<int> recv_neig;

        for (int n = 0; n < global_dim[curr_dim]; ++n) {
          /*
          if (localRank == 13)
              std::cout << "DEBUG: " << n << " " <<
                  old_border.at(0) << " " <<
                  old_border.at(1) << " / " <<
                  borders.at(2*n) << " " <<
                  borders.at(2*n+1) << " / " <<
                  (borders.at(2*n) <= old_border.at(0)) << " " <<
                  (borders.at(2*n+1) > old_border.at(0)) << " " <<
                  (borders.at(2*n) <= old_border.at(1)) << " " <<
                  (borders.at(2*n+1) >= old_border.at(1)) << std::endl;
          */
          if ((borders.at(2 * n) <= old_border.at(0) &&
               borders.at(2 * n + 1) >= old_border.at(0)) ||
              (borders.at(2 * n) <= old_border.at(1) &&
               borders.at(2 * n + 1) >= old_border.at(1)) ||
              (borders.at(2 * n) > old_border.at(0) &&
               borders.at(2 * n + 1) < old_border.at(1)))
            send_neig.push_back(n);
          if ((old_borders.at(2 * n) <= local_borders.at(0) &&
               old_borders.at(2 * n + 1) >= local_borders.at(0)) ||
              (old_borders.at(2 * n) <= local_borders.at(1) &&
               old_borders.at(2 * n + 1) >= local_borders.at(1)) ||
              (old_borders.at(2 * n) > local_borders.at(0) &&
               old_borders.at(2 * n + 1) < local_borders.at(1)))
            recv_neig.push_back(n);
        }

        // vectors to transfer and received points
        std::vector<std::vector<double>> send_vec(send_neig.size());

        // all points are sorted into the correct vectors
        for (auto p : points) {
          for (int n = 0; n < (int)send_neig.size(); ++n) {
            int neig_id = send_neig.at(n);
            /*
            if (localRank == 1)
            {
                std::cout << p[curr_dim] << " "
                          << borders.at(2*neig_id) << " "
                          << borders.at(2*neig_id+1) << std::endl;
            }
            */
            if ((borders.at(2 * neig_id) <= p[curr_dim] &&
                 borders.at(2 * neig_id + 1) > p[curr_dim])) {
              send_vec.at(n).push_back(p[0]);
              send_vec.at(n).push_back(p[1]);
              send_vec.at(n).push_back(p[2]);
              send_vec.at(n).push_back(p.getWeight());
              send_vec.at(n).push_back((double)(p.get_id()));
            }
          }
        }
        // clear old point vector
        points.clear();
        n_points = 0;

        std::vector<int> n_send(send_neig.size());
        for (int n = 0; n < (int)send_neig.size(); ++n)
          n_send.at(n) = send_vec.at(n).size();

        // communicate number of particles to be send to neighbors
        std::vector<int> n_recv(recv_neig.size());

        std::vector<MPI_Request> sreq(send_neig.size());
        std::vector<MPI_Request> rreq(recv_neig.size());
        std::vector<MPI_Status> ssta(send_neig.size());
        std::vector<MPI_Status> rsta(recv_neig.size());

        for (int n = 0; n < (int)send_neig.size(); ++n) {
          MPI_Isend(n_send.data() + n, 1, MPI_INT, send_neig.at(n), 2000,
                    comm_col, sreq.data() + n);
        }
        for (int n = 0; n < (int)recv_neig.size(); ++n) {
          MPI_Irecv(n_recv.data() + n, 1, MPI_INT, recv_neig.at(n), 2000,
                    comm_col, rreq.data() + n);
        }
        MPI_Waitall(send_neig.size(), sreq.data(), ssta.data());

        int recvs = 0;

        /*
        std::cout << localRank << " "  << send_neig.at(0) << " "
                                        << ((send_neig.size() >=
        2)?send_neig.at(1):-1) << " "
                                        << ((send_neig.size() >=
        3)?send_neig.at(2):-1) << " "
                                        << send_neig.size() << " | "
                                        << recv_neig.at(0) << " "
                                        << ((recv_neig.size() >=
        2)?recv_neig.at(1):-1) << " "
                                        << ((recv_neig.size() >=
        3)?recv_neig.at(2):-1) << " "
                                        << recv_neig.size() << " | "
                                        << old_border.at(0) << " " <<
        old_border.at(1) << " | "
                                        << local_borders.at(0) << " " <<
        local_borders.at(1)
                                        << std::endl;
        */

        std::vector<std::vector<double>> recv_vec(recv_neig.size());
        while (recvs < (int)recv_neig.size()) {
          int idx;
          MPI_Waitany(recv_neig.size(), rreq.data(), &idx, rsta.data());
          recv_vec.at(idx).resize(n_recv.at(idx));
          recvs++;
        }

        // transfer points from old domains to new domains
        for (int n = 0; n < (int)send_neig.size(); ++n) {
          MPI_Isend(send_vec.at(n).data(), send_vec.at(n).size(), MPI_DOUBLE,
                    send_neig.at(n), 3000, comm_col, sreq.data() + n);
        }
        for (int n = 0; n < (int)recv_neig.size(); ++n) {
          MPI_Irecv(recv_vec.at(n).data(), recv_vec.at(n).size(), MPI_DOUBLE,
                    recv_neig.at(n), 3000, comm_col, rreq.data() + n);
        }
        MPI_Waitall(send_neig.size(), sreq.data(), ssta.data());

        recvs = 0;
        while (recvs < (int)recv_neig.size()) {
          int idx;
          MPI_Waitany(recv_neig.size(), rreq.data(), &idx, rsta.data());
          for (int p = 0; p < (int)recv_vec.at(idx).size(); p += 5) {
            ALL::Point<double> tmp(3, recv_vec.at(idx).data() + p,
                                  recv_vec.at(idx).at(p + 3),
                                  (long)(recv_vec.at(idx).at(p + 4)));
            points.push_back(tmp);
          }
          recvs++;
        }
        n_points = points.size();

        MPI_Comm_free(&comm_col);
        break;
      }
      default:
        break;
      }

      if (i_loop <= 100 || i_loop % OUTPUT_INTV == 0) {
#ifdef ALL_VORONOI_ACTIVE
        if (chosen_method != ALL::LB_t::VORONOI) {
#endif
#ifdef ALL_VTK_OUTPUT
          // if (localRank == 0)
          //  std::cout << "creating vtk outlines output" << std::endl;
          if (chosen_method != ALL::LB_t::FORCEBASED)
            lb_obj.printVTKoutlines(output_step);
          // if (localRank == 0)
          //  std::cout << "creating vtk vertices output" << std::endl;
          if (chosen_method == ALL::LB_t::FORCEBASED)
            lb_obj.printVTKvertices(output_step);
          // if (localRank == 0)
          //  std::cout << "creating points output" << std::endl;
          // if (i_loop == 0)
          print_points(points, i_loop / OUTPUT_INTV, chosen_method, cart_comm);
#endif
          output_step++;
#ifdef ALL_VORONOI_ACTIVE
        } else {
#ifdef ALL_VTK_OUTPUT
          print_points(points, i_loop / OUTPUT_INTV, chosen_method, cart_comm);
#endif
        }
#endif
      }

      // calculate quality parameters
      if (!weighted_points) {
        n_local = (double)n_points;
      } else {
        n_local = 0.0;
        for (auto p = points.begin(); p != points.end(); ++p)
          n_local += p->getWeight();
      }

      double n_total_points;
      MPI_Allreduce(&n_local, &n_total_points, 1, MPI_DOUBLE, MPI_SUM,
                    cart_comm);

      MPI_Allreduce(&n_local, &n_total, 1, MPI_DOUBLE, MPI_SUM, cart_comm);
      avg_work = n_total / (double)n_ranks;
      MPI_Allreduce(&n_local, &n_min, 1, MPI_DOUBLE, MPI_MIN, cart_comm);
      MPI_Allreduce(&n_local, &n_max, 1, MPI_DOUBLE, MPI_MAX, cart_comm);
      d_min = n_min / avg_work;
      d_max = n_max / avg_work;
      d_ratio = (d_max - d_min) / (d_max + d_min);

      double loc_diff, tot_diff;
      loc_diff = (n_local - avg_work) * (n_local - avg_work);
      MPI_Reduce(&loc_diff, &tot_diff, 1, MPI_DOUBLE, MPI_SUM, 0, cart_comm);

      d_min = n_min / avg_work;
      d_max = n_max / avg_work;
      d_ratio = (d_max - d_min) / (d_max + d_min);

      if (localRank == 0) {
        tot_diff = sqrt(tot_diff);
        std::ofstream of;
        if (!weighted_points)
          of.open("stddev.dat", std::ios::out | std::ios::app);
        else
          of.open("stddev_w.dat", std::ios::out | std::ios::app);
        of << i_loop + 1 << " " << tot_diff << std::endl;
        of << std::flush;
        of.close();

        if (!weighted_points)
          of.open("minmax.dat", std::ios::out | std::ios::app);
        else
          of.open("minmax_w.dat", std::ios::out | std::ios::app);
        of << i_loop + 1 << " " << d_min << " " << d_max << " " << d_ratio
           << " " << max_neighbors << std::endl;
        of << std::flush;
        of.close();
        if (i_loop % OUTPUT_INTV == 0) {
          std::cout << "d_min: " << d_min << std::endl;
          std::cout << "d_max: " << d_max << std::endl;
          std::cout << "d_ratio: " << d_ratio << std::endl;
          std::cout << std::flush;
        }
        if (d_ratio < min_ratio) {
          min_ratio = d_ratio;
          min_step = i_loop;
        }
      }

      // output of borders / contents
      for (int i = 0; i < n_ranks; ++i) {
        if (localRank == i) {
          std::ofstream of;
          if (!weighted_points)
            of.open("domain_data.dat", std::ios::out | std::ios::app);
          else
            of.open("domain_data_w.dat", std::ios::out | std::ios::app);
          of << (i_loop + 1) << " " << localRank << " ";

          if (chosen_method == ALL::LB_t::HISTOGRAM) {
            int *nNeighbors;
            std::vector<int> neighbors = lb_obj.getNeighbors();
            nNeighbors = neighbors.data();
            int n_neig = 0;
            for (int j = 0; j < 3; ++j) {
              n_neig += nNeighbors[j];
            }
            of << n_neig << " ";
          }
          /*
          for (int j = 0; j < sys_dim; ++j)
          {
              of << " " << local_coords[j]
                  << " " << lp[j] << " " << up[j] << " ";
          }
          */

          of << n_local << " ";

          of << std::endl;
          if (i == n_ranks - 1)
            of << std::endl;
          of.close();
          MPI_Barrier(cart_comm);
        } else
          MPI_Barrier(cart_comm);
      }

      if (std::abs(n_total_points - total_points) > 1e-6) {
        std::cout << std::setprecision(14) << n_total_points
                  << " != " << total_points << std::endl;
        return (-1);
      }

      if (i_loop == max_loop - 1) {
        // output of the last configuration after lb-step!
        lb_obj.setVertices(vertices);
        double vis_work = 0.0;
        for (int p = 0; p < n_points; ++p) {
          vis_work += points.at(p).getWeight();
        }
        lb_obj.setWork(vis_work);
#ifdef ALL_VORONOI_ACTIVE
        if (chosen_method != ALL::LB_t::VORONOI) {
#endif
#ifdef ALL_VTK_OUTPUT
          if (chosen_method != ALL::LB_t::FORCEBASED) {
            lb_obj.printVTKoutlines(output_step);
            double width[3];
            double volume = 1.0;

            for (int i = 0; i < 3; ++i) {
              width[i] = vertices.at(1)[i] - vertices.at(0)[i];
              volume *= width[i];
            }

            double surface = 2.0 * width[0] * width[1] +
                             2.0 * width[1] * width[2] +
                             2.0 * width[0] * width[2];

            for (int r = 0; r < n_ranks; ++r) {
              if (r == localRank) {
                std::ofstream of;
                of.open("domain_data.dat", std::ios::out | std::ios::app);
                of << localRank << " " << width[0] << " " << width[1] << " "
                   << width[2] << " " << vis_work << " " << volume << " "
                   << surface << " " << (surface / volume) << " " << std::endl;
                of.close();
              }
              MPI_Barrier(cart_comm);
            }
          }
          if (chosen_method == ALL::LB_t::FORCEBASED)
            lb_obj.printVTKvertices(output_step);
#endif
          output_step++;
#ifdef ALL_VORONOI_ACTIVE
        } else {
#ifdef ALL_VTK_OUTPUT
          print_points(points, i_loop / OUTPUT_INTV, chosen_method, cart_comm);
#endif
        }
#endif
      }
    }

    // create binary output (e.g. for HimeLB)
    // format:
    // n_ranks integers (number of particles per domain)
    // n_ranks integers (offset in file)
    // <n_particles[0]> * (long + 3*double) double values (particle positions)

    // open file
    int err;
    MPI_File outfile;
    err =
        MPI_File_open(MPI_COMM_WORLD, "particles_output.bin",
                      MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &outfile);

    // write out the number of particles on each process
    MPI_File_write_at(outfile, (MPI_Offset)(localRank * sizeof(int)), &n_points,
                      1, MPI_INT, MPI_STATUS_IGNORE);

    // compute the offset for each process
    long n_p = (long)n_points;
    long offset;

    MPI_Scan(&n_p, &offset, 1, MPI_LONG, MPI_SUM, cart_comm);
    offset -= n_p;
    // size of data for a single block
    long block_size = sizeof(long); // + 3 * sizeof(double);

    offset *= block_size;
    offset += n_ranks * (sizeof(int) + sizeof(long));

    MPI_File_write_at(
        outfile, (MPI_Offset)(n_ranks * sizeof(int) + localRank * sizeof(long)),
        &offset, 1, MPI_LONG, MPI_STATUS_IGNORE);

    for (int n = 0; n < n_points; ++n) {
      long id = points.at(n).get_id();
      double loc[3];
      loc[0] = points.at(n)[0];
      loc[1] = points.at(n)[1];
      loc[2] = points.at(n)[2];
      MPI_File_write_at(outfile, (MPI_Offset)((offset + n * block_size)), &id,
                        1, MPI_LONG, MPI_STATUS_IGNORE);
      /*
      MPI_File_write_at(outfile,
                        (MPI_Offset)((offset+n*block_size)+sizeof(long)),
                        loc,
                        3,
                        MPI_DOUBLE,
                        MPI_STATUS_IGNORE);
      */
    }

    MPI_File_close(&outfile);

    /*
    // output of borders / contents
    if (n_ranks < 216)
    {
    for (int i = 0; i < n_ranks; ++i)
    {
    if (localRank == i)
    {
    std::ofstream of;
    if (!weighted_points)
    of.open("domain_data.dat", std::ios::out | std::ios::app);
    else
    of.open("domain_data_w.dat", std::ios::out | std::ios::app);
    of << 0 << " " << localRank << " ";

    for (int j = 0; j < sys_dim; ++j)
    {
    of << " " << local_coords[j] << " "
    << lp[j] << " "
    << up[j] << " ";
    }

    of << " " << n_local << " ";

    of << std::endl;
    if (i == n_ranks - 1) of << std::endl;
    of.close();
    MPI_Barrier(cart_comm);
    }
    else
    MPI_Barrier(cart_comm);
    }
    }
     */

    delete transfer;

    if (localRank == 0) {
      std::cout << std::endl;
      std::cout << "min_ratio: " << min_ratio << std::endl;
      std::cout << "min_step: " << min_step << std::endl;
      std::cout << std::endl;
    }
    MPI_Barrier(cart_comm);
    MPI_Finalize();
    return EXIT_SUCCESS;
  } catch (ALL::CustomException &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
}
// vim: sw=2 et ts=2
