/*
Copyright 2018-2020 Rene Halver, Forschungszentrum Juelich GmbH, Germany
Copyright 2018-2020 Godehard Sutmann, Forschungszentrum Juelich GmbH, Germany
Copyright 2019-2020 Stephan Schulz, Ruhr-Universität Bochum, Germany

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

#include "../include/ALL.hpp"
#include "../include/ALL_CustomExceptions.hpp"
#include "../include/ALL_Defines.h"
#include <cassert>
#include <cstring>
#include <vector>

static int ALL_errno;
static const char* ALL_errdesc;

typedef ALL::ALL<double, double> ALL_t;

#ifdef ALL_FORTRAN_ERROR_ABORT
#define ALL_try
#define ALL_catch
#else
#define ALL_try try {
#define ALL_catch } catch (ALL::CustomException &e) {  \
  ALL_errno = e.get_error_id();                   \
  ALL_errdesc = e.what();                         \
}
#endif

/*
This interface is unstable and subject to change! It is only used for the
Fortran interface. Several safety checks are only done on the Fortran side. Use
the C++ API or Fortran interfaces
*/

extern "C" {
// wrapper to create a ALL::ALL<double,double> object (should be the
// most commonly used one)
ALL_t *all_init_c(ALL::LB_t method, const int dim, double gamma) {
  ALL_try
  return new ALL::ALL<double, double>(method, dim, gamma);
  ALL_catch
  return NULL;
}

// delete ALL object
void all_finalize_c(ALL_t *all_obj) {
  delete all_obj;
  // todo(s.schulz): does this throw if all_obj is already deleted? can we detect it?
}

void all_set_gamma_c(ALL_t *all_obj, double gamma) {
  ALL_try
  all_obj->setGamma(gamma);
  ALL_catch
}

// set grid parameters
void all_set_proc_grid_params_c(ALL_t *all_obj, int nloc,
                                int *loc, int nsize, int *size) {
  ALL_try
  int dim = all_obj->getDimension();
  if (dim != nloc)
    throw ALL::InvalidArgumentException(
        __FILE__, __func__, __LINE__,
        "Length of index array does not match dimension");
  if (dim != nsize)
    throw ALL::InvalidArgumentException(
        __FILE__, __func__, __LINE__,
        "Length of size array does not match dimension");
  std::vector<int> vloc(loc, loc+nloc);
  std::vector<int> vsize(size, size+nsize);
  all_obj->setProcGridParams(vloc, vsize);
  ALL_catch
}

// wrapper to set the minimum domain size
void all_set_min_domain_size_c(ALL_t *all_obj, int dim,
                               double *domain_size) {
  ALL_try
  if (all_obj->getDimension() != dim)
    throw ALL::InvalidArgumentException(
        __FILE__, __func__, __LINE__,
        "Length of array does not match dimension");
  std::vector<double> t_domain_size(domain_size, domain_size+dim);
  all_obj->setMinDomainSize(t_domain_size);
  ALL_catch
}

// wrapper to set the work (scalar only for the moment)
void all_set_work_c(ALL_t *all_obj, double work) {
  ALL_try
  all_obj->setWork(work);
  ALL_catch
}

// set multidimensional work
void all_set_work_multi_c(ALL_t *all_obj, double *work, int dim) {
  ALL_try
  std::vector<double> t_work(work, work+dim);
  all_obj->setWork(t_work);
  ALL_catch
}

// wrapper to set the vertices (using an array of double values and dimension)
void all_set_vertices_c(ALL_t *all_obj, const int n,
                        const int dim, const double *vertices) {
  ALL_try
  std::vector<ALL::Point<double>> points(n, ALL::Point<double>(dim));
  if (dim != all_obj->getDimension())
    throw ALL::PointDimensionMissmatchException(
        __FILE__, __func__, __LINE__,
        "Dimension of ALL::Points in input vector do not match dimension of "
        "ALL "
        "object.");
  for (int i = 0; i < n; ++i) {
    for (int d = 0; d < dim; ++d) {
      points.at(i)[d] = vertices[i * dim + d];
    }
  }
  all_obj->setVertices(points);
  ALL_catch
}

// wrapper to set the communicator
void all_set_communicator_c(ALL_t *all_obj, MPI_Fint fcomm) {
  ALL_try
  MPI_Comm cComm = MPI_Comm_f2c(fcomm);
  all_obj->setCommunicator(cComm);
  ALL_catch
}

void all_set_sys_size_c(ALL_t *all_obj, double *size, int dim) {
  ALL_try
  std::vector<double> t_size(size, size+dim);
  all_obj->setSysSize(t_size);
  ALL_catch
}

// wrapper to setup routine
void all_setup_c(ALL_t *all_obj) {
  ALL_try
  all_obj->setup();
  ALL_catch
}

void all_set_method_data_histogram_c(ALL_t *all_obj, int *nbins) { all_obj->setMethodData(nbins); }

// wrapper to call balance routine
void all_balance_c(ALL_t *all_obj) {
  ALL_try
  all_obj->balance();
  ALL_catch
}

void all_get_gamma_c(ALL_t *all_obj, double *gamma) { *gamma =  all_obj->getGamma(); }

// wrapper to get number of new vertices
void all_get_number_of_vertices_c(ALL_t *all_obj, int *n_vertices) {
  ALL_try
  *n_vertices = (all_obj->getVertices()).size();
  ALL_catch
}

// wrapper to return new vertices
void all_get_vertices_c(ALL_t *all_obj, int n_vertices,
                               double *vertices) {
  ALL_try
  std::vector<ALL::Point<double>> tmp_vertices = all_obj->getVertices();
  int dimension = all_obj->getDimension();
  assert(n_vertices = tmp_vertices.size());
  for (int i = 0; i < n_vertices; ++i) {
    for (int j = 0; j < dimension; ++j) {
      vertices[i * dimension + j] = tmp_vertices.at(i)[j];
    }
  }
  ALL_catch
}

// set process tag
void all_set_proc_tag_c(ALL_t *all_obj, int tag) {
  ALL_try
  all_obj->setProcTag(tag);
  ALL_catch
}

// wrapper to return new vertices
void all_get_prev_vertices_c(ALL_t *all_obj, int n_vertices,
                             double *prevVertices) {
  ALL_try
  std::vector<ALL::Point<double>> tmp_vertices = all_obj->getPrevVertices();
  int dimension = all_obj->getDimension();
  assert(n_vertices = tmp_vertices.size());
  for (int i = 0; i < n_vertices; ++i) {
    for (int j = 0; j < dimension; ++j) {
      prevVertices[i * dimension + j] = tmp_vertices.at(i)[j];
    }
  }
  ALL_catch
}

// get currently set dimension
void all_get_dimension_c(ALL_t *all_obj, int *dim) {
  ALL_try
  *dim = all_obj->getDimension();
  ALL_catch
}

void all_get_length_of_work_c(ALL_t *all_obj, int *length) {
  ALL_try
  std::vector<double> work;
  all_obj->getWork(work);
  *length = work.size();
  ALL_catch
}

void all_get_work_c(ALL_t *all_obj, double *work) {
  ALL_try
  *work = all_obj->getWork();
  ALL_catch
}

void all_get_work_array_c(ALL_t *all_obj, double *work, int length) {
  ALL_try
  std::vector<double> all_work;
  all_obj->getWork(all_work);
  assert((int)all_work.size() == length);
  memcpy(work,&all_work[0],length*sizeof(*work));
  ALL_catch
}

void all_get_number_of_neighbors_c(ALL_t *all_obj, int *count) {
  ALL_try
  std::vector<int> neighbors = all_obj->getNeighbors();
  *count = neighbors.size();
  ALL_catch
}

void all_get_neighbors_c(ALL_t *all_obj, int *neighbors, int count) {
  ALL_try
  std::vector<int> all_neighbors = all_obj->getNeighbors();
  assert((int)all_neighbors.size() == count);
  memcpy(neighbors,&all_neighbors[0],count*sizeof(*neighbors));
  ALL_catch
}

#ifdef ALL_VTK_OUTPUT
// print VTK outlines
void all_print_vtk_outlines_c(ALL_t * all_obj known_unused, int known_unused step) {
  // todo(s.schulz): Should this function even be declared if ALL_VTK_OUTPUT isn't?
  ALL_try
  all_obj->printVTKoutlines(step);
  ALL_catch
}

void all_print_vtk_vertices_c(ALL_t * all_obj known_unused, int known_unused step) {
  // todo(s.schulz): Should this function even be declared if ALL_VTK_OUTPUT isn't?
  ALL_try
  all_obj->printVTKvertices(step);
  ALL_catch
}
#endif

// retrieve error number
int all_errno_c(void) { return ALL_errno; }

// set error number
void all_reset_errno_c() {
  ALL_errno = 0;
  ALL_errdesc = NULL;
}

// retrieve error description
void all_errdesc_c(char *description, size_t len) {
  if(ALL_errdesc) {
    strncpy(description, ALL_errdesc, len);
  } else {
    strncpy(description, "No error", len);
  }
  // Use space padding instead of zero padding
  // WARNING: The result is no longer a zero terminated string!
  size_t msg_len = strlen(description);
  memset(description+msg_len, ' ', len-msg_len);
}
}//extern "C"

// VIM modeline for indentation
// vim: sw=2 et ts=2
